
import numpy as np
import pandas as pd

from scipy.interpolate import RectBivariateSpline
import scipy.constants

from ..pvt.pvt import get_hydrostatic_P, get_pvt

'''Some global parameters'''
G       = scipy.constants.g   #9.81 m/s2 gravity acceleration
BAR2PA  = scipy.constants.bar #10**5 Going from bars to Pascal: 1 bar = 10**5 Pascal
REGR_A  = -0.000116           #Intercept term in regression equation for the proxy. Consider as a input
REGR_B  = 0.000002725         #Inclination-term in regression equation for the proxy

def compute_CO2_pressures(well_header: dict, p_init: dict, base_co2: float, *, pvt_path: str) -> pd.DataFrame:
    '''The pressure and density for H2O and CO2  along the columns are calculated using an approximate integration
        Hydrostatic pressure - caculatong downwards from msl
        Pressure and density assuming a water column - starting at top reservoir and the given overpressure RP
        Pressure and density assuming a CO2   column - starting at CO2-reference level (CO2-column could start below top reservoir) and the given overpressure RP
        Shmin

        Columns are
            depth_msl        temp          hs_p         RHOH20            Shmin   RPx_h20        RPx_h20_rho     RPx_co2                            RPx_co2_rho     RPx_h2o_rho_in_co2_column
            depth below msl   temperature   hydrostatic  water density             hs_p+          densities at    Pressure given a CO2 column        Corresponding   water densitites if we are
            depth ref for                   pressure     at hydrostatic            overpressure   RPx_h20         and overpressure RP                CO2 densities   in a CO2-column.
            all other values                column       pressure                  RPx
        
        

        ---> The RP are all RP-input + hydrostatic_pressure. Hence if e.g. RP1 is hydrostatic pressure RP1-columns and hydrostatic_pressure-columns are identical.
    '''
    # well_header = self.header
    # p_init      = self.reservoir_P
    # base_co2    = self.co2_datum

    t_vec, p_vec, rho_co2_vec, rho_h2o_vec = get_pvt(pvt_path)
    p_msl = scipy.constants.atm/scipy.constants.bar  #1.01325 bar Pressure at MSL

    get_rho_h2o = RectBivariateSpline(p_vec, t_vec, rho_h2o_vec)
    get_rho_co2 = RectBivariateSpline(p_vec, t_vec, rho_co2_vec)

    # TODO(hzh): no return value?????
    def get_rho(phase, p, t):
        if phase=='co2':
            rho  = get_rho_co2(p,t)[0, 0]
        elif phase == 'h2o':
            rho  = get_rho_h2o(p,t)[0, 0]
        return rho

    #Get hydrostatic pressure from msl and downwards
    pt_df = get_hydrostatic_P(well_header, pvt_path=pvt_path)  #Calculate depth, temp(depth) hydrostatic_pressure(depth), H2ORHO(depth for hydrostatic pressure)

    #Compute Shmin
    wd = well_header['sf_depth_msl']                           #Depth to sea floor
    wp_ml = np.interp(wd, pt_df['depth_msl'], pt_df['hs_p'])   #Pressure at sea floor
    pt_df['Shmin'] = pt_df['hs_p']                             #Initialize with hydrostatic pressure
    shmin_query = pt_df.query('depth_msl>=@wd')                #Take only values below sea floor
    pt_df.loc[pt_df['depth_msl']>=wd, 'Shmin'] = wp_ml + (shmin_query['depth_msl']-wd)*0.1695

    #Find hydrostatic pressure and temperature at reference depth - typically at top reservoir
    ref_z = p_init['depth_msl']
    ref_z_hsp = np.interp(ref_z, pt_df['depth_msl'], pt_df['hs_p'])
    ref_z_temp = np.interp(ref_z, pt_df['depth_msl'], pt_df['temp'])

    #Integrate from reference depth and upwards given different starting pressures
    for key, value in p_init.items():                 #depth_msl, RP1, RP2, hydrostatic_pressure
        if key == 'depth_msl':                        #No calculations - just report
            print(f"Reference depth: {value}")
        else:
            phases_zval = [('h2o', ref_z), ('co2', base_co2)]
            rp = key                                  #pressure name
            p0 = value                                #initial pressure value

            for phase, ref_depth in phases_zval:
                    #crate column names in dataframe
                    colname_p   = rp+'_'+phase
                    colname_rho = rp+'_'+phase+'_rho'

                    #create empty columns for each phase
                    pt_df[colname_p] = np.nan         #pressure
                    pt_df[colname_rho] = np.nan       #density

                    #Assign pressure-value at reference point if depth value exists
                    pt_df.loc[pt_df.depth_msl == ref_depth, colname_p] = p0
                    
                    #Create tuple lists for sections for integration
                    segments_list = []

                    #Split vector in two sections for iterations: upwards(-1) and downwards(+1)
                    above_ref_query = pt_df.query('depth_msl<=@ref_depth')
                    segments_list.append((above_ref_query, -1))
                    
                    #below segment only evaluated for water
                    if phase == 'h2o':
                            below_ref_query = pt_df.query('depth_msl>@ref_depth')
                            segments_list.append((below_ref_query, 1))

                    #if CO2 - then p0 is updated to H2O pressure at CO2 datum
                    if phase == 'co2':
                            p0 = np.interp(ref_depth, pt_df['depth_msl'], pt_df[rp+'_h2o'])
                            #Update pressure-value at reference point if depth value exists (is it really needed?)
                            pt_df.loc[pt_df.depth_msl == ref_depth, colname_p] = p0

                    #interpolate temperature to given value
                    t0 = np.interp(ref_z, pt_df['depth_msl'], pt_df['temp'])

                    #Update temperature-value and density-value at reference point if depth value exists (is it really needed?)
                    pt_df.loc[pt_df.depth_msl == ref_depth, colname_rho] = get_rho(phase, p0, t0)    # TODO(hzh): no return value?????
                    pt_df.loc[pt_df.depth_msl == ref_depth, colname_p]   = p0

                    for segment, sign in segments_list:
                            p  = p0
                            z0 = ref_depth

                            for z_idx, row in segment[::sign].iterrows():
                                    t = row['temp']
                                    z = row['depth_msl']

                                    if phase == 'h2o':
                                            rho = get_rho_h2o(p, t)[0,0]
                                    else:
                                            rho = get_rho_co2(p, t)[0,0]

                                    dz = z - z0
                                    p += (rho*G*dz)/scipy.constants.bar    #1e-5
                                    if p<p_msl:
                                            p=p_msl
                                    
                                    pt_df.loc[z_idx, colname_p] = p
                                    pt_df.loc[z_idx, colname_rho] = rho
                                    z0 = z

            #We need the density for water given the CO2-pressures, too
            temps = pt_df['temp']
            pressures = pt_df[f'{rp}_co2']
            rhos = []
            for t, p in zip(temps, pressures):
                rhos.append(get_rho_h2o(p,t)[0,0])

            pt_df[f"{rp}_h2o_rho_in_co2_column"] = rhos[:]

    #pt_df.to_csv("pt.csv")
    return pt_df
