
import numpy as np
import pandas as pd

from scipy.interpolate import RectBivariateSpline
import scipy.constants
from typing import Union

from ..pvt.pvt import get_hydrostatic_P, get_pvt

'''Some global parameters'''
G       = scipy.constants.g   #9.81 m/s2 gravity acceleration
BAR2PA  = scipy.constants.bar #10**5 Going from bars to Pascal: 1 bar = 10**5 Pascal
SHMIN_FAC = 0.1695

'''Global names. Should be replaced by a global class containing these names? TODO'''
SF_DEPTH_NAME = 'sf_depth_msl'
DEPTH_NAME    = 'depth_msl'
SHMIN_NAME    = 'Shmin'
TEMP_NAME     = 'temp'
RHO_NAME      = 'rho'
MAX_PRESSURE_NAME = 'max_pressure'


def _get_shmin(well_header: dict, pt_df_in: pd.DataFrame) -> pd.DataFrame:
    '''Compute Shmin based on an empirical formula
    '''
    pt_df = pt_df_in.copy()
    sf_depth = well_header['sf_depth_msl']                                 #Depth to sea floor
    sf_pressure = np.interp(sf_depth, pt_df['depth_msl'], pt_df['hs_p'])   #Pressure at sea floor
    pt_df[SHMIN_NAME] = pt_df['hs_p']                                         #Initialize with hydrostatic pressure
    shmin_query = pt_df.query('depth_msl>=@sf_depth')                      #Take only values below sea floor
    pt_df.loc[pt_df['depth_msl']>=sf_depth, SHMIN_NAME] = sf_pressure + (shmin_query['depth_msl']-sf_depth)*SHMIN_FAC

    return pt_df


def _get_rho_in_pressure_column(pt_df_in: pd.DataFrame, colname_p: str, colname_rho: str, get_rho: callable) -> pd.DataFrame:
    ''' Calculate the water density
        in the CO2 column with its calculated pressure and temperatur:
    '''
    pt_df = pt_df_in.copy()
    temps = pt_df['temp']
    pressures = pt_df[colname_p]
    rhos = []
    for t, p in zip(temps, pressures):
        rhos.append(get_rho(p,t)[0,0])

    pt_df[colname_rho] = rhos[:]

    return pt_df

def _get_max_pressure(pt_df_in: pd.DataFrame, max_pressure_pos: Union[dict, list, float, int], get_rho: callable) -> pd.DataFrame:
    '''
    Calculates downwards maximum pressure at a depth deeper than the input depth.
    I.e. the pressure at deeper locations such that witha CO2-column the pressure will be equal to Shmin at the "start/input"-depth
    Input:
      max_pressure_pos is a depth from where max pressure (wrt Shmin) is calculated. It can be a dict (my_well.barriers) a list of numbers or scalars.
      If my_well.barriers is given then it is calculated from the base of each barrier
    '''

    pt_df = pt_df_in.copy()
    if isinstance(max_pressure_pos, dict):   #Then max_pressure_pos is the same as barriers - and max pressure is calcualted for each barrier
        for idx, key in max_pressure_pos['barrier_name'].items():
            barr_depth = max_pressure_pos['bottom_msl'][idx]
            colname_p = f"{MAX_PRESSURE_NAME}_{key}"
            print(f"Calculating max pressure below barrier {key} from depth {barr_depth}")            
            p0 = np.interp(barr_depth, pt_df['depth_msl'], pt_df[SHMIN_NAME])
            pt_df = _integrate_pressure(pt_df, get_rho, barr_depth, p0, 'down', colname_p)
    elif isinstance(max_pressure_pos, (list, float, int)):
        if isinstance(max_pressure_pos, (float, int)): #Make it a list with one element
            barriers = [max_pressure_pos]
        for depth in max_pressure_pos:
            colname_p = f"{MAX_PRESSURE_NAME}_at_{int(depth)}" 
            print(f"Calculating max pressure from depth {depth}")            
            p0 = np.interp(float(depth), pt_df['depth_msl'], pt_df[SHMIN_NAME])
            pt_df = _integrate_pressure(pt_df, get_rho, float(depth), p0, 'down', colname_p)

    return pt_df


def _integrate_pressure(pt_df_in: pd.DataFrame, get_rho: callable, reference_depth: float, reference_pressure: float, direction: str, colname_p: str) -> pd.DataFrame:
    '''
    Simple integration to find pressure
    Starting point: reference_depth at reference pressure.  Reference temperature is found in pt_df
                    Then iterate upwards (up) or downwards (down) to subtract or add pressure.
                    Recalculates rho at each step.
    '''
    pt_df = pt_df_in.copy()

    ## Initialization
    #Presure at MSL
    p_msl = scipy.constants.atm/BAR2PA  #1.01325 bar Pressure at MSL

    #New columns needed in DataFrame
    if colname_p not in pt_df.columns:
        pt_df[colname_p] = np.nan

    colname_rho = colname_p+'_rho'
    if colname_rho not in pt_df.columns:
        pt_df[colname_rho] = np.nan

    #Take out the part of the dataframe that is either below or above the reference depth.
    if direction.lower() == "up":
        query = pt_df.query('depth_msl<=@reference_depth')
        sign  = -1
    elif direction.lower() == "down":
        query = pt_df.query('depth_msl>=@reference_depth')
        sign = 1
    else:
        print(f"ERROR: Not a valid direction. It should be 'up' or 'down'. You wrote {direction.lower()}")

    ###Do the calculations
    #Starting pressure and depth
    p  = reference_pressure
    z0 = reference_depth

    #Loop through the rows in the dataframe - upwards from the bottom of downwards from the top
    for z_idx, row in query[::sign].iterrows():
        t = row['temp']
        z = row['depth_msl']

        rho = get_rho(p,t)[0,0]        
        dz = z - z0
        p += (rho*G*dz)/BAR2PA   #/1e-5

        #Ensure you do not get negative pressures
        if p<p_msl:
            p=p_msl
                
        pt_df.loc[z_idx, colname_p] = p
        pt_df.loc[z_idx, colname_rho] = rho
        z0 = z

    return pt_df
##################################################################################   




def compute_CO2_pressures(well_header: dict, p_init: dict, base_co2: float, *, pvt_path: str, max_pressure_pos: Union[dict, list, float, int] = None) -> pd.DataFrame:
    '''The pressure and density for H2O and CO2  along the columns are calculated using an approximate integration
        Hydrostatic pressure - caculatong downwards from msl
        Pressure and density assuming a water column - starting at top reservoir and the given overpressure RP
        Pressure and density assuming a CO2   column - starting at CO2-reference level (CO2-column could start below top reservoir) and the given overpressure RP
        Shmin

        Input:
        max_pressure_pos is a depth from where max pressure (wrt Shmin) is calculated. It can be a dict (my_well.barriers) a list of numbers or scalars.
                If my_well.barriers is given then it is calculated from the base of each barrier

        Columns are
            depth_msl        temp          hs_p         RHOH20            Shmin   RPx_h20        RPx_h20_rho     RPx_co2                            RPx_co2_rho     RPx_h2o_rho_in_co2_column
            depth below msl   temperature   hydrostatic  water density             hs_p+          densities at    Pressure given a CO2 column        Corresponding   water densitites if we are
            depth ref for                   pressure     at hydrostatic            overpressure   RPx_h20         and overpressure RP                CO2 densities   in a CO2-column.
            all other values                column       pressure                  RPx
        
        

        ---> The RP are all RP-input + hydrostatic_pressure. Hence if e.g. RP1 is hydrostatic pressure RP1-columns and hydrostatic_pressure-columns are identical.
    '''

    #Get PVT-input
    t_vec, p_vec, rho_co2_vec, rho_h2o_vec = get_pvt(pvt_path)
    p_msl = scipy.constants.atm/BAR2PA  #1.01325 bar Pressure at MSL

    #An interpolator. Used later to look-up densities given pressure p and temperature t
    get_rho_h2o = RectBivariateSpline(p_vec, t_vec, rho_h2o_vec)
    get_rho_co2 = RectBivariateSpline(p_vec, t_vec, rho_co2_vec)


    #Get hydrostatic pressure from msl and downwards
    #A DataFame with the columns depth, temp(depth) hydrostatic_pressure(depth), RHOH2O
    #Going all the way to (well_td_rkb - well_rkb) + 200 m
    pt_df = get_hydrostatic_P(well_header, pvt_path=pvt_path)

    #Get the shmin -values
    pt_df = _get_shmin(well_header, pt_df)

    #Get max pressure - to not get pressures above Shmin below barriers or at a given depth.
    if max_pressure_pos is not None:
        pt_df =  _get_max_pressure(pt_df, max_pressure_pos, get_rho_co2)
           

    #Find hydrostatic pressure and temperature at reference depth - typically at top reservoir
    #This is the same depth the over-pressure cases are set to
    ref_z = p_init['depth_msl']
    
    print(f"Reference depth and pressure scenarios to be calculated there: \n     p_init {p_init}")
    print(f"Top reservoir {ref_z}")
    print(f"From where there is CO2 - hence affecting pressure upwards: base co2: {base_co2}")
    print(f"Hence CO2 column in the reservoir is then {base_co2 - ref_z} m")

    #Integrate from reference depth and upwards given different starting pressures
    for key, value in p_init.items():                 #depth_msl, RP1, RP2, hydrostatic_pressure
        if key == 'depth_msl':                        #No calculations - just report
            print(f"Reference depth: {value}")
        else:
            rp = key                                  #pressure name: RP1, RP2 etc
            p0 = value                                #Initial pressure value set for RP1, RP2 etc

            #Water
            water_p_colname = rp+'_h2o'
            pt_df = _integrate_pressure(pt_df, get_rho_h2o, ref_z, p0, 'up', water_p_colname)
            pt_df = _integrate_pressure(pt_df, get_rho_h2o, ref_z, p0, 'down', water_p_colname)

            #CO2
            #Find water pressure at base_co2
            p0 = np.interp(base_co2, pt_df['depth_msl'], pt_df[water_p_colname])
            co2_p_colname  = rp+'_co2'
            pt_df = _integrate_pressure(pt_df, get_rho_co2, base_co2, p0, 'up', co2_p_colname)

            #We need the density for water given the CO2-pressures, too
            pt_df = _get_rho_in_pressure_column(pt_df, co2_p_colname, f"{rp}_h2o_rho_in_co2_column", get_rho_h2o)


    return pt_df

