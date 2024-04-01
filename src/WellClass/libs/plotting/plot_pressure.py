
import pandas as pd
import matplotlib.pyplot as plt 
import numpy as np

from typing import Union
from ..well_pressure import Pressure

'''Global names. Should be replaced by a global class containing these names? TODO'''
SHMIN_NAME        = 'Shmin'
SF_DEPTH_NAME     = 'sf_depth_msl'
DEPTH_NAME        = 'depth_msl'
MAX_PRESSURE_NAME = 'max_pressure'
RHO_NAME          = 'rho'

def plot_pressure(my_pressure: Pressure, geology_dict: dict, barriers: dict, ax=None, 
                  plot_HSP:bool = True, 
                  plot_RP:bool = True, 
                  plot_MSAD:bool = True, 
                  plot_maxP:bool = True):
    """
    pressure vs depth
    Takes the Pressure class, geology and barriers tables
    plot_HSP to plot brine hydrostatci pressure
    plot_RP to plot reservoir pressure scenarios
    plot_MSAD top plot Minimum Safatey Abandonement pressure for each reservoir scenario
    plot_maxP to plot max pressurization at a given deth
    """

    if ax is None:
            fig, ax = plt.subplots()
        
#     if not hasattr(my_pressure, "pressure_CO2"):
#         my_pressure._compute_CO2_pressures()

    pt_df = my_pressure.pressure_CO2

    well_header = my_pressure.header
    sf_depth    = well_header['sf_depth_msl']
    
    barriers_df = pd.DataFrame(barriers)
    geology_df  = pd.DataFrame(geology_dict)
    
    #define plot spatial references
    base_deepest_rsrv = geology_df[geology_df.reservoir_flag]['base_msl'].max()
    ymax = max([base_deepest_rsrv,my_pressure.co2_datum])+100
    xmax = pt_df.query('depth_msl>@ymax')['Shmin'].iloc[0]
    xmin = 0
    
    # Draw cement plugs
    for idx, row in barriers_df.iterrows():
        barrier = ax.axhspan(row['top_msl'], row['bottom_msl'], color='lightgray', zorder=-20, label = 'cement plug')

    #Plot hydrostatic pressure gradient
    if plot_HSP:
        pt_df.plot(x='hs_p', y='depth_msl', ax=ax, label='$p_{hs}$', color='steelblue', lw = 0.75)

    #Plot minimum horizontal stress
    pt_df.plot(x='Shmin', y='depth_msl', ax=ax, label='$\sigma_{h min}$', color='k', lw = 0.75)

    #Plot fluid pressure scenarios
    ls_list = ['solid','dashed','dashdot', 'dotted']


    if plot_RP:
        counter = 0


        for key in my_pressure.reservoir_P:
                if key != 'depth_msl':
                        RP_label = f'$CO_2$ at ({key})'

                        if plot_MSAD:
                            msad = my_pressure.MSAD[key]
                            z_MSAD = msad['z_MSAD']
                            p_MSAD = msad['p_MSAD']

                            if ~np.isnan(z_MSAD):
                                RP_label = f'{RP_label}\nMSAD = {z_MSAD:.2f} mTVDMSL'

                            ax.scatter(p_MSAD, z_MSAD, color='firebrick')

                             
                        pt_df.query('depth_msl>=@sf_depth').plot(x=key+'_h2o', y='depth_msl', ax=ax, label = '_nolegend_', color='steelblue', legend=False, lw = 0.75, ls=ls_list[counter])
                        pt_df.query('depth_msl>=@sf_depth').plot(x=key+'_co2', y='depth_msl', ax=ax, label = RP_label, color='firebrick', legend=True, lw = 0.75, ls=ls_list[counter])

                        counter+=1
                        counter = counter%(len(ls_list))  #If more cases than in ls_list then restart counter  


        #Plot MSAD
        if plot_MSAD:
            for key, value in my_pressure.MSAD.items():
                z_MSAD = value['z_MSAD']
                p_MSAD = value['p_MSAD']
        

    #Plot max pressure cases
    if plot_maxP:
        counter = 0
        for key in pt_df.columns: 
           if key.startswith(MAX_PRESSURE_NAME) and RHO_NAME not in key:
               colname = f"{key}"
               pt_df.query('depth_msl>=@sf_depth').plot(x=colname, y='depth_msl', ax=ax, label = colname, color='firebrick', legend=False, lw = 0.75, ls=ls_list[counter])
    
               counter+=1
               counter = counter%(len(ls_list))  #If more cases than in ls_list then restart counter  



    #Optimize legend
    ax.legend()
    handles, labels = ax.get_legend_handles_labels()  
    lgd = dict(zip(labels, handles))
    ax.legend(lgd.values(), lgd.keys())
    
    ax.set_xlim(xmin, xmax)
    ax.set_xlabel('pressure [bar]')
    ax.set_ylim(0, ymax)
    ax.invert_yaxis()

    if 'fig' in locals():
            return fig, ax
    else:
            return ax
