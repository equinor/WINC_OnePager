
import pandas as pd
import matplotlib.pyplot as plt 
import numpy as np
from itertools import cycle

from typing import Union
from ..well_pressure import Pressure


def plot_pressure(my_pressure: Pressure, 
                  geology_dict: dict = None,
                  barriers: dict = None, 
                  ax = None, 
                  plot_HSP: bool = True, 
                  plot_MSAD: bool = True, 
                  legend: bool = True,
                  plot_selected_scenarios: list = None,
                  plot_resrv: bool = True,  # Option to plot (p_resrv, z_resrv)
                  plot_fluid_contact: bool = True,
                  plot_delta_p: bool = True  ):  # Option to plot (p_fluid_contact, z_fluid_contact)
    
    """
    pressure vs depth
    Takes the Pressure class, geology and barriers tables
    plot_HSP to plot brine hydrostatic pressure
    plot_RP to plot reservoir pressure scenarios
    plot_MSAD to plot Minimum Safety Abandonment pressure for each reservoir scenario
    plot_maxP to plot max pressurization at a given depth
    plot_selected_scenarios to plot only a selection of scenarios
    """

    if ax is None:
        fig, ax = plt.subplots()

    # Use the collated_profiles DataFrame from the Pressure class
    # pt_df = my_pressure.collate_all_profiles()

    #List to store the reference detoth values to define the spatial references of the plot
    depth_values = []
    if geology_dict is not None:
        geology_df  = pd.DataFrame(geology_dict)
        base_deepest_rsrv = float(geology_df[geology_df.reservoir_flag]['base_msl'].max())
        depth_values.append(base_deepest_rsrv)
    
    #define plot spatial references
    if isinstance(my_pressure.z_fluid_contact, float):
        depth_values.append(my_pressure.z_fluid_contact)
    print(depth_values)
    # Define plot spatial references
    ymax = np.ceil(max(depth_values))
    xmax = my_pressure.init_curves.query('depth<=@ymax')['min_horizontal_stress'].max()
    xmin = 0

    # Draw cement plugs
    if barriers:
        barriers_df = pd.DataFrame(barriers)
        for idx, row in barriers_df.iterrows():
            barrier = ax.axhspan(row['top_msl'], row['bottom_msl'], color='lightgray', zorder=-20, label = 'cement plug')

    # Plot hydrostatic pressure gradient
    if plot_HSP:
        ax.plot(my_pressure.init_curves['hydrostatic_pressure'], my_pressure.init_curves['depth'], label='Hydrostatic Pressure', color='steelblue', lw=0.75)

    # Plot minimum horizontal stress
    ax.plot(my_pressure.init_curves['min_horizontal_stress'], my_pressure.init_curves['depth'], label='Min Horizontal Stress', color='black', lw=0.75)

    # Plot fluid pressure scenarios
    scenarios_summary = my_pressure.scenario_manager.get_scenarios_summary()

    # Define line styles for multiple scenarios
    line_styles = ['-', '--', '-.', ':']
    line_style_cycle = cycle(line_styles)

    
    if len(scenarios_summary) > 0:
        
        for scenario_name, scenario in scenarios_summary.iterrows():
            if plot_selected_scenarios and scenario_name not in plot_selected_scenarios:
                continue
            
            sc_type         = scenario['from_resrvr']
            sc_name         = scenario['name']
            sc_msad_p       = scenario['p_MSAD']
            sc_msad_z       = scenario['z_MSAD']
            sc_z_resrv         = scenario['z_resrv']    
            sc_p_resrv         = scenario['p_resrv']
            sc_z_fluid_contact = scenario['z_fluid_contact']
            sc_p_fluid_contact = scenario['p_fluid_contact']
            sc_delta_p      = scenario['p_delta']

            sc_curves = my_pressure.scenario_manager.scenarios[sc_name].init_curves

            #define legend if it is a reservoir pressure scenario
            if sc_type:
                sc_label = f'$CO_2$ P ($\\Delta$P = {sc_delta_p:.0f} bar)'

            #define legend if it is a maximum pressure scenario
            else:
                sc_label = f'max $CO_2$ P to ($\\Delta$P = {sc_delta_p:.0f} bar)'

            #include MSAD
            if plot_MSAD:
                if ~np.isnan(sc_msad_z):
                    sc_label = f'{sc_label}\nMSAD = {sc_msad_z:.0f} mTVDMSL'
                
                ax.scatter(sc_msad_p, sc_msad_z, color='firebrick')

            # Plot brine pressure profile if different from hydrostatic
            linestyle = next(line_style_cycle)
            ax.plot(sc_curves['brine_pressure'], sc_curves['depth'], label='brine pressure', color='steelblue', lw=0.75, ls = linestyle)
            ax.plot(sc_curves['fluid_pressure'], sc_curves['depth'], label='fluid pressure', color='firebrick', lw=0.75, ls = linestyle)

            if plot_fluid_contact:
                ax.scatter(sc_p_fluid_contact, sc_z_fluid_contact, color='blue', label=f'fluid contact')

            # Plot (p_resrv, z_resrv) points if option is turned on

            if plot_resrv:
                ax.scatter(sc_p_resrv, sc_z_resrv, color='green', label=f'reservoir pressure')

            if plot_delta_p and sc_delta_p != 0:
                delta_p_lims = [sc_p_fluid_contact, sc_p_fluid_contact - sc_delta_p]
                ax.hlines(sc_z_fluid_contact, min(delta_p_lims), max(delta_p_lims), color='black', linestyle=':', label=f'$\\Delta$P = {scenario["p_delta"]:.0f} bar')


    #Optimize legend
    if legend:
        ax.legend()
        handles, labels = ax.get_legend_handles_labels()  
        lgd = dict(zip(labels, handles))
        ax.legend(lgd.values(), lgd.keys())

    ax.set_xlim(xmin, xmax)
    ax.set_xlabel('Pressure [bar]')
    ax.set_ylim(0, ymax)
    ax.invert_yaxis()
    ax.grid(visible=True, linewidth=0.5)

    if 'fig' in locals():
        return fig, ax
    else:
        return ax
