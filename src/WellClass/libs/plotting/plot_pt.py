import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from itertools import cycle
from ..well_pressure import Pressure


def plot_pt(my_pressure: Pressure, 
            fig : mpl.figure.Figure = None, 
            ax : plt.Axes =None, 
            file_only=False, 
            file_name="pt",
            legend: bool = True,
            plot_selected_scenarios: list = None,
            plot_HSP:bool = False, 
            plot_Shmin:bool = False, 
            plot_MSAD:bool = False,                    
            plot_resrv: bool = False,  # Option to plot (p_resrv, z_resrv)
            plot_fluid_contact: bool = False,
            plot_delta_p: bool = False  ):  # Option to plot (p_fluid_contact, z_fluid_contact    
     
    """ 
    pressure vs temperature
    Takes the presssure class and plots pressure vs temperature.
    Background of the image is colored by the fluid density

    my_pressure: pressure class
    fig: matplotlib fig object
    ax: matploltib ax object
    file_only: boolean to print the image to file instead of screen
    file_name: preffix of file name
    plot_HSP to plot brine hydrostatci pressure
    plot_RP to plot reservoir pressure scenarios
    plot_MSAD top plot Minimum Safatey Abandonement pressure for each reservoir scenario
    plot_maxP to plot max pressurization at a given deth
    """

    if fig is None:
        fig, ax = plt.subplots()

    # Retrieve curves from PVT data
    t = my_pressure.pvt_data['temperature']
    p = my_pressure.pvt_data['pressure']
    fluid = my_pressure.fluid_type
    rho_fluid = my_pressure.pvt_data[fluid]['rho']

    #Plot density colormap
    rho_pcm = ax.pcolormesh(t, p, rho_fluid, alpha=0.5)

    #Plot phase boundary and critical point
    if fluid == 'pure_co2':
        t_co2 = np.array([-50,-48.35,-46.69,-45.04,-43.38,-41.73,-40.08,-38.42,-36.77,-35.11,-33.46,-31.81,-30.15,-28.5,-26.84,-25.19,-23.53,-21.88,-20.23,-18.57,-16.92,-15.26,-13.61,-11.96,-10.3,-8.65,-6.99,-5.34,-3.69,-2.03,-0.38,1.28,2.93,4.58,6.24,7.89,9.55,11.2,12.86,14.51,16.16,17.82,19.47,21.13,22.78,24.43,31.05])
        p_co2 = np.array([6.8,7.27,7.77,8.29,8.83,9.4,10,10.63,11.28,11.97,12.68,13.43,14.21,15.02,15.87,16.75,17.66,18.62,19.61,20.64,21.7,22.81,23.96,25.15,26.38,27.66,28.98,30.34,31.76,33.21,34.72,36.28,37.89,39.54,41.25,43.01,44.83,46.7,48.63,50.61,52.65,54.75,56.91,59.12,61.4,63.75,73.76])
        ax.plot(t_co2, p_co2, color='k', lw=1.5, label = r'$CO_2$ phase env.')
        ax.scatter(t_co2.max(), p_co2.max(), c='k')



    
    # if not hasattr(my_pressure, "pressure_CO2"):
    #     my_pressure._compute_CO2_pressures()
    # pt_df = my_pressure.pressure_CO2

    # co2_datum = my_pressure.co2_datum  # noqa: F841
    
    #Plot fluid pressure scenarios
    ls_list = ['solid','dashed','dashdot', 'dotted']
    counter = 0

    if plot_HSP:
        my_pressure.init_curves.plot(y='hydrostatic_pressure', x='temperature', ax=ax, label='$p_{hs}$', color='steelblue', lw = 0.75)

    if plot_Shmin:
        my_pressure.init_curves.plot(x='temperature', y='min_horizontal_stress', ax=ax, label='$\\sigma_{min}$', color='k', lw = 0.75)

    ymax = 0
  

    #Retrieve pressures as well
    scenarios_summary = my_pressure.scenario_manager.get_scenarios_summary()

    # Define line styles for multiple scenarios
    line_styles = ['-', '--', '-.', ':']
    line_style_cycle = cycle(line_styles)

    #Read pressure scenarios
    # scenarios = pd.DataFrame(my_pressure.pressure_scenarios).T

    xmax = 0


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


            sc_curves =  my_pressure.scenario_manager.scenarios[sc_name].init_curves

            sc_t_temp = np.interp(float(sc_msad_z), sc_curves['depth'], sc_curves['temperature'])
            sc_t_resrv = np.interp(float(sc_z_resrv), sc_curves['depth'], sc_curves['temperature'])
            sc_t_fluid_contact = np.interp(float(sc_z_fluid_contact), sc_curves['depth'], sc_curves['temperature'])

            #plot co2 gradient
            ax.plot(sc_curves['temperature'], sc_curves['fluid_pressure'], label = f'{fluid}{sc_name}', color='firebrick', lw = 0.75, ls=next(line_style_cycle))

            #plot water gradient
            ax.plot(sc_curves['temperature'], sc_curves['brine_pressure'], label = f'brine {sc_name}', color='steelblue', lw = 0.75, ls=next(line_style_cycle))

            if plot_MSAD:
                ax.scatter(sc_t_temp, sc_msad_p, color='firebrick')

            # find max temp and pressure to use as xmax and ymax
            base_msl = sc_curves[sc_curves['depth']>sc_z_fluid_contact+100]['brine_pressure'].iloc[0]
            max_temp = sc_curves[sc_curves['depth']>sc_z_fluid_contact+100]['temperature'].iloc[0]

            if max_temp > xmax:
                xmax = max_temp

            if base_msl > ymax:
                ymax = base_msl

            if plot_fluid_contact:
                ax.scatter(sc_t_fluid_contact, sc_p_fluid_contact, color='blue', label=f'fluid contact')

            if plot_resrv:
                ax.scatter(sc_t_resrv, sc_p_resrv, color='green', label=f'reservoir pressure')

            if plot_delta_p and sc_delta_p != 0:
                delta_p_lims = [sc_p_fluid_contact, sc_p_fluid_contact - sc_delta_p]
                ax.vlines(sc_t_fluid_contact, min(delta_p_lims), max(delta_p_lims), color='black', linestyle=':', label=f'$\\Delta$P = {scenario["p_delta"]:.0f} bar')
            
  


    # xmax = pt_df['init'].query('depth_msl>(@co2_datum)+50')['temp'].iloc[0]
    
    #xmax = pt_df['temp'].max()

    #Optimize legend
    if legend:
        handles, labels = ax.get_legend_handles_labels()  
        unique_labels = dict(zip(labels, handles))
        ax.legend(unique_labels.values(), unique_labels.keys(), loc='upper left', bbox_to_anchor=(1.3,1), fontsize='small')


    # Adjust plot layout to fit colorbar and make room for the legend
    plt.subplots_adjust(right=0.8)  # Adjust this value to fit the colorbar and leave space for the legend


    ax.set_ylabel('p [bar]')
    ax.set_xlabel('T [$\\degree$C]')
    ax.set_xlim(1, xmax)
    ax.set_ylim(1, ymax)
    fig.colorbar(rho_pcm, label=r'$\rho_{CO_2}$ [$kg/m^3$]')

    # fig.tight_layout()

    if file_only:
        plt.savefig(f'{file_name}.png')
    else:
        plt.show()