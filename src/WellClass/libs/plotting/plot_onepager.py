import matplotlib.pyplot as plt 
from .plot_pressure import plot_pressure
from .plot_sketch import plot_sketch
from ..well_class.well_class import Well
from ..well_pressure.Pressure import Pressure

from IPython.display import display

#Plot sketch, pressures

def plot_onepager(well:Well, pressure:Pressure,
                  width: float = 7.68,
                  height: float = 5.5,
                  save_fig: bool = False,
                  fname: str = None,
                  #plot_sktech arguments
                  draw_drillings=True,
                  draw_casings=True,
                  draw_casing_shoes=True,
                  draw_open_hole=False,    # Note: default is False
                  draw_welded=True,
                  draw_cement_bond=True,
                  draw_barriers=True,
                  draw_geology=True,
                  draw_annotation=True,
                  #plot_pressure arguments
                  geology_dict: dict = None,
                  barriers: dict = None, 
                  plot_HSP: bool = False, 
                  plot_MSAD: bool = False, 
                  legend: bool = True,
                  plot_selected_scenarios: list = None,
                  plot_resrv: bool = False,  # Option to plot (p_resrv, z_resrv)
                  plot_fluid_contact: bool = False,
                  plot_fluid_pressure: bool = True,
                  plot_delta_p: bool = False  ):
    """
    Plot a one-pager of the well and pressure profile.
    
    Args:
        my_well (Well): Well object containing well data.
        my_pressure (Pressure): Pressure object containing pressure data.
        ax (matplotlib.axes.Axes, optional): Axes to plot on. If None, a new figure is created.
    """
    
    


    # Plot pressure profile
    display(pressure.scenarios_summary())
    fig, (ax1, ax2) = plt.subplots(1,2, sharey=True, figsize=(width, height))
    plot_sketch(well,  ax=ax1, 
                draw_drillings=draw_drillings,
                draw_casings=draw_casings,
                draw_casing_shoes=draw_casing_shoes,
                draw_open_hole=draw_open_hole,
                draw_welded=draw_welded,
                draw_cement_bond=draw_cement_bond,
                draw_barriers=draw_barriers,
                draw_geology=draw_geology,
                draw_annotation=draw_annotation)
    plot_pressure(pressure, ax=ax2, 
                  geology_dict=geology_dict,
                  barriers=barriers, 
                  plot_HSP=plot_HSP, 
                  plot_MSAD=plot_MSAD, 
                  legend=legend,
                  plot_selected_scenarios=plot_selected_scenarios,
                  plot_resrv=plot_resrv,  # Option to plot (p_resrv, z_resrv)
                  plot_fluid_contact=plot_fluid_contact,
                  plot_fluid_pressure=plot_fluid_pressure,
                  plot_delta_p=plot_delta_p)

    fig.suptitle(f'{well.header['well_name']}')
    fig.tight_layout()
    fig.subplots_adjust(wspace=0)

    if save_fig:
        if fname is None:
            wname = well.header['well_name']
            wname = wname.replace("/", "_")
            wname = wname.replace(" ", "_")

            fname = f"{wname}_onepager.png"
        fig.savefig(fname, dpi=300)
        print(f"Figure saved as {fname}")

    return fig, (ax1, ax2)
