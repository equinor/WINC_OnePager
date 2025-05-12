import matplotlib.pyplot as plt 
from .plot_pressure import plot_pressure
from .plot_sketch import plot_sketch
from ..well_class.well_class import Well
from ..well_pressure.Pressure import Pressure
from .plot_wellbore_elements import geology_plotter
from IPython.display import display
import pandas as pd

def plot_onepager(well: Well, pressure: Pressure,
                  include_sketch: bool = True,
                  include_pressure: bool = True,
                  filename: str = None,
                  width: float = 10,
                  height: float = 8,
                  # sketch arguments
                  sketch_drillings: bool = True,
                  sketch_casings: bool = True,
                  sketch_casing_shoes: bool = True,
                  sketch_open_hole: bool = False,
                  sketch_welded: bool = True,
                  sketch_cement_bond: bool = True,
                  sketch_barriers: bool = True,
                  sketch_geology: bool = True,
                  sketch_annotation: bool = True,
                  # plot_pressure arguments
                  pressure_geology: bool = True,
                  pressure_barriers: bool = True,
                  pressure_HSP: bool = False,
                  pressure_MSAD: bool = True,
                  pressure_legend: bool = True,
                  pressure_selected_scenarios: list = None,
                  pressure_resrv: bool = False,
                  pressure_fluid_contact: bool = False,
                  pressure_fluid_pressure: bool = True,
                  pressure_delta_p: bool = False):
    """
    Plot a one-pager of the well and pressure profile.
    
    Args:
        well (Well): Well object containing well data.
        pressure (Pressure): Pressure object containing pressure data.
    """

    display(pressure.scenarios_summary())

    if not include_sketch:
        fig, (ax_sketch, ax_pressure) = plt.subplots(1,2, sharey=True, figsize=(width, height), width_ratios=(1, 5))
        ax = (ax_sketch, ax_pressure)

    elif not include_pressure:
        fig, ax_sketch = plt.subplots(figsize=(width, height))
        ax = ax_sketch
        
    else:
        fig, (ax_sketch, ax_pressure) = plt.subplots(1,2, sharey=True, figsize=(width, height))
        ax = (ax_sketch, ax_pressure)

    if include_sketch:
        plot_sketch(well, ax=ax_sketch, 
                    draw_drillings=sketch_drillings,
                    draw_casings=sketch_casings,
                    draw_casing_shoes=sketch_casing_shoes,
                    draw_open_hole=sketch_open_hole,
                    draw_welded=sketch_welded,
                    draw_cement_bond=sketch_cement_bond,
                    draw_barriers=sketch_barriers,
                    draw_geology=sketch_geology,
                    draw_annotation=sketch_annotation)
    else:
        geology_df = pd.DataFrame(well.geology)
        geology_plotter(axis=ax_sketch, df_geol=geology_df, w_header=well.header, 
                        geol_bool=True, annot_bool=sketch_annotation, width=2, 
                        x_txt_pos=-1.8, txt_size=6)
        ax_sketch.set_xlim(-2, 2)
        ax_sketch.xaxis.set_visible(False)

    if include_pressure:
        if pressure_geology:
            geology_dict = well.geology
        else:
            geology_dict = None
            
        if pressure_barriers:
            barriers_dict = well.barriers
        else:
            barriers_dict = None

        plot_pressure(pressure, ax=ax_pressure, 
                geology_dict=geology_dict,
                barriers=barriers_dict, 
                plot_HSP=pressure_HSP, 
                plot_MSAD=pressure_MSAD, 
                legend=pressure_legend,
                plot_selected_scenarios=pressure_selected_scenarios,
                plot_resrv=pressure_resrv,  # Option to plot (p_resrv, z_resrv)
                plot_fluid_contact=pressure_fluid_contact,
                plot_fluid_pressure=pressure_fluid_pressure,
                plot_delta_p=pressure_delta_p)

    fig.suptitle(f"{well.header['well_name']}")
    fig.tight_layout()
    fig.subplots_adjust(wspace=0)

    if filename is not None:
        wname = well.header['well_name'].replace("/", "_").replace(" ", "_")
        fname = f"{wname}_onepager.png"
        fig.savefig(fname, dpi=300)
        print(f"Figure saved as {fname}")

    return fig, (ax_sketch, ax_pressure) if include_pressure else ax_sketch
