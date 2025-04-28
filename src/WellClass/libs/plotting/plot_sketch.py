
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as plticker

from ..well_class.well_class import Well
from .plot_wellbore_elements import *


def plot_sketch(mywell: Well, ax=None, 
                *, 
                draw_drillings=True,
                draw_casings=True,
                draw_casing_shoes=True,
                draw_open_hole=False,    # Note: default is False
                draw_welded=True,
                draw_cement_bond=True,
                draw_barriers=True,
                draw_geology=True,
                draw_annotation=True,
                save_file=None,
                ):

    """ 
    plot well sketch
    
    """

    if ax is None:
        fig, ax = plt.subplots()

    #Automatically define the y scale in 250 m intervals
    loc = plticker.MultipleLocator(base=250) # this locator puts ticks at regular intervals
    ax.yaxis.set_major_locator(loc)
            

    drilling_df =     pd.DataFrame(mywell.drilling)
    casings_df =      pd.DataFrame(mywell.casings)
    cb_df =           pd.DataFrame(mywell.cement_bond)
    borehole_df =     pd.DataFrame(mywell.borehole)
    barriers_df =     pd.DataFrame(mywell.barriers)
    barriers_fmt_df = pd.DataFrame(mywell.barriers_mod)
    geology_df =      pd.DataFrame(mywell.geology)
    well_header =     mywell.header



    #define plot spatial references
    AX_WIDTH = 2 * drilling_df['diameter_m'].max()/2 #plot width
    XCOORD_LEFT = -drilling_df['diameter_m'].max()/2 #well construction text
    XCOORD_RIGHT = drilling_df['diameter_m'].max()/2 #geology text
    TXT_FS_LEFT = 7
    TXT_FS_RIGHT = 6
    STEELCOLOR = '#702F00'
    base_deepest_rsrv = geology_df[geology_df.reservoir_flag]['base_msl'].max()
    ymax = max([base_deepest_rsrv,mywell.co2_datum])+100

    # Draw drilling (Bit size)
    hole_plotter(axis = ax, df = drilling_df,  hole_bool=draw_drillings)    
   
    #Draw casings
    casings_plotter(axis = ax, df = casings_df, color_tone=STEELCOLOR, txt_size=TXT_FS_LEFT,
                 x_txt_pos=XCOORD_LEFT, annot_bool=draw_annotation, casings_bool=draw_casings,
                 c_shoe_bool = draw_casing_shoes, c_weld_bool=draw_welded)


    #Draw cement bond
    cement_bond_plotter(axis = ax, df = cb_df, cement_bond_bool=draw_cement_bond)


    #draw barriers
    cement_plug_plotter(axis = ax, df_barriers=barriers_df, df_barriers_mod=barriers_fmt_df, plug_bool = draw_barriers,
                        annot_bool = draw_annotation, txt_size= TXT_FS_LEFT)

    #Draw open hole (borehole/pipe) for testing only
    hole_plotter(axis = ax, df = borehole_df,  hole_bool=draw_open_hole, fill_bool=False, z_order = 100)    


    #Draw geological information
    geology_plotter(axis = ax, df_geol = geology_df, w_header = well_header, geol_bool=draw_geology, 
                    annot_bool=draw_annotation, width = AX_WIDTH, x_txt_pos=XCOORD_RIGHT, txt_size=TXT_FS_RIGHT)


    ax.set_xlim(-AX_WIDTH, AX_WIDTH)
    ax.set_ylim(0, ymax)
    ax.invert_yaxis()
    ax.set_ylabel('depth [mMSL]')
    ax.set_xlabel('radius [m]')
    
    # Hide x-axis ticks and labels for the left subplot (ax1)
    ax.set_xticks([])  # Remove x-axis ticks
    ax.set_xticklabels([])  # Remove x-axis tick labels
    ax.set_xlabel("")  # Remove x-axis label



    # save figure to the disk
    if save_file:
          plt.savefig(save_file)

    if 'fig' in locals():
            return fig, ax
    else:
            return ax
    

