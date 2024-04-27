
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.ticker as plticker

from ..well_class.well_class import Well
from ..utils.fraction_float import float_to_fraction_inches

def hole_plotter(axis, df, hole_bool, fill_bool = True, z_order = 0):

    """
    Draws all open hole elements. Applies for both drilling and borehole dataframes
    """

    if hole_bool:
        for idx, row in df.iterrows():
                xy = (-row['diameter_m']/2, row['top_msl'])
                width = row['diameter_m']
                height = row['bottom_msl'] - row['top_msl']
                axis.add_patch(Rectangle(xy, width, height, zorder=z_order, fill = fill_bool, facecolor=r'#CB8A58'))


 

def casings_plotter(axis, df, color_tone, txt_size, x_txt_pos,
                 annot_bool, casings_bool, c_shoe_bool, c_weld_bool):
    
    
    """
    Draws all components linked to the tubular assembly:
    - Casings
    - Casing shoes/welded transitions
    - Annotations
    """
    
    y_base  =  df['bottom_msl']
    y_top   =  df['top_msl']
    x_left  = -df['diameter_m']/2
    x_right =  df['diameter_m']/2
    shoe_size = 3



    # Create marker for casing shoe
    left_shoe = [[0, 0], [-shoe_size, 0], [0, shoe_size], [0, 0]]
    right_shoe = [[0, 0], [shoe_size, 0], [0, shoe_size], [0, 0]]

    
    # query dataframe for shoe items
    shoe_query = df[df['shoe']]

    # define x and y positions
    x_pos_shoe = shoe_query['diameter_m']/2
    y_pos_shoe = shoe_query['bottom_msl']

    # query dataframe for welded sections
    weld_query = df[~df['shoe'].astype(bool)]
    
    # Draw casings
    if casings_bool:
        #draw right hand casing
        axis.vlines(x =  x_right, ymin = y_top, ymax = y_base,  color = color_tone, lw= 1.5, zorder = 10)

        #draw left hand casing
        axis.vlines(x =  x_left,  ymin = y_top, ymax = y_base,  color = color_tone, lw= 1.5, zorder = 10)

    # Draw casing shoes
    if c_shoe_bool:
        #draw left casing shoe
        axis.scatter( x_pos_shoe, y_pos_shoe, marker = right_shoe, c=color_tone, zorder=10)

        #draw right casing shoe
        axis.scatter(-x_pos_shoe, y_pos_shoe, marker = left_shoe,  c=color_tone, zorder=10)

    # Draw welded
    if c_weld_bool:
        
        for idx, row in weld_query.iterrows():
                max_D = row['diameter_m']
                max_Z = row['bottom_msl']

                query = casings_df.query('diameter_m<@max_D & top_msl==@max_Z')
                min_D = query.iloc[0]['diameter_m']

                axis.plot([ max_D/2,  min_D/2], [max_Z]*2, c=color_tone, zorder=10)
                axis.plot([-max_D/2, -min_D/2], [max_Z]*2, c=color_tone, zorder=10)
          

    # Draw annotations
    if annot_bool:
        for idx, row in shoe_query.iterrows():
                ycoord = row['bottom_msl']
                d_in =   row['diameter_in']
                shoe_label = float_to_fraction_inches(d_in)+' shoe'
                
                axis.annotate(shoe_label, xy = (x_txt_pos, ycoord), fontsize = txt_size, va = 'center', ha='right')


def cement_bond_plotter(axis, df, cement_bond_bool):
    if cement_bond_bool:
        for idx, row in df.iterrows():
            width = (row['od_m'] - row['id_m'])/2
            height = row['bottom_msl'] - row['top_msl']

            right_xy = (row['id_m']/2, row['top_msl'])
            left_xy = (-row['od_m']/2, row['top_msl'])

            axis.add_patch(Rectangle(right_xy, width, height, facecolor='lightgray', zorder=5, hatch='\\\\\\'))
            axis.add_patch(Rectangle(left_xy, width, height, facecolor='lightgray', zorder=5 , hatch='///'))
            
      

def cement_plug_plotter(axis, df_barriers, df_barriers_mod, plug_bool, annot_bool, txt_size):

    """
    axis: Matplotlib object where items will be plotted
    df_barriers: Dataframe listing the barriers 
    df_barriers_mod: Dataframe that describes the visual display of barriers.
    plug_bool: Boolean if barriers are to be displayed
    annot_bool: Boolean if annotations are to be included

    """

    if plug_bool:
        for idx, row in df_barriers_mod.iterrows():

                xy = (-row['diameter_m']/2, row['top_msl'])
                width = row['diameter_m']
                height = row['bottom_msl'] - row['top_msl']
                axis.add_patch(Rectangle(xy, width, height, facecolor='gray', zorder=1))

    if annot_bool:
        for idx, row in df_barriers.iterrows():
                ycoord = (row['top_msl'] + row['bottom_msl'])/2
                axis.annotate(text = row['barrier_name'], xy = (0, ycoord), fontsize = txt_size, va = 'center', ha='center')

def geology_plotter(axis, df_geol, w_header, geol_bool, annot_bool, width, x_txt_pos, txt_size):
    if geol_bool:
        axis.hlines(y=df_geol['top_msl'], xmin=-width, xmax=width, zorder=-4, lw=.25, color='k')
        axis.axhspan(0, w_header['sf_depth_msl'], color='lightblue', alpha=0.5, zorder=-20)
        axis.axhspan(w_header['sf_depth_msl'], w_header['well_td_rkb'], color='tan', alpha=0.5, zorder=-20)

    if annot_bool:
        for index, row in df_geol.iterrows():
                if row['reservoir_flag']:
                        axis.axhspan(row['top_msl'], row['base_msl'], color='yellow', zorder=-10)
        
                ycoord = (row['top_msl'] + row['base_msl'])/2
                axis.annotate(text = row['geol_unit'], xy = (x_txt_pos, ycoord), fontsize = txt_size, va = 'center')


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
                save_file=None):

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
    
    # save figure to the disk
    if save_file:
          plt.savefig(save_file)

    if 'fig' in locals():
            return fig, ax
    else:
            return ax
    

