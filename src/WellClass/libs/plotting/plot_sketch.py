
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

from ..well_class.well_class import Well
from ..utils.fraction_float import float_to_fraction_inches

def plot_casings(axis, df, color_tone, txt_size, x_txt_pos,
                 annot_bool, casings_bool, c_shoe_bool, c_weld_bool):
    
    steelcolor = '#702F00'
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
                draw_annotation=True):
    """ plot well sketch
    """

    if ax is None:
        fig, ax = plt.subplots()
            

    drilling_df =     pd.DataFrame(mywell.drilling)
    casings_df =      pd.DataFrame(mywell.casings)
    cb_df =           pd.DataFrame(mywell.cement_bond)
    borehole_df =     pd.DataFrame(mywell.borehole)
    barriers_df =     pd.DataFrame(mywell.barriers)
    barriers_fmt_df = pd.DataFrame(mywell.barriers_mod)
    geology_df =      pd.DataFrame(mywell.geology)
    well_header =     mywell.header



    #define plot spatial references
    ax_width = 2 * drilling_df['diameter_m'].max()/2 #plot width
    XCOORD_LEFT = -drilling_df['diameter_m'].max()/2 #well construction text
    XCOORD_RIGHT = drilling_df['diameter_m'].max()/2 #geology text
    TXT_FS_LEFT = 7
    TXT_FS_RIGHT = 6
    STEELCOLOR = '#702F00'
    base_deepest_rsrv = geology_df[geology_df.reservoir_flag]['base_msl'].max()
    ymax = max([base_deepest_rsrv,mywell.co2_datum])+100

    # Draw drilling (Bit size)    
    if draw_drillings:
        for idx, row in drilling_df.iterrows():
                xy = (-row['diameter_m']/2, row['top_msl'])
                width = row['diameter_m']
                height = row['bottom_msl'] - row['top_msl']
                ax.add_patch(Rectangle(xy, width, height, zorder=0, facecolor=r'#CB8A58'))

    #Draw casings
    plot_casings(axis = ax, df = casings_df, color_tone=STEELCOLOR, txt_size=TXT_FS_LEFT,
                 x_txt_pos=XCOORD_LEFT, annot_bool=draw_annotation, casings_bool=draw_casings,
                 c_shoe_bool = draw_casing_shoes, c_weld_bool=draw_welded)


    #Draw cement bond
    if draw_cement_bond:
        for idx, row in cb_df.iterrows():
                width = (row['od_m'] - row['id_m'])/2
                height = row['bottom_msl'] - row['top_msl']

                right_xy = (row['id_m']/2, row['top_msl'])
                left_xy = (-row['od_m']/2, row['top_msl'])

                ax.add_patch(Rectangle(right_xy, width, height, facecolor='lightgray', zorder=5, hatch='\\\\\\'))
                ax.add_patch(Rectangle(left_xy, width, height, facecolor='lightgray', zorder=5 , hatch='///'))

    #draw barriers
    if draw_barriers:
        for idx, row in barriers_fmt_df.iterrows():

                xy = (-row['diameter_m']/2, row['top_msl'])
                width = row['diameter_m']
                height = row['bottom_msl'] - row['top_msl']
                ax.add_patch(Rectangle(xy, width, height, facecolor='gray', zorder=1))

    if draw_annotation:
        for idx, row in barriers_df.iterrows():
                ycoord = (row['top_msl'] + row['bottom_msl'])/2
                ax.annotate(text = row['barrier_name'], xy = (0, ycoord), fontsize = TXT_FS_LEFT, va = 'center', ha='center')

    #Draw open hole (borehole/pipe) for testing only
    if draw_open_hole:
        for idx, row in borehole_df.iterrows():
                xy = (-row['id_m']/2, row['top_msl'])
                width = row['id_m']
                height = row['bottom_msl'] - row['top_msl']
                ax.add_patch(Rectangle(xy, width, height, zorder=9, fill=False, edgecolor='k', lw=2, ))


    #Draw geological information
    if draw_geology:
        ax.hlines(y=geology_df['top_msl'], xmin=-ax_width, xmax=ax_width, zorder=-4, lw=.25, color='k')
        ax.axhspan(0, well_header['sf_depth_msl'], color='lightblue', alpha=0.5, zorder=-20)
        ax.axhspan(well_header['sf_depth_msl'], well_header['well_td_rkb'], color='tan', alpha=0.5, zorder=-20)

    if draw_annotation:
        for index, row in geology_df.iterrows():
                if row['reservoir_flag']:
                        ax.axhspan(row['top_msl'], row['base_msl'], color='yellow', zorder=-10)
        
                ycoord = (row['top_msl'] + row['base_msl'])/2
                ax.annotate(text = row['geol_unit'], xy = (XCOORD_RIGHT, ycoord), fontsize = TXT_FS_RIGHT, va = 'center')

    ax.set_xlim(-ax_width, ax_width)
    ax.set_ylim(0, ymax)
    ax.invert_yaxis()
    ax.set_ylabel('depth [mMSL]')
    ax.set_xlabel('radius [m]')
    
    if 'fig' in locals():
            return fig, ax
    else:
            return ax
    

