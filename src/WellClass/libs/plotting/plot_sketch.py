
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

from ..well_class.well_class import Well
from ..utils.fraction_float import float_to_fraction_inches

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
    xcoord_left = -drilling_df['diameter_m'].max()/2 #well construction text
    xcoord_right = drilling_df['diameter_m'].max()/2 #geology text
    txt_fs_left = 7
    txt_fs_right = 6
    steelcolor = '#702F00'
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
    if draw_casings:
        ax.vlines(x= casings_df['diameter_m']/2, ymin=casings_df['top_msl'], ymax=casings_df['bottom_msl'],  color=steelcolor, lw=1.5, zorder=10)
        ax.vlines(x=-casings_df['diameter_m']/2, ymin=casings_df['top_msl'], ymax=casings_df['bottom_msl'],  color=steelcolor, lw=1.5, zorder=10)

    #Draw casing shoes
    if draw_casing_shoes:
        
        shoe_size = 3

        left_shoe = [[0, 0], [-shoe_size, 0], [0, shoe_size], [0, 0]]
        right_shoe = [[0, 0], [shoe_size, 0], [0, shoe_size], [0, 0]]

        shoe_query = casings_df[casings_df['shoe']]

        ax.scatter( shoe_query['diameter_m']/2, shoe_query['bottom_msl'], marker = right_shoe, c=steelcolor, zorder=10)
        ax.scatter(-shoe_query['diameter_m']/2, shoe_query['bottom_msl'], marker = left_shoe, c=steelcolor, zorder=10)

    if draw_annotation:
        for idx, row in shoe_query.iterrows():
                ycoord = row['bottom_msl']
                d_in = row['diameter_in']
                shoe_label = float_to_fraction_inches(d_in)+' shoe'
                ax.annotate(shoe_label, xy = (xcoord_left, ycoord), fontsize = txt_fs_left, va = 'center', ha='right')

    #Draw welded
    if draw_welded:
        weld_query = casings_df[~casings_df['shoe'].astype(bool)]

        for idx, row in weld_query.iterrows():
                max_D = row['diameter_m']
                max_Z = row['bottom_msl']

                query = casings_df.query('diameter_m<@max_D & top_msl==@max_Z')
                min_D = query.iloc[0]['diameter_m']

                ax.plot([ max_D/2,  min_D/2], [max_Z]*2, c=steelcolor, zorder=10)
                ax.plot([-max_D/2, -min_D/2], [max_Z]*2, c=steelcolor, zorder=10)

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
                ax.annotate(text = row['barrier_name'], xy = (0, ycoord), fontsize = txt_fs_left, va = 'center', ha='center')

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
                ax.annotate(text = row['geol_unit'], xy = (xcoord_right, ycoord), fontsize = txt_fs_right, va = 'center')

    ax.set_xlim(-ax_width, ax_width)
    ax.set_ylim(0, ymax)
    ax.invert_yaxis()
    ax.set_ylabel('depth [mMSL]')
    ax.set_xlabel('radius [m]')
    
    if 'fig' in locals():
            return fig, ax
    else:
            return ax
    

