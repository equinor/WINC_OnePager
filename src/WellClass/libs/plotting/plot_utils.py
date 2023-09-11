
import matplotlib.pyplot as plt

import matplotlib.colors as colors

# plots
from .plot_sketch import plot_sketch

def plot_well_perm(my_well, *, x, y, Z, on_coarse=True):
    """ well sketch and permeability, x-z slice
    """

    fig, (ax2, ax) = plt.subplots(1,2, sharey=True)

    # 1. well sketch
    plot_sketch(my_well, ax2)

    # 2. property corresponding to grid coordinates
    plot_perm(x, y, Z, ax=ax, on_coarse=on_coarse)

    fig.tight_layout(h_pad = 0)

    plt.show()

def plot_perm(x, y, Z, *, ax, on_coarse=True):
    """ plot permeability, x-z slice
    """

    # min/max permeability
    zmin = Z.min()
    zmax = Z.max()
    if zmin == 0:
        zmin = 1e-4

    # print(f' zmin = {zmin}, zmax = {zmax}')
    # print(f' x shape = {x.shape}, y shape = {y.shape}')
    # print(f' Z shape = {Z.shape}')

    if on_coarse:    
        # ax.pcolor(Z,  norm=colors.LogNorm(vmin=Z.min(), vmax=Z.max()))
        ax.pcolormesh(x, y, Z,  norm=colors.LogNorm(vmin=zmin, vmax=zmax))
        # ax.pcolormesh(x, y, Z,  norm=colors.LogNorm(vmin=zmin, vmax=zmax), edgecolor='gray', lw = 0.01)
        ax.set_xlim(-1e3, 1e3)

    else:
         
        # ax.pcolor(Z,  norm=colors.LogNorm(vmin=Z.min(), vmax=Z.max()))
        ax.pcolormesh(x, y, Z,  norm=colors.LogNorm(vmin=zmin, vmax=zmax), edgecolor='gray', alpha=0.5, lw = 0.01)

        # for idx, row in drilling_df.iterrows():
        #         xy = (-row['cart_dist']/2, row['top_msl'])
        #         width = row['cart_dist']
        #         height = row['bottom_msl'] - row['top_msl']
        #         ax.add_patch(Rectangle(xy, width, height, zorder=10, facecolor=r'#CB8A58', alpha=0.5))


        # #Draw casings
        # ax.vlines(x= casings_df['cart_dist']/2, ymin=casings_df['top_msl'], ymax=casings_df['bottom_msl'],  color='k', lw=1.5, zorder=10)
        # ax.vlines(x=-casings_df['cart_dist']/2, ymin=casings_df['top_msl'], ymax=casings_df['bottom_msl'],  color='k', lw=1.5, zorder=10)

        ax.set_xlim(-.8, .8)
        # ax.set_ylim(710, 670)
