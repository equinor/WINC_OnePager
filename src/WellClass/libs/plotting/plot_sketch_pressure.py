import matplotlib.pyplot as plt

from ..well_class.well_class import Well
from ..well_pressure import Pressure
from .plot_pressure import plot_pressure
from .plot_sketch import plot_sketch


def plot_sketch_pressure(my_well: Well, my_pressure: Pressure, *, save_file=None):
    fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(10, 8))

    plot_sketch(my_well, ax=ax1)

    plot_pressure(my_pressure, my_well.geology, my_well.barriers, ax=ax2)

    title = f"{my_well.header['well_name']}\n{my_pressure.fluid_type}: {my_pressure.fluid_composition}"

    fig.suptitle(title, fontsize=10)

    fig.tight_layout()
    fig.subplots_adjust(wspace=0)

    if save_file:
        fig.savefig(save_file)
