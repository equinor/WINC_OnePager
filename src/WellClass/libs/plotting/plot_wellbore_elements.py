import pandas as pd
from matplotlib.axes import Axes
from matplotlib.patches import Rectangle

from ..utils.fraction_float import float_to_fraction_inches


def hole_plotter(axis: Axes, df: pd.DataFrame, hole_bool: bool, fill_bool: bool = True, z_order: int = 0) -> None:
    """
    Draws all open hole elements. Applies for both drilling and borehole dataframes
    """
    if hole_bool:
        for idx, row in df.iterrows():
            xy = (-row["diameter_m"] / 2, row["top_msl"])
            width = row["diameter_m"]
            height = row["bottom_msl"] - row["top_msl"]
            axis.add_patch(Rectangle(xy, width, height, zorder=z_order, fill=fill_bool, facecolor=r"#CB8A58"))


def casings_plotter(
    axis: Axes,
    df: pd.DataFrame,
    color_tone: str,
    txt_size: int,
    x_txt_pos: float,
    annot_bool: bool,
    casings_bool: bool,
    c_shoe_bool: bool,
    c_weld_bool: bool,
) -> None:
    """
    Draws all components linked to the tubular assembly:
    - Casings
    - Casing shoes/welded transitions
    - Annotations
    """
    y_base = df["bottom_msl"]
    y_top = df["top_msl"]
    x_left = -df["diameter_m"] / 2
    x_right = df["diameter_m"] / 2
    shoe_size = 3

    # Create marker for casing shoe
    left_shoe = [[0, 0], [-shoe_size, 0], [0, shoe_size], [0, 0]]
    right_shoe = [[0, 0], [shoe_size, 0], [0, shoe_size], [0, 0]]

    # query dataframe for shoe items
    shoe_query = df[df["shoe"]]

    # define x and y positions
    x_pos_shoe = shoe_query["diameter_m"] / 2
    y_pos_shoe = shoe_query["bottom_msl"]

    # query dataframe for welded sections
    weld_query = df[~df["shoe"].astype(bool)]

    # Draw casings
    if casings_bool:
        # draw right hand casing
        axis.vlines(x=x_right, ymin=y_top, ymax=y_base, color=color_tone, lw=1.5, zorder=10)

        # draw left hand casing
        axis.vlines(x=x_left, ymin=y_top, ymax=y_base, color=color_tone, lw=1.5, zorder=10)

    # Draw casing shoes
    if c_shoe_bool:
        # draw left casing shoe
        axis.scatter(x_pos_shoe, y_pos_shoe, marker=right_shoe, c=color_tone, zorder=10)

        # draw right casing shoe
        axis.scatter(-x_pos_shoe, y_pos_shoe, marker=left_shoe, c=color_tone, zorder=10)

    # Draw welded
    if c_weld_bool:
        for idx, row in weld_query.iterrows():
            max_D = row["diameter_m"]
            max_Z = row["bottom_msl"]

            query = df.query("diameter_m<@max_D & top_msl==@max_Z")
            min_D = query.iloc[0]["diameter_m"]

            axis.plot([max_D / 2, min_D / 2], [max_Z] * 2, c=color_tone, zorder=10)
            axis.plot([-max_D / 2, -min_D / 2], [max_Z] * 2, c=color_tone, zorder=10)

    # Draw annotations
    if annot_bool:
        for idx, row in shoe_query.iterrows():
            ycoord = row["bottom_msl"]
            d_in = row["diameter_in"]
            shoe_label = float_to_fraction_inches(d_in) + " shoe"

            axis.annotate(shoe_label, xy=(x_txt_pos, ycoord), fontsize=txt_size, va="center", ha="right")


def cement_bond_plotter(axis: Axes, df: pd.DataFrame, cement_bond_bool: bool) -> None:
    if cement_bond_bool:
        for idx, row in df.iterrows():
            width = (row["od_m"] - row["id_m"]) / 2
            height = row["bottom_msl"] - row["top_msl"]

            right_xy = (row["id_m"] / 2, row["top_msl"])
            left_xy = (-row["od_m"] / 2, row["top_msl"])

            axis.add_patch(Rectangle(right_xy, width, height, facecolor="lightgray", zorder=5, hatch="\\\\\\"))
            axis.add_patch(Rectangle(left_xy, width, height, facecolor="lightgray", zorder=5, hatch="///"))


def cement_plug_plotter(
    axis: Axes, df_barriers: pd.DataFrame, df_barriers_mod: pd.DataFrame, plug_bool: bool, annot_bool: bool, txt_size: int
) -> None:
    """
    axis: Matplotlib object where items will be plotted
    df_barriers: Dataframe listing the barriers
    df_barriers_mod: Dataframe that describes the visual display of barriers.
    plug_bool: Boolean if barriers are to be displayed
    annot_bool: Boolean if annotations are to be included

    """
    if plug_bool:
        for idx, row in df_barriers_mod.iterrows():
            xy = (-row["diameter_m"] / 2, row["top_msl"])
            width = row["diameter_m"]
            height = row["bottom_msl"] - row["top_msl"]
            axis.add_patch(Rectangle(xy, width, height, facecolor="gray", zorder=1))

    if annot_bool:
        for idx, row in df_barriers.iterrows():
            ycoord = (row["top_msl"] + row["bottom_msl"]) / 2
            axis.annotate(text=row["barrier_name"], xy=(0, ycoord), fontsize=txt_size, va="center", ha="center")


def geology_plotter(
    axis: Axes, df_geol: pd.DataFrame, w_header: dict, geol_bool: bool, annot_bool: bool, width: float, x_txt_pos: float, txt_size: int
) -> None:
    if geol_bool:
        axis.hlines(y=df_geol["top_msl"], xmin=-width, xmax=width, zorder=-4, lw=0.25, color="k")
        axis.axhspan(0, w_header["sf_depth_msl"], color="lightblue", alpha=0.5, zorder=-20)
        axis.axhspan(w_header["sf_depth_msl"], w_header["well_td_rkb"], color="tan", alpha=0.5, zorder=-20)

    if annot_bool:
        for index, row in df_geol.iterrows():
            if row["reservoir_flag"]:
                axis.axhspan(row["top_msl"], row["base_msl"], color="yellow", zorder=-10)

            ycoord = (row["top_msl"] + row["base_msl"]) / 2
            axis.annotate(text=row["geol_unit"], xy=(x_txt_pos, ycoord), fontsize=txt_size, va="center")
