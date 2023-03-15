# GaP
A tool to build legacy well representations on reservoir simulation grids

<h2>Introduction</h2>

To be able to model a legacy well in reservoir scale, we need to make sure all of the elements including  multiple casings with different OD, cement bonds, barriers, the open area between barriers inside the casings, etc. are considered.
Then, the existing simple well models available in the commercial simulators (Eclipse, PFT, IX) are not able to include those details. They just introduce a node where the flow will be discharged from/to the grid to the node.

The way around it is to define the well as a part of the reservoir by manipulating the local grids (LGR) and properties of the grid, in that setting:

- **The pipes** can be mimicked by a very narrow (in the size of ID of a pipe-20 to 50 cm) and high permeability grids. Those grids could be surrounded by grids with zero transmissibility.
- **The open hole sections** can be modelled the same as pipe, but without zero transmissibility around the high perm. area to allow moving of fluids side ways from the piper. This is particularly important in the cases where there is a drilled but uncased hole in the setting.
- **The barriers** can be mimicked by very low permeability with the same size of pipe ID and have their own start and end depth.
- **The cement bond** can be mimicked by low perm vertical layer adjustment to the casing with given depth and interval.

Each unit above should have specified start and end depths.

In this context, the GaP script generates the LGR, changes the properties of the LGR to reflect all requested elements, and isolates reservoir from the overburden, opens the communication between reservroir and overburden only through the well. 
