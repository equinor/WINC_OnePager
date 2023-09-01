# GaP code

A tool to build legacy well representations on reservoir simulation grids

## Introduction ##

To be able to model a legacy well in reservoir scale, we need to make sure all of the elements including  multiple casings with different OD, cement bonds, barriers, the open area between barriers inside the casings, etc. are considered.
Then, the existing simple well models available in the commercial simulators (Eclipse, PFT, IX) are not able to include those details. They just introduce a node where the flow will be discharged from/to the grid to the node.

The way around it is to define the well as a part of the reservoir by manipulating the local grids (LGR) and properties of the grid, in that setting:

- **The pipes** can be mimicked by a very narrow (in the size of ID of a pipe-20 to 50 cm) and high permeability grids. Those grids could be surrounded by grids with zero transmissibility.
- **The open hole sections** can be modelled the same as pipe, but without zero transmissibility around the high perm. area to allow moving of fluids side ways from the piper. This is particularly important in the cases where there is a drilled but uncased hole in the setting.
- **The barriers** can be mimicked by very low permeability with the same size of pipe ID and have their own start and end depth.
- **The cement bond** can be mimicked by low perm vertical layer adjustment to the casing with given depth and interval.

Each unit above should have specified start and end depths.

In this context, the GaP script generates the LGR, changes the properties of the LGR to reflect all requested elements, and isolates reservoir from the overburden, opens the communication between reservroir and overburden only through the well. 

## Experiments
There are a couple of examples for verifying the codes.

### - Configuration format

We define GaP input `yaml` format similar to `kubernetes`. Here is the sample for `smeaheia` dataset:
```yaml
apiVersion: 'gap/v0.1'

metadata:
  name: 'smeaheia'
  
spec:

  sim_case:
    folder: "data/smeaheia"
    filename: "GEN_NOLGR_PH2"

  lgr_out:
    folder: "data/smeaheia"
    filename: "LEG_HIRES"

  defaults:
    mindz_ob: 10.0

  casings:

    - type: "conductor"
      ID: 0.762
      pipe:
        strt_depth: 312
        end_depth: 371
        perm: 10000      # permeability of tube
      oph:
        strt_depth: 4
        end_depth: 376
      cement:
        strt_depth: 312
        end_depth: 371
        perm: 5          # permeability of cement bond

    - type: "production"
      ID: 0.244
      pipe:
        strt_depth: 599
        end_depth: 1114
        perm: 10000      # permeability of tube
      oph:
        strt_depth: 599
        end_depth: 1622
      cement:
        strt_depth: 683
        end_depth: 1114
        perm: 5          # permeability of cement bond

    - type: "surface"
      ID: 0.3397
      pipe:
        strt_depth: 312
        end_depth: 686
        perm: 10000      # permeability of tube
      oph:
        strt_depth: 312
        end_depth: 697
      cement:
        strt_depth: 312
        end_depth: 686
        perm: 5          # permeability of cement bond

  barriers:

    - type: "barrier"
      ID: 0.3397
      pipe:
        strt_depth: 352
        end_depth: 560
        perm: 0.5

    - type: "barrier"
      ID: 0.244
      pipe:
        strt_depth: 1020
        end_depth: 1046
        perm: 100
```
Note the input order of casings and barriers doesn't matter. And we intentionally mess up the oder of ```production casing``` and ```sufface casing```. The codes know how to sort them according to their `ID`s. And some of the fields are optional.

### - Run the examples
There are a couple of ways to test the codes: one is to run the test codes in commandline and the other is to run notebooks in the browser with jupyter app.

To run it in **commandline**, type the following inside the `GaP` directiory:
```
$ python -m experiments.GaP_HIRES --config-file data/smeaheia/config.yaml --display
```
This is an example for two-barriers. The output by default will be in the same directory as input. To check the differences between the current run with original output, try the following:
```
$ diff data/smeaheia/LEG_HIRES.grdecl data/smeaheia/LEG_HIRES.grdecl.original
```
Nothing will show up if they are exactly the same. Otherwise, there will be bunch of differences.

For one-barrier example, please try the following:
```
$ python -m experiments.GaP_HIRES --config-file data/smeaheia_onepipe/config.yaml --display
```
The following will be used for comparisons:
```
$ diff data/smeaheia_onepipe/LEG_HIRES_BD_V1.grdecl data/smeaheia_onepipe/LEG_HIRES_BD_V1.grdecl.original
```

To run the **jupyter** notebooks, change from current directory to the `notebooks` directory, and launch jupyter lab app:
```
$ jupyter-lab
```
and play with the `demo.ipynb` for details.
