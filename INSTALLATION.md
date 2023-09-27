# SCREEN

This repository contains codes for SCREEN project.

## Clone the repository
Locate a folder at your local machine that you intend to investigate the codes, and then clone the repository
```
$ git clone https://github.com/equinor/SCREEN
```
Note here `$` is a linux command prompt and therefore there is no need to type it. By this time you should have a folder named `SCREEN` at your local machine. Now change the directory with linux command:
```
$ cd SCREEN
```

It's normal for us to make a new branch if we indend to make some changes of the codes. This can be done with the `-b` option, for example:
```
$ git checkout -b hzh/cleanup
```
This would generate a new branch, named `hzh/cleanup`. Here I use my `equinor` account short name `hzh` concatenated with a feature description `cleanup`. There is no need to follow my convention. However, please note branch names have limitations. You probably should pick a branch name that makes sense. 

## Virtual environment

It's a common practice to work on a project within a python virtual environment. I have been using python's builtin module `venv` for a long while. So I am going to stick to it here as an example to set up the virtual environment. But you are free to use any other virtual environment setups that you feel comfortable with, such as *poetry*, *conda*, etc. 
```
$ python -m venv venv_screen
```
This will build a virtual environment `venv_screen`. You only need to creat it once.

To activate this `venv_screen`, run the following command: 
```
$ source venv_screen/bin/activate.csh
```
when your linux Shell is `csh`. 

If you are using `bash` or plain `sh`, you can activate it with the following command:
```
$ source venv_screen/bin/activate
```

We pack needed python packages into a file, such as `requirements.txt`. And install those python packages to this virtual environment by running the following command:
```
$ pip install -r requirements.txt
```
You should now be ready to play with the source codes.


## The code structures

The following represents the current code structures:

```
├── experiments
│   ├── gap_pflotran.py
│   ├── gap_wellclass.py
│   ├── __init__.py
│   └── LEG_HIRES.grdecl
├── INSTALLATION.md
├── mkdocs.yml
├── notebooks
│   ├── GaP-WellClass2.ipynb
│   ├── GaP-WellClass.ipynb
│   ├── LEG_HIRES.grdecl
│   ├── LEG_HIRES.grdecl.smeaheia
│   ├── pflotran_gap.ipynb
│   ├── Pressure-WellClass.ipynb
│   ├── Pressure-WellClass_test.ipynb
│   └── WellClass_csv_yaml.ipynb
├── README.md
├── requirements.txt
├── src
│   ├── GaP
│   │   ├── data
│   │   ├── experiments
│   │   ├── __init__.py
│   │   ├── libs
│   │   ├── notebooks
│   │   └── README.md
│   ├── __init__.py
│   ├── PressCalc
│   │   ├── 1D_PresCalc.ipynb
│   │   ├── __init__.py
│   │   ├── phase_envelope.png
│   │   ├── Pressure_plot.png
│   │   ├── PT_01012996
│   │   ├── PT_010153
│   │   └── Readme.md
│   ├── WellClass
│   │   ├── __init__.py
│   │   ├── libs
│   │   ├── notebooks
│   │   └── README.md
│   └── WellViz
│       ├── __init__.py
│       ├── Readme.md
│       └── WellViz_Jan23_Dash_v4.py
└── test_data
    ├── examples
    │   ├── cosmo
    │   ├── cosmo-pflotran
    │   ├── frigg
    │   ├── simple_well
    │   ├── smeaheia_v1
    │   └── smeaheia_v2
    └── pvt_constants
        ├── pressure.txt
        ├── rho_co2.txt
        ├── rho_h2o.txt
        └── temperature.txt
```
It was generated with the linux command `tree`:
```
$ tree -I 'docs|site|venv_screen|*pycache*|Equinor*|originals' -L 3
```
## Experiments
There are at least two ways to test the codes. One is to run the experiments with Jupyter lab, and the other is commnadline option.

### Jupyter notebooks
Jupyter notebooks are located in directory `notebooks`. To test its functionaries, change current directory to `notebooks` and launch jupyter notebooks at the commandline:
```
$ jupyter-lab
```
There are currently six notebooks in the directory:
- Notebook **GaP-WellClass.ipynb** is used to test the integration of `GaP` and `WellClass`. It involves testing of many low-level functions. We don't recommend that you start with this one.
- Notebook **GaP-WellClass2.ipynb** is another test example for the integration of `GaP` and `WellClass`. It has a better structure and hides many details. Both this and previous notebook require eclipse `.EGRID` and `.INIT` files.
- Notebook **pflotran-gap.ipynb** integrates `GaP` and `WellClass` too. But instead of the user-provided `.EGRID` and `.INIT` files, both files are generates by calling pflotran scripts.
- Notebook **Pressure-WellClass_test.ipynb** is Alejandro's tests on deviated wells
- Notebook **Pressure-WellClass.ipynb** is used to test pressure. 
- Notebook **WellClass_csv_yaml.ipynb** is used to test pressure and loading `.csv` and `.yaml` input files.

### Commandline option
Two python scripts for commandline option are available in directory `experiments`. One script, **gap_plotran.py**, generates Eclipse `.EGRID` and `.INIT` on the fly, while the other script, **gap_wellclass.py**, requires the user to provide these two grid files.  

The followings are some of  the sample runs. In either way, you should run the python script inside the ```SCREEN``` directory. 

1. To test **gap_wellclass.py**, run either of the followings:
```python
# 1. for smeaheia_v1

$ python -m experiments.gap_wellclass --sim-path ./test_data/examples/smeaheia_v1 --well smeaheia.yaml --sim-case GEN_NOLGR_PH2 --plot 

# 2. for smeaheia_v2

$ python -m experiments.gap_wellclass --sim-path ./test_data/examples/smeaheia_v2 --well smeaheia.yaml --sim-case TEMP-0 --plot

# 3. for cosmo

$ python -m experiments.gap_wellclass --sim-path ./test_data/examples/cosmo --well cosmo.yaml --sim-case TEMP-0 --plot

```
This will generate an output file `LEG_HIRES.grdecl` in `experiments` directory.

2. To test **gap_plotran.py**, run the following commnad at the directory ``SCREEN``:
```python
$ python -m experiments.gap_pflotran \
    --sim-path ./test_data/examples/cosmo-pflotran \
    --well cosmo.yaml \
    --sim-case1 TEMP-0_NOSIM \
    --sim-case2 TEMP-0 \
    --plot
```
## Test data
The test data are located in the folder `test_data/examples`. Three of its sub-directories, i.el, `cosmo`, `smeaheia_v1` and `smeaheia_v2`, contain the data for testing **gap_wellclass.py**. 

One sub-directory, `cosmo-plotran`, contains configuration parameters for testing **gap_pflotran.py**, i.e., use pfloatran to generate `.EGRID` and `.INIT`. 

Another sub-directory `frigg` contains information for testing deviated wells.

In addition, the PVT values are included in the directory `pvt_contants` for self-consistent testing of pressure-related computes.

## Documentation

To deploy the document to github pages, type the following at the command line:
```
$ mkdocs gh-deploy
```
It may take some minutes until the documentation goes live. And the documentation page can be found at [SCREEN docs](https://redesigned-dollop-m5l6pme.pages.github.io/).
