# SCREEN

This repository contains codes for SCREEN project.

## Clone the repository
Locate a folder at your local machine that you intend to investigate the codes, and then clone the repository
```
$ git clone https://github.com/equinor/SCREEN
```
Note here `$` is a linux command prompt. There is no need to type it. Now you should have a folder named `SCREEN` at your local machine. then `cd` to it.

It's normal for us to make a new branch if we indend to make some changes of the codes. This can be done with the `-b` option, for example:
```
$ git checkout -b hzh/cleanup
```
This would generate a new branch, named `hzh/cleanup`.

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
├── INSTALLATION.md
├── mkdocs.yml
├── notebooks
│   ├── GaP-WellClass.ipynb
│   ├── LEG_HIRES.grdecl
│   ├── LEG_HIRES.grdecl.smeaheia
│   ├── Pressure-WellClass.ipynb
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
│       └── Readme.md
└── test_data
    ├── examples
    │   ├── cosmo
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
To run the experiments, change directory to `notebooks` and launch jupyter notebooks:
```
$ jupyter-lab
```
You can now play with the notebook `GaP-WellClass.ipynb` to test the integration of `GaP` and `WellClass`. And notebook `Pressure-WellClass.ipynb` to test pressure. And notebook `WellClass_csv_yaml.ipynb` to test loading `.csv` and `.yaml` input files.

The test data is located at folder `test_data`. Three sub-directories, such as `cosmo`, `smeaheia_v1` and `smeaheia_v2`, contain the data for testing. In addition, the PVT values are included in the directory `pvt_contants` for self-consistent testing of pressure-related computes.

## Documentation

To deploy the document to github pages, type the following command:
```
$ mkdocs gh-deploy
```
It may take some minutes until the documentation goes live. And the documentation page can be found at [SCREEN docs](https://redesigned-dollop-m5l6pme.pages.github.io/).
