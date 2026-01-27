# WINC_OnePager

This repository contains source codes and documentation for the WINC_OnePager project - an analytical wellbore pressure estimation tool for CO2 leakage risk assessment through abandoned wells.

> **Note:** This project was initially developed as part of the [SCREEN project](https://colab.equinor.com/technologies/4FAAF5BD-19C3-46A3-ACB6-5D38DD2C7A89/summary), but the analytical pressure estimation module has been developed and maintained as a separate, standalone entity.

[![WINC_OnePager-unittest](https://github.com/equinor/WINC_OnePager/actions/workflows/pytest.yaml/badge.svg)](https://github.com/equinor/WINC_OnePager/actions/workflows/pytest.yaml)
[![WINC_OnePager-docs](https://github.com/equinor/WINC_OnePager/actions/workflows/mkdocs.yaml/badge.svg)](https://github.com/equinor/WINC_OnePager/actions/workflows/mkdocs.yaml)
[![WINC_OnePager-lint](https://github.com/equinor/WINC_OnePager/actions/workflows/ruff.yaml/badge.svg)](https://github.com/equinor/WINC_OnePager/actions/workflows/ruff.yaml)

## Installation Instructions

### Prerequisites: Python Installation

This code has been tested with Python versions 3.9 through 3.12. The recommended way to install and manage Python is using [uv](https://docs.astral.sh/uv/).

> **Equinor users:** Please follow the internal guidelines available at: https://wiki.equinor.com/wiki/Using_Python_on_Windows_11_with_uv

**For Windows users**, install uv using winget:

```shell
winget install --id=astral-sh.uv -e
```

**For Linux/macOS users**:

```shell
curl -LsSf https://astral.sh/uv/install.sh | sh
```

Once uv is installed, reload your terminal and install Python:

```shell
uv python install 3.12
```

For other installation methods, see the [uv installation docs](https://docs.astral.sh/uv/getting-started/installation/).

### Clone and Install

```shell
git clone https://github.com/equinor/WINC_OnePager
cd WINC_OnePager
uv sync
```

This will create a virtual environment and install all dependencies including dev, docs, test, and dash groups.

### Run Commands

Use `uv run` to execute commands in the project environment:

```shell
uv run python -m experiments.well_sketch_pressure --config-file ./test_data/examples/wildcat/wildcat.yaml
```

Or activate the environment:

```shell
source .venv/bin/activate  # Linux/macOS
# or
.venv\Scripts\activate     # Windows
```


### With pip

You can also install using pip in a virtual environment:

#### 1. Create and Activate a Virtual Environment

```shell
python -m venv .venv
source .venv/bin/activate  # Linux/macOS
# or
.venv\Scripts\activate     # Windows
```

#### 2. Install the Project

```shell
pip install .
```

Or install with development dependencies:

```shell
pip install -e ".[dev,docs,test]"
```


## Experiments
There are at least two ways to make experimenal runs of the codes. One is to run the experiments with Jupyter lab (notbeooks folder), and the other is commandline option. While Jupyter notebooks are mainly for QC tests and research purposes, the commandline option is aiming for production run.

### 1. Jupyter notebooks
Jupyter notebooks are located in directory `notebooks`. To test its functionaries, change current directory to `notebooks` and launch jupyter notebooks at the commandline:
```
jupyter-lab
```
Or if you prefer, you can run these Jupyter notebooks from Microsoft's VS code.

There exist several Jupyter notebooks in the directory:

- Notebook **Pressure-WellClass.ipynb** is used to test pressure. 
- Notebook **WellClass_csv_yaml.ipynb** is used to test pressure and loading `.csv` and `.yaml` input files.

### 2. Commandline option
Two python scripts for commandline option are available in directory `experiments`. One script, **well_sketch.py**, can be used for generating a well sketch, **well_sketch_pressure.py** can be used for generating both a well sketch and a pressure plot. 

The followings are some of  the sample runs. In either way, you should run the python script inside the ```WINC_OnePager``` directory. 

1. To test **well_sketch_pressure.py**, run either of the followings:
```shell
# for smeaheia_v1
uv run python -m experiments.well_sketch_pressure --config-file ./test_data/examples/smeaheia_v1/smeaheia.yaml

# for wildcat
uv run python -m experiments.well_sketch_pressure --config-file ./test_data/examples/wildcat/wildcat.yaml
```

### 3. Test data
In order for a quick test of the codes, we include some test dataset in the folder `test_data/examples`. 

```
├── frigg
│   ├── GaP_input_Frigg_v3.csv
│   └── X_5687dev.txt
├── simple_well
│   ├── Simple_well.csv
│   └── Simple_well.yaml
├── smeaheia_v1
│   ├── GaP_input_Smeaheia_v3.csv
│   └── smeaheia.yaml
├── wildcat
    ├── GaP_input_Wildcat_v3.csv
    └── wildcat.yaml


## Unit testing and code coverage
We are using `pytest` for unit testing and code coverage. The unit testing utilizes `wildcat` as the testing example. So please make sure the saved .pkl files in ```test_data/examples/wildcat/pytest``` exists and is updated. Here is a commandline example:
```shell
uv run pytest tests
```
This will report the unit testing results. And the following will report not only unit testing but also code coverage:
```shell
uv run pytest --cov tests
```
or a little bit more complex command:
```shell
uv run pytest --cov --cov-branch --cov-report term-missing tests
```

## Documentation

The documentation can be automatically generated and deployed to GitHub Pages. To do that, type the following at the command line:
```shell
uv run mkdocs gh-deploy
```
It may take some minutes until the documentation goes live. And the generated documentation page can be found at [WINC_OnePager docs](https://equinor.github.io/WINC_OnePager/).

## The code structures

The following represents the current code structures:

```
.
├── CITATION.cff
├── experiments
│   ├── __init__.py
│   ├── well_pressure_tables.py
│   ├── well_sketch_pressure.py
│   └── well_sketch.py
├── LICENSE
├── mkdocs.yml
├── notebooks
│   ├── Pressure-WellClass.ipynb
│   ├── PVT_data.ipynb
│   ├── WellClass_csv_yaml.ipynb
│   └── WellClass-onepager.ipynb
├── pyproject.toml
├── README.md
├── requirements.txt
├── requirements.txt.frozen
├── src
│   ├── __init__.py
│   └── WellClass
│       ├── __init__.py
│       ├── libs
│       ├── notebooks
│       ├── README.md
│       └── tools
├── test_data
│   └── examples
│       ├── frigg
│       ├── simple_well
│       ├── smeaheia_v1
│       ├── wildcat
│       ├── wildcat-pflotran
│       └── wildcat-pflotran-2
└── tests
    ├── conftest.py
    └── well_class
        └── test_well_class.py
```

It was generated with the linux command `tree`:
```shell
tree -I 'docs|site|venv_screen|*pycache*|Equinor*|originals' -L 3
```

