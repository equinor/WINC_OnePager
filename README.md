# WINC_OnePager

This repository contains source codes and documentation for WINC_OnePager project.

[![WINC_OnePager-unittest](https://github.com/equinor/WINC_OnePager/actions/workflows/pytest.yaml/badge.svg)](https://github.com/equinor/WINC_OnePager/actions/workflows/pytest.yaml)
[![WINC_OnePager-docs](https://github.com/equinor/WINC_OnePager/actions/workflows/mkdocs.yaml/badge.svg)](https://github.com/equinor/WINC_OnePager/actions/workflows/mkdocs.yaml)
[![WINC_OnePager-lint](https://github.com/equinor/WINC_OnePager/actions/workflows/ruff.yaml/badge.svg)](https://github.com/equinor/WINC_OnePager/actions/workflows/ruff.yaml)

## Installation Instructions
There are two methods to install and run the project - with cloning and without cloning. You can choose depending on whether you need to work with the repository directly or just want to use the package.


### With Poetry (Recommended for Dependency Management in Equinor)
Poetry is Equinor's recommended tool complying with IT policy for Python dependency management.

#### 1. Prerequisites
    
Ensure you have Python `^3.9` installed and accessible in your path.

#### 2. Installing Poetry
    
If you don't have Poetry installed, you can do so with the following command compatible with Python 3.8:

```shell
curl -sSL https://install.python-poetry.org | python3 - --version 1.2.0
```

After installation, verify that Poetry is correctly installed:
```shell
poetry --version
```
#### 3A. Install the Project Using Poetry Without Cloning

To install the project without cloning the repository:

```shell
mkdir my-project
cd my-project
# Create a pyproject.toml file with the content described in the original README, then execute:
poetry install
```

#### 3B. Install the Project Using Poetry With Cloning
To install the project after cloning the repository:

```shell
git clone https://github.com/equinor/WINC_OnePager
cd WINC_OnePager
# Optionally, create a new branch
# Then execute:
poetry install
```

#### 4. Activate the Poetry environment:
```shell
poetry shell
```



### Installation Using pip and a Virtual Environment

The installation of the WINC_OnePager project can be done using `pip`, which is a straightforward approach regardless of cloning. First, ensure that you are in a Python virtual environment to isolate the project dependencies.

#### 1. Creating a Virtual Environment

If you haven't already set up a virtual environment, you can create one using Python's built-in `venv`:

```shell
python -m venv venv_screen
source venv_screen/bin/activate
```

For Windows users, activate the virtual environment with:

```shell
.\venv_screen\Scripts\activate.bat
```
#### 2. Installing the Project

-   **Using `pip install .`**: This method works both when the repository has been cloned and when you have a project directory set up with a pyproject.toml or setup.py file. It installs the current directory as a package along with its dependencies:

    ```shell
    pip install .
    ```

    This command tells pip to install the current directory (i.e., the project) as a package.

-   **Using `pip install -r requirements.txt`**: This method is specific to situations where the repository has been cloned. It will install the dependencies specified in the requirements.txt file:
    ```shell
    pip install -r requirements.txt
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

1. To test **well_sketch.py**, run either of the followings:
```
# 1. for smeaheia_v1

python -m experiments.well_sketch_pressure --config-file ./test_data/examples/smeaheia_v1/smeaheia.yaml -pvt ./test_data/pvt_constants 

# 3. for wildcat

python -m experiments.well_sketch_pressure --config-file ./test_data/examples/wildcat/wildcat.yaml -pvt ./test_data/pvt_constants 
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
│   ├── GaP_input_Wildcat_v3.csv
│   └── wildcat.yaml
├── wildcat-pflotran
│   └── wildcat.yaml
└── wildcat-pflotran-2
    └── wildcat.yaml
```

## Unit testing and code coverage
We are using `pytest` for unit testing and code coverage. The unit testing utilizes `wildcat` as the testing example. So please make sure the saved .pkl files in ```test_data/examples/wildcat/pytest``` exists and is updated. Here is a commandline example:
```python
python -m pytest tests
```
This will report the unit testing results. And the following will report not only unit testing but also code coverage:
```python
python -m pytest --cov tests
```
or a litle bit more complex command:
```python
python -m pytest --cov --cov-branch --cov-report term-missing tests
```

## Documentation

The document can be automatically generated and deployed to github pages. To do that, type the following at the command line:
```
mkdocs gh-deploy
```
It may take some minutes until the documentation goes live. And the generated documentation page can be found at [WINC_OnePager docs](https://redesigned-dollop-m5l6pme.pages.github.io/).

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
├── poetry.lock
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

