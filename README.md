# SCREEN

This repository contains source codes and documentation for SCREEN_WellView project.

[![SCREEN-unittest](https://github.com/equinor/SCREEN/actions/workflows/pytest.yaml/badge.svg)](https://github.com/equinor/SCREEN/actions/workflows/pytest.yaml)
[![SCREEN-docs](https://github.com/equinor/SCREEN/actions/workflows/mkdocs.yaml/badge.svg)](https://github.com/equinor/SCREEN/actions/workflows/mkdocs.yaml)
[![SCREEN-lint](https://github.com/equinor/SCREEN/actions/workflows/ruff.yaml/badge.svg)](https://github.com/equinor/SCREEN/actions/workflows/ruff.yaml)

## Clone the repository
Locate a folder at your local machine that you intend to investigate the codes, and then clone the repository
```
git clone https://github.com/equinor/SCREEN_WellView
```
By this time you should have a folder named `SCREEN_WellView` at your local machine. Now change the directory with linux command:
```
cd SCREEN
```

It's normal for us to make a new branch if we indend to make some changes of the codes. This can be done with the `-b` option, for example:
```
git checkout -b xyz/cleanup
```
This would generate a new branch, named `xyz/cleanup`. Here the branch name is created by concatenating a short name, such as `xyz`,  of `equinor` account with a feature description `cleanup`. There is no need to follow this convention. You could simply pick any branch name as long as it makes sense. However, please note branch names have limitations.

## Virtual environment

It's a common practice to work on a project within a python virtual environment. I have been using python's builtin module `venv` for a long while. So I am going to stick to it here as an example to set up the virtual environment. But you are free to use any other virtual environment setups that you feel comfortable with, such as `conda`, etc. 
```
python -m venv venv_screen
```
This will build a virtual environment `venv_screen`. You only need to creat it once.

To activate this `venv_screen`, run the following command: 
```
source venv_screen/bin/activate.csh
```
when your linux Shell is `csh`. 

If you are using `bash` or plain `sh`, you can activate it with the following command:
```
source venv_screen/bin/activate
```

We pack needed python packages into a file, such as `requirements.txt`. And install those python packages to this virtual environment by running the following command:
```
pip install -r requirements.txt
```
You should now be ready to play with the source codes.

## Poetry for dependency management
`poetry`, as an alternative tool for dependency management, is required in Equinor to be compliant with IT policy.

### 1. Installation of poetry
In case that `poetry` is not available in the system, you may have to download and install it yourself. Somehow the newer version of `poetry` doesn't work with python `3.8`. So we have to get an earlier of `poetry`, such as `1.2.0` with the following command: 
```
curl -sSL https://install.python-poetry.org | python3 - --version 1.2.0
```
This should install `poetry` in your home directory. For example, in my case, it is installed at `/private/hzh/.local/bin`. Make sure to set this path as an environment variable in your `.cshrc` or `.bashrc` file so that you can access `poetry` from anywhere.

You can test that everything is set up by executing:

```
poetry --version
```

### 2. Installation of python dependencies
To install python dependencies, run the following command<strong>*</strong>:
```
poetry install
```

To check what packages have been installed, try the following command:
```
poetry show --tree
```

To show where the virtual environment is located, run
```
poetry env info
```
This will display the python path for the activated virtual environment. The executable path can be used in VS Code to set up python interpreter path for builtin jupyter notebooks.

The virtual envrionment can be activated with the command<strong>*</strong>:
```
poetry shell
```
This generated shell will be used to run standalone python scripts.

> Note: For whatever reason, in case you decide to uninstall `poetry` package, you can do it with
> ```
> curl -sSL https://install.python-poetry.org | python3 - --uninstall
> ```

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

The followings are some of  the sample runs. In either way, you should run the python script inside the ```SCREEN``` directory. 

1. To test **well_sketch.py**, run either of the followings:
```
# 1. for smeaheia_v1

python -m experiments.well_sketch_pressure --config-file ./test_data/examples/smeaheia_v1/smeaheia.yaml -pvt ./test_data/pvt_constants 

# 3. for wildcat

python -m experiments.well_sketch_pressure --config-file ./test_data/examples/wildcat/wildcat.yaml -pvt ./test_data/pvt_constants 


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
It may take some minutes until the documentation goes live. And the generated documentation page can be found at [SCREEN docs](https://redesigned-dollop-m5l6pme.pages.github.io/).

## The code structures

The following represents the current code structures:

```
├── CITATION.cff
├── LICENSE
├── README.md
├── experiments
│   ├── __init__.py
│   ├── well_sketch.py
│   └── well_sketch_pressure.py
├── mkdocs.yml
├── notebooks
│   ├── Pressure-WellClass.ipynb
│   └── WellClass_csv_yaml.ipynb
├── poetry.lock
├── pyproject.toml
├── requirements.txt
├── src
│   ├── WellClass
│   │   ├── README.md
│   │   ├── __init__.py
│   │   ├── libs
│   │   └── notebooks
│   └── __init__.py
├── test_data
│   ├── examples
│   │   ├── frigg
│   │   ├── simple_well
│   │   ├── smeaheia_v1
│   │   ├── wildcat
│   │   ├── wildcat-pflotran
│   │   └── wildcat-pflotran-2
│   └── pvt_constants
│       ├── pressure.txt
│       ├── rho_co2.txt
│       ├── rho_h2o.txt
│       └── temperature.txt
└── tests
    ├── conftest.py
    └── well_class
        ├── conftest.py
        └── test_well_class.py
```
It was generated with the linux command `tree`:
```shell
tree -I 'docs|site|venv_screen|*pycache*|Equinor*|originals' -L 3
```

