# Installation

Porechop_ABI can be installed using several methods, that can be separated in two categories:
* Using package managers such as Conda or Pip.
* Installing from sources with pip or maturin.

Except for conda installation, all requirements described in [README](README.md) must be satisfied before install.


# Table of Content
* [Install with Conda](#install-with-conda)
* [Fetching source code](#fetching-source-code)
* [Install with PIP](#install-with-pip)
* [Build a wheel](#build-a-wheel)
* [Build and run without installation](#build-and-run-without-installation)
* [Troubleshooting](#troubleshooting)


## Install with Conda
[Conda](https://docs.conda.io/en/latest/) is a cross platform package / environment management system.
Porechop_ABI was added as a package to the [bioconda](https://bioconda.github.io/) channel to make the installation easier.

If you do not have conda installed, you can get started [here](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).

With Miniconda3 or Anaconda  installed, it is recommended to keep conda up to date.

```bash
conda update conda
```

The different package required to install Porechop_ABI are scattered between the conda-forge and bioconda channels.
You will need to add them both to your available channels:
```bash
conda config --add channels conda-forge
conda config --add channels bioconda
```
An alternative is to specify the channels during installation with *-c* flag after the package name.

With a fresh environment, install Porechop_ABI using either the Anaconda GUI, or the following command lines:

```bash
conda create -n my_env
conda activate my_env
conda install porechop_abi
porechop_abi -h
```


## Fetching source code

All the following methods require you to either download and unzip this repository, or simply clone it:

```bash
git clone https://github.com/bonsai-team/Porechop_ABI.git
cd Porechop_ABI
````


## Install with pip
If you are not using Conda but have pip installed, we recommand using pip to install Porechop_ABI, as it is the second easiest way.
Building from source requires a Rust toolchain (>= 1.85, see [rustup](https://rustup.rs)); pip compiles the extension module automatically.

In the Porechop_ABI folder, use pip as a python module.
```bash
python3.[X] -m pip install .
porechop_abi -h
```
Be sure the python version match the requirement (>= 3.8 and networkx installed)


## Build a wheel
To install on machines without Rust, build a wheel once with [maturin](https://www.maturin.rs) and install it with pip:
```bash
python3.[X] -m pip install maturin
maturin build --release        # wheel written to rust/target/wheels/
python3.[X] -m pip install rust/target/wheels/porechop_abi-*.whl
```
The wheel uses the stable Python ABI, so one wheel per platform works with every Python version >= 3.8.


## Build and run without installation

`maturin develop` compiles the extension module into the `porechop_abi` folder (it needs an active virtual environment). The program can then be executed by calling the `porechop-abi-runner.py` script.

```bash
cd Porechop_ABI
python3.[X] -m venv .venv && . .venv/bin/activate
pip install maturin networkx
maturin develop --release
./porechop-abi-runner.py -h
```


## Troubleshooting
* If pip reports that `cargo` or `rustc` cannot be found, install Rust with [rustup](https://rustup.rs) and open a new shell, or install a prebuilt wheel.

* If you have already the original Porechop installed, Porechop_ABI *will not* replace it. Despite having the same abilities, Porechop_ABI is a different software. It must be run using the `porechop_abi` command, and not just `porechop`.
