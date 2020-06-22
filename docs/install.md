# Install

GMGC-mapper runs on Python 3.6-3.8 and requires
[prodigal](https://github.com/hyattpd/Prodigal) to be available for genome
mode.

## Conda install

The easiest way to install GMGC-mapper is through bioconda, which will ensure
all dependencies (including `prodigal`) are installed automatically:

```bash
conda install -c bioconda gmgc-mapper
```

## pip install

Alternatively, `GMGC-mapper` is available from PyPI, so can be installed
through pip:

```bash
pip install GMGC-mapper
```

Note that this does not install `prodigal` (which is necessary for the
genome-based workflow).

## Install from source

Finally, especially if you are retrieving the cutting edge version from Github,
you can install with the standard

```bash
python setup.py install
```
