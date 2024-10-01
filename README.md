# mici_princ_char_alignment
Character level alignment of the Mići Princ dataset


# Log:
## 2024-10-01T15:03:47
Implemented kaldi grapheme alignment pipeline.

# Running the pipeline

## Set up the environment

Prerequisites: You should have `mamba` (or `conda`) installed.

All `mamba` commands in this section can be adapted to be run with `conda` instead by simply replacing `mamba` with `conda`.

The environments needed here can be recreated exactly with

```mamba create -f p37.yml```
and
```mamba create -f miciprincalign.yml```

This will install all the necessary software. Next we activate it by running

```mamba activate miciprincalign```

Check where your python 3.7 environment lives with `mamba activate p37; which python; mamba deactivate` and edit the location of the python in file Snakemake, line 67.

## Downloading data

I use a script to download data from CLARIN repo. Run it with:

```bash 0_data_download.sh```

## Running the pipeline

With `mamba`:

```mamba activate miciprincalign; snakemake -j 1 --use-conda```

With `conda`:

```conda activate miciprincalign; snakemake -j 1 --use-conda --conda-frontend conda```

