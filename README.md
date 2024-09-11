# mici_princ_char_alignment
Character level alignment of the Mići Princ dataset


# First steps:
Chapter 13 to be used as a proof of concept.

As of 2024-09-11T08:49:35 the character alignment is implemented, testing revealed that in the presence of music, the alignment will likely fail.

![alignment plot](imgs/chars.png)


# Running the pipeline

## Set up the environment

Prerequisites: You should have mamba (or conda) installed. All mamba commands in this section can be adapted to be run with conda instead by simply replacing `mamba` with `conda`.

The environment can be recreated exactly with
```mamba create -f transformers.yml```

This will install all the necessary software. Next we activate it by running

```mamba activate transformers```

## Running the pipeline

With mamba :

```snakemake -j 1 --use-conda```

With conda:

```snakemake -j 1 --use-conda --conda-frontend conda```

