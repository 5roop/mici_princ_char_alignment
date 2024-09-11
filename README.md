# mici_princ_char_alignment
Character level alignment of the Mići Princ dataset


# Log:
Chapter 13 to be used as a proof of concept.

As of 2024-09-11T08:49:35 the character alignment is implemented, testing revealed that in the presence of music, the alignment will likely be suboptimal.

With this in mind, I opted to go for a targeted attack. I take the wav files and MP_NN.json files. For all the entries (=sentences) in the JSON, I trim the audio and align text to it. This is simpler, removes the problems with the music, and should  be far more accurate as the search space is magnitudes shorter.

In addition, it has been observed that Danijel's alignment is unreliable in cases where speech is followed by silence. In this case, the last element is prolonged into the silence, until there is speech again. This was removed with post-processing using word boundaries obtained with Kaldi. After producing character alignments, I check if a character starts within a word, but ends outside of it. In such cases, I trim the character alignment so that it ends when the word does.

![alignment plot](imgs/chars.png)


# Running the pipeline

## Set up the environment

Prerequisites: You should have `mamba` (or `conda`) installed. All `mamba` commands in this section can be adapted to be run with `conda` instead by simply replacing `mamba` with `conda`.

The environment can be recreated exactly with
```mamba create -f transformers.yml```

This will install all the necessary software. Next we activate it by running

```mamba activate transformers```

## Downloading data

I use a script to download data from CLARIN repo. Run it with:
```bash 0_data_download.sh```

## Running the pipeline

With `mamba` :

```snakemake -j 1 --use-conda```

With `conda`:

```snakemake -j 1 --use-conda --conda-frontend conda```

