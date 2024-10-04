try:
    alignmentfile = snakemake.input.alignment
    wordalignmentfile = snakemake.input.wordalignment
    # wavfile = snakemake.input.wav
    outfile = snakemake.output[0]
except NameError:
    alignmentfile = "kaldidebug/trans_phones_with_symbols.ctm"
    wordalignmentfile = "kaldidebug/ali.ctm"
    # wavfile = "data/MPmp3_wav/MP_01_0.33-12.69.wav"
    outfile = "brisi.TextGrid"

import polars as pl
import textgrids

df = (
    pl.read_csv(
        alignmentfile,
        separator=" ",
        has_header=False,
        new_columns="sampleid channels start duration phoneme idk idk2".split(),
    )
    .with_columns((pl.col("start") + pl.col("duration")).alias("end"))
    .select(pl.exclude("^idk.*$"))
)


_, mins, maxs = (
    df.group_by("sampleid")
    .agg(pl.col("start").min().alias("min"), pl.col("end").max().alias("max"))
    .to_numpy()
    .reshape(-1)
)

tg = textgrids.TextGrid()
tg.xmax = maxs
tg.xmin = mins

intervals = []
for row in df.iter_rows(named=True):
    intervals.append(
        textgrids.Interval(
            text=row["phoneme"], xmin=round(row["start"], 2), xmax=round(row["end"], 2)
        )
    )
tier = textgrids.Tier(data=intervals, xmax=maxs, xmin=mins)
tg["CharAlign"] = tier

df = pl.read_csv(
    wordalignmentfile,
    separator=" ",
    has_header=False,
    new_columns="sampleid channels start duration word idk idk2".split(),
).with_columns(end=pl.col("start")+pl.col("duration"))

intervals = []
for row in df.iter_rows(named=True):
    intervals.append(
        textgrids.Interval(
            text=row["word"], xmin=round(row["start"], 2), xmax=round(row["end"], 2)
        )
    )
tier = textgrids.Tier(data=intervals, xmax=maxs, xmin=mins)
tg["WordAlign"] = tier


tg.write(outfile)
