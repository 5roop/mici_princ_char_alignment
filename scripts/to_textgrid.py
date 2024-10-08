try:
    alignmentfile = snakemake.input.alignment
    wordalignmentfile = snakemake.input.wordalignment
    graphemealignment = snakemake.input.graphemealignment
    outfile = snakemake.output[0]
except NameError:
    alignmentfile = "kaldidebug/trans_phones_with_symbols.ctm"
    wordalignmentfile = "kaldidebug/ali_words.ctm"
    graphemealignment = "kaldidebug/ali_graphemes.ctm"
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
tg["PhonAlign"] = tier

# Words:
df = pl.read_csv(
    wordalignmentfile,
    separator=" ",
    has_header=False,
    new_columns="sampleid channels start duration word idk idk2".split(),
).with_columns(end=pl.col("start")+pl.col("duration"))

intervals = []
for row in df.iter_rows(named=True):
    try:
        if (previous_end := intervals[-1].xmax) < row["start"]:
                intervals.append(
                    textgrids.Interval(
                        text="", xmin=round(previous_end, 2), xmax=round(row["start"], 2)
                    )
                )
    except IndexError:
         pass
    intervals.append(
        textgrids.Interval(
            text=row["word"], xmin=round(row["start"], 2), xmax=round(row["end"], 2)
        )
    )
tier = textgrids.Tier(data=intervals, xmax=maxs, xmin=mins)
tg["WordAlign"] = tier


# # Graphemes:
# df = pl.read_csv(
#     graphemealignment,
#     separator=" ",
#     has_header=False,
#     new_columns="sampleid channels start duration word idk idk2".split(),
# ).with_columns(end=pl.col("start")+pl.col("duration"))

# intervals = []
# for row in df.iter_rows(named=True):
#     try:
#         if (previous_end := intervals[-1].xmax) < row["start"]:
#                 intervals.append(
#                     textgrids.Interval(
#                         text="", xmin=round(previous_end, 2), xmax=round(row["start"], 2)
#                     )
#                 )
#     except IndexError:
#          pass
#     intervals.append(
#         textgrids.Interval(
#             text=row["word"], xmin=round(row["start"], 2), xmax=round(row["end"], 2)
#         )
#     )
# tier = textgrids.Tier(data=intervals, xmax=maxs, xmin=mins)
# tg["CharAlign"] = tier


tg.write(outfile)
