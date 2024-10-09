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
# import textgrids
from praatio import textgrid
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

tg = textgrid.Textgrid()
intervals = []
for row in df.iter_rows(named=True):
    intervals.append(
        (
             round(row["start"], 2),
             round(row["end"], 2),
             row["phoneme"]
         )
    )
tier = textgrid.IntervalTier("PhoneAlign", intervals, mins, maxs)
tg.addTier(tier)

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
        if (previous_end := intervals[-1][1]) < row["start"]:
                intervals.append(
                    (
                         round(previous_end, 2),
                        round(row["start"], 2),
                        "( )",
                     )
                )
    except IndexError:
         pass
    intervals.append(
        (
             round(row["start"], 2),
             round(row["end"], 2),
             row["word"],
        )
    )
tier = textgrid.IntervalTier("WordAlign", intervals, mins, maxs)
tg.addTier(tier)

tg.save(outfile, format="long_textgrid", includeBlankSpaces=True)
