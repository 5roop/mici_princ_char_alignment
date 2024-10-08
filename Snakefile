from pathlib import Path
import json
import polars as pl


def preprocess_text(s: str) -> str:
    from string import punctuation

    for p in punctuation + "…":
        s = s.replace(p, "")
    s = s.casefold()
    s = s.replace("x", "ks").replace("y", "i").replace("werthu", "vertu")
    s = " ".join(s.split())
    return s


jsons = list(Path("data/MPasrjson").glob("MP_*.asr.json"))
data = []
for i in jsons:
    data += json.loads(i.read_text())
expected_output = [Path("output", i["audio"]).with_suffix(".wav") for i in data]

wildcard_constraints:
    file="MP_\d+_\d+\.*\d*-\d+\.*\d*"
rule gather_output:
    input:
        wav=expected_output,
        tg=[i.with_suffix(".TextGrid") for i in expected_output]
rule prep_files:
    input: "output/{file}.wav"
    output:
        s = temp("tmp/{file}.spk2utt"),
        t = temp("tmp/{file}.text"),
        w = temp("tmp/{file}.wav.scp"),
        g = temp("tmp/{file}.text_for_graphemes")
    run:
        for entry in data:
            if str(Path("output", entry["audio"]).with_suffix(".wav")) == input[0]:
                Path(output.s).write_text("speaker1 sampleid")
                Path(output.t).write_text("sampleid "+preprocess_text(entry["normalized_text"]))
                Path(output.g).write_text("sampleid "+ " ".join(preprocess_text(entry["normalized_text"])).replace("   ", " "))
                Path(output.w).write_text("sampleid "+input[0])
                break

rule do_kaldi:
    input:
        spk2utt="tmp/{file}.spk2utt",
        text="tmp/{file}.text",
        wavscp="tmp/{file}.wav.scp",
        text_for_graphemes="tmp/{file}.text_for_graphemes",
    output:
        temp(directory("tmp/{file}")),
        temp("tmp/{file}/trans_phones_with_symbols.ctm"),
        temp("tmp/{file}/ali_words.ctm"),
        temp("tmp/{file}/ali_graphemes.ctm")
    shell:
        """
        wavscp={input.wavscp}
        spk2utt={input.spk2utt}
        text={input.text}
        tmpdir={output[0]}

        models=models
        kaldi=/opt/kaldi

        export LC_ALL=C

        mkdir -p $tmpdir

        cut -f2- -d' ' $text | tr ' ' '\n' | sort -u > $tmpdir/wlist

        /home/rupnik/anaconda3/envs/p37/bin/python lexicon.py $tmpdir/wlist $models/phonetisaurus-hr/model.fst $tmpdir
        $kaldi/src/featbin/compute-mfcc-feats  --config=$models/nnet3/conf/mfcc.conf scp:$wavscp ark:$tmpdir/mfcc.ark
        $kaldi/src/online2bin/ivector-extract-online2 --config=$models/nnet3/conf/ivector.conf ark:$spk2utt ark:$tmpdir/mfcc.ark ark:$tmpdir/ivec.ark
        $kaldi/egs/wsj/s5/utils/sym2int.pl -f 2- $tmpdir/words.txt $text  > $tmpdir/text.int

        $kaldi/src/bin/compile-train-graphs --read-disambig-syms=$tmpdir/disambig.int $models/nnet3/tdnn1a_sp/tree $models/nnet3/tdnn1a_sp/final.mdl $tmpdir/L.fst ark:$tmpdir/text.int ark:$tmpdir/graphs.fsts
        $kaldi/src/nnet3bin/nnet3-latgen-faster --allow_partial=True --online-ivectors=ark:$tmpdir/ivec.ark --online-ivector-period=10 $models/nnet3/tdnn1a_sp/final.mdl ark:$tmpdir/graphs.fsts ark:$tmpdir/mfcc.ark ark:$tmpdir/ali.lat
        $kaldi/src/latbin/lattice-align-words $tmpdir/word_boundary.int $models/nnet3/tdnn1a_sp/final.mdl ark:$tmpdir/ali.lat ark:- | $kaldi/src/latbin/lattice-to-ctm-conf ark:- - | $kaldi/egs/wsj/s5/utils/int2sym.pl -f 5 $tmpdir/words.txt - > $tmpdir/ali.ctm

        # # # My additions:
        cp $tmpdir/ali.ctm $tmpdir/ali_words.ctm
        $kaldi/src/latbin/lattice-align-phones --replace-output-symbols=true $models/nnet3/tdnn1a_sp/final.mdl ark:$tmpdir/ali.lat ark:$tmpdir/phone_aligned.lats
        $kaldi/src/latbin/lattice-to-ctm-conf --inv-acoustic-scale=10 --decode-mbr ark:$tmpdir/phone_aligned.lats $tmpdir/trans_phones.ctm
        $kaldi/egs/wsj/s5/utils/int2sym.pl -f 5 $tmpdir/phones.txt $tmpdir/trans_phones.ctm > $tmpdir/trans_phones_with_symbols.ctm

        touch $tmpdir/ali_graphemes.ctm
        exit
        # # # # Grapheme alignment:
        # text={input.text_for_graphemes}
        # cut -f2- -d' ' $text | tr ' ' '\n' | sort -u > $tmpdir/wlist
        # /home/rupnik/anaconda3/envs/p37/bin/python lexicon.py $tmpdir/wlist $models/phonetisaurus-hr/model.fst $tmpdir
        # $kaldi/src/featbin/compute-mfcc-feats  --config=$models/nnet3/conf/mfcc.conf scp:$wavscp ark:$tmpdir/mfcc.ark
        # $kaldi/src/online2bin/ivector-extract-online2 --config=$models/nnet3/conf/ivector.conf ark:$spk2utt ark:$tmpdir/mfcc.ark ark:$tmpdir/ivec.ark
        # $kaldi/egs/wsj/s5/utils/sym2int.pl -f 2- $tmpdir/words.txt $text  > $tmpdir/text.int

        # $kaldi/src/bin/compile-train-graphs --read-disambig-syms=$tmpdir/disambig.int $models/nnet3/tdnn1a_sp/tree $models/nnet3/tdnn1a_sp/final.mdl $tmpdir/L.fst ark:$tmpdir/text.int ark:$tmpdir/graphs.fsts
        # $kaldi/src/nnet3bin/nnet3-latgen-faster --allow_partial=True --online-ivectors=ark:$tmpdir/ivec.ark --online-ivector-period=10 $models/nnet3/tdnn1a_sp/final.mdl ark:$tmpdir/graphs.fsts ark:$tmpdir/mfcc.ark ark:$tmpdir/ali.lat
        # $kaldi/src/latbin/lattice-align-words $tmpdir/word_boundary.int $models/nnet3/tdnn1a_sp/final.mdl ark:$tmpdir/ali.lat ark:- | $kaldi/src/latbin/lattice-to-ctm-conf ark:- - | $kaldi/egs/wsj/s5/utils/int2sym.pl -f 5 $tmpdir/words.txt - > $tmpdir/ali.ctm
        # cp $tmpdir/ali.ctm $tmpdir/ali_graphemes.ctm
        """

rule do_tg_compilation:
    input:
        alignment = "tmp/{file}/trans_phones_with_symbols.ctm",
        wordalignment = "tmp/{file}/ali_words.ctm",
        graphemealignment = "tmp/{file}/ali_graphemes.ctm",
        thedir = "tmp/{file}"
    output:
        "output/{file}.TextGrid"
    conda:
        "miciprincalign.yml"
    script:
        "scripts/to_textgrid.py"
rule cp_wav:
    input:
        "data/MPmp3_wav/{file}.wav"
    output:
        "output/{file}.wav"
    shell:
        "cp {input} {output}"
