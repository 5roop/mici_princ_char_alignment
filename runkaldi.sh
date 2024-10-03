
tmpdir=kaldidebug/
mkdir -p "$tmpdir"
echo "speaker1 sampleid">kaldidebug/spk2utt
echo "sampleid data/MPmp3_wav/MP_27_121.62-132.06.wav" > kaldidebug/wavscp
echo "sampleid jenu večer je pozabiv stakleno zvono ili je morda ovca po noći tiho odašla ta put se si zvončići store va suze to je kako bi reć vela tajna " > kaldidebug/text


wavscp=kaldidebug/wavscp
spk2utt=kaldidebug/spk2utt
text=kaldidebug/text

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

# # # Phoneme alignment:
$kaldi/src/latbin/lattice-align-phones --replace-output-symbols=true $models/nnet3/tdnn1a_sp/final.mdl ark:$tmpdir/ali.lat ark:$tmpdir/phone_aligned.lats
$kaldi/src/latbin/lattice-to-ctm-conf --inv-acoustic-scale=10 --decode-mbr ark:$tmpdir/phone_aligned.lats $tmpdir/trans_phones.ctm
$kaldi/egs/wsj/s5/utils/int2sym.pl -f 5 $tmpdir/phones.txt $tmpdir/trans_phones.ctm > $tmpdir/trans_phones_with_symbols.ctm