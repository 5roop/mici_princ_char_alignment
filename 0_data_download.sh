mkdir -p data
cd data
mkdir MPwav MPasrjson MPmp3 MPjson MPexb

cd MPwav
curl --remote-name-all https://www.clarin.si/repository/xmlui/bitstream/handle/11356/1765/MP.wav.tgz
tar xzvf *.tgz
rm *.tgz .*.wav

cd ../MPasrjson
curl --remote-name-all https://www.clarin.si/repository/xmlui/bitstream/handle/11356/1765/MP.asr.json.tgz
tar xzvf *.tgz
rm *.tgz

cd ../MPmp3
curl --remote-name-all https://www.clarin.si/repository/xmlui/bitstream/handle/11356/1765/MP.mp3.tgz
tar xzvf *.tgz
rm *.tgz

cd ../MPjson
curl --remote-name-all https://www.clarin.si/repository/xmlui/bitstream/handle/11356/1765/MP.json.tgz
tar xzvf *.tgz
rm *.tgz

cd ../MPexb
curl --remote-name-all https://www.clarin.si/repository/xmlui/bitstream/handle/11356/1765/MP.exb.tgz
tar xzvf *.tgz
rm *.tgz