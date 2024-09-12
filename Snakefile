
rule GatherExbs:
    input:
        expand("data/MPinspection/MP_{n:02}.{ext}", n=range(3), ext=["exb", "wav"]),
    output:
        "MP_wav+exb.zip",
    shell:
        "zip {output} {input}"


rule GatherAligned:
    input:
        expand("data/MPcharjson/MP_{n:02}.json", n=range(29)),
    output:
        "MP_char_aligned.zip",
    shell:
        "zip {output} {input}"


rule ProduceEXB:
    input:
        json="data/MPcharjson/MP_{n}.json",
        exb="data/MPexb/MP_{n}.exb",
    output:
        exb="data/MPinspection/MP_{n}.exb",
    conda:
        "transformers.yml"
    script:
        "scripts/produce_exbs.py"


rule CopyAudio:
    input:
        "data/MPwav/MP_{n}.wav",
    output:
        "data/MPinspection/MP_{n}.wav",
    shell:
        "ffmpeg -i {input} -ar 16000 -ac 1 {output}"


rule DoOne:
    input:
        wav="data/MPwav/MP_{n}.wav",
        json="data/MPjson/MP_{n}.json",
    output:
        json="data/MPcharjson/MP_{n}.json",
    params:
        cuda="4",
    conda:
        "transformers.yml"
    script:
        "scripts/get_char_alignment.py"
