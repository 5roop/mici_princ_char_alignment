
rule GatherAligned:
    input:
        expand("data/MPcharjson/MP_{n:02}.json", n=range(29))
    output:
        "MP_char_aligned.zip"
    shell:
        "zip {output} {input}"


rule DoOne:
    input:
        wav="data/MPwav/MP_{n}.wav",
        json="data/MPjson/MP_{n}.json"
    output:
        json="data/MPcharjson/MP_{n}.json"
    params:
        cuda="4",
    conda:
        "transformers.yml"
    script:
        "scripts/get_char_alignment.py"
