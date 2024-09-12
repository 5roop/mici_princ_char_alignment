try:
    # wavpath = snakemake.input.wav
    jsonpath = snakemake.input.json
    exbinpath = snakemake.input.exb
    outputpath = snakemake.output.exb
except NameError:
    # wavpath = "data/MPwav/MP_00.wav"
    jsonpath = "data/MPcharjson/MP_00.json"
    exbinpath = "data/MPexb/MP_00.exb"
    outputpath = "brisi.exb"

from lxml import etree as ET
from pathlib import Path
import json

charaligned = json.loads(Path(jsonpath).read_text())
doc = ET.fromstring(Path(exbinpath).read_bytes())
timeline = doc.find(".//common-timeline")
last_tier = doc.findall(".//tier")[-1]
new_tier = ET.Element(
    "tier",
    attrib={
        "id": "characters",
        "category": "a",
        "type": "a",
        "display-name": "Characters",
    },
)
newevent_counter = 0
for i in charaligned:
    chars = i["chars"]
    for char in chars:
        ts = char["time_s"]
        te = char["time_e"]
        c = char["char"]
        tis = f"Tchar{newevent_counter}"
        tie = f"Tchar{newevent_counter + 1}"
        timeline.append(ET.Element("tli", attrib={"id": tis, "time": str(ts)}))
        timeline.append(ET.Element("tli", attrib={"id": tie, "time": str(te)}))
        newevent = ET.Element("event", attrib={"start": tis, "end": tie})
        newevent.text = c
        new_tier.append(newevent)
        newevent_counter += 2
last_tier.getparent().insert(
    last_tier.getparent().index(last_tier) + 1,
    new_tier,
)
timeline[:] = sorted(timeline, key=lambda e: float(e.get("time", "0.0")))
ET.indent(doc)

Path(outputpath).write_bytes(
    ET.tostring(
        doc,
        pretty_print=True,
        encoding="utf8",
        xml_declaration=True,
    ).replace(b"utf8", b"utf-8")
)
