#!/usr/bin/env python3
"""Set figure and table captions in a pandoc-produced docx to a given font size.

Pandoc writes our captions as ordinary body paragraphs, so there is no Word
caption style to restyle. This walks document.xml and shrinks:
  - paragraphs starting with "Table N." or "Figure N."  (the caption proper)
  - the explanatory paragraph directly following each table (the table footnote)
Run as: caption_fontsize.py <docx> [pt]
"""
import re, shutil, sys, zipfile
from xml.etree import ElementTree as ET

W = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
ET.register_namespace("w", W)
def q(t): return f"{{{W}}}{t}"

# w:sz must precede these siblings inside w:rPr; insert before the first one present.
AFTER = [q(x) for x in ("u", "effect", "bdr", "shd", "fitText", "vertAlign",
                        "rtl", "cs", "em", "lang", "eastAsianLayout")]

def set_size(rpr, half_pt):
    for tag in ("sz", "szCs"):
        for old in rpr.findall(q(tag)):
            rpr.remove(old)
    pos = len(rpr)
    for i, child in enumerate(rpr):
        if child.tag in AFTER:
            pos = i
            break
    for off, tag in enumerate(("sz", "szCs")):
        el = ET.Element(q(tag)); el.set(q("val"), str(half_pt))
        rpr.insert(pos + off, el)

def shrink(par, half_pt):
    ppr = par.find(q("pPr"))
    if ppr is None:
        ppr = ET.Element(q("pPr")); par.insert(0, ppr)
    mark = ppr.find(q("rPr"))
    if mark is None:
        mark = ET.SubElement(ppr, q("rPr"))
    set_size(mark, half_pt)
    for run in par.findall(q("r")):
        rpr = run.find(q("rPr"))
        if rpr is None:
            rpr = ET.Element(q("rPr")); run.insert(0, rpr)
        set_size(rpr, half_pt)

def text_of(par):
    return "".join(t.text or "" for t in par.iter(q("t")))

def main(path, pt=10):
    half = int(round(pt * 2))
    with zipfile.ZipFile(path) as z:
        items = {n: z.read(n) for n in z.namelist()}
    doc = items["word/document.xml"]
    # Preserve the original namespace prefixes. Without this ElementTree invents
    # ns1, ns2, ... on output, renaming the math (m:) and drawing (wp:, a:, pic:)
    # prefixes throughout the document.
    for prefix, uri in re.findall(rb'xmlns:([A-Za-z0-9]+)="([^"]+)"', doc[:4000]):
        ET.register_namespace(prefix.decode(), uri.decode())
    root = ET.fromstring(doc)
    body = root.find(q("body"))
    caption = re.compile(r"^\s*(Table|Figure)\s+\d+[.:]")
    n_cap = n_note = n_tbl = 0
    # A table's explanatory text sits before the table for Tables 1 and 3 and
    # after it for the rest, so everything from a "Table N." caption up to and
    # including the first paragraph after its table counts as caption text.
    # Any heading, image or other caption ends the run. Documents whose tables
    # carry no caption (the German progress report) are left untouched.
    state = None   # None | "pre" (caption seen, table not yet) | "post" (table passed)
    for child in list(body):
        if child.tag == q("tbl"):
            if state == "pre":
                # Table body too, but only for captioned tables, so documents
                # without table captions keep their own formatting.
                for par in child.iter(q("p")):
                    shrink(par, half)
                n_tbl += 1
                state = "post"
            else:
                state = None
            continue
        if child.tag != q("p"):
            continue
        txt = text_of(child)
        if not txt.strip():
            continue
        if caption.match(txt):
            shrink(child, half); n_cap += 1
            state = "pre" if txt.lstrip().startswith("Table") else None
        elif state == "pre":
            shrink(child, half); n_note += 1          # explanation before the table
        elif state == "post":
            shrink(child, half); n_note += 1; state = None
        else:
            state = None
    items["word/document.xml"] = ET.tostring(root, encoding="UTF-8", xml_declaration=True)
    shutil.copy(path, path + ".bak")
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as z:
        for name, data in items.items():
            z.writestr(name, data)
    print(f"  {path}: {n_cap} captions + {n_note} notes + {n_tbl} table bodies set to {pt}pt")

if __name__ == "__main__":
    main(sys.argv[1], float(sys.argv[2]) if len(sys.argv) > 2 else 10)
