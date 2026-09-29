#!/usr/bin/env bash
# Rebuild the manuscript docx from its single markdown source.
# All tables and figures are in the main body (no separate supplement).
# progress_report.md is a separate German progress report for the supervisors.
# Run from anywhere: bash manuscript/build_docx.sh

set -euo pipefail

# Resolve script dir so we can cd into manuscript/ regardless of where it's called from
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$script_dir"

FILTER="--lua-filter=color_span.lua"
# reference_thesis.docx: Pandoc-Standardvorlage, in der jede Farbangabe auf 000000
# steht. Ohne sie faerbt Word Ueberschriften, Titel und Links blau (accent1).
REF="--reference-doc=reference_thesis.docx --no-highlight"

# Vorlage reference_thesis.docx setzt die Canvas-Vorgaben: Times New Roman 12pt,
# doppelter Zeilenabstand. Danach in Word nur noch EndNote "Update Citations and
# Bibliography" laufen lassen, das Literaturverzeichnis entsteht dort.
# Dieses Skript schreibt deshalb nach masterarbeit.docx und NICHT
# mehr nach masterarbeit.docx -- sonst wuerde ein Build die Abgabefassung
# ueberschreiben.
echo "Building masterarbeit.docx from masterarbeit.md ..."
pandoc $FILTER $REF masterarbeit.md -o masterarbeit.docx

# progress_report.md ist der interne deutsche Bericht und nicht Teil des
# oeffentlichen Exports. Ohne diese Pruefung bricht der Build dort ab.
if [ -f progress_report.md ]; then
  echo "Building progress_report.docx from progress_report.md ..."
  pandoc $FILTER $REF progress_report.md -o progress_report.docx
else
  echo "progress_report.md not present, skipping."
fi

# Figure and table captions to 10 pt. Pandoc emits them as ordinary body
# paragraphs, so there is no caption style to restyle; this rewrites the sizes
# in the finished docx. Re-run automatically on every build.
echo "Setting figure/table captions to 12 pt ..."
# Canvas: Legenden zaehlen zum Text und bleiben 12pt. Nur der Tabellenkoerper
# darf kleiner sein, den setzt caption_fontsize.py separat.
python3 caption_fontsize.py masterarbeit.docx 12
[ -f progress_report.docx ] && python3 caption_fontsize.py progress_report.docx 12
rm -f masterarbeit.docx.bak progress_report.docx.bak

echo ""
echo "Done. Files:"
ls -la masterarbeit.docx progress_report.docx 2>/dev/null || ls -la masterarbeit.docx
