#!/bin/bash
# Telluric correction for WASP-18 / CRIRES+ Y1029
# 3-pass iterative workflow, no FTS gas cell (Y-band)

set -e
cd "$(dirname "$0")/../.."
source .venv/bin/activate

DATA=data/WASP18
FILES="$DATA/cr2res*.fits"
OSET="1:28"  # 9 orders x 3 detectors for Y-band

echo "=== Pass 1: initial fit ==="
viper "$FILES" -inst CRIRES \
  -createtpl -nocell -fts None -telluric add -tsig 10 -tpl_wave tell \
  -deg_norm 0 -deg_wave 2 -deg_bkg 0 -oversampling 1 -kapsig 15 \
  -oset $OSET -tag $DATA/tpl1

echo ""
echo "=== Pass 2: refine with template ==="
viper "$FILES" $DATA/tpl1_tpl.fits -inst CRIRES \
  -createtpl -nocell -fts None -telluric add -tsig 1 -tpl_wave tell \
  -deg_norm 2 -deg_wave 2 -deg_bkg 0 -oversampling 1 -kapsig 15 4.5 \
  -oset $OSET -lookctpl -tag $DATA/tpl2



