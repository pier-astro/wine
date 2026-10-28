#!/usr/bin/env bash
# Run after initializing HEASoft; products are written in example/.
set -euo pipefail
cd "$(dirname "$0")/.."

exposure=100000
seed=456
min_counts=10

mkdir -p results
rm -f simulated.fak simulated.grp
xspec - > results/simulation.log 2>&1 <<EOF
query yes
lmod wine ../xspec
@model.xcm
energies 0.01 40 6000 log
xset seed $seed
fakeit none
response.rmf
response.arf
y

simulated.fak
$exposure
exit
EOF

ftgrouppha infile=simulated.fak outfile=simulated.grp \
    grouptype=optmin groupscale="$min_counts" respfile=response.rmf \
    backfile=none clobber=yes
fthedit 'simulated.grp[1]' BACKFILE add value=NONE
