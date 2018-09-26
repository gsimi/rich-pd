#!/bin/bash

work_dir="$PWD" 

# Helper
usage() { echo "Usage: $0 <data_dir> <ped_tag> <run_tag1> <run_tag2> ..." 1>&2; exit 1; }

# External arguments
data_dir="$1"
ped_tag="$2"
runs_tag=${@:1,2,3}

if [ -z "${data_dir}" ] || [ -z "${runs_tag}" ] || [ -z "${ped_tag}" ]; then usage; fi

# Loop over runs
for run_tag in ${runs_tag[@]}; do
  echo Analyzing run $run_tag with pedestal $ped_tag ...
  root -b -l -e 'gROOT->LoadMacro("'$work_dir'/waveforms.C")' -x ''$work_dir'/macros/analyze.C("'$data_dir'", "'$run_tag'", "'$ped_tag'")' &
done

wait
echo Done!
