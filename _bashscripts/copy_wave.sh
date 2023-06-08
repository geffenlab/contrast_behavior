#!/bin/bash

mice=("/cygdrive/e/kilosort_backup/CA070/" \
	  "/cygdrive/e/kilosort_backup/CA073/" \
	  "/cygdrive/e/kilosort_backup/CA074/" \
	  "/cygdrive/d/kilosort_backup/CA102/" \
	  "/cygdrive/d/kilosort_backup/CA104/" \
	  "/cygdrive/d/kilosort_backup/CA105/" \
	  "/cygdrive/d/kilosort_backup/CA106/" \
	  "/cygdrive/d/kilosort_backup/CA107/" \
	  "/cygdrive/f/CA118/" \
	  "/cygdrive/f/CA119/" \
	  "/cygdrive/f/CA121/" \
	  "/cygdrive/f/CA122/")

for f in ${mice[@]}; do
    #--include "*/" --include "*waveforms.mat" --exclude "*" \
    rsync -avhP -e ssh --exclude-from '/Users/chris/OEPs_exclude_list.txt' \
	  OKcomputer@130.91.96.237:$f ~/data/kilosort/${f: -6}
    #echo ${f: -6}

done
