#!bin/bash
for f in $1*.NLR.txt ; do sort -k1,1 $f > $(basename $f .txt).sort ;done
for f in *.sort ; do python /tgac/workarea/group-tg/projects/erinphd/script/processing_scripts/primary_transcript.py $f > $(basename $f .sort).primaryT ; done
for f in *primaryT ; do cat $f | sed 's/\]//g' | sed 's/\[//g' | sed "s/'//g" > $(basename $f ).rm ; done 
