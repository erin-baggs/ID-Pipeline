#!bin/bash
for f in $1*.out ; do perl /tgac/workarea/group-tg/projects/erinphd/script/plant_rgenes/processing_scripts/K-parse_Pfam_domains_v3.1.pl -p $1 -o $(basename $1 .out).parsed.verbose -v T ; done 
