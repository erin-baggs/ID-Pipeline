#!/bin/bash -e
#SBATCH -p  ei-long # partition (queue)
#SBATCH -c 8 
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem 32G # memory pool for all cores
#SBATCH -J pfamscan
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=baggse@nbi.ac.uk # send-to address

#Updated HMMER and Pfam 
# Usage: bash slurm_run_pfamscan_updated.sh dir
#Pfam31.0 compatible 

# Dependencies
# 1. HMMER software (http://hmmer.janelia.org/)
# including pfam_scan.pl (part of HMMER) Move in same directory as this script or set path at command string.
# 2. Pfam database (http://pfam.xfam.org/)
# 3. File names should be consistent with Phytozome and include Species_*.protein.fa

#adaptation to HPC
source perl_activeperl-5.18.2.1804
source hmmer-3.1b2 #UPDATE

# specify directory in which you want to search for all *protein.fa files and annotate them with pfamscan
IN_DIR="$1" 

# set location of Pfam db
Pfam=/tgac/workarea/group-tg/projects/erinphd/data/Pfam31.0 #UPDATE

# record current date and time and create a commands file.
TSTAMP=$(date +%Y%m%d-%H%M%S)
CMDFILE=pfamscan_commands.${TSTAMP}.txt

protfiles=$(find $IN_DIR -name '*protein.fa') #this line searches for protein fasta files with 'find' command
for FILE in $protfiles
do
    echo "FILE: $FILE"
    basename=$(basename $FILE .protein.fa)
    echo "BASENAME: $basename"
    dirname=$(dirname $FILE)
    echo "DIR: $dirname"
    species=$(echo $basename | cut -f1 -d"_")
    echo "SPECIES: $species"
    if [ ! -e "$dirname/pfam/" ]; then
	mkdir $dirname/pfam
        echo "mkdir $dirname/pfam"
    fi
    export cmd="time perl /tgac/workarea/group-tg/projects/erinphd/script/plant_rgenes/processing_scripts/pfam_scan_EMBL_update.pl -e_seq 1 -e_dom 1 -as -outfile $dirname/pfam/${basename}_pfamscan-$(date +%m-%d-%Y).out -cpu 8 -fasta $FILE -dir $Pfam";

    #check -pfamB still correct
    srun bash -c "$cmd";
    sacct -j ${SLURM_JOB_ID};

    CMDSTR+="'$cmd"
    CMDSTR+="'"$'\n'
done
echo "$CMDSTR" > $CMDFILE
