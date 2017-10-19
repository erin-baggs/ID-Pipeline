
#!/bin/bash -e
source t_coffee-9.03.r1318
source blast-2.6.0
source lmod-6.1
ml raxml/8.2.9
# Use name of python mother script as argument 1

python $1 ../../data/NLR_ID_txt/test_set.txt ../../data/parsed_verbose/Athaliana_167_TAIR10_pfamscan-12-13-2016.parsed.verbose ../../data/proteome/Athaliana_167_TAIR10.protein.fa ../../data/GFF/TAIR10_GFF3_genes.gff ./db.fa > test_run1_output
