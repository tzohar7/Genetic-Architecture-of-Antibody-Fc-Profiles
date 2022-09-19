########################################################
# SETUP
########################################################

cd ~/data/software/annovar/

input_file="/net/bmc-pub14/data/dallgglab/users/tzohar/results/GWAS/fast-lmm/inv/results_table_mini.avinput"
output_file="/net/bmc-pub14/data/dallgglab/users/tzohar/results/GWAS/fast-lmm/inv/anno"

########################################################
# RUN ANNOTATION
########################################################

perl table_annovar.pl $input_file humandb/ -buildver hg19 -out $output_file -remove -protocol refGene,tfbsConsSites,wgRna,gwasCatalog,dbnsfp42a -operation g,r,r,r,f -arg '--separate','','','','' -nastring . -csvout -polish
