# Set directory
cd data/results

# Load bcftools
module load bcftools

########################################################
# GET RSIDs
########################################################

awk -v OFS="\t" '(NR>1) {print $2,$3}' results_table.txt > get_rsid.tmp
bcftools query -f '%CHROM:%POS\t%ID\t%REF\t%ALT\n' -i 'INFO/VC="SNV"' -R get_rsid.tmp ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/All_20180423.vcf.gz > convert_rsID.txt
rm get_rsid.tmp
