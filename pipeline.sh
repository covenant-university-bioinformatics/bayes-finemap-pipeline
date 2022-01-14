#!/usr/bin/bash
#Functionally informed fine mapping
#Specifically, the fine-mapping script takes two types of inputs:
#1. A summary statistics file with the following columns: SNP, CHR, BP, A1, A2, Z (Z-score),
#and (optionally) SNPVAR (per-SNP heritability, which is proportional to prior causal probability)
#2. LD information, taken from one of three possible data sources: 
#i. A plink file with genotypes from a reference panel ii. A bgen file with genotypes from a reference panel
#iii. A pre-computed LD matrix
#Two ways to estimate prior causal probabilities (SNPVAR):
#1("One"). Computing  SNPVAR via an L2-regularized extension of stratified LD-score regression (S-LDSC)
#This is a relatively simple approach, but the prior causal probabilities may not be robust to modeling misspecification
#2("Two"). SNPVAR non-parametrically.
#This is the most robust approach, but it is computationally intensive  
#***********************************************************************************************************************

#path to the sumstat. given by user
sumstats=$1 #~/pGWAS/data/AFR/boltlmm_sumstats.gz for testing
#Fine-mapping methods: susie, FINEMAP
method=$2 #susie or finemap. default method is susie
#Null, One, Two
#If Null, i.e no Functionally informed fine mapping
PriorType=$3 
# n is the sample size used to generate summary statistics $3
n=$4 #327209  from testing sumstat file

pop=$5 #Population folders: AFR, EUR, AMR, EASIA, SASIA,
Mainfolder=pGWAS # e.g This could be changed
dir=~/pGWAS/data #directory to the input files i.e population files and python files
#input directory contains, e.g summary statistics, annotation, reference genome, LD score
#will combine pop folder with input directory to point to the population specific directory
#will combine input directory with binary_dir to point to a specific Python script
chr=$6 # the specific chromosome we intend to finemap e.g chr1, chr2,....
#the start and end positions of the target locus to finemap on a given chromosome
start=$7 #46000001
end=$8 #4900000
max_num_causal=$( if [ ! $9 = "" ]; then echo $9;  else echo 5;fi) 
#The two methods allow to specify the expected causal variant e.g defualt value is 5
output=$dir/output # Output directory
#mkdir -p $output if not exists
if [ ! -d "$output" ]; then
    mkdir -p $output 
    echo "Creating the $output directory"
fi
#munge summary statistics file in a PolyFun-friendly parquet format.
#excluding SNPs with INFO score<0.6, with MAF<0.001 or in the MHC region
python $dir/binary_dir/munge_polyfun_sumstats.py \
    --sumstats $sumstats \
    --n $n \
    --out $output/sumstats.parquet \
    --min-info 0.6 \
    --min-maf 0.001
#df = pd.read_parquet("$output/sumstats.parquet")
if [ "$PriorType" = 'One' ]; then
    echo "Computing prior causal probabilities via an L2-regularized extension of stratified LD-score regression (S-LDSC)"
    python $dir/binary_dir/polyfun.py \
        --compute-h2-L2 \
        --no-partitions \
        --output-prefix $output/testrun \
        --sumstats $output/sumstats.parquet \
        --ref-ld-chr $dir/$pop/annotations. \
        --w-ld-chr $dir/$pop/weights.\
        --allow-missing
    cat $output/testrun.$chr.snpvar_ridge_constrained.gz | zcat | head
fi

if [ "$PriorType" = 'Two' ]; then
    echo "Computing prior causal probabilities non-parametrically"
    #Step 1. will partition the SNPs into bins
    echo "1. partition the SNPs into bins"
    python $dir/binary_dir/polyfun.py \
        --compute-h2-L2 \
        --output-prefix $output/testrun \
        --sumstats $output/sumstats.parquet \
        --ref-ld-chr $dir/$pop/annotations. \
        --w-ld-chr $dir/$pop/weights.\
        --allow-missing\
        --skip-Ckmedian\
        --num-bins 20\
    
        #Step 2. Compute LD-scores for each SNP bin
        echo "2. Compute LD-scores for each SNP bin"
        echo "Creating BASH commands for multiprocess"
        #remove the commands from path if already exist
        rm -f $dir/$pop/ldscorebin.txt
        for chr in {1..22}; do echo python $dir/binary_dir/polyfun.py --compute-ldscores --output-prefix $output/testrun --bfile-chr $dir/$pop/reference. --chr $chr >> $dir/$pop/ldscorebin.txt; done
        #Call multiproces
        ./multiprocess.sh $dir/$pop/ldscorebin.txt
        #3. Re-estimate per-SNP heritabilities via S-LDSC
        echo "3. Re-estimate per-SNP heritabilities via S-LDSC"
        python $dir/binary_dir/polyfun.py \
            --compute-h2-bins \
            --output-prefix $output/testrun \
            --sumstats $output/sumstats.parquet \
            --w-ld-chr $dir/$pop/weights.
        cat $output/testrun.$chr.snpvar_constrained.gz | zcat | head
fi

echo "FINE-MAPPING STARTS FROM HERE"
echo "Checking readig method..."
if [[ -f $output/testrun.$chr.snpvar_ridge_constrained.gz && $PriorType = "One" ]];then
    sumstatfile=testrun.$chr.snpvar_ridge_constrained.gz
    read=read_csv
elif [[ -f $output/testrun.$chr.snpvar_constrained.gz && $PriorType = "Two" ]];then
    sumstatfile=testrun.$chr.snpvar_constrained.gz
    read=read_csv
elif [[ -f $output/sumstats.parquet && $PriorType = 'Null' ]];then
    sumstatfile=sumstats.parquet
    read=read_parquet
else
  echo "Summary statistics file could not be found in the path"
fi

if [ -f $output/$sumstatfile ];then
    #extracting the summary statistics for the intended fine-mapping CHR and save to the output directory
    sample_size=$(python $dir/binary_dir/popsize.py $read $chr $output/$sumstatfile $output/chr$chr.sumstats.txt.gz)
    echo "updated sample size: $sample_size"
    echo "Fine-mapping with $method, using genotypes from a plink file"
    sumstats=$output/chr$chr.sumstats.txt.gz
    mkdir -p $output/LD_cache
    path_to_FINEMAP=$dir/package/finemap_v1.4_x86_64/finemap_v1.4_x86_64
    python $dir/binary_dir/finemapper.py \
        --geno $dir/$pop/reference.$chr \
        --sumstats $sumstats \
        --n $sample_size \
        --chr $chr \
        --start $start \
        --end $end \
        --method $method \
        $(if [ $method = 'finemap' ]; then echo '--finemap-exe '$path_to_FINEMAP' ';fi)\
        --max-num-causal $max_num_causal \
        --cache-dir $output/LD_cache \
        --out $output/finemap.$start.$end.gz\
        $(if [ $PriorType = 'Null' ]; then echo '--non-funct';fi)\
        --allow-missing
    cat $output/finemap.$start.$end.gz | zcat | head
fi

#To execute the pipeline ...
#./pipeline.sh ~/pGWAS/data/AFR/boltlmm_sumstats.gz susie Null 327209 ~/data AFR 1 46000001 49000001 5
#Teting with small file
#~/pGWAS/./pipeline.sh ~/pGWAS/data/test_data/boltlmm_sumstats.gz finemap Null 327209 test_data 1 46000001 49000001 5













