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

#path to the sumstat. Supplied by user
sumstats=$1 #~data/EUR/boltlmm_sumstats.gz for testing
#Fine-mapping methods: susie, FINEMAP
method=$2 #susie i.e default
#One, Two, Notspecified 
#If Notspecified, i.e no Functionally informed fine mapping
PriorType=$3 
# n is the sample size used to generate summary statistics $3
n=$4 #327209  from testing sumstat file
#Populations: AFR,EUR,AMR,ASIA,
pop=$5 #AFR
dir=~/data 
#binary_dir
#input directory e.g summary stat annotation reference genome LD score
#merge pop with input to point to the population specific directory
#merge input with binary_dir to point to a specific Python script
start=$6 #46000001
end=$7 #49000001
# the start and end positions of the target locus to finemap
chr=$8 #1 
max_num_causal=$9 #defualt value is 5

#mkdir -p $output
#d summary statistics file in a PolyFun-friendly parquet format.
#excluding SNPs with INFO score<0.6, with MAF<0.001 or in the MHC region
python $dir/binary_dir/munge_polyfun_sumstats.py \
    --sumstats $sumstats \
    --n $n \
    --out $dir/output/sumstats.parquet \
    --min-info 0.6 \
    --min-maf 0.001
#df = pd.read_parquet("$dir/output/sumstats.parquet")
if [ "$PriorType" = 'One' ]; then
    echo "Computing prior causal probabilities via an L2-regularized extension of stratified LD-score regression (S-LDSC)"
    python $dir/binary_dir/polyfun.py \
        --compute-h2-L2 \
        --no-partitions \
        --output-prefix $dir/output/testrun \
        --sumstats $dir/output/sumstats.parquet \
        --ref-ld-chr $dir/$pop/annotations. \
        --w-ld-chr $dir/$pop/weights.\
        --allow-missing
    cat $dir/output/testrun.$chr.snpvar_ridge_constrained.gz | zcat | head
fi

if [ "$PriorType" = 'Two' ]; then
    echo "Computing prior causal probabilities non-parametrically"
    #1. will then partition the SNPs into bins
    echo "1. partition the SNPs into bins"
    python $dir/binary_dir/polyfun.py \
        --compute-h2-L2 \
        --output-prefix $dir/output/testrun \
        --sumstats $dir/output/sumstats.parquet \
        --ref-ld-chr $dir/$pop/annotations. \
        --w-ld-chr $dir/$pop/weights.\
        --allow-missing\
        --skip-Ckmedian\
        --num-bins 20\
    
        #2. Compute LD-scores for each SNP bin
        echo "2. Compute LD-scores for each SNP bin"
        echo "Creating BASH commands for multiprocess"
        #remove the commands from path if already exist
        rm -f $dir/$pop/ldscorebin.txt
        for chr in {1..22}; do echo python $dir/binary_dir/polyfun.py --compute-ldscores --output-prefix $dir/output/testrun --bfile-chr $dir/$pop/reference. --chr $chr >> $dir/$pop/ldscorebin.txt; done
        #Call multiproces
        ./multiprocess.sh $dir/$pop/ldscorebin.txt
        #3. Re-estimate per-SNP heritabilities via S-LDSC
        echo "3. Re-estimate per-SNP heritabilities via S-LDSC"
        python $dir/binary_dir/polyfun.py \
            --compute-h2-bins \
            --output-prefix $dir/output/testrun \
            --sumstats $dir/output/sumstats.parquet \
            --w-ld-chr $dir/$pop/weights.
        cat $dir/output/testrun.$chr.snpvar_constrained.gz | zcat | head
fi

echo "FINE-MAPPING STARTS FROM HERE"
if [[ -f $dir/output/testrun.$chr.snpvar_ridge_constrained.gz && $PriorType = "One" ]];then
    sumstatfile=testrun.$chr.snpvar_ridge_constrained.gz
    read=read_csv
elif [[ -f $dir/output/testrun.$chr.snpvar_constrained.gz && $PriorType = "Two" ]];then
    sumstatfile=testrun.$chr.snpvar_constrained.gz
    read=read_csv
elif [[ -f $dir/output/sumstats.parquet && $PriorType = 'Notspecified' ]];then
    sumstatfile=sumstats.parquet
    read=read_parquet
else
  echo "Summary statistics file could not be found in the path"
fi

if [ -f $dir/output/$sumstatfile ];then
    #extract the summary statistics for the intended fine-mapping CHR and save to the user's directory
    snippet="import pandas as pd;
    df= pd.$read('$dir/output/$sumstatfile', $(if [ $read = read_csv ];then echo "sep='\t'";fi));
    df=df[df['CHR']==$chr];
    df.to_csv('$dir/output/chr$chr.sumstats.txt.gz', index=False,  sep='\t',  compression='gzip');
    print(df['N'].unique()[0])"
    sample_size=$(python <( echo $snippet))
    
    echo "Fine-mapping with $method, using genotypes from a plink file"
    sumstats=$dir/output/chr$chr.sumstats.txt.gz
    mkdir -p $dir/output/LD_cache
    path_to_FINEMAP=$dir/package/finemap_v1.4_x86_64/finemap_v1.4_x86_64
    python $dir/binary_dir/finemapper.py \
        --geno $dir/$pop/reference.$chr \
        --sumstats $sumstats \
        --n $sample_size \
        --chr $chr \
        --start $start \
        --end $end \
        --method $method \
        $(if [ $method = 'FINEMAP' ]; then echo '--finemap-exe '$path_to_FINEMAP' ';fi)\
        --max-num-causal $max_num_causal \
        --cache-dir $dir/output/LD_cache \
        --out $dir/output/finemap.$start.$end.gz\
        $(if [ $PriorType = 'Notspecified' ]; then echo '--non-funct';fi)\
        --allow-missing
    cat $dir/output/finemap.$start.$end.gz | zcat | head
fi
    













