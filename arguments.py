import argparse
def Parser(script_name):

    parser = argparse.ArgumentParser()

    if script_name == 'partition.py' or script_name == 'ld_score.py':

        #partitioning-related parameters
        parser.add_argument('--num-bins', type=int, default=None, help='Number of bins to partition SNPs into. If not specified, PolyFun will automatically select this number based on a BIC criterion')
        parser.add_argument('--anno', default=None, help='Comma-delimited list of annotations to use (if not specified, will use all annotations)')
        parser.add_argument('--skip-Ckmedian', default=False, action='store_true', help='If specified, use a regular K-means algorithm instead of the R Ckmeans.1d.dp package')
        
        #mode related parameters
        parser.add_argument('--compute-ldscores', default=False, action='store_true', help='If specified, PolyFun will compute LD-scores of SNP bins')
        parser.add_argument('--compute-h2-L2', default=False, action='store_true', help='If specified, PolyFun will compute per-SNP h2 using L2-regularized S-LDSC')
        parser.add_argument('--compute-h2-bins', default=False, action='store_true', help='If specified, PolyFun will robustly compute per-SNP h2 based on SNP bins')
        parser.add_argument('--no-partitions', default=False, action='store_true', help='If specified, PolyFun will not partition SNPs into bins (do this only if you want to use the per-SNP h^2 from L2-regularized S-LDSC as prior causal probabilities, which is typically not recommended)')
        
        #ld-score related parameters
        parser.add_argument('--chr', type=int, default=None, help='Chromosome number (only applicable when only specifying --ldscores). If not set, PolyFun will compute LD-scores for all chromosomes')
        parser.add_argument('--ld-wind-cm', type=float, default=None, help='window size to be used for estimating LD-scores in units of centiMorgans (cM).')
        parser.add_argument('--ld-wind-kb', type=int, default=None, help='window size to be used for estimating LD-scores in units of Kb.')
        parser.add_argument('--ld-wind-snps', type=int, default=None, help='window size to be used for estimating LD-scores in units of SNPs.')
        parser.add_argument('--chunk-size',  type=int, default=50, help='chunk size for LD-scores calculation')
        parser.add_argument('--keep',  default=None, help='File with ids of individuals to use when computing LD-scores')
        
        #per-SNP h2 related parameters
        parser.add_argument('--q', type=float, default=100, help='The maximum ratio between the largest and smallest truncated per-SNP heritabilites')

        #data input/output parameters
        parser.add_argument('--sumstats', help='Input summary statistics file')
        parser.add_argument('--ref-ld-chr', help='Suffix of LD-score files (as in ldsc)')
        parser.add_argument('--w-ld-chr', help='Suffix of LD-score weights files (as in ldsc)')
        parser.add_argument('--bfile-chr', default=None, help='Prefix of plink files (used to compute LD-scores)')
        parser.add_argument('--ld-ukb', default=False, action='store_true', help='If specified, PolyFun will use UKB LD matrices to compute LD-scores')
        parser.add_argument('--ld-dir', default=None, help='The path of a directory with UKB LD files (if not specified PolyFun will create a temporary directory)')
        parser.add_argument('--output-prefix', required=True, help='Prefix of all PolyFun output file names')    
        parser.add_argument('--allow-missing', default=False, action='store_true', help='If specified, PolyFun will not terminate if some SNPs with sumstats are not found in the annotations files')
        
        #LDSC parameters
        parser.add_argument('--nnls-exact', default=False, action='store_true', help='If specified, S-LDSC will estimate non-negative taus using an exact instead of an approximate solver (this will be slower but slightly more accurate)')
    

    if script_name == 'susie.py' or script_name == 'finemap.py':
        #Fine-Map
        #general parameters
        parser.add_argument('--method', required=True, help='Fine-mapping method (currently susie and finemap are supported)')
        parser.add_argument('--sumstats', required=True, help='Name of sumstats file')
        parser.add_argument('--chr', required=True, type=int, help='Target chromosome')
        parser.add_argument('--start', required=True, type=int, help='First base-pair in the region to finemap')
        parser.add_argument('--end', required=True, type=int, help='Last base-pair in the region to finemap')
        parser.add_argument('--n', required=True, type=int, help='Sample size')
        parser.add_argument('--geno', default=None, help='Genotypes file (plink or bgen format)')
        parser.add_argument('--ld', default=None, help='prefix or fill name of an LD matrix file')
        
        #LDstore related parameters
        parser.add_argument('--ldstore2', default=None, help='Path to an LDstore 2.0 executable file')
        parser.add_argument('--finemap-exe', default=None, help='Path to FINEMAP v1.4 executable file')
        parser.add_argument('--memory', type=int, default=1, help='Maximum amount of memory in GB to allocate to LDStore')
        parser.add_argument('--threads', type=int, default=None, help='The number of CPU cores LDstore will use (if not specified, LDstore will use the max number of CPU cores available')
        parser.add_argument('--cache-dir', default=None, help='If specified, this is a path of a directory that will cache LD matrices that have already been computed')
        parser.add_argument('--debug-dir', default=None, help='If specified, this is a path of a directory that will include files for debugging problems')    
        parser.add_argument('--susie-outfile', default=None, help='If specified, the SuSiE object will be saved to an output file')
        
        
        
        parser.add_argument('--max-num-causal', required=True, type=int, help='Number of causal SNPs')
        parser.add_argument('--non-funct', action='store_true', default=False, help='Perform non-functionally informed fine-mapping')
        parser.add_argument('--hess', action='store_true', default=False, help='If specified, estimate causal effect variance via HESS')
        parser.add_argument('--verbose', action='store_true', default=False, help='If specified, show verbose output')
        parser.add_argument('--allow-missing', default=False, action='store_true', help='If specified, SNPs with sumstats that are not \
                                found in the LD panel will be omitted. This is not recommended, because the omitted SNPs may be causal,\
                                which could lead to false positive results')
        
        parser.add_argument('--sample-file', default=None, help='BGEN files must be used together with a sample file')
        parser.add_argument('--incl-samples', default=None, help='A single-column text file specifying the ids of individuals to include in fine-mapping')
        parser.add_argument('--out', required=True, help='name of the output file')
    
    

    

    return parser

