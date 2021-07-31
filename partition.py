
import os

from polyfun_utils import configure_logger, get_file_name
from compute import PolyFun
from arguments import Parser


def Verify_args(args): 
    #verify that the requested computations are valid
    if args.compute_h2_L2:
        if args.sumstats is None:
            raise ValueError('--sumstats must be specified when using --compute-h2-L2')
        if args.ref_ld_chr is None:
            raise ValueError('--ref-ld-chr must be specified when using --compute-h2-L2')
        if args.w_ld_chr is None:
            raise ValueError('--w-ld-chr must be specified when using --compute-h2-L2')

        #verify partitioning parameters
    if args.skip_Ckmedian and (args.num_bins is None or args.num_bins<=0):
        raise ValueError('You must specify --num-bins when using --skip-Ckmedian')

    return args

    
def Verify_files(args):        
    #check that required input files exist
    if args.compute_h2_L2:
        if not os.path.exists(args.sumstats):
            raise IOError('Cannot find sumstats file %s'%(args.sumstats))
        for chr_num in range(1,23):
            get_file_name(args, 'ref-ld', chr_num, verify_exists=True, allow_multiple=True)
            get_file_name(args, 'w-ld', chr_num, verify_exists=True)
            get_file_name(args, 'annot', chr_num, verify_exists=True, allow_multiple=True)
            
    '''if args.compute_ldscores:
        if args.chr is None: chr_range = range(1,23)            
        else: chr_range = range(args.chr, args.chr+1)
        
        for chr_num in chr_range:
            if args.bfile_chr is not None:
                get_file_name(args, 'bim', chr_num, verify_exists=True)
                get_file_name(args, 'fam', chr_num, verify_exists=True)
                get_file_name(args, 'bed', chr_num, verify_exists=True)
            if not args.compute_h2_L2:
                get_file_name(args, 'snpvar_ridge', chr_num, verify_exists=True)
                get_file_name(args, 'bins', chr_num, verify_exists=True)
                
    if args.compute_h2_bins and not args.compute_ldscores:
        for chr_num in range(1,23):
            get_file_name(args, 'w-ld', chr_num, verify_exists=True)
            if not args.compute_h2_L2:
                get_file_name(args, 'bins', chr_num, verify_exists=True) 
        
  
    '''


if __name__ == '__main__':


    #This command will estimate per-SNP heritabilities for SNPs on odd (resp. even) chromosomes 
    #by applying L2-regularized S-LDSC to even (resp. odd) chromosomes, and will then partition the SNPs into bins.
    
    #configure logger
    #extract args
    args = Parser().parse_args()
    #check and fix args
    args = Verify_args(args) 

    Verify_files(args)

    #check that the output directory exists
    if len(os.path.dirname(args.output_prefix))>0 and not os.path.exists(os.path.dirname(args.output_prefix)):
        raise ValueError('output directory %s doesn\'t exist'%(os.path.dirname(args.output_prefix)))    

    #configure logger
    configure_logger(args.output_prefix)

    #create and run PolyFun object
    polyfun_obj = PolyFun()
    polyfun_obj.polyfun_main(args)