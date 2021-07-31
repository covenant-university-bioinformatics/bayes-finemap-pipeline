
import os

from polyfun_utils import configure_logger, get_file_name
from compute import PolyFun
from arguments import Parser


def Verify_args(args): 
    #verify that the requested computations are valid
    if args.compute_h2_bins:
        if args.sumstats is None:
            raise ValueError('--sumstats must be specified when using --compute-h2-bins')
        if args.w_ld_chr is None:
            raise ValueError('--w-ld-chr must be specified when using --compute-h2-bins')
        if args.ref_ld_chr is not None and not args.compute_ldscores:
            raise ValueError('--ref-ld-chr should not be specified when using --compute-h2-bins, unless you also use --compute-ldscores')

    return args

    
def Verify_files(args):        
    #check that required input files exist
                
    if args.compute_h2_bins and not args.compute_ldscores:
        for chr_num in range(1,23):
            get_file_name(args, 'w-ld', chr_num, verify_exists=True)
            if not args.compute_h2_L2:
                get_file_name(args, 'bins', chr_num, verify_exists=True) 


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