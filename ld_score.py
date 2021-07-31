
import os
import logging
from polyfun_utils import configure_logger, get_file_name
from compute import PolyFun
from arguments import Parser


def Verify_args(args): 
    #verify LD-score related parameters
    if args.chr is not None:
        if args.compute_h2_L2 or args.compute_h2_bins:
            raise ValueError('--chr can only be specified when using only --compute-ldscores')
    if args.bfile_chr is not None:
        if not args.compute_ldscores:
            raise ValueError('--bfile-chr can only be specified when using --compute-ldscores')
    if args.ld_ukb:
        if not args.compute_ldscores:
            raise ValueError('--ld-ukb can only be specified when using --compute-ldscores')
    if args.no_partitions:
        if not args.compute_h2_L2:
            raise ValueError('cannot specify --no-partitions without specifying --compute-h2-L2')
        if args.compute_ldscores:
            raise ValueError('cannot specify both --no-partitions and --compute-ldscores')    
        if args.compute_h2_bins:
            raise ValueError('cannot specify both --no-partitions and --compute-h2-bins')
            
    if args.compute_ldscores and args.compute_h2_bins and not args.compute_h2_L2:
        raise ValueError('cannot use both --compute-ldscores and --compute_h2_bins without also specifying --compute-h2-L2') 

    if args.ld_dir is not None and not args.ld_ukb:
        raise ValueError('You cannot specify --ld-dir without also specifying --ld-ukb')
    if args.bfile_chr is not None and args.ld_ukb:
        raise ValueError('You can specify only one of --bfile-chr and --ld-ukb')
    if not args.compute_ldscores:
        if not (args.ld_wind_cm is None and args.ld_wind_kb is None and args.ld_wind_snps is None):
            raise ValueError('--ld-wind parameters can only be specified together with --compute-ldscores')
        if args.keep is not None:
            raise ValueError('--keep can only be specified together with --compute-ldscores')
        if args.chr is not None:
            raise ValueError('--chr can only be specified together with --compute-ldscores')

    if args.compute_ldscores:
        if args.bfile_chr is None and not args.ld_ukb:
            raise ValueError('You must specify either --bfile-chr or --ld-ukb when you specify --compute-ldscores')    
        if not args.ld_ukb and (args.ld_wind_cm is None and args.ld_wind_kb is None and args.ld_wind_snps is None):
            args.ld_wind_cm = 1.0
            logging.warning('no ld-wind argument specified.  PolyFun will use --ld-cm 1.0')

    return args

    
def Verify_files(args):        
    #check that required input files exist            
    if args.compute_ldscores:
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

if __name__ == '__main__':

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