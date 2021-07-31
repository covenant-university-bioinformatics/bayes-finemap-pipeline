import numpy as np; np.set_printoptions(precision=4, linewidth=200)
import pandas as pd; pd.set_option('display.width', 200)
import os
import time
import logging
import tempfile
import shutil
import glob
from importlib import reload
from polyfun_utils import TqdmUpTo, configure_logger
import urllib.request
from urllib.parse import urlparse
from ldstore.bcor import bcor
import scipy.sparse as sparse

from fine_map import Fine_Mapping
from arguments import Parser
import sys

def uri_validator(x):
    '''
    code taken from: https://stackoverflow.com/questions/7160737/python-how-to-validate-a-url-in-python-malformed-or-not
    '''
    try:
        result = urlparse(x)
        return all([result.scheme, result.netloc, result.path])
    except:
        return False

def load_ld_bcor(ld_prefix):
    bcor_file = ld_prefix+'.bcor'
    if not os.path.exists(bcor_file):
        raise IOError('%s not found'%(bcor_file))
    logging.info('Loading LD file %s'%(bcor_file))
    t0 = time.time()
    bcor_obj = bcor(bcor_file)
    df_ld_snps = get_bcor_meta(bcor_obj)
    ld_arr = bcor_obj.readCorr([])
    assert np.all(~np.isnan(ld_arr))
    logging.info('Done in %0.2f seconds'%(time.time() - t0))
    return ld_arr, df_ld_snps

def load_ld_npz(ld_prefix):

    logging.info('Loading LD file %s'%(ld_prefix))
    t0 = time.time()
        
    #load SNPs info
    snps_filename_parquet = ld_prefix+'.parquet'
    snps_filename_gz = ld_prefix+'.gz'
    if os.path.exists(snps_filename_parquet):
        df_ld_snps = pd.read_parquet(snps_filename_parquet)
    elif os.path.exists(snps_filename_gz):
        df_ld_snps = pd.read_table(snps_filename_gz, sep='\s+')
        df_ld_snps.rename(columns={'allele1':'A1', 'allele2':'A2', 'position':'BP', 'chromosome':'CHR', 'rsid':'SNP'}, inplace=True, errors='ignore')
    else:
        raise ValueError('couldn\'t find SNPs file %s or %s'%(snps_filename_parquet, snps_filename_gz))
        
    #load LD matrix
    R_filename = ld_prefix+'.npz'
    if not os.path.exists(R_filename):
        raise IOError('%s not found'%(R_filename))
    ld_arr = sparse.load_npz(R_filename).toarray()
    ld_arr = ld_arr+ld_arr.T
    assert np.allclose(np.diag(ld_arr), 1.0)
    assert np.all(~np.isnan(ld_arr))
    
    #sanity checks
    assert ld_arr.shape[0] == ld_arr.shape[1]
    if ld_arr.shape[0] != df_ld_snps.shape[0]:
        raise ValueError('LD matrix has a different number of SNPs than the SNPs file')
    
    logging.info('Done in %0.2f seconds'%(time.time() - t0))
    return ld_arr, df_ld_snps


def get_bcor_meta(bcor_obj):
    df_ld_snps = bcor_obj.getMeta()
    df_ld_snps.rename(columns={'rsid':'SNP', 'position':'BP', 'chromosome':'CHR', 'allele1':'A1', 'allele2':'A2'}, inplace=True, errors='raise')
    ###df_ld_snps['CHR'] = df_ld_snps['CHR'].astype(np.int)
    df_ld_snps['BP'] = df_ld_snps['BP'].astype(np.int)
    return df_ld_snps
        
    
def read_ld_from_file(ld_file):
    #if ld_file is a prefix, make it into a full file name
    if not ld_file.endswith('.bcor') and not ld_file.endswith('.npz'):
        if os.path.exists(ld_file+'.npz'):
            ld_file = ld_file + '.npz'
        elif os.path.exists(ld_file+'.bcor'):
            ld_file = ld_file + '.bcor'
        else:
            raise IOError('No suitable LD file found')

    #read the LD file
    if ld_file.endswith('.bcor'):
        ld_arr, df_ld_snps = load_ld_bcor(ld_file[:-5])
    elif ld_file.endswith('.npz'):
        ld_arr, df_ld_snps = load_ld_npz(ld_file[:-4])
    else:
        raise ValueError('unknown LD format')
    assert np.all(~np.isnan(ld_arr))
    return ld_arr, df_ld_snps


        
def download_ld_file(url_prefix):
    temp_dir = tempfile.mkdtemp()
    filename_prefix = os.path.join(temp_dir, 'ld')
    for suffix in ['npz', 'gz']:
        url = url_prefix + '.' + suffix
        suffix_file = filename_prefix + '.' + suffix
        with TqdmUpTo(unit='B', unit_scale=True, unit_divisor=1024, miniters=1, desc='downloading %s'%(url)) as t:                  
            try:
                urllib.request.urlretrieve(url, filename=suffix_file, reporthook=t.update_to)
            except urllib.error.HTTPError as e:
                if e.code == 404:
                    raise ValueError('URL %s wasn\'t found'%(url))
                else:
                    raise
                
    return filename_prefix


    
def Verify_args(args): 
    #verify that the requested computations are valid
    #check params
    if args.max_num_causal==1:
        if args.geno is not None or args.ld is not None:
            raise ValueError('When max_num_causal=1 fine-mapping please omit the flags --geno and --ld (we cannot use LD information in that setting)')
    else:
        if args.geno is None:
            if args.ld is None:
                raise ValueError('must specify either --geno or --ld')
            if args.ldstore2 is not None:
                raise ValueError('cannot specify both --ld and --ldstore2')
        if args.geno is not None:
            if args.ld is not None:
                raise ValueError('cannot specify both --geno and --ld')
            if args.geno.endswith('.bgen') and args.ldstore2 is None:
                raise ValueError('You must specify --ldstore2 when --geno that points to a bgen file')

    return args

    

class SUSIE(Fine_Mapping):

    def __init__(self, genotypes_file, sumstats_file, n, chr_num, ldstore_exe, sample_file=None,
                 incl_samples=None, cache_dir=None, n_threads=None, memory=None):

        super(SUSIE, self).__init__(genotypes_file, sumstats_file, n, chr_num, ldstore_exe=ldstore_exe, sample_file=sample_file, incl_samples=incl_samples, cache_dir=cache_dir, n_threads=n_threads, memory=memory)
                       
        #load SuSiE R package
        import rpy2
        import rpy2.robjects.numpy2ri as numpy2ri
        import rpy2.robjects as ro
        ro.conversion.py2ri = numpy2ri
        numpy2ri.activate()
        from rpy2.robjects.packages import importr
        self.susieR = importr('susieR')
        self.R_null = ro.rinterface.NULL
        #self.RNULLType = rpy2.rinterface.RNULLType
        
        
    
        
    def finemap(self, locus_start, locus_end, num_causal_snps, use_prior_causal_prob=True, prior_var=None, residual_var=None, hess=False, verbose=False, ld_file=None, debug_dir=None, allow_missing=False, susie_outfile=None):
    
        #check params
        if use_prior_causal_prob and 'SNPVAR' not in self.df_sumstats.columns:
            raise ValueError('SNPVAR column not found in sumstats file')
            
        #set locus
        self.set_locus(locus_start, locus_end)
        
        #download LD file if it's a url
        if uri_validator(ld_file):
            ld_file = download_ld_file(ld_file)
            delete_ld_files_on_exit = True
        else:
            delete_ld_files_on_exit = False
    
        #Load LD data into memory if num_causal_snps>1
        if num_causal_snps==1:
            if hess:
                raise ValueError('Cannot use HESS-based variance estimator when assuming a single causal SNP per locus')
            self.df_ld = pd.DataFrame(np.eye(self.df_sumstats_locus.shape[0]), index=self.df_sumstats_locus.index, columns=self.df_sumstats_locus)
            self.df_ld_snps = self.df_sumstats_locus
        else:
            if ld_file is None:
                ld_arr, df_ld_snps = self.get_ld_data(locus_start, locus_end, need_bcor=False, verbose=verbose)
            else:
                ld_arr, df_ld_snps = read_ld_from_file(ld_file)
            assert np.all(~np.isnan(ld_arr))
            self.sync_ld_sumstats(ld_arr, df_ld_snps, allow_missing=allow_missing)
            del ld_arr
            del df_ld_snps
        
        #define prior causal probabilities
        if use_prior_causal_prob:
            prior_weights = self.df_sumstats_locus['SNPVAR'].copy().values
            prior_weights /= prior_weights.sum()
            assert np.isclose(prior_weights.sum(), 1)
            
        #flip effect sizes if needed
        assert np.all(self.df_ld_snps['BP'] == self.df_sumstats_locus['BP'])
        is_flipped = self.df_ld_snps['A1'] == self.df_sumstats_locus['A2']
        is_not_flipped = self.df_ld_snps['A1'] == self.df_sumstats_locus['A1']
        assert np.all(is_flipped | is_not_flipped)
        bhat = self.df_sumstats_locus['Z'].values.copy()
        if np.any(is_flipped):
            bhat[is_flipped.values] *= -1
            logging.info('Flipping the effect-sign of %d SNPs that are flipped compared to the LD panel'%(is_flipped.sum()))
            
            
        #Use HESS to estimate causal effect sizes
        if hess:
            if prior_var is not None:
                raise ValueError('cannot specify both hess and a custom prior_var')
            prior_var = self.estimate_h2_hess() / num_causal_snps
            if prior_var <= 0:
                raise ValueError('HESS estimates that the locus causally explains zero heritability')
            logging.info('HESS estimated causal effect size variance: %0.4e'%(prior_var))
    
        #rpy2 bug fix
        import rpy2.robjects.numpy2ri as numpy2ri
        reload(numpy2ri)
        numpy2ri.activate()
        
        #run SuSiE
        t0 = time.time()
        m = self.df_sumstats_locus.shape[0]        
        logging.info('Starting %s SuSiE fine-mapping for chromosome %d BP %d-%d (%d SNPs)'%(
            ('functionally-informed' if use_prior_causal_prob else 'non-functionally informed'),
            self.chr,
            locus_start,
            locus_end,
            self.df_ld.shape[0]
            ))
            
        #save variables to debug dir if needed
        if debug_dir is not None:
            os.makedirs(debug_dir, exist_ok=True)
            logging.info('Saving debug info to: %s'%(debug_dir))
            self.df_sumstats_locus.to_csv(os.path.join(debug_dir, 'df_sumstats_locus.txt'), index=False, sep='\t')
            np.savetxt(os.path.join(debug_dir, 'bhat.txt'), bhat)
            #np.savez_compressed(os.path.join(debug_dir, 'R.npz'), R=self.df_ld.values)
            np.savetxt(os.path.join(debug_dir, 'n.txt'), [self.n])
            np.savetxt(os.path.join(debug_dir, 'L.txt'), [num_causal_snps])
            np.savetxt(os.path.join(debug_dir, 'residual_var.txt'), [np.nan] if (residual_var is None) else [residual_var])
            np.savetxt(os.path.join(debug_dir, 'prior_var.txt'), [np.nan] if (prior_var is None) else [prior_var])
            np.savetxt(os.path.join(debug_dir, 'prior_weights.txt'), prior_weights if use_prior_causal_prob else [np.nan])
            
            #create a zipped debug file
            import zipfile
            debug_files = glob.glob(os.path.join(debug_dir, '*.txt'))
            zip_file = os.path.join(debug_dir, 'debug.zip')
            zf = zipfile.ZipFile(zip_file, mode='w')
            for debug_file in debug_files:
                zf.write(debug_file, os.path.basename(debug_file), compress_type=zipfile.ZIP_DEFLATED)
                

        assert self.df_ld.notnull().all().all()
            
        # susie_obj = self.susieR.susie_z(
                # z=self.df_sumstats_locus['Z'].values.reshape((m,1)),
                # R=self.df_ld.values,
                # n=self.n,
                # L=num_causal_snps,
                # prior_variance=(0.0001 if (prior_var is None) else prior_var),
                # estimate_prior_variance=(prior_var is None),
                # residual_variance=(self.R_null if (residual_var is None) else residual_var),
                # estimate_residual_variance=(residual_var is None),
                # verbose=verbose,
                # prior_weights=(prior_weights.reshape((m,1)) if use_prior_causal_prob else self.R_null)
            # )
        try:
            susie_obj = self.susieR.susie_suff_stat(
                    bhat=bhat.reshape((m,1)),
                    shat=np.ones((m,1)),
                    R=self.df_ld.values,
                    n=self.n,
                    L=num_causal_snps,
                    scaled_prior_variance=(0.0001 if (prior_var is None) else prior_var),
                    estimate_prior_variance=(prior_var is None),
                    residual_variance=(self.R_null if (residual_var is None) else residual_var),
                    estimate_residual_variance=(residual_var is None),
                    verbose=verbose,
                    prior_weights=(prior_weights.reshape((m,1)) if use_prior_causal_prob else self.R_null)
                )
        except:
            susie_obj = self.susieR.susie_bhat(
                    bhat=bhat.reshape((m,1)),
                    shat=np.ones((m,1)),
                    R=self.df_ld.values,
                    n=self.n,
                    L=num_causal_snps,
                    scaled_prior_variance=(0.0001 if (prior_var is None) else prior_var),
                    estimate_prior_variance=(prior_var is None),
                    residual_variance=(self.R_null if (residual_var is None) else residual_var),
                    estimate_residual_variance=(residual_var is None),
                    verbose=verbose,
                    prior_weights=(prior_weights.reshape((m,1)) if use_prior_causal_prob else self.R_null)
                )
        susie_time = time.time()-t0        
        logging.info('Done in %0.2f seconds'%(susie_time))
        
        #extract pip and beta_mean
        pip = np.array(self.susieR.susie_get_pip(susie_obj))
        beta_mean = np.array(self.susieR.coef_susie(susie_obj)[1:])
        assert np.allclose(beta_mean, np.sum(np.array(susie_obj.rx2('mu')) * np.array(susie_obj.rx2('alpha')), axis=0) / np.array(susie_obj.rx2('X_column_scale_factors')))

        #compute the posterior mean of beta^2
        s_alpha = np.array(susie_obj.rx2('alpha'))
        s_mu = np.array(susie_obj.rx2('mu'))
        s_mu2 = np.array(susie_obj.rx2('mu2'))
        s_X_column_scale_factors = np.array(susie_obj.rx2('X_column_scale_factors'))
        beta_var = np.sum(s_alpha*s_mu2 - (s_alpha*s_mu)**2, axis=0) / (s_X_column_scale_factors**2)
        assert np.all(beta_var>=0)
        
        #create output df
        df_susie = self.df_sumstats_locus.copy()
        df_susie['PIP'] = pip
        df_susie['BETA_MEAN'] = beta_mean
        df_susie['BETA_SD'] = np.sqrt(beta_var)
        
        #add distance from center
        start = df_susie['BP'].min()
        end = df_susie['BP'].max()
        middle = (start+end)//2
        df_susie['DISTANCE_FROM_CENTER'] = np.abs(df_susie['BP'] - middle)        
        
        #mark causal sets
        self.susie_dict = {key:np.array(susie_obj.rx2(key)) for key in list(susie_obj.names)}
        df_susie['CREDIBLE_SET'] = 0
        susie_sets = self.susie_dict['sets'][0]
        #if type(susie_sets) != self.RNULLType:
        try:
            for set_i, susie_set in enumerate(susie_sets):
                is_in_set = np.zeros(df_susie.shape[0], dtype=np.bool)
                is_in_set[np.array(susie_set)-1] = True
                is_in_set[df_susie['CREDIBLE_SET']>0] = False
                df_susie.loc[is_in_set, 'CREDIBLE_SET'] = set_i+1
        except TypeError:
            pass
            
            
        #save SuSiE object if requested
        if susie_outfile is not None:
            from rpy2.robjects.packages import importr
            R_base = importr('base', robject_translations = {'print.me': 'print_dot_me', 'print_me': 'print_uscore_me'})
            R_base.saveRDS(susie_obj, file=susie_outfile)
            logging.info('Saved SuSiE object to RDS file: %s'%(susie_outfile))

        
        #delete the LD file if needed
        if delete_ld_files_on_exit:
            ld_file_dir = os.path.dirname(ld_file)
            if os.path.exists(ld_file_dir): shutil.rmtree(ld_file_dir)
        
        return df_susie

 

if __name__ == '__main__':


    #This command will estimate per-SNP heritabilities for SNPs on odd (resp. even) chromosomes 
    #by applying L2-regularized S-LDSC to even (resp. odd) chromosomes, and will then partition the SNPs into bins.
    
    #configure logger
    #extract args
    script_name = os.path.basename(__file__)
    
    args = Parser(script_name).parse_args()
    #check and fix args
    args = Verify_args(args) 

    

    #check that the output directory exists
    if len(os.path.dirname(args.out))>0 and not os.path.exists(os.path.dirname(args.out)):
        raise ValueError('output directory %s doesn\'t exist'%(os.path.dirname(args.out)))    
    if args.susie_outfile is not None and len(os.path.dirname(args.susie_outfile))>0 and not os.path.exists(os.path.dirname(args.susie_outfile)):
        raise ValueError('output directory %s doesn\'t exist'%(os.path.dirname(args.susie_outfile)))        

    #configure logger
    configure_logger(args.out)

    #Create a fine-mapping with susie
    if args.method == 'susie':
        finemap_obj = SUSIE(genotypes_file=args.geno, sumstats_file=args.sumstats, n=args.n, chr_num=args.chr, 
                                    sample_file=args.sample_file, incl_samples=args.incl_samples,
                                    ldstore_exe=args.ldstore2, n_threads=args.threads,
                                    cache_dir=args.cache_dir, memory=args.memory)


    #run fine-mapping
    df_finemap = finemap_obj.finemap(locus_start=args.start, locus_end=args.end, num_causal_snps=args.max_num_causal,
                 use_prior_causal_prob=not args.non_funct, prior_var=None, residual_var=None, hess=args.hess,
                 verbose=args.verbose, ld_file=args.ld, debug_dir=args.debug_dir, allow_missing=args.allow_missing, susie_outfile=args.susie_outfile)
    logging.info('Writing fine-mapping results to %s'%(args.out))
    df_finemap.sort_values('PIP', ascending=False, inplace=True)
    if args.out.endswith('.parquet'):
        df_finemap.to_parquet(args.out, index=False)
    else:
        df_finemap.to_csv(args.out, sep='\t', index=False, float_format='%0.5e')