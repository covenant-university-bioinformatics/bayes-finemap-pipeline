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
import subprocess
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

def run_executable(cmd, description, good_returncode=0, measure_time=True, check_errors=True, show_output=False, show_command=False):
    proc = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    logging.info('Running %s...'%(description))
    if show_command:
        logging.info('Command: %s'%(' '.join(cmd)))
    t0 = time.time()
    stdout = []
    if show_output:
        for line in proc.stdout:
            if len(line.strip()) > 0:
                line_str = line.strip().decode("utf-8")
                stdout.append(line_str)
                print(line_str)
        print()
        stdout = '\n'.join(stdout)
        _, stderr = proc.communicate()
    else:
        stdout, stderr = proc.communicate()
        if stdout is not None:
            stdout = stdout.decode('ascii')
            if len(stdout)==0: stdout=None
    if stderr is not None:
        stderr = stderr.decode('ascii')
        if len(stderr)==0: stderr=None        
        
    #if (stderr is not None or proc.returncode != good_returncode):
    if proc.returncode != good_returncode:
        if stderr is not None:            
            logging.error('stderr:\n%s'%(stderr))
        if stdout is not None and not show_output:            
            logging.error('stdout:\n%s'%(stdout))
        raise RuntimeError('%s error'%(description))
    if measure_time:
        logging.info('done in %0.2f seconds'%(time.time() - t0))
        
    if check_errors and stdout is not None:        
        for l in stdout.split('\n'):
            if 'error' in l.lower():
                logging.error(l)
                raise RuntimeError('%s reported an error'%(description))
    if check_errors and stderr is not None:
            for l in stderr.split('\n'):
                if 'error' in l.lower():
                    logging.error(l)
                    raise RuntimeError('%s reported an error'%(description))
        
    return stdout, stderr
    
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

    

class FINEMAP(Fine_Mapping):

    def __init__(self, genotypes_file, sumstats_file, n, chr_num, finemap_exe, ldstore_exe, sample_file=None,
                 incl_samples=None, cache_dir=None, n_threads=None, memory=None):

        super(FINEMAP, self).__init__(genotypes_file, sumstats_file, n, chr_num, ldstore_exe=ldstore_exe, sample_file=sample_file, incl_samples=incl_samples, cache_dir=cache_dir, n_threads=n_threads, memory=memory)
        self.finemap_exe = finemap_exe

    
        
    def finemap(self, locus_start, locus_end, num_causal_snps, use_prior_causal_prob=True, prior_var=None, residual_var=None, hess=False, verbose=False, ld_file=None, debug_dir=None, allow_missing=False, susie_outfile=None):
    
        #check params
        if use_prior_causal_prob and 'SNPVAR' not in self.df_sumstats.columns:
            raise ValueError('SNPVAR column not found in sumstats file')
        if hess:
            raise ValueError('FINEMAP cannot be used with a HESS-based variance estimator')
        if residual_var is not None:
            raise ValueError('cannot specify residual_var for FINEMAP')
        if debug_dir is not None:
            raise NotImplementedError('FINEMAP object does not support --debug-dir')
        # if allow_missing:
            # raise ValueError('FINEMAP object does not support --allow-missing')
            
        #download LD file if it's a url
        if uri_validator(ld_file):
            ld_file = download_ld_file(ld_file)            
            
        #create prefix of output files
        finemap_dir = tempfile.mkdtemp()
        assert os.path.isdir(finemap_dir)        
        finemap_output_prefix = os.path.join(finemap_dir, 'finemap')
                    
    
        #set locus
        self.set_locus(locus_start, locus_end)
        
        #find or create a suitable ld_file
        if num_causal_snps==1:
            if ld_file is not None:
                raise ValueError('cannot specify an ld file when assuming a single causal SNP per locus')
            ld_file = finemap_output_prefix+'.ld'
            np.savetxt(ld_file, np.eye(self.df_sumstats_locus.shape[0], dtype=np.int), fmt='%s')
        else:
            if ld_file is None:
                ld_data = self.get_ld_data(locus_start, locus_end, need_bcor=True, verbose=verbose)
                if isinstance(ld_data, str):
                    ld_file = ld_data
                    assert ld_file.endswith('.bcor')
                    assert os.path.exists(ld_file)
                elif isinstance(ld_data, tuple):
                    assert len(ld_data)==2
                    ld_arr, df_ld_snps = ld_data[0], ld_data[1]
                    self.sync_ld_sumstats(ld_arr, df_ld_snps, allow_missing=allow_missing)
                    del ld_arr, df_ld_snps
                    ld_file = finemap_output_prefix + '.ld'
                    np.savetxt(ld_file, self.df_ld.values, fmt='%0.5f')
            elif ld_file.endswith('.bcor'):
                pass
            elif ld_file.endswith('.npz') or os.path.exists(ld_file+'.npz'):
                    ld_arr, df_ld_snps = read_ld_from_file(ld_file)
                    self.sync_ld_sumstats(ld_arr, df_ld_snps, allow_missing=allow_missing)
                    del ld_arr, df_ld_snps
                    ld_file = finemap_output_prefix + '.ld'
                    np.savetxt(ld_file, self.df_ld.values, fmt='%0.5f')
            else:
                raise ValueError('unknown LD file format for file: %s'%(ld_file))
        

        #define file names
        master_file = finemap_output_prefix+'.master'
        snp_filename = finemap_output_prefix+'.snp'
        config_filename = finemap_output_prefix+'.config'
        cred_filename = finemap_output_prefix+'.cred'
        log_filename = finemap_output_prefix+'.log'        
        z_filename = finemap_output_prefix+'.z'
        
        #flip some of the alleles
        if num_causal_snps == 1:
            is_flipped = np.zeros(self.df_sumstats_locus.shape[0], dtype=np.bool)
        else:
            if ld_file.endswith('.bcor'):
                bcor_obj = bcor(ld_file)
                df_ld_snps = get_bcor_meta(bcor_obj)
                self.sync_ld_sumstats(None, df_ld_snps, allow_missing=allow_missing)
                del df_ld_snps
            assert np.all(self.df_ld_snps['BP'] == self.df_sumstats_locus['BP'])
            is_flipped = self.df_ld_snps['A1'] == self.df_sumstats_locus['A2']
            is_not_flipped = self.df_ld_snps['A1'] == self.df_sumstats_locus['A1']
            assert np.all(is_flipped | is_not_flipped)
            if np.any(is_flipped):
                logging.info('Flipping the effect-sign of %d SNPs that are flipped compared to the LD panel'%(is_flipped.sum()))
            
        #create df_z and save it to disk
        df_z = self.df_sumstats_locus[['SNP', 'CHR', 'BP', 'A1', 'A2', 'Z']].copy()
        df_z.loc[is_flipped, 'A1'] = self.df_sumstats_locus.loc[is_flipped, 'A2']
        df_z.loc[is_flipped, 'A2'] = self.df_sumstats_locus.loc[is_flipped, 'A1']
        df_z.loc[is_flipped, 'Z'] *= (-1)
        
        df_z.rename(columns={'SNP':'rsid', 'CHR':'chromosome', 'BP':'position', 'A1':'allele1', 'A2':'allele2', 'Z':'beta'}, inplace=True, errors='ignore')
        df_z['se'] = 1
        df_z['maf'] = 0.05
        df_z = df_z[['rsid', 'chromosome', 'position', 'allele1', 'allele2', 'maf', 'beta', 'se']]
        if use_prior_causal_prob:
            df_z['prob'] = self.df_sumstats_locus['SNPVAR'] / self.df_sumstats_locus['SNPVAR'].sum()
        df_z.to_csv(z_filename, header=True, index=False, sep=' ', float_format='%0.5f')

        #create the master file
        df_master = pd.DataFrame()
        df_master['z'] = [z_filename]
        df_master['snp'] = [snp_filename]
        df_master['config'] = [config_filename]
        df_master['cred'] = [cred_filename]
        df_master['log'] = [log_filename]
        df_master['n_samples'] = [self.n]
        if ld_file.endswith('.bcor'):
            df_master['bcor'] = [ld_file]
        elif ld_file.endswith('.ld'):
            df_master['ld'] = [ld_file]
        else:
            raise ValueError('Illegal LD file format')
        df_master.to_csv(master_file, sep=';', header=True, index=False)
    
        #prepare the FINEMAP command
        finemap_cmd = [self.finemap_exe]        
        finemap_cmd += ['--in-files', master_file, '--sss']
        finemap_cmd += ['--force-n-samples']
        finemap_cmd += ['--log', log_filename]
        finemap_cmd += ['--n-causal-snps', str(num_causal_snps)]
        finemap_cmd += ['--std-effects']
        ###finemap_cmd += ['--flip-beta']
        if self.n_threads is not None:
            finemap_cmd += ['--n-threads', str(self.n_threads)]
        if prior_var is not None:
            finemap_cmd += ['--prior-std', str(np.sqrt(prior_var))]
        if use_prior_causal_prob:
            finemap_cmd += ['--prior-snps']

        #run FINEMAP
        t0 = time.time()
        m = self.df_sumstats_locus.shape[0]        
        logging.info('Starting %s FINEMAP fine-mapping for chromosome %d BP %d-%d (%d SNPs)'%(
            ('functionally-informed' if use_prior_causal_prob else 'non-functionally informed'),
            self.chr,
            locus_start,
            locus_end,
            self.df_sumstats_locus.shape[0]
            ))
        run_executable(finemap_cmd, 'FINEMAP', measure_time=True, show_output=verbose, show_command=verbose)
        if not os.path.exists(log_filename+'_sss'):
            raise IOError('FINEMAP output files not found')
            
        #load log file
        found_post_csnps = False
        with open(log_filename+'_sss') as f:
            for line in f:
                if line.startswith('- Post-expected # of causal SNPs'):
                    post_mean_num_csnps = float(line[line.index(': ')+2:-1])
                    found_post_csnps = True
                    break
        if not found_post_csnps:
            raise IOError('corrupt log file found: %s'%(log_filename+'_sss'))
                
        #load results
        df_finemap = pd.read_table(snp_filename, sep=' ', usecols=['rsid', 'chromosome', 'position', 'allele1', 'allele2', 'prob', 'mean', 'sd'])
        df_finemap.rename(columns={'rsid':'SNP', 'position':'BP', 'chromosome':'CHR', 'prob':'PIP', 'mean':'BETA_MEAN', 'sd':'BETA_SD', 'allele1':'A1', 'allele2':'A2'}, inplace=True, errors='raise')
        df_finemap.sort_values('PIP', inplace=True, ascending=False)
        
        #read log10bf
        log10bf = None
        with open(log_filename+'_sss', 'r') as f:
            for line in f:
                if line.startswith('- Log10-BF'):
                    log10bf = float(line.split()[-1])
                    break
        if log10bf is None:
            raise ValueError('FINEMP did not report Log10-BF')
            
        #add distance from center
        start = df_finemap['BP'].min()
        end = df_finemap['BP'].max()
        middle = (start+end)//2
        df_finemap['DISTANCE_FROM_CENTER'] = np.abs(df_finemap['BP'] - middle)
        
        #add causal set info
        df_finemap['CREDIBLE_SET'] = 0
        cred_file = None
        for m in range(num_causal_snps, 0, -1):
            if os.path.exists(cred_filename+str(m)):
                cred_file = cred_filename+str(m)
                break
        if cred_file is None:
            raise IOError('cred file not found')
        df_cred = pd.read_table(cred_file, sep=' ', usecols=(lambda c: c.startswith('cred')), comment='#')
        df_finemap.set_index('SNP', inplace=True, drop=False)
        for c_i, c in enumerate(df_cred.columns):
            df_finemap.loc[df_cred[c].dropna(), 'CREDIBLE_SET'] = c_i+1
        df_finemap.reset_index(inplace=True, drop=True)
            
        finemap_time = time.time()-t0
        logging.info('Done in %0.2f seconds'%(finemap_time))
        
        
        return df_finemap

 

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
    if args.method == 'finemap':
        finemap_obj = FINEMAP(genotypes_file=args.geno, sumstats_file=args.sumstats, n=args.n, chr_num=args.chr, 
                                    sample_file=args.sample_file, incl_samples=args.incl_samples,
                                    ldstore_exe=args.ldstore2, finemap_exe=args.finemap_exe, n_threads=args.threads,
                                    cache_dir=args.cache_dir, memory=args.memory)
    else:
        raise ValueError('unknown method specified in --method')

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