import numpy as np; np.set_printoptions(precision=4, linewidth=200)
import pandas as pd; pd.set_option('display.width', 200)
import os
import time
import scipy.stats as stats
import logging
import gzip
from tqdm import tqdm
import tempfile
import shutil
import glob
import subprocess
from importlib import reload
from polyfun_utils import set_snpid_index, TqdmUpTo
from pyarrow import ArrowIOError
from pyarrow.lib import ArrowInvalid
from ldstore.bcor import bcor
import scipy.sparse as sparse
from pandas_plink import read_plink
from sklearn.impute import SimpleImputer
from polyfun_utils import set_snpid_index, TqdmUpTo, configure_logger, check_package_versions
import urllib.request
from urllib.parse import urlparse


def save_ld_to_npz(ld_arr, df_ld_snps, npz_file):
    
    assert npz_file.endswith('.npz')
    logging.info('Saving LD file %s'%(npz_file))
    t0 = time.time()

    #save meta file
    meta_file = npz_file[:-4] + '.gz'
    df_ld_snps.to_csv(meta_file, sep='\t', index=False)
    
    #save .npz file
    R = np.tril(ld_arr).astype(np.float32)
    np.fill_diagonal(R, np.diag(R)/2.0)    
    R = sparse.coo_matrix(R)
    sparse.save_npz(npz_file, R, compressed=True)
    logging.info('Done in %0.2f seconds'%(time.time() - t0))

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

class Fine_Mapping(object):
    def __init__(self, genotypes_file, sumstats_file, n, chr_num, ldstore_exe, 
                    sample_file=None, incl_samples=None, cache_dir=None, n_threads=None, memory=None):
    
        #check that data is valid
        if genotypes_file is not None:        
            if genotypes_file.endswith('.bgen'):
                if sample_file is None:
                    raise IOError('sample-file must be provided with a bgen file')
    
    
        #read sumstats and filter to target chromosome only
        logging.info('Loading sumstats file...')
        t0 = time.time()
        try:
            df_sumstats = pd.read_parquet(sumstats_file)
        except (ArrowIOError, ArrowInvalid):
            df_sumstats = pd.read_table(sumstats_file, sep='\s+')
        if not np.any(df_sumstats['CHR'] == chr_num):
            raise IOError('sumstats file does not include any SNPs in chromosome %s'%(chr_num))
        if np.any(df_sumstats['CHR'] != chr_num):
            df_sumstats = df_sumstats.query('CHR==%s'%(chr_num)).copy()
        df_sumstats = set_snpid_index(df_sumstats)
        if 'P' not in df_sumstats.columns:
            df_sumstats['P'] = stats.chi2(1).sf(df_sumstats['Z']**2)
        logging.info('Loaded sumstats for %d SNPs in %0.2f seconds'%(df_sumstats.shape[0], time.time()-t0))

        
        #save class members
        self.genotypes_file = genotypes_file
        self.n = n
        self.sample_file = sample_file
        self.df_sumstats = df_sumstats
        self.incl_samples = incl_samples
        self.ldstore_exe = ldstore_exe
        self.cache_dir = cache_dir
        self.n_threads = n_threads
        self.chr = chr_num
        self.memory = memory
        
        
    def sync_ld_sumstats(self, ld_arr, df_ld_snps, allow_missing=False):
        df_ld_snps = set_snpid_index(df_ld_snps)
        
        if ld_arr is None:
            df_ld = pd.DataFrame(np.zeros(len(df_ld_snps.index), dtype=np.int), index=df_ld_snps.index, columns=['dummy'])
        else:
            assert ld_arr.shape[0] == df_ld_snps.shape[0]
            assert ld_arr.shape[0] == ld_arr.shape[1]
            df_ld = pd.DataFrame(ld_arr, index=df_ld_snps.index, columns=df_ld_snps.index)
        
        #make sure that all SNPs in the sumstats file are in the LD file
        if not np.all(self.df_sumstats_locus.index.isin(df_ld.index)):
            if allow_missing:
                num_missing = np.sum(~self.df_sumstats_locus.index.isin(df_ld.index))
                logging.warning('%d variants with sumstats were not found in the LD file and will be omitted (please note that this may lead to false positives if the omitted SNPs are causal!)'%(num_missing))            
                self.df_sumstats_locus = self.df_sumstats_locus.loc[self.df_sumstats_locus.index.isin(df_ld.index)]
                assert np.all(self.df_sumstats_locus.index.isin(df_ld.index))
            else:
                error_msg = ('not all SNPs in the sumstats file were found in the LD matrix!'
                            ' You could drop the missing SNPs with the flag --allow-missing, but please note that'
                            ' these omitted SNPs may be causal, in which case you may get false positive results...')
                raise ValueError(error_msg)
        
        #filter LD to only SNPs found in the sumstats file
        assert not np.any(self.df_sumstats_locus.index.duplicated())
        if df_ld.shape[0] != self.df_sumstats_locus.shape[0] or np.any(df_ld.index != self.df_sumstats_locus.index):
            if ld_arr is None:
                df_ld = df_ld.loc[self.df_sumstats_locus.index]
            else:
                df_ld = df_ld.loc[self.df_sumstats_locus.index, self.df_sumstats_locus.index]
            df_ld_snps = df_ld_snps.loc[df_ld.index]

        #do a final verification that we're synced
        assert np.all(df_ld.index == self.df_sumstats_locus.index)
        assert np.all(df_ld_snps.index == self.df_sumstats_locus.index)
        
        #add leading zero to sumstats CHR column if needed
        if np.any(df_ld_snps['CHR'].astype(str).str.startswith('0')):
            self.df_sumstats_locus = self.df_sumstats_locus.copy()
            self.df_sumstats_locus['CHR'] = self.df_sumstats_locus['CHR'].astype(str)
            is_1digit = self.df_sumstats_locus['CHR'].str.len()==1
            self.df_sumstats_locus.loc[is_1digit, 'CHR'] = '0' + self.df_sumstats_locus.loc[is_1digit, 'CHR']
        
        #update self.df_ld
        self.df_ld = df_ld
        self.df_ld_snps = df_ld_snps
        assert self.df_ld.notnull().all().all()
        


    def find_cached_ld_file(self, locus_start, locus_end, need_bcor=False):
    
        #if there's no cache dir, return None
        if self.cache_dir is None:
            return None
    
        if self.incl_samples is None:
            fname_pattern = '%s.%d'%(os.path.basename(self.genotypes_file), self.chr)
        else:
            fname_pattern = '%s.%s.%d'%(os.path.basename(self.genotypes_file), os.path.basename(self.incl_samples), self.chr)
            
        #search for suitable LD files
        bcor_files = glob.glob(os.path.join(self.cache_dir, fname_pattern+'*.bcor'))
        npz_files = glob.glob(os.path.join(self.cache_dir, fname_pattern+'*.npz'))
        if need_bcor:
            ld_files = bcor_files + npz_files
        else:
            ld_files = npz_files + bcor_files
            
        for ld_file in ld_files:
            if os.stat(ld_file).st_size==0:
                os.remove(ld_file)
                continue
            ld_basename = os.path.basename(ld_file)
            bp1 = int(ld_basename.split('.')[-3])
            bp2 = int(ld_basename.split('.')[-2])
            assert bp1 < bp2
            if (bp1 > locus_start) or (bp2 < locus_end): continue
            
            #get the list of SNPs in the LD file
            if ld_file.endswith('.npz'):
                meta_file = ld_file[:-4] + '.gz'
                if not os.path.exists(meta_file): continue
                df_ld_snps = pd.read_table(meta_file)
            elif ld_file.endswith('.bcor'):
                bcor_obj = bcor(ld_file)
                df_ld_snps = bcor_obj.getMeta()
                del bcor_obj
                df_ld_snps.rename(columns={'rsid':'SNP', 'position':'BP', 'chromosome':'CHR', 'allele1':'A1', 'allele2':'A2'}, inplace=True, errors='raise')
                ###df_ld_snps['CHR'] = df_ld_snps['CHR'].astype(np.int)
                df_ld_snps['BP'] = df_ld_snps['BP'].astype(np.int)
            else:
                raise IOError('unknown file extension')
            df_ld_snps = set_snpid_index(df_ld_snps)
            
            #make sure that the LD file includes data for all the SNPs in the locus
            if not np.all(self.df_sumstats_locus.index.isin(df_ld_snps.index)):
                continue
    
            #if we got here than we found a suitable d file
            logging.info('Found a cached LD file containing all SNPs with sumstats in chromosome %d BP %d-%d: %s'%(self.chr, locus_start, locus_end, ld_file))
            return ld_file
            
            
    def get_ld_output_file_prefix(self, locus_start, locus_end, output_dir=None):
        if self.cache_dir is None:
            if output_dir is None:
                output_dir = tempfile.mkdtemp()    
            output_prefix = os.path.join(output_dir, 'ld')
        else:
            if self.incl_samples is None:
                output_prefix = os.path.join(self.cache_dir, '%s.%d.%d.%d'%(os.path.basename(self.genotypes_file), self.chr, locus_start, locus_end))
            else:
                output_prefix = os.path.join(self.cache_dir, '%s.%s.%d.%d.%d'%(os.path.basename(self.genotypes_file), os.path.basename(self.incl_samples), self.chr, locus_start, locus_end))
                
        return output_prefix
        
            
    def compute_ld_bgen(self, locus_start, locus_end, verbose=False):
    
        #create df_z
        df_z = self.df_sumstats_locus[['SNP', 'CHR', 'BP', 'A1', 'A2']].copy()
        if df_z['CHR'].iloc[0]<10: df_z['CHR'] = '0' + df_z['CHR'].astype(str)
        df_z.rename(columns={'SNP':'rsid', 'CHR':'chromosome', 'BP':'position', 'A1':'allele1', 'A2':'allele2'}, inplace=True)
                
        #Create LDstore input files
        temp_dir = tempfile.mkdtemp()    
        incl_file = os.path.join(temp_dir, 'incl.incl')
        master_file = os.path.join(temp_dir, 'master.master')
        z_file = os.path.join(temp_dir, 'chr%s.%s_%s.z'%(self.chr, locus_start, locus_end))
        dose_file = os.path.join(temp_dir, 'dosages.bdose')
        df_z.to_csv(z_file, sep=' ', index=False)
        
        #find number of samples
        if self.incl_samples is None:
            num_samples = pd.read_table(self.sample_file).shape[0]-1
        else:
            num_samples = pd.read_table(self.incl_samples, header=None).shape[0]
            
        #get output file name
        bcor_file = self.get_ld_output_file_prefix(locus_start, locus_end, temp_dir) + '.bcor'
        
        #Create LDstore master file
        df_master = pd.DataFrame(columns=['z','bgen','bgi','bcor','dose','sample','n_samples'])
        df_master['z'] = [z_file]
        df_master['bgen'] = [self.genotypes_file]
        df_master['bgi'] = [self.genotypes_file+'.bgi']
        df_master['bcor'] = [bcor_file]
        df_master['bdose'] = [dose_file]
        df_master['sample'] = [self.sample_file]
        df_master['n_samples'] = num_samples
        if self.incl_samples is not None:
            df_master['incl'] = self.incl_samples
        df_master.to_csv(master_file, sep=';', header=True, index=False)    
        
        #run LDstore
        ldstore_cmd = [self.ldstore_exe, '--in-files', master_file, '--write-bcor', '--write-bdose', '--bdose-version', '1.0']
        if self.memory is not None:
            ldstore_cmd += ['--memory', str(self.memory)]
        if self.n_threads is not None:
            ldstore_cmd += ['--n-threads', str(self.n_threads)]
        run_executable(ldstore_cmd, 'LDStore', measure_time=True, show_output=verbose, show_command=verbose)
        
        if not os.path.exists(bcor_file):
            raise IOError('Could not find output BCOR file')
            
        return bcor_file
    
    
    def read_plink_genotypes(self, bed):
        X = bed.compute().astype(np.float32)
        if np.any(np.isnan(X)):
            imp = SimpleImputer(missing_values=np.nan, strategy='mean', copy=False)
            imp.fit(X)
            X = imp.transform(X)
        X -= X.mean(axis=0)
        assert not np.any(np.isnan(X))
        is_polymorphic = X.std(axis=0)>0
        X[:, is_polymorphic] /= X[:, is_polymorphic].std(axis=0)
        return X
    
    
    def compute_ld_plink(self, locus_start, locus_end, verbose):
        logging.info('Computing LD from plink fileset %s region %s-%d'%(self.genotypes_file, locus_start, locus_end))
        t0 = time.time()
        
        #read the plink file
        df_bim, df_fam, bed = read_plink(self.genotypes_file)
        df_bim.rename(columns={'snp':'SNP', 'pos':'BP', 'chrom':'CHR', 'a0':'A2', 'a1':'A1'}, inplace=True)
        df_bim['A1'] = df_bim['A1'].astype('str')
        df_bim['A2'] = df_bim['A2'].astype('str')
        df_bim['CHR'] = df_bim['CHR'].astype(np.int)
        del df_bim['i']
        del df_bim['cm']
        bed = bed.T
        
        #zoom in on target locus
        is_snp_in_region = df_bim['BP'].between(locus_start, locus_end)
        df_bim = df_bim.loc[is_snp_in_region]
        df_ld_snps = df_bim
        bed = bed[:, is_snp_in_region.values]
        
        #compute chunk size, using the formula MEM = bed.shape[0] * chunk_size * 4 / 2**30
        if self.memory is None:
            mem_limit = 1
        else:
            mem_limit = self.memory
        chunk_size = np.int((np.float(mem_limit) * 0.8) / bed.shape[0] / 4 * (2**30))
        if chunk_size==0: chunk_size=1
        if chunk_size > bed.shape[1]: chunk_size = bed.shape[1]
        num_chunks = np.int(np.ceil(bed.shape[1] / chunk_size))
        if num_chunks>1:
            assert chunk_size * (num_chunks-2) < bed.shape[1]-1
        if chunk_size * (num_chunks-1) >= bed.shape[1]:
            num_chunks-=1
        
        #compute LD in chunks
        logging.info('Found %d SNPs in target region. Computing LD in %d chunks...'%(bed.shape[1], num_chunks))
        ld_arr = np.empty((bed.shape[1], bed.shape[1]), dtype=np.float32)
        for chunk_i in tqdm(range(num_chunks)):
            chunk_i_start = chunk_i*chunk_size
            chunk_i_end = np.minimum(chunk_i_start+chunk_size, bed.shape[1])
            X_i = self.read_plink_genotypes(bed[:, chunk_i_start:chunk_i_end])
            ld_arr[chunk_i_start:chunk_i_end, chunk_i_start:chunk_i_end] = X_i.T.dot(X_i) / X_i.shape[0]
            for chunk_j in range(chunk_i+1, num_chunks):
                chunk_j_start = chunk_j*chunk_size
                chunk_j_end = np.minimum(chunk_j_start+chunk_size, bed.shape[1])
                X_j = self.read_plink_genotypes(bed[:, chunk_j_start:chunk_j_end])
                ld_arr[chunk_i_start:chunk_i_end, chunk_j_start:chunk_j_end] = X_i.T.dot(X_j) / X_i.shape[0]
                ld_arr[chunk_j_start:chunk_j_end, chunk_i_start:chunk_i_end] = ld_arr[chunk_i_start:chunk_i_end, chunk_j_start:chunk_j_end].T
        ld_arr = np.nan_to_num(ld_arr, copy=False)
        ld_diag = np.diag(ld_arr).copy()
        if np.any(np.isclose(ld_diag, 0.0)):
            ld_diag[np.isclose(ld_diag, 0.0)] = 1.0
            np.fill_diagonal(ld_arr, ld_diag)
                            
        logging.info('Done in %0.2f seconds'%(time.time() - t0))
        return ld_arr, df_ld_snps
        
    
            
    def set_locus(self, locus_start, locus_end):
    
        #update self.df_sumstats_locus
        self.df_sumstats_locus = self.df_sumstats.query('%d <= BP <= %d'%(locus_start, locus_end))
        if self.df_sumstats_locus.shape[0] == 0:
            raise ValueError('No SNPs found in sumstats file in the BP range %d-%d'%(locus_start, locus_end))
            
            
            
    def get_ld_data(self, locus_start, locus_end, need_bcor=False, verbose=False):
    
        ld_arr, df_ld_snps, ld_file = None, None, None
        
        #check if we already have a suitable LD file in the cache dir
        ld_file = self.find_cached_ld_file(locus_start, locus_end, need_bcor=need_bcor)
            
        #compute LD if we couldn't find a suitable LD file
        if ld_file is None:
            if self.genotypes_file.endswith('.bgen'):
                if not os.path.exists(self.genotypes_file):
                    raise IOError('%s doesn\'t exist'%(self.genotypes_file))
                ld_file = self.compute_ld_bgen(locus_start, locus_end, verbose=verbose)
            elif os.path.exists(self.genotypes_file+'.bed'):
                ld_arr, df_ld_snps = self.compute_ld_plink(locus_start, locus_end, verbose=verbose)
            else:
                raise ValueError('no suitable file found for: %s'%(self.genotypes_file))
                
        #arrange the LD data
        assert ld_file is None or (ld_arr is None and df_ld_snps is None)
        
        #if there is no LD file, return the LD data directly
        if ld_file is None:
            #cache output if possible
            if self.cache_dir is not None:
                npz_file = self.get_ld_output_file_prefix(locus_start, locus_end) + '.npz'
                save_ld_to_npz(ld_arr, df_ld_snps, npz_file)
            return ld_arr, df_ld_snps
            
        #if we have an LD file, return it if it's a bcor and we want a bcor, or return the LD data directly otherwise
        else:
            if ld_file.endswith('.bcor'):
                if need_bcor:
                    return ld_file
                else:
                    ld_arr, df_ld_snps = read_ld_from_file(ld_file)
                    
                    #cache output if possible
                    if self.cache_dir is not None and ld_file.endswith('.bcor'):
                        npz_file = self.get_ld_output_file_prefix(locus_start, locus_end) + '.npz'
                        save_ld_to_npz(ld_arr, df_ld_snps, npz_file)
                    return ld_arr, df_ld_snps
            
            else:
                ld_arr, df_ld_snps = read_ld_from_file(ld_file)
                return ld_arr, df_ld_snps
    


    def finemap(self):
        raise NotImplementedError()
        
    def estimate_h2_hess(self, prop_keep=0.005, R_cutoff=0.99, pvalue_bound=None):
        '''
            prop_keep:  Proprtion of SNPs to use in the estimation (only the ones with the smallest p-values)
            R_cutoff: Exclude one of each pair of SNPs with with magnitude of correlation greater than this value
            pvalue_bound: An upper bound on the p-value cutoff (i.e., SNPs with P greater than this cutoff will never be used in the estimation)
        '''
        
        #keep only potential causal SNPs
        pvalue_cutoff = self.df_sumstats_locus['P'].quantile(prop_keep)
        if pvalue_cutoff==0:
            pvalue_cutoff = np.min(self.df_sumstats_locus['P'].loc[lambda p:p>0])
        if pvalue_bound is not None and pvalue_cutoff>pvalue_bound:
            pvalue_cutoff = pvalue_bound
        is_potential_csnp = self.df_sumstats_locus['P'].values<pvalue_cutoff
        if np.any(is_potential_csnp):
            R_pot_csnp = self.df_ld.loc[is_potential_csnp, is_potential_csnp].values
        else:
            return 0

        #take a maximally independent subset
        np.fill_diagonal(R_pot_csnp,0)
        import networkx as nx
        G = nx.from_numpy_matrix(np.abs(R_pot_csnp)>R_cutoff)
        np.fill_diagonal(R_pot_csnp,1)
        inds = np.sort(nx.maximal_independent_set(G))
        
        #estimate h2 using HESS
        R_subset = R_pot_csnp[np.ix_(inds, inds)]
        alpha_subset = self.df_sumstats_locus.loc[is_potential_csnp, 'Z'].iloc[inds].values / np.sqrt(self.n)
        h2_hess = alpha_subset.dot(np.linalg.solve(R_subset, alpha_subset)) - R_subset.shape[0]/self.n
        
        return h2_hess
        
        
    def estimate_h2_hess_wrapper(self, prop_keep=0.005, R_cutoff=0.99, min_h2=1e-4, num_samples=100):
        '''
            prop_keep:  Proprtion of SNPs to use in the estimation (only the ones with the smallest p-values)
            R_cutoff: Exclude one of each pair of SNPs with with magnitude of correlation greater than this value
            min_h2: Exclude SNPs that tag less than this amount of heritability
            num_samples: Number of random samples of indepdendent SNPs to draw        
        '''

        if min_h2 is None:
            pvalue_bound = None
        else:
            pvalue_bound = stats.chi2(1).sf(min_h2 * self.n)
        
        h2_hess_list = [self.estimate_h2_hess(prop_keep=prop_keep, R_cutoff=R_cutoff, pvalue_bound=pvalue_bound) \
                        for try_num in range(num_samples)]
        h2_hess = np.mean(h2_hess_list)
        return h2_hess        