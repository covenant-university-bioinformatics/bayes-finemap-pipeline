import numpy as np
import pandas as pd
import os
import sys
import time
from ldsc_polyfun import ldscore, parse
import logging
from tqdm import tqdm
from pandas.api.types import is_numeric_dtype
from polyfun_utils import configure_logger, set_snpid_index, SNP_COLUMNS
from pyarrow import ArrowIOError
from pyarrow.lib import ArrowInvalid

SNP_COLUMNS = ['CHR', 'SNP', 'BP', 'A1', 'A2']

annots=['Coding_UCSC', 'Coding_UCSC.extend.500', 'Promoter_UCSC',
         'Promoter_UCSC.extend.500','Conserved_LindbladToh', 
         'Conserved_LindbladToh.extend.500', 'Enhancer_Hoffman',
         'Enhancer_Hoffman.extend.500', 'Enhancer_Andersson', 
         'Enhancer_Andersson.extend.500'
        ]



#read bim/snp 
for Chr in tqdm(range(1,23)): 
    array_snps = parse.PlinkBIMFile(f'/content/gdrive/MyDrive/Fine-mapping/project1/g1000_afrcm/g1000.afr.{Chr}.bim')
    df_bim = array_snps.df
  #read annotation files
    for annot in annots:
        df_annot = pd.read_csv(f'/content/gdrive/MyDrive/Fine-mapping/project1/annot/{annot}.{Chr}.gz')
        df_bim[annot] = df_annot['ANNOT'].values
  
    df_bim.drop('CM', inplace=True, axis=1)
    df_bim.to_parquet(f'/content/gdrive/MyDrive/Fine-mapping/project1/annotations/annotations.{Chr}.annot.parquet')



