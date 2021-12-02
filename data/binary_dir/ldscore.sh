#!/bin/bash
chr=$1
python compute_ldscores.py  \
--bfile /content/gdrive/MyDrive/Fine-mapping/project1/g1000_afrcm/g1000.afr.${chr} \
--annot /content/gdrive/MyDrive/Fine-mapping/project1/annotations/annotations.${chr}.annot.parquet \
--out /content/gdrive/MyDrive/Fine-mapping/project1/ldscore/annotations.${chr}.ldscore.parquet
