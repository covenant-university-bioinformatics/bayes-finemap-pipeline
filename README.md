# PolyFun, PolyLoc and PolyPred
**PolyFun** (POLYgenic FUNctionally-informed fine-mapping)


This page contains the code of the methods **PolyFun** for functionally-informed fine-mapping, **PolyLoc** for polygenic localization of complex trait heritability, and **PolyPred** for complex trait prediction. **PolyFun** and **PolyLoc** are described in [Weissbrod et al. 2020 Nat Genet](https://www.nature.com/articles/s41588-020-00735-5). 
<br><br>
**PolyFun** estimates prior causal probabilities for SNPs, which can then be used by fine-mapping methods like [SuSiE](https://github.com/stephenslab/susieR) or [FINEMAP](http://www.christianbenner.com/). Unlike previous methods for functionally-informed fine-mapping, **PolyFun** can aggregate polygenic data from across the entire genome and hundreds of functional annotations.
<br><br>




The files in the `ldstore` directory are an adaptation of the [ldstore package](https://pypi.org/project/ldstore/) (written by [Christian Benner](http://www.christianbenner.com)) to Python 3.

<br>
<br>

# Manual
The detailed manual of PolyFun, PolyLoc and PolyPred in the [Wiki page](https://github.com/omerwe/polyfun/wiki). If you run into any issues, please check [the FAQ](https://github.com/omerwe/polyfun/wiki/7.-FAQ) first.



<br><br>
# Installation

## Create an Anaconda environment
The easiest way to install polyfun is by creating a dedicated environment through the [Anaconda Python distribution](https://www.anaconda.com/download). To do this, please install Anaconda on your machine and then type the following commands:
```

```
This will install all the dependencies except for [SuSiE](https://github.com/stephenslab/susieR), [FINEMAP](http://www.christianbenner.com), and [LDstore](http://www.christianbenner.com)
You can use PolyFun without these packages to compute prior causal probabilities, but you won't be able to apply the actual fine-mapping. Please see installation instructions for these three packages below.

After the installation, you can always invoke the PolyFun environment with the command `conda activate polyfun`.
We recommend that you frequently make sure you have the latest version of polyfun installed by going to the polyfun directory and typing `git pull`.


## Install option 2: Manually install packages
PolyFun and PolyLoc are designed for Python >=3.6 and require the following freely available Python packages:
* [numpy](http://www.numpy.org/) and [scipy](http://www.scipy.org/)
* [scikit-learn](http://scikit-learn.org/stable/)
* [pandas](https://pandas.pydata.org/getpandas.html) (version >=0.25.0)
* [tqdm](https://github.com/tqdm/tqdm)
* [pyarrow](https://arrow.apache.org/docs/python/install.html)
* [bitarray](https://github.com/ilanschnell/bitarray)
* [networkx](https://github.com/networkx/networkx) (only required for HESS-based estimation of effect size variance)
* [pandas-plink](https://github.com/limix/pandas-plink)

It is recommended (but not required) to also install the following:
* [rpy2](https://rpy2.bitbucket.io/)  (a Python package)
* [R version 3.5.1 or higher](https://www.r-project.org/)
* [Ckmeans.1d.dp](https://cran.r-project.org/web/packages/Ckmeans.1d.dp/index.html) (a package for R, that will be invoked from python via the rpy2 package).

If rpy2 or Ckmeans.1d.dp are not installed, PolyFun will fallback to suboptimal clustering via scikit-learn.

The `finemapper` script also requires the following:
1. A fine-mapping package you'd like to use. At the moment we support [susieR](https://github.com/stephenslab/susieR) and [FINEMAP v1.4](http://www.christianbenner.com). Please see installation instructions for these packages below.
2. (optional) The program [LDstore 2.0](http://www.christianbenner.com) for computing LD directly from .bgen files (imputed genotypes)

We recommend running PolyFun/PolyLoc via the [Anaconda Python distribution](https://www.anaconda.com/download/). In Anaconda, you can install all the Python packages with the command "conda install \<package_name\>". Alternatively, the Python packages can be installed with the command "pip install --user \<package_name\>".


## Installing SuSiE
To install SuSiE, please start an R shell (usually by typing `R`) and then type: <br>
```
devtools::install_github("stephenslab/susieR@0.8.0",build_vignettes=FALSE)
```
If this doesn't work, please refer to the [SuSiE website](https://github.com/stephenslab/susieR) for more information, or contact the SuSiE authors through the [SuSiE Github page](https://github.com/stephenslab/susieR).

## Installing FINEMAP v1.4
To install FINEMAP v1.4, please type one of the following two commands:
<br>
If you use Linux:
```
wget http://www.christianbenner.com/finemap_v1.4_x86_64.tgz
tar xvf finemap_v1.4_x86_64.tgz
```
If you use Mac OS X :
```
wget http://www.christianbenner.com/finemap_v1.4_MacOSX.tgz
tar xvf finemap_v1.4_MacOSX.tgz
```

## Installing LDstore 2.0
To install LDstore, please type one of the following two commands:
<br>
If you use Linux:
```
wget http://www.christianbenner.com/ldstore_v2.0_x86_64.tgz
tar xzvf ldstore_v2.0_x86_64.tgz
```
If you use Mac OS X :
```
wget http://www.christianbenner.com/ldstore_v2.0_MacOSX.tgz
tar xzvf ldstore_v2.0_MacOSX.tgz
```


<br><br>




