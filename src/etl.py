#!/usr/bin/env python
# coding: utf-8

# In[17]:


import pandas as pd


# In[18]:


# !pip install --user pysam
# !pip install --user scikit-allel
# !pip install --user biopython
# !pip install --user html5lib
# !pip install --user pytabix
# !pip install --user wget
# !pip install --user numcodecs
# !pip install --user zarr
# !pip install --user ipytree
import pysam # module to read BAM alignment files and API to samtools
import allel # module to read VCF files
import Bio # module to read FASTQ files (Use SeqIO)
import wget
import numcodecs
import zarr
import ipytree
import urllib

print('zarr', zarr.__version__, 'numcodecs', numcodecs.__version__)
print(pysam.__version__)
print(allel.__version__)
print(Bio.__version__)


# In[3]:


# TENTATIVE: configs for the future??? population, number samples, filetype


# base_url = 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/'
# outpath = './biodata'


# ## Test VCF

# In[4]:


import os
# os.getcwd()


# In[5]:


import shutil
import urllib.request as request
from contextlib import closing
from urllib.error import URLError

import sys
import time
import urllib
import subprocess

url = [("http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz")]
outpath = "./biodata"

### FROM BLOG: https://blog.shichao.io/2012/10/04/progress_speed_indicator_for_urlretrieve_in_python.html
# adds download progress bar to urlretrieve
def reporthook(count, block_size, total_size):
    global start_time
    if count == 0:
        start_time = time.time()
        return
    duration = time.time() - start_time
    progress_size = int(count * block_size)
    speed = int(progress_size / (1024 * duration))
    percent = int(count * block_size * 100 / total_size)
    sys.stdout.write("\r...%d%%, %d MB, %d KB/s, %d seconds passed" %
                    (percent, progress_size / (1024 * 1024), speed, duration))
    sys.stdout.flush()
###

def get_biodata(urls, outdir):
    """ Feed config file with array of urls, and directory to save at.
    """
    files_downloaded = []
    for url in urls:
        files_downloaded.append(get_biofile(url, outdir))
    return files_downloaded

        
# in the future maybe specify number, chromosome, population, etc.
def get_biofile(url, outdir):
    """ Helper to specifically download one file from url to out directory.
    """
    filename = url.split('/')[-1]
    
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    fullpath = os.path.join(outdir, filename) 

    print("Downloading...")
    urllib.request.urlretrieve(url, fullpath, reporthook) # downloads to fullpath
    
    for_check_zip = '.gz'
    if url.endswith(for_check_zip):
        print()
        print("Unzipping")
        subprocess.run(["unzip",fullpath])
        
    return fullpath        


# In[6]:


# can only do vcf for now
def process_file(fullpath):
    """ Converts VCF file to zarr format and file.
        Fullpath is the path to the file. Outpath is where to put the zarr file.
    """
#     callset = allel.read_vcf(fullpath, fields=['numalt'], log=sys.stdout)
    filename = fullpath.split('/')[-1]

     # choose a name for the zarr storage file
    zarrpath = "".join(filename.split('.')[:-2])
    print(zarrpath)

    # CONVERTING TO ZARR STORAGE
    allel.vcf_to_zarr(fullpath, zarrpath, #group
        fields='*', alt_number=8, log=sys.stdout,
        compressor=numcodecs.Blosc(cname='zstd', clevel=1, shuffle=False))



# In[10]:


filepaths = get_biodata(url, outpath)


# In[12]:


filepaths


# In[14]:


# for file in filepaths:    
#     process_file(file)


# In[16]:


process_file('/datasets/dsc180a-wi20-public/Genome/vcf/ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz')


# In[22]:


# allel.vcf_to_zarr(vcf_path, zarr_path, group='22', fields='*', log=sys.stdout, overwrite=True)
zarr.open_group('chr22TEST.genotypes.zarr', mode='r')

