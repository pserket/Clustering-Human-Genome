import subprocess
import pandas as pd

import shutil
import urllib.request as request
from contextlib import closing
from urllib.error import URLError
import os


import sys
import time
import urllib
import subprocess


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



def run_process(command, print_out=1):
    p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    outputs = []
    if print_out:
        for outline in p.stdout.readlines():
            print(outline.strip())
            outputs.append(outline.strip().decode("utf-8"))
        
    retval = p.wait()
    return outputs

def prune_filter_vcf(fp, out_path, out_prefix):
    """ ex: fp = '/datasets/dsc180a-wi20-public/Genome/vcf/ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz'
    out_path = '../data/vcf/filtered/'
    out_prefix = 'chr22filtered'
    MAF = '0.1'
    """
    from pathlib import Path
    Path(out_path).mkdir(parents=True, exist_ok=True)
    
    filtered_vcf = out_path + 'filt_' +  out_prefix
    
    # plink2 filter based on maf, missing genotype, sample missingness
    run_process('plink2 --vcf ' + fp + ' --geno 0.1 --mind 0.1 --maf 0.1 --make-bed --recode vcf -out ' + filtered_vcf)
    
    # plink2 filters maf and creates prune.in for pca
    cmd = ["plink2 --vcf " + filtered_vcf + '.vcf' + " --indep-pairwise 50 10 0.1 --out " + out_path + out_prefix]
    
    run_process(cmd)
    return [out_path + out_prefix + '.prune.in', filtered_vcf+ '.vcf']
    
    


    
def make_pca(fp, prune_path, pca_out_path, out_prefix):
    """ fp = '/datasets/dsc180a-wi20-public/Genome/vcf/ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz'
    prune_path = '../data/vcf/filtered/chr22filtered.prune.in'
    
    """
    from pathlib import Path
    Path(pca_out_path).mkdir(parents=True, exist_ok=True)
    
#     if os.path.exists(prune_path.replace('~', str(Path.home())).replace('\'', '')): # check if my computer
#         print('hi')
        
    pca_file_name = pca_out_path + out_prefix + '_pca'
    
    cmd = 'plink2 --vcf ' + fp + ' --extract ' + prune_path + ' --make-bed --pca --out ' + pca_file_name
    
    run_process(cmd)
    
    return pca_file_name
    
def plot_from_pca(pca_file_name, population_df):

    import matplotlib.pyplot as plt


    eigvec = pd.read_table(pca_file_name + ".eigenvec", delimiter=' ', header = None)
    eigval = pd.read_table(pca_file_name + ".eigenval", delimiter=' ', header = None)

    eigvec_wPop = pop_df.merge(eigvec, left_on='sample', right_on=0)

    to_plot = eigvec_wPop.copy()
    to_plot[0] = to_plot[0].apply(lambda x: x[:2])
    to_plot[1] = to_plot[1].apply(lambda x: x[:2])

    # ax = to_plot.plot.scatter(x=2, y=3, label='pop',legend=False)

    import colorsys

    N = 26
    # HSV_tuples = [(x*1.0/N, 0.5, 0.5) for x in range(N)]
    # RGB_tuples = map(lambda x: colorsys.hsv_to_rgb(*x), HSV_tuples)

    # colors = list(RGB_tuples)
    from random import randint

    colors=[]
    for i in range(N):
        colors.append('#%06X' % randint(0, 0xFFFFFF))

    colors = dict(zip(list(to_plot['super_pop'].unique()), colors))
    ###
    _, ax = plt.subplots()
    for key,group in to_plot.groupby('super_pop'):
        group.plot.scatter(ax=ax, x=2, y=3, label=key, color = colors[key]);
        ax.set_title('Principal Components 1 and 2')
        ax.set_ylabel('PC 2')
        ax.set_xlabel('PC 1')
        plt.savefig('pca1_2.png')

#     plt.show()
    
    ax.legend()

    _, ax = plt.subplots()
    for key,group in to_plot.groupby('super_pop'):
        group.plot.scatter(ax=ax, x=2, y=4, label=key, color = colors[key]);
        ax.set_title('Principal Components 1 and 3')
        ax.set_ylabel('PC 3')
        ax.set_xlabel('PC 1')
        plt.savefig('pca1_3.png')

    ###
#     plt.show()
    

    ax.legend()
    perc_var = (eigval[:16])/eigval[0].sum()*100

    ax2 = perc_var.plot(kind='bar', legend=False)
    ax2.set_title('Percent Variance Explained by PC')
    
    plt.savefig('pca_var.png')
    
def index_bwa(ref, name_idx):
    ''' name_idx = 'refhg38'
    '''
    run_process('bwa index -p ' +name_idx+ ' -a bwtsw ' + ref)
    print('done indexing reference')

def fa_to_sam(fp1, fp2, ref, refname, outpath, index=0,header=0):
    """ http://bioinformatics-core-shared-training.github.io/cruk-bioinf-sschool/Day1/Sequence%20Alignment_July2015_ShamithSamarajiwa.pdf
        refname: nickname after indexing the ref file
    """
    if index == 1: # set refrence index
        index_bwa(ref, name_idx)
    space = ' '
    
    # 1000g is paired ends, use sampe
    cmd = 'bwa aln -t 4 ' + refname + ' ' + fp1 +'> ' +outpath + fp1 + ".sai"
    sai1 = outpath + fp1 + ".sai"
    
    cmd = 'bwa aln -t 4 ' + refname + ' ' + fp1 +'> '+ outpath + fp2 + ".sai"
    sai2 =outpath + fp2 + ".sai"
    
    cmd = 'bwa sampe '+ refname + space + sai1 + space + sai2 + space +fp1+space+fp2 +' > '+ outpath + fp1 + '.sam'
    
    return 1

def sam_to_bam(sam, ref,header=0):
    """if sam has header use header=1"""
    if header==0:
        cmd= ['samtools view -bT '+ ref+ ' ' +sam+ ' > '+ sam +'.bam'] # when no header
    else:
        cmd= ['samtools view -bS '+ sam + ' > '+ sam +'.bam']
        
    run_process(cmd)
    run_process('samtools sort ' + sam +'.bam ' + sam +'_sorted.bam')
    return 1

def bam_to_vcf(bam, ref):
    """"TO DO"""
    run_process('samtools faidx ' + ref)


