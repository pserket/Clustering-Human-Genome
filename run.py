#!/usr/bin/env python
# coding: utf-8

# In[4]:


import subprocess
import pandas as pd
import sys

def run_process(command, print_out=1):
    p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    outputs = []
    if print_out:
        for outline in p.stdout.readlines():
            print(outline.strip())
            outputs.append(outline.strip().decode("utf-8"))
        
    retval = p.wait()
    return outputs

def filter_vcf(fp, out_path, out_prefix, MAF='0.1'):
    """ ex: fp = '/datasets/dsc180a-wi20-public/Genome/vcf/ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz'
    out_path = '../data/vcf/filtered/'
    out_prefix = 'chr22filtered'
    MAF = '0.1'
    """
    from pathlib import Path
    Path(out_path).mkdir(parents=True, exist_ok=True)
    
    cmd = ["cd "+ out_path + "; plink2 --vcf " + fp + " --maf " + str(MAF) + " --indep-pairwise 50 10 0.1 --out "+ out_prefix]
    run_process(cmd)
    return out_path + out_prefix + '.prune.in'
    
def make_pca(fp, prune_path, pca_out_path, out_prefix):
    """ fp = '/datasets/dsc180a-wi20-public/Genome/vcf/ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz'
    prune_path = '../data/vcf/filtered/chr22filtered.prune.in'
    
    """
    from pathlib import Path
    Path(pca_out_path).mkdir(parents=True, exist_ok=True)
    
#     if os.path.exists(prune_path.replace('~', str(Path.home())).replace('\'', '')): # check if my computer
#         print('hi')
        
    pca_file_name = pca_out_path + out_prefix + '_pca'
#     print(pca_file_name)

    cmd = 'plink2 --vcf ' + fp + ' --extract ' + prune_path + ' --make-bed --pca --out ' + pca_file_name
    run_process(cmd)
    
    return pca_file_name
    
def plot_from_pca(pca_file_name):
#     to_cd = '/'.join(pca_file_name.split('/')[:-1])
#     def cd(path):
#         os.chdir(os.path.expanduser(path))
#     cd(to_cd)

    eigvec = pd.read_table(pca_file_name + ".eigenvec", delimiter=' ', header = None)
    eigval = pd.read_table(pca_file_name + ".eigenval", delimiter=' ', header = None)

    to_plot = eigvec.copy()
    to_plot[0] = to_plot[0].apply(lambda x: x[:2])
    to_plot[1] = to_plot[1].apply(lambda x: x[:2])

    ax = to_plot.plot.scatter(x=2, y=3, label='0', legend=False)

    ax.set_title('Principal Components 1 and 2')
    ax.set_ylabel('PC 2')
    ax.set_xlabel('PC 1')

    perc_var = (eigval[:16])/eigval[0].sum()*100

    ax2 = perc_var.plot(kind='bar', legend=False)
    ax2.set_title('Percent Variance Explained by PC')
    
    
def main(targets):
    # grab small amount of data
    if 'data-test' in targets:
        fp = './test/chr22_test.vcf.gz'
        
    # make the data target, NOT IMPLEMENTED
#     if 'data' in targets:
#         cfg = load_params(DATA_PARAMS)
#         get_data(**cfg)

    # process and run data analysis
    if 'process' in targets:
        prunefp = filter_vcf(fp, "./data/vcf/filtered/", fp.split('/')[-1])
        pca_out_path = './data/pca/'
        pca_file_name = make_pca(fp, prunefp, pca_out_path, fp.split('/')[-1])
        plot_from_pca(pca_file_name)
        
    return


if __name__ == '__main__':
    targets = sys.argv[1:]
    main(targets)


# In[ ]:




