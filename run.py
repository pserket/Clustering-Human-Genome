
import sys
import json

sys.path.append('./src')
from etl import *
    
def main(targets):
    # grab small amount of data
    
    if 'test-project' in targets:
        cfg = json.load(open('./config/test-params.json'))
        
        # Download test data (Chromosome 22) from 1000 genomes site (if up).
        fp = get_biodata(cfg['test_url'], '.' + cfg['outpath'])
        
        pop_file = get_biodata(cfg['pop_url'], cfg['pop_outpath']) 
        
#         fp = './test/chr22_test.vcf.gz'
        
    # Prune/filter test data and save in filtered directory
        filt_out_path = "./data/vcf/filtered/"
        
        fp = fp[0]
        prunefp, filtered_vcf = prune_filter_vcf(fp, filt_out_path, fp.split('/')[-1])
        
        pca_out_path = './data/pca/'
        
        # since it's test sample, filtering will get rid of too much of the variance, so input the regular data to plot 
        # usually input filtered_vcf in place of fp
        pca_file_name = make_pca(fp, prunefp, pca_out_path, fp.split('/')[-1])
        
        pop_df = pd.read_csv(pop_file[0], delimiter='\t').drop(columns=['Unnamed: 4', 'Unnamed: 5'])
        
        plot_from_pca(pca_file_name, pop_df)

    return


if __name__ == '__main__':
    targets = sys.argv[1:]
    main(targets)




