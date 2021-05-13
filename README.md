# iswitch

## installation

not yet published

to install it locally, clone it from github: 

    git clone https://github.com/strasserle/iswitch
    
change into the outer folder 'iswitch' and install the package

    cd iswtich
    pip install .
    
in python, use the following line to try the example date 
    
    from iswitch import isoform_switch_analyzer as isa
    isa.detect_iswitches()

## input format

currently the input needs to be exactly like the TCGA example data
sources: 
- [pheno](https://xenabrowser.net/datapages/?dataset=TCGA_phenotype_denseDataOnlyDownload.tsv&host=https%3A%2F%2Fpancanatlas.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443)
- [Data_tsv and annotation_file](https://xenabrowser.net/datapages/?dataset=tcga_Kallisto_tpm&host=https%3A%2F%2Ftoil.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443)
