# iswitch

## Installation

not yet published

to install it locally, clone it from github: 

    git clone strasserle/iswitch
    
change into the outer folder 'iswitch' and install the package:

    cd iswtich
    pip install .
    
## Run example data

Download the following data into the data directory: 
 
- [pheno](https://xenabrowser.net/datapages/?dataset=TCGA_phenotype_denseDataOnlyDownload.tsv&host=https%3A%2F%2Fpancanatlas.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443)
- [Data_tsv and annotation_file](https://xenabrowser.net/datapages/?dataset=tcga_Kallisto_tpm&host=https%3A%2F%2Ftoil.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443)

In python, use the following line to try the example data:
    
    from iswitch import isoform_switch_analyzer as isa
    switches, switch_pair_df = isa.detect_iswitches()
    
The last line is equivalent to:

    switches, switch_pair_df = sa.detect_iswitches(pheno="iswitch\\data\\TCGA_phenotype_denseDataOnlyDownload.tsv.gz",
                     data="iswitch\\data\\tcga_Kallisto_tpm.gz",
                     disease="iswitch\\kidney-clear-cell-carcinoma",
                     anno="iswitch\\data\\gencode.v23.annotation.transcript.probemap")

## Input format

To detect isoform switches in custom data, specify the following params as follows: 

- **pheno**: Tab separated gz file with a column "sample" and the according phenotype in a column called "primary_disease" as well as a column "sample_type_id" where a 1 denotes case samples and a 11 control samples. 
Note that this notation is inferred from TCGA. 
- **data**: Tab separated table where the columns represent the samples (same sample ids as in pheno) and the rows represent the isoforms (same ids as in annotation_file).
- **disease**: The primary disease you are interested in.
- **anno**: Tab separated file with the columns "id" and "gene" which give the gene id for each transcript.

## Visualize output

To visualize all isoform switches detected for a certain gene, try: 
    
    example_plot_switches(gene_id, switch_pair_df, switches, casecolor, controlcolor) # or 
    example_plot_switches_reordered(gene_id, switch_pair_df, switches, casecolor, controlcolor)
    
To do an enrichment analysis, install [the CORUM core complexes file](https://mips.helmholtz-muenchen.de/corum/download/coreComplexes.txt.zip) into the data directory and run: 

    prot_complex_table, prot_complex_dict = enrichment(switch_pairs)
    
To draw a graph of a certain protein complex, use afterwards: 
    
    draw_protein_graph(complex_name, switch_df, prot_complex_table, prot_complex_dict)
