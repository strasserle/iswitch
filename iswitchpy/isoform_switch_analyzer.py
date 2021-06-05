import csv
import math
import os
from collections import Counter
from zipfile import ZipFile

import gseapy
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sb
from matplotlib import colors
from matplotlib import rcParams
from scipy.stats import t
from statsmodels.stats.multitest import fdrcorrection


# from pyvis.network import Network
# from matplotlib_venn import venn2, venn2_circles, venn2_unweighted


## Isoform Switch Analyzer
#################################################

# Decection of isoform switches


def detect_iswitches(method="isa",
                     pheno="data/TCGA_phenotype_denseDataOnlyDownload.tsv.gz",
                     data="data/tcga_Kallisto_tpm.gz",
                     disease="kidney-clear-cell-carcinoma",
                     anno="data//gencode.v23.annotation.transcript.probemap",
                     quiet=False,
                     *args):
    """ wrapper for isoforms swtich detection and analysis

    Parameters
    -----------
    method: String (default: "isa"
        "isa" or "spada"
    pheno: String
        Tab separated gz file with a column "sample" and the according phenotype in a column called "primary_disease" as well as a column "sample_type_id" where a 1 denotes case samples and a 11 control samples. Note that this notation is inferred from TCGA.
    data: String
        Tab separated table where the columns represent the samples (same sample ids as in pheno) and the rows represent the isoforms (same ids as in annotation_file).
    disease: String
        The primary disease you are interested in.
    anno: String
        Tab separated file with the columns "id" and "gene" which give the gene id for each transcript.
    quiet: boold (default: False)
        whether to suppress status messages
    *args:
        parameters of the according isoform switch detection methods, eg signif_threshold for isa
    """

    if method == "isa":
        return detect_iswitches_isa(pheno, data, disease, anno, quiet)

    elif method == "spada":
        return detect_iswitches_spada(pheno, data, disease, anno, quiet, args)


def detect_iswitches_spada(pheno, data, disease, anno, signif_threshold=0.3, quiet=False):
    """ detect isoform switches with SPADA

    Parameters
    ------------
    pheno: String
        Tab separated gz file with a column "sample" and the according phenotype in a column called "primary_disease" as well as a column "sample_type_id" where a 1 denotes case samples and a 11 control samples. Note that this notation is inferred from TCGA.
    data: String
        Tab separated table where the columns represent the samples (same sample ids as in pheno) and the rows represent the isoforms (same ids as in annotation_file).
    disease: String
        The primary disease you are interested in.
    anno: String
        Tab separated file with the columns "id" and "gene" which give the gene id for each transcript.
    signif_threshold: float
        minimum fraction of samples
    quiet: boold (default: False)
        whether to suppress status messages
    """

    ## get case and control sample ids for the given cancer type

    # read phenotype file
    phenotype = pd.read_csv(pheno, sep="\t", compression="gzip")

    # select samples for given cancer type
    disease = disease.replace("-", " ")
    phenotype = phenotype[phenotype._primary_disease == disease]

    # select ids for tumor samples
    cases = pd.DataFrame({"sampleID": phenotype["sample"][phenotype.sample_type_id == 1], "condition": "case"})

    # select ids for solid tissue samples
    controls = pd.DataFrame({"sampleID": phenotype["sample"][phenotype.sample_type_id == 11], "condition": "control"})

    # create condition matrix for later usage
    conditionMatrix = cases.append(controls)

    ## check which samples are also in tpm file

    # read header from tpm file
    all_samples = pd.read_csv(data, sep="\t", nrows=0)

    # extract ids
    all_samples = all_samples.columns.to_list()[1:]

    # select ids from the given cancer type (case and control)
    needed_samples = [x for x in all_samples if x in conditionMatrix.sampleID.tolist()]
    needed_cases = [x for x in cases.sampleID if x in needed_samples]
    needed_controls = [x for x in controls.sampleID if x in needed_samples]

    # check wheter spada expression input needs to be created from tpm matrix
    dis_ease = disease.replace(" ", "_")
    if not os.path.isfile(dis_ease + 'Case.tsv') or not os.path.isfile(dis_ease + 'Control.tsv'):

        # read needed columns (samples) from tpm file
        abundance = pd.read_csv(data, sep="\t", compression="gzip", usecols=lambda x: x in set(["sample"] + needed_samples),index_col=0, nrows=100)

        # create input files for spada
        case_expr = 2 ** abundance.loc[:, needed_cases] - 0.001
        case_expr.to_csv(dis_ease + 'Case.tsv', sep='\t')
        control_expr = 2 ** abundance.loc[:, needed_controls] - 0.001
        control_expr.to_csv(dis_ease + 'Control.tsv', sep='\t')

    elif not os.path.isfile("switches" + dis_ease + ".csv"):
        # read needed columns (samples) from tpm file
        abundance = pd.read_csv(data, sep="\t", compression="gzip", usecols=lambda x: x in set(["sample"] + needed_samples), index_col=0)

        switches = isoform_report("spada", abundance, anno, conditionMatrix, dis_ease)

    else:
        switches = pd.read_csv("switches" + dis_ease + ".csv")

    if not quiet:
        print("Starting SPADA run ...")

    # run spada (script spada is in cwd, package spada needs to be installed
    os.system('python spada switches --expression-control ' + dis_ease + 'Control.tsv '
                                                                         '--expression-case ' + dis_ease + 'Case.tsv')

    # read spada ouput, reformat for further usage:
    spada_out = pd.read_csv("switches_spada.tsv", sep="\t")

    # calculate whether significant
    spada_out["frac"] = None
    switch_dict = {}
    for index, row in spada_out.iterrows():
        sams = spada_out.at[index, 'Samples'].split(',')
        symb = spada_out.at[index, 'Symbol']
        ca = spada_out.at[index, 'Case_transcript']
        co = spada_out.at[index, 'Control_transcript']
        nbsams = len(sams)
        frac = nbsams / len(needed_cases)
        # fill frac column with correct fraction of supporting samples
        spada_out['frac'][index] = frac
        key = symb + '/' + co + '-' + ca
        switch_dict[key] = frac

    # filter df to only significant switches
    # and write to tsv file again
    significant = spada_out['frac'] >= signif_threshold
    signSwitches = spada_out[significant]
    out = signSwitches.to_csv(dis_ease + '_significantSwitches_spada.tsv', sep='\t', index=False)

    switch_pairs = signSwitches[["GeneId", "Symbol", "Case_transcript", "Control_transcript"]]

    # CREATE GENE REPORT
    geneDf = signSwitches[['GeneId', 'Symbol']]
    tupelList = geneDf.values.tolist()

    geneList = []
    # convert tupel list to glist: ['geneId1,symbol1', 'geneId2,symbol2',...]
    for ge in tupelList:
        g = ge[0]
        s = ge[1]
        ge = g + ',' + s
        geneList.append(ge)

    # print(geneList)

    # create counter to get counts & retrieve desc. gene order
    geneCounter = Counter(geneList)
    orderedList = geneCounter.most_common(len(geneCounter.keys()))

    # for el in orderedList:
    #    print(el)

    with open('geneReport.tsv', 'w+') as geneReport:
        writer = csv.writer(geneReport, delimiter='\t')
        h1 = 'Gene'
        h2 = 'Symbol'
        h3 = 'Number-of-switches'
        writer.writerow([h1] + [h2] + [h3])
        for g in orderedList:
            writer.writerow([g[0].split(',')[0]] + [g[0].split(',')[1]] + [g[1]])

    # case_trs = spada_out[["GeneId", "Case_transcript"]]
    # control_trs = spada_out[["GeneId", "Control_transcript"]]
    # switchlist = case_trs.append(control_trs, ignore_index=True)
    # switches = pd.DataFrame({"gene_id": spada_out["GeneId"].append(spada_out["GeneId"]), "isoform_id": spada_out["Case_transcript"].append(spada_out["Control_transcript"])})

    return switches, switch_pairs



def detect_iswitches_isa(pheno, data, disease, anno, dIFcutoff=0.1, alpha=0.05, foldChangePseudoCount=0.01,
                         quiet=False):
    """ detect isoform switches with a python reimplementation of the IsoformSwitchAnalyzeR

    Parameters
    -------------
    pheno: String
        Tab separated gz file with a column "sample" and the according phenotype in a column called "primary_disease" as well as a column "sample_type_id" where a 1 denotes case samples and a 11 control samples. Note that this notation is inferred from TCGA.
    data. String
        Tab separated table where the columns represent the samples (same sample ids as in pheno) and the rows represent the isoforms (same ids as in annotation_file).
    disease: String
        The primary disease you are interested in.
    anno: String
        Tab separated file with the columns "id" and "gene" which give the gene id for each transcript.
    dIFcutoff: float (default: 0.1)
        cutoff for the dIF
    alpha: float (default: 0.05)
        cutoff for isoform switch q-value
    foldChangePseudoCount: float (default: 0.01)
        fold change pseudo count
    quiet: boold (default: False)
        whether to suppress status messages
    """

    if not quiet:
        print("Starting ISA run ...")

    ### Part 0: read files
    ######################################

    ## get case and control sample ids for the given cancer type

    # read phenotype file
    phenotype = pd.read_csv(pheno, sep="\t", compression="gzip")

    # select samples for given cancer type
    disease = disease.replace("-", " ")
    phenotype = phenotype[phenotype._primary_disease == disease]
    disease = disease.replace(" ", "_")

    # select ids for tumor samples
    cases = pd.DataFrame({"sampleID": phenotype["sample"][phenotype.sample_type_id == 1], "condition": "case"})

    # select ids for solid tissue samples
    controls = pd.DataFrame({"sampleID": phenotype["sample"][phenotype.sample_type_id == 11], "condition": "control"})

    # create condition matrix for later usage
    conditionMatrix = cases.append(controls)

    ## check which samples are also in tpm file

    # read header from tpm file
    all_samples = pd.read_csv(data, sep="\t", nrows=0)

    # extract ids
    all_samples = all_samples.columns.to_list()[1:]

    # select ids from the given cancer type (case and control)
    needed_samples = [x for x in all_samples if x in conditionMatrix.sampleID.tolist()]

    ### Part 1: from abundances/counts to dIF [importRdata / importIsoformExpression (tximport: estimates counts from abundances, incorporating bias correction, inter-library normalization)]
    ##########################################################################################################################################################################################

    ## select all needed data

    # TCGA data: TPM
    # https://xenabrowser.net/datapages/?dataset=tcga_Kallisto_tpm&host=https%3A%2F%2Ftoil.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443

    # read needed columns (samples) from tpm file
    abundance = pd.read_csv(data, sep="\t", compression="gzip", usecols=lambda x: x in set(["sample"] + needed_samples))

    # rename columns
    abundance.columns = ['isoform_id'] + [x for x in abundance.columns[1:]]

    # remove ids from condition matrix if not in tpm matrix
    conditionMatrix = conditionMatrix.merge(pd.DataFrame({"sampleID": abundance.columns}), how="inner", on="sampleID")

    ##  read annotation file to get the gene id for each isoform

    if not quiet:
        print("CALCULATING ISOFORM FRACTIONS")
        print('Step 1 of 6: Obtaining annotation...')

    # read file
    annotation = pd.read_csv(anno, sep="\t")

    # select columns
    annotation = annotation[["id", "gene"]]

    # rename columns
    annotation.columns = ["isoform_id", "gene_id"]

    # add gene id to abundance matrix
    abundance = abundance.merge(annotation, how='left', on="isoform_id")
    # attention ! NANs !!!

    # print information about the given data
    initial_number_of_genes = len(abundance.gene_id.drop_duplicates())
    initial_number_of_isoforms = len(abundance.isoform_id)
    print("Number of genes in input: " + str(initial_number_of_genes))
    print("Number of isoforms in input: " + str(initial_number_of_isoforms))

    ## sum up abundance values from all isoforms belonging to a gene to get the gene expression (method isoformToGeneExp) (needed for calculating IF)

    if not quiet:
        print('Step 2 of 6: Calculating gene expression...')

    geneRepExpression = abundance.groupby(by="gene_id", axis=0).sum()

    geneRepExpression.reset_index(inplace=True)

    ## which conditions need to be compared, how many replicates are present for each condition

    conditions = conditionMatrix.groupby(by="condition", axis=0).count()

    conditions.reset_index(inplace=True)

    conditions.columns = ["condition", "nrReplicates"]

    conditions1 = list(conditions.condition)

    conditions2 = conditions1.copy()

    comparisonsToMake = pd.DataFrame()

    for t1 in conditions1:

        conditions2.remove(t1)

        for t2 in conditions2:
            comparisonsToMake = comparisonsToMake.append(
                pd.DataFrame([[t1, t2]], columns=["condition_1", "condition_2"]), ignore_index=False)

    ## for each gene / isoform in each condition the mean expression and standard error (of mean (measurement), s.e.m.) expression are calculated

    conditionList = [conditionMatrix[conditionMatrix.condition == x] for x in conditions.condition]

    conditionSummary = dict()

    for sampleVec in conditionList:
        isoSummary = pd.DataFrame({"gene_id": abundance.gene_id, "isoform_id": abundance.isoform_id,
                                   "iso_value": abundance[sampleVec.sampleID].mean(axis=1),
                                   "iso_stderr": abundance[sampleVec.sampleID].std(axis=1).divide(
                                       math.sqrt(len(sampleVec.condition)))})

        geneSummary = pd.DataFrame(
            {"gene_id": geneRepExpression.gene_id, "gene_value": geneRepExpression[sampleVec.sampleID].mean(axis=1),
             "gene_stderr": geneRepExpression[sampleVec.sampleID].std(axis=1).divide(
                 math.sqrt(len(sampleVec.condition)))})

        combinedData = pd.merge(isoSummary, geneSummary, how="outer", on="gene_id")

        conditionSummary[sampleVec.condition.unique()[0]] = combinedData

    ## for each pairwise comparison of condition (or as controlled by the comparisonsToMake argument) the mean gene and isoform expression values are then used to calculate log2fc and IF

    if not quiet:
        print('Step 3 of 6: Making comparisons...')

    isoAnnot = pd.DataFrame()

    for comp in comparisonsToMake.iterrows():
        cond1data = conditionSummary[comp[1][0]]

        cond2data = conditionSummary[comp[1][1]]

        verticalCombined = pd.merge(cond1data, cond2data, how="inner", on=["gene_id", "isoform_id"],
                                    suffixes=["_1", "_2"])

        verticalCombined["condition_1"] = comp[1][0]

        verticalCombined["condition_2"] = comp[1][1]

        isoAnnot = isoAnnot.append(verticalCombined)

    ## add comparison data

    if not quiet:
        print('Step 4 of 6: Calculating DIF ...')

    # log2FC

    ps = foldChangePseudoCount

    isoAnnot["iso_log2_fc"] = np.log2((isoAnnot["iso_value_2"] + ps) / (isoAnnot["iso_value_1"] + ps))

    isoAnnot["gene_log2_fc"] = np.log2((isoAnnot["gene_value_2"] + ps) / (isoAnnot["gene_value_1"] + ps))

    # qvalues

    isoAnnot["iso_q_value"] = None

    isoAnnot["gene_q_value"] = None

    # isoform fraction values

    isoAnnot["IF1"] = np.round(isoAnnot["iso_value_1"] / isoAnnot["gene_value_1"], decimals=4)

    isoAnnot["IF2"] = np.round(isoAnnot["iso_value_2"] / isoAnnot["gene_value_2"], decimals=4)

    isoAnnot["dIF"] = isoAnnot["IF2"] - isoAnnot["IF1"]

    # switch values

    isoAnnot["isoform_switch_q_value"] = None

    isoAnnot["gene_switch_q_value"] = None

    ### Part 2: from dIF to p-value to q-value
    #################################################################################################################################

    # Adapted from Isoform Switch Analyzer version 0.99 (function IsoformSwitchTest)
    # https://github.com/kvittingseerup/IsoformSwitchAnalyzeR/blob/a31694e8be21602a559d14f3b392e6235ccb7e45/R/test_isoform_switches.R

    ## filtering for eligible data

    myDiffData = isoAnnot

    if not quiet:
        print('Step 5 of 6: Filtering for eligible data...')

    ## extract genes with multiple isoforms

    # get unique gene - isoform pairs
    geneIsoOverview = myDiffData[["isoform_id", "gene_id"]].drop_duplicates()

    # count isoforms per gene
    iso_counter = Counter(geneIsoOverview["gene_id"])

    # select genes with at least 2 isoforms
    genes_of_interest = [x for x in myDiffData["gene_id"] if iso_counter[x] > 1]

    # reduce data to selected isoforms
    myDiffData = myDiffData.loc[myDiffData["gene_id"].isin(genes_of_interest)]

    ## calculations

    # add replicate number

    myDiffData["nrReplicates_1"] = [conditions["nrReplicates"][conditions["condition"] == x].values[0] for x in
                                    myDiffData["condition_1"]]

    myDiffData["nrReplicates_2"] = [conditions["nrReplicates"][conditions["condition"] == x].values[0] for x in
                                    myDiffData["condition_2"]]

    ## statistical calculations

    # critical value
    myDiffData["gene_cv1"] = myDiffData["gene_stderr_1"] * myDiffData["nrReplicates_1"] ** (1 / 2) / myDiffData[
        "gene_value_1"]

    # confidence level
    ci = 0.95  # hardcoded since this is not subject to change

    # confidence interval
    myDiffData["gene_lower_CI_1"] = myDiffData["gene_value_1"] - myDiffData["gene_stderr_1"] * [
        t.ppf(ci / 2 + 0.5, myDiffData["nrReplicates_1"].values[x] - 1) for x in range(0, len(myDiffData))]

    myDiffData["gene_lower_CI_2"] = myDiffData["gene_value_2"] - myDiffData["gene_stderr_2"] * [
        t.ppf(ci / 2 + 0.5, myDiffData["nrReplicates_2"].values[x] - 1) for x in range(0, len(myDiffData))]

    # filter on CV and CI - this is nesseary for the implemnted version of Fieller's theorem

    myDiffData2 = myDiffData[
        [((myDiffData["gene_lower_CI_1"].values[x] > 0) & (myDiffData["gene_lower_CI_2"].values[x] > 0)) for x in
         range(0, len(myDiffData))]]

    myDiffData2 = myDiffData2[[((myDiffData2["gene_stderr_1"].values[x] * myDiffData2["nrReplicates_1"].values[x] ** (
                1 / 2) < myDiffData2["gene_value_1"].values[x] / 2) & (myDiffData2["gene_stderr_2"].values[x] *
                                                                       myDiffData2["nrReplicates_2"].values[x] ** (
                                                                                   1 / 2) <
                                                                       myDiffData2["gene_value_2"].values[x] / 2)) for x
                               in range(0, len(myDiffData2))]]

    # calculate expression variance

    myDiffData2["gene_var_1"] = (myDiffData2["gene_stderr_1"] * myDiffData2["nrReplicates_1"] ** (1 / 2)) ** 2

    myDiffData2["gene_var_2"] = (myDiffData2["gene_stderr_2"] * myDiffData2["nrReplicates_2"] ** (1 / 2)) ** 2

    myDiffData2["iso_var_1"] = (myDiffData2["iso_stderr_1"] * myDiffData2["nrReplicates_1"] ** (1 / 2)) ** 2

    myDiffData2["iso_var_2"] = (myDiffData2["iso_stderr_2"] * myDiffData2["nrReplicates_2"] ** (1 / 2)) ** 2

    # calculate var of isoform fraction

    myDiffData2["IF_var_1"] = (1 / myDiffData2["gene_value_1"] ** 2) * (
                myDiffData2["iso_var_1"] + (myDiffData2["IF1"] ** 2 * myDiffData2["gene_var_1"]) - (
                    2 * myDiffData2["IF1"] * myDiffData2["iso_var_1"]))

    myDiffData2["IF_var_2"] = (1 / myDiffData2["gene_value_2"] ** 2) * (
                myDiffData2["iso_var_2"] + (myDiffData2["IF2"] ** 2 * myDiffData2["gene_var_2"]) - (
                    2 * myDiffData2["IF2"] * myDiffData2["iso_var_2"]))

    # filter by var of IF

    myDiffData2 = myDiffData2[
        [((myDiffData2["IF_var_1"].values[x] > 0) & (myDiffData2["IF_var_2"].values[x] > 0)) for x in
         range(0, len(myDiffData2))]]

    ## Do statistical test of IF differences

    if not quiet:
        print('Step 6 of 6: Testing isoform usage')

    # split in comparisonss (since different comparisons might have different numbers of replicates)

    comparisons = myDiffData2[["condition_1", "condition_2"]].drop_duplicates()

    # wrange data to get all comparisons to make
    myDiffData2List = {}

    comparisons = myDiffData2[["condition_1", "condition_2"]].drop_duplicates().transpose()

    n = 0
    for i in comparisons:
        myDiffData2List[n] = myDiffData2[(myDiffData2["condition_1"] == comparisons[i]['condition_1']) & (
                    myDiffData2["condition_2"] == comparisons[i]['condition_2'])]
        n = n + 1

    # for each condition, do test

    for x in myDiffData2List:
        # standard errror of dIF
        myDiffData2List[x]["dIF_std_err"] = ((myDiffData2List[x]["IF_var_1"] / myDiffData2List[x]["nrReplicates_1"]) + (
                    myDiffData2List[x]["IF_var_2"] / myDiffData2List[x]["nrReplicates_1"])) ** (1 / 2)

        # test statistic
        myDiffData2List[x]["t_statistics"] = myDiffData2List[x]["dIF"] / myDiffData2List[x]["dIF_std_err"]

        # degrees of freedom
        myDiffData2List[x]["deg_free"] = ((myDiffData2List[x]["IF_var_1"] / myDiffData2List[x]["nrReplicates_1"]) + (
                    myDiffData2List[x]["IF_var_2"] / myDiffData2List[x]["nrReplicates_2"])) ** 2 / (((myDiffData2List[
                                                                                                          x][
                                                                                                          "IF_var_1"] ** 2) / (
                                                                                                                 myDiffData2List[
                                                                                                                     x][
                                                                                                                     "nrReplicates_1"] ** 2 * (
                                                                                                                             myDiffData2List[
                                                                                                                                 x][
                                                                                                                                 "nrReplicates_1"] - 1))) + (
                                                                                                                (
                                                                                                                            myDiffData2List[
                                                                                                                                x][
                                                                                                                                "IF_var_2"] ** 2) / (
                                                                                                                            myDiffData2List[
                                                                                                                                x][
                                                                                                                                "nrReplicates_2"] ** 2 * (
                                                                                                                                        myDiffData2List[
                                                                                                                                            x][
                                                                                                                                            "nrReplicates_2"] - 1))))

        # p-value
        myDiffData2List[x]["p_value"] = 2 * (1 - t.cdf(abs(myDiffData2List[x]["t_statistics"]), df=myDiffData2List[x][
            "deg_free"]))  # the 2* to make it two tailed

        # p-value correction
        myDiffData2List[x]["isoform_switch_q_value"] = fdrcorrection(myDiffData2List[x]["p_value"], method="indep")[1]

    ## write output files

    for i in range(len(myDiffData2List)):
        myDiffData2List[i].to_csv(disease + "_significantSwitches" + str(i) + "_isa.csv", index=False)

    ### Part 3: Prepare the data for further usage
    ############################################

    # read results from steps before into one df

    switches = pd.DataFrame()
    for i in range(len(myDiffData2List)):
        other = pd.read_csv(disease + "_significantSwitches" + str(i) + "_isa.csv")
        switches = switches.append(other)

    # significant switches

    signif_switches = switches[(switches.isoform_switch_q_value < alpha) & (abs(switches.dIF) > dIFcutoff)].isoform_id

    # describe the direction of the switch

    switches["switchDirection"] = ["up" if x > 0 else "down" for x in switches.dIF]

    # get pairs of isoforms

    gene_ref = switches[["gene_id", "condition_1", "condition_2"]].drop_duplicates()

    switch_pairs = pd.DataFrame()

    for idx in gene_ref.index:

        # for each set of isoforms from one gene and the same condition-comparison
        g_id = gene_ref.loc[idx].gene_id
        c1 = gene_ref.loc[idx].condition_1
        c2 = gene_ref.loc[idx].condition_2

        # get the statistics
        app = switches[(switches.gene_id == g_id) & (switches.condition_1 == c1) & (switches.condition_2 == c2)]

        # divide into up- and downregulated isoforms
        ups = app[app.switchDirection == "up"].isoform_id
        downs = app[app.switchDirection == "down"].isoform_id

        # get pairs of one up- and one downregulated isoform, where at least one of them is significant
        paired_ids = [(u, d) for u in ups for d in downs if
                      ((u in set(signif_switches.isoform_id)) | (d in set(signif_switches.isoform_id)))]

        # hardcoded conditions: condition_1 = case, condition_2 = control
        for pair in paired_ids:
            # produce list in the format needed for visualization
            iso_pair = pd.DataFrame({"GeneId": [g_id], "Control_transcipt": [pair[0]], "Case_transcript": [pair[1]]})
            switch_pairs = switch_pairs.append(iso_pair)

    # write file
    switch_pairs.to_csv(disease.replace(" ", "-") + "_" + "isoformSwitchAnalyzePy.tsv", sep="\t", index=False)

    # # Part 4: Summarize and visualize output
    summarize_iswitches(switches, switch_pairs, alpha, dIFcutoff)

    return switches, switch_pairs


### Part 4: Summarize and visualize output 
############################################

# colors
#    controlcolor = "darkgrey"
#    casecolor = "dimgrey"
#    spada_light = "#9acd9a"
#    spada_dark = "#438943"
#    isa_light = '#9a9aff'
#    isa_dark = "#4d4dff"

def summarize_iswitches(switches, switch_pairs, alpha, dIFcutoff):
    # all isoforms with a significant change in isoform usage
    signif_switches = switches[(switches.isoform_switch_q_value < alpha) & (abs(switches.dIF) > dIFcutoff)].isoform_id

    ## print some basic summary statistics:
    print("Isoform Switch Summary\n------------------------------\n")
    print("Number of switches: " + str(len(switch_pairs)))
    print("Number of genes: " + str(len(switch_pairs.GeneId.drop_duplicates())) + " (" + str(
        len(switches.gene_id.drop_duplicates())) + ")")
    print("Number of significant isoforms: " + str(len(signif_switches.drop_duplicates())))  # 197044
    tabl = switch_pairs.GeneId.value_counts().value_counts().sort_index()
    print("Number of switches per gene:")
    print(pd.DataFrame({"#switches": tabl.index, "count": tabl.values}))

    ## plot the number of isoforms per gene
    plt.rcParams['figure.figsize'] = [6, 2]
    fig, ax = plt.subplots()
    signif_isoform_count = switches[switches.isoform_id.isin(signif_switches)].gene_id.value_counts().value_counts()
    plt.bar(signif_isoform_count.index, signif_isoform_count.values, color="red", label="significant")  # , alpha = 0.5)
    plt.bar(switches.gene_id.value_counts().value_counts().index, switches.gene_id.value_counts().value_counts().values,
            color="black", alpha=0.5, label="not significant")
    my_xticks = [1, 2, 3, 4] + list(range(10, 61, 10))
    plt.xticks(my_xticks)  # color = "red")
    plt.xlabel("isoform count")
    plt.ylabel("gene count")
    plt.title("Number of Isoforms per Gene")
    plt.legend()
    plt.savefig("ISA_isoforms_per_gene.png")

    ## plot the distibution of the corrected p-value
    fig, ax = plt.subplots()
    plt.plot(range(len(switches)), switches.isoform_switch_q_value.sort_values())
    plt.axhline(alpha, color="red")
    plt.ylabel("p_adj")
    plt.title("Isoform Switch adjusted p-value Distribution")
    plt.xlabel("transcripts")
    plt.xticks([], [])
    plt.savefig("ISA_pvalue_distribution.png")

    ## plot the distribution of the dIF
    fig, ax = plt.subplots()
    plt.plot(range(len(switches)), switches.dIF.sort_values())
    plt.axhline(dIFcutoff, color="red")
    plt.ylabel("dIF")
    plt.title("dIF Distribution")
    plt.xlabel("transcripts")
    plt.xticks([], [])
    plt.savefig("ISA_dIF_distribution.png")


## plot example genes

def example_plot_switches(gene_id, switch_pair_df, switches, casecolor="dimgrey", controlcolor="darkgrey"):
    pairs = switch_pair_df[switch_pair_df.Symbol == gene_id]
    isos = set(pairs.Control_transcript)
    isos.update(set(pairs.Case_transcript))
    example = switches[switches.isoform_id.isin(isos)]
    plt.rcParams['figure.figsize'] = [5, 2]
    fig, ax = plt.subplots()
    x = np.arange(len(example.isoform_id))
    width = 0.4
    ax.bar(x - width / 2, height=example.IF2, label="control", width=width, color=controlcolor, capsize=2,
           yerr=example.IF_var_2 ** (1 / 2))
    ax.bar(x + width / 2, height=example.IF1, label="case", width=width, color=casecolor, capsize=2,
           yerr=example.IF_var_1 ** (1 / 2))
    ax.set_xticks(x)
    # plt.xticks(rotation = 45)
    xticklabels = [x.split(".")[0] for x in example.isoform_id]
    ax.set_xticklabels(xticklabels, ha="right", rotation=45)
    # ax.set_xticklabels([], [])
    plt.xlabel("isoforms")
    plt.ylabel("IF")
    plt.legend()
    plt.title("Gene " + gene_id)
    plt.savefig(f'ISA_example_plot_{gene_id}.png')


def example_plot_switches_reordered(gene_id, switch_pair_df, switches, casecolor="dimgrey", controlcolor="darkgrey"):
    pairs = switch_pair_df[switch_pair_df.Symbol == gene_id]
    isos_case = pairs.Case_transcript
    isos_control = pairs.Control_transcript
    example_case = [switches[switches.isoform_id == x].IF1.values[0] for x in isos_case]
    example_control = [switches[switches.isoform_id == x].IF2.values[0] for x in isos_control]
    case_yerr = [switches[switches.isoform_id == x].IF_var_1.values[0] ** (1 / 2) for x in isos_case]
    control_yerr = [switches[switches.isoform_id == x].IF_var_2.values[0] ** (1 / 2) for x in isos_control]
    # plt.rcParams['figure.figsize'] = [10, 2]
    fig, ax = plt.subplots()
    x = np.arange(len(pairs))
    width = 0.45
    ax.bar(x - width / 2, height=example_control, label="control", width=width, color=controlcolor, yerr=control_yerr,
           capsize=2, )
    ax.bar(x + width / 2, height=example_case, label="case", width=width, color=casecolor, yerr=case_yerr, capsize=2, )
    xtickmove = [-width / 1, -0.5 - width / 2] * int((len(pairs)) + 1)
    xtickpos = [(x + y) / 2 for x, y in zip(xtickmove, np.arange(len(pairs) * 2))]
    ax.set_xticks(xtickpos)
    xticklabels = [x.split(".")[0] for y in zip(isos_control, isos_case) for x in y]
    ax.set_xticklabels(xticklabels, ha="right", rotation=45)
    # ax.set_xticklabels([], [])
    plt.xlabel("switches")
    plt.ylabel("IF")
    dev_line = mlines.Line2D([], [], color='black', marker='_', label="std dev")
    handles, old_labels = ax.get_legend_handles_labels()
    handles.append(dev_line)
    plt.legend(bbox_to_anchor=(1, 1), handles=handles, loc='upper left')
    plt.title("Gene " + gene_id)
    plt.ylim(0, None)
    plt.savefig(f'ISA_example_plot_reordered_{gene_id}.png')


## example call:

# example_genes = switch_pairs.Symbol.value_counts().index[0:10].values
# for example_gene in example_genes:
#     example_plot_switches_reordered(example_gene, switch_pairs)
#     example_plot_switches(example_gene, switch_pairs)

## read corum library

def enrichment(switch_pairs):
    prot_complex_table = pd.read_table("data\coreComplexes.txt.zip", compression="zip")
    prot_complex_dict_entrez = dict(
        zip(prot_complex_table.ComplexName, prot_complex_table["subunits(Entrez IDs)"].str.split(pat=";")))
    # list(prot_complex_dict_entrez.items())[:10]
    prot_complex_dict_pfam = dict(
        zip(prot_complex_table.ComplexName, prot_complex_table["subunits(Entrez IDs)"].str.split(pat=";")))
    # list(prot_complex_dict_pfam.items())[:10]
    prot_complex_dict = dict(
        zip(prot_complex_table.ComplexName, prot_complex_table["subunits(Gene name)"].str.split(pat=";")))
    # list(prot_complex_dict.items())[:10]
    prot_complex_ids_dict = dict(
        zip(prot_complex_table.ComplexID, prot_complex_table["subunits(Gene name)"].str.split(pat=";")))
    # list(prot_complex_ids_dict.items())[:10]

    # run enrichment on CORUM protein complexes
    enrichr_complex = gseapy.enrichr(gene_list=switch_pairs.Symbol.drop_duplicates(), gene_sets=prot_complex_dict,
                                     outdir="Enrichr_ISA_complex")

    # prepare data for plotting
    o_list_isa = enrichr_complex.results.Overlap.str.split("/")
    enrichr_complex.results.Overlap = [int(x[0]) for x in o_list_isa]
    enrichr_complex.results["Protein_ratio"] = [int(x[0]) / int(x[1]) for x in o_list_isa]
    enrichr_complex.results["Adjusted P-value"] = enrichr_complex.results["Adjusted P-value"].round(decimals=4)

    ## plot gsea results, function adapted from gprofiler

    def scale_data_5_75(data):
        mind = np.min(data)
        maxd = np.max(data)

        if maxd == mind:
            maxd = maxd + 1
            mind = mind - 1

        drange = maxd - mind
        return ((((data - mind) / drange * 0.70) + 0.05) * 100)

    def plot_enrich(data, n_terms=10):
        # Test data input
        if not isinstance(data, pd.DataFrame):
            raise ValueError('Please input a Pandas Dataframe output by gprofiler.')

        if not np.all([term in data.columns for term in ['Adjusted P-value', 'Term', 'Overlap']]):
            raise TypeError('The data frame {} does not contain enrichment results from gprofiler.'.format(data))

        data_to_plot = data.iloc[:n_terms, :].copy()
        data_to_plot['go.id'] = data_to_plot.Term

        min_pval = data_to_plot['Adjusted P-value'].min()
        max_pval = data_to_plot['Adjusted P-value'].max()

        # Scale intersection_size to be between 5 and 75 for plotting
        # Note: this is done as calibration was done for values between 5 and 75
        data_to_plot['scaled.overlap'] = scale_data_5_75(data_to_plot['Overlap'])
        # data_to_plot['scaled.overlap'] = data_to_plot['Overlap']

        norm = colors.LogNorm(min_pval, max_pval)
        sm = plt.cm.ScalarMappable(cmap="viridis", norm=norm)
        sm.set_array([])

        rcParams.update({'font.size': 34, 'font.weight': 'normal'})

        sb.set(style="whitegrid")

        path = plt.scatter(x='Protein_ratio', y="Term", c='Adjusted P-value', cmap='viridis',
                           norm=colors.LogNorm(min_pval, max_pval),
                           data=data_to_plot, linewidth=1, edgecolor="grey",
                           s=[(i + 9) ** 1.4 for i in data_to_plot['scaled.overlap']])
        ax = plt.gca()
        ax.invert_yaxis()

        ax.set_ylabel('')
        ax.set_xlabel('protein ratio', fontweight='normal', fontsize=14)
        plt.yticks(size=18)
        ax.xaxis.grid(False)
        ax.yaxis.grid(True)

        # Shrink current axis by 20%
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

        # Get tick marks for this plot
        # Note: 6 ticks maximum
        min_tick = np.floor(np.log10(min_pval)).astype(int)
        max_tick = np.ceil(np.log10(max_pval)).astype(int)
        tick_step = np.ceil((max_tick - min_tick) / 6).astype(int)

        # Ensure no 0 values
        if tick_step == 0:
            tick_step = 1
            min_tick = max_tick - 1

        ticks_vals = [10 ** i for i in range(max_tick, min_tick - 1, -tick_step)]
        ticks_labs = ['$10^{' + str(i) + '}$' for i in range(max_tick, min_tick - 1, -tick_step)]

        # Colorbar
        fig = plt.gcf()
        cbaxes = fig.add_axes([0.8, 0.03, 0.03, 0.35])
        cbaxes.xaxis.grid(False)
        cbaxes.yaxis.grid(False)
        cbaxes.set_yticklabels([])
        cbaxes.set_xticklabels([])
        cbar = ax.figure.colorbar(sm, ticks=ticks_vals, anchor=(0, 0), cax=cbaxes)
        cbar.ax.set_yticklabels(ticks_labs)
        cbar.set_label("$p_{adj}$", fontsize=14, fontweight="normal")  # fontsize=14, , fontweight='bold'

        # Size legend
        min_olap = data_to_plot['Overlap'].min()
        max_olap = data_to_plot['Overlap'].max()
        olap_range = max_olap - min_olap

        # Note: approximate scaled 5, 25, 50, 75 values are calculated
        #  and then rounded to nearest number divisible by 5
        # size_leg_vals = [np.round(i / 5) * 5 for i in
        #                  [min_olap, min_olap + (20 / 70) * olap_range, min_olap + (45 / 70) * olap_range, max_olap]]
        step = int(np.round(olap_range / 4))
        step = int(np.round((olap_range + step - 1) / 3))
        size_leg_vals = [x for x in range(min_olap, max_olap + step, step)]
        size_leg_vals = [1, 15, 30, 45]
        size_leg_scaled_vals = scale_data_5_75(size_leg_vals)
        # step = int(np.round(olap_range/4))
        # step = int(np.round((olap_range+step-1)/3))
        # size_leg_scaled_vals = [x for x in range(min_olap, max_olap + step, step)] *10
        # size_leg_vals = size_leg_scaled_vals

        cbaxes2 = fig.add_axes([0.795, 0.2, 0.03, 0.31])
        cbaxes2.xaxis.grid(False)
        cbaxes2.yaxis.grid(False)
        cbaxes2.set_frame_on(False)
        cbaxes2.set_yticklabels([])
        cbaxes2.set_xticklabels([])
        l1 = plt.scatter([], [], s=(size_leg_scaled_vals[0] + 9) ** 1.4, edgecolors='none', color='black')
        l2 = plt.scatter([], [], s=(size_leg_scaled_vals[1] + 9) ** 1.4, edgecolors='none', color='black')
        l3 = plt.scatter([], [], s=(size_leg_scaled_vals[2] + 9) ** 1.4, edgecolors='none', color='black')
        l4 = plt.scatter([], [], s=(size_leg_scaled_vals[3] + 9) ** 1.4, edgecolors='none', color='black')

        labels = [str(int(i)) for i in size_leg_vals]
        # fontsize=12,
        leg = plt.legend([l1, l2, l3, l4], labels, ncol=1, frameon=False,
                         handlelength=1, loc='center left', borderpad=1.5, labelspacing=0.6,
                         handletextpad=0.8, title='protein overlap', scatterpoints=1, bbox_to_anchor=(-2, 1.5),
                         facecolor='black', fontsize=14, title_fontsize=14)

        # fig.tight_layout()
        # plt.gcf().subplots_adjust(right = -3)
        # plt.gcf().subplots_adjust(left = -3)
        # plt.autoscale()
        print(plt.margins())
        plt.margins(10, 10, tight=None)
        print(plt.xlim())
        plt.ylim(-2, 5)
        plt.xlim(-2, 5)
        plt.margins(15, 15)

        plt.savefig(f'ISA_corum_enrichment_top10.png')

    plt.rcParams['figure.figsize'] = [4, 3]
    plot_enrich(enrichr_complex.results.sort_values(by="Adjusted P-value"))

    return prot_complex_table, prot_complex_dict


def draw_protein_graph(complex_name, switch_df, prot_complex_table, prot_complex_dict):
    ## plot example complexes as networks

    # prepare biogrid file (human)
    biogrid_path = "data/BIOGRID-ORGANISM-4.2.193.tab3.zip"

    # ZipFile(biogrid_path).namelist()
    biogrid_human_path = ZipFile(biogrid_path).extract(member="BIOGRID-ORGANISM-Homo_sapiens-4.2.193.tab3.txt",
                                                       path="data")

    # read file
    biogrid_human = pd.read_csv(biogrid_human_path, sep="\t",
                                usecols=["#BioGRID Interaction ID", "Official Symbol Interactor A",
                                         "Official Symbol Interactor B"])
    biogrid_human.columns = ["interaction_id", "gene_a", "gene_b"]

    # needs:
    #   protein_complex_table,
    #   prot_complex_dict,
    #   biogrid_human
    #   switch_pairs

    complex_id = prot_complex_table.ComplexID[prot_complex_table.ComplexName == complex_name].values[0]

    # all proteins in this complex: 
    proteins_of_complex = prot_complex_dict[complex_name]

    # the interactions between these proteins: 
    interactions_of_complex = biogrid_human[
        (biogrid_human["gene_a"].isin(proteins_of_complex)) & (biogrid_human["gene_b"].isin(proteins_of_complex))]

    # the proteins with switches
    proteins_with_switch = set(proteins_of_complex).intersection(switch_df.Symbol.drop_duplicates())

    # interactions between proteins with switches
    interactions_with_switches = interactions_of_complex[
        interactions_of_complex.gene_a.isin(proteins_with_switch) | interactions_of_complex.gene_b.isin(
            proteins_with_switch)]

    edgelist = list(zip(proteins_of_complex, proteins_of_complex))

    edgelist_complex = list(zip(interactions_of_complex.gene_a, interactions_of_complex.gene_b))

    edgelist_switches = list(zip(interactions_with_switches.gene_a, interactions_with_switches.gene_b))

    labels = dict(zip(proteins_of_complex, proteins_of_complex))

    # draw graph 
    # fig, ax = plt.subplots()
    # plt.rcParams['figure.figsize'] = [4, 3]
    G = nx.Graph(edgelist)
    pos = nx.spring_layout(G, seed=0)
    # pos = nx.random_layout(G, seed = 0)
    node_colors = ['red' if x in edgelist_switches else "#97b3f0" for x in edgelist]
    nx.draw_networkx_nodes(G, pos, nodelist=proteins_of_complex, node_size=1000, node_color="darkgrey")
    nx.draw_networkx_nodes(G, pos, nodelist=proteins_with_switch, node_size=1000, node_color="red")
    nx.draw_networkx_edges(G, pos, edgelist=edgelist_complex, width=0.8, edge_color="grey")
    nx.draw_networkx_edges(G, pos, edgelist=edgelist_switches, width=0.8, edge_color="red")
    nx.draw_networkx_labels(G, pos, labels, font_size=8)
    plt.tight_layout()
    plt.axis("off")
    # fig = plt.figure()
    # fig.canvas.mpl_connect('draw_event', on_draw)
    # plt.subplots_adjust = 0.5
    # plt.savefig("test.pdf")
    plt.margins(0.1, 0.1)
    # plt.rcParams.update({'figure.autolayout': True})
    plt.title(complex_name)  # , size = 14)
    grey_patch = mpatches.Patch(color='grey', label=f'unaffected')
    red_patch = mpatches.Patch(color='red', label='affected')
    interaction_line = mlines.Line2D([], [], color='black', marker='_', label="interaction")
    protein_circle = mlines.Line2D([], [], color='black', marker='o', label="protein", markersize=22)
    plt.legend(handles=[protein_circle, interaction_line, red_patch, grey_patch], loc='lower center',
               bbox_to_anchor=(1.2, 0.3), borderpad=1, labelspacing=0.7)  # below: (0.5, -0.4)
    plt.savefig(f'ISA_example_complex_{complex_name}.png')


## example invocation:

# example_complexes = enrichr_complex.results.sort_values(by = "Adjusted P-value")[["Term"]][0:10]
# plt.rcParams['figure.figsize'] = [4, 3]
# for example_complex in example_complexes:
#     draw_protein_graph("MIB complex", switch_pairs)


# detect_iswitches(method="spada")


def isoform_report(method, abundance, anno, conditionMatrix, disease="kidney_clear_cell_carcinoma", foldChangePseudoCount=0.001, quiet=True):

    if not quiet:
        print("CREATING ISOFORM REPORT")
        print('Step 1 of 6: Obtaining annotation...')

    # read file
    annotation = pd.read_csv(anno, sep="\t")

    # select columns
    annotation = annotation[["id", "gene"]]

    # rename columns
    annotation.columns = ["isoform_id", "gene_id"]

    # add gene id to abundance matrix
    abundance = abundance.merge(annotation, how='left', on="isoform_id")
    # attention ! NANs !!!

    # add gene id to abundance matrix
    abundance = abundance.merge(annotation, how='left', on="isoform_id")
    # attention ! NANs !!!

    # print information about the given data
    initial_number_of_genes = len(abundance.gene_id.drop_duplicates())
    initial_number_of_isoforms = len(abundance.isoform_id)
    print("Number of genes in input: " + str(initial_number_of_genes))
    print("Number of isoforms in input: " + str(initial_number_of_isoforms))

    ## sum up abundance values from all isoforms belonging to a gene to get the gene expression (method isoformToGeneExp) (needed for calculating IF)

    if not quiet:
        print('Step 2 of 6: Calculating gene expression...')

    geneRepExpression = abundance.groupby(by="gene_id", axis=0).sum()

    geneRepExpression.reset_index(inplace=True)

    ## which conditions need to be compared, how many replicates are present for each condition

    conditions = conditionMatrix.groupby(by="condition", axis=0).count()

    conditions.reset_index(inplace=True)

    conditions.columns = ["condition", "nrReplicates"]

    conditions1 = list(conditions.condition)

    conditions2 = conditions1.copy()

    comparisonsToMake = pd.DataFrame()

    for t1 in conditions1:

        conditions2.remove(t1)

        for t2 in conditions2:
            comparisonsToMake = comparisonsToMake.append(
                pd.DataFrame([[t1, t2]], columns=["condition_1", "condition_2"]), ignore_index=False)

    ## for each gene / isoform in each condition the mean expression and standard error (of mean (measurement), s.e.m.) expression are calculated

    conditionList = [conditionMatrix[conditionMatrix.condition == x] for x in conditions.condition]

    conditionSummary = dict()

    for sampleVec in conditionList:
        isoSummary = pd.DataFrame({"gene_id": abundance.gene_id, "isoform_id": abundance.isoform_id,
                                   "iso_value": abundance[sampleVec.sampleID].mean(axis=1),
                                   "iso_stderr": abundance[sampleVec.sampleID].std(axis=1).divide(
                                       math.sqrt(len(sampleVec.condition)))})

        geneSummary = pd.DataFrame(
            {"gene_id": geneRepExpression.gene_id, "gene_value": geneRepExpression[sampleVec.sampleID].mean(axis=1),
             "gene_stderr": geneRepExpression[sampleVec.sampleID].std(axis=1).divide(
                 math.sqrt(len(sampleVec.condition)))})

        combinedData = pd.merge(isoSummary, geneSummary, how="outer", on="gene_id")

        conditionSummary[sampleVec.condition.unique()[0]] = combinedData

    ## for each pairwise comparison of condition (or as controlled by the comparisonsToMake argument) the mean gene and isoform expression values are then used to calculate log2fc and IF

    if not quiet:
        print('Step 3 of 6: Making comparisons...')

    isoAnnot = pd.DataFrame()

    for comp in comparisonsToMake.iterrows():
        cond1data = conditionSummary[comp[1][0]]

        cond2data = conditionSummary[comp[1][1]]

        verticalCombined = pd.merge(cond1data, cond2data, how="inner", on=["gene_id", "isoform_id"],
                                    suffixes=["_1", "_2"])

        verticalCombined["condition_1"] = comp[1][0]

        verticalCombined["condition_2"] = comp[1][1]

        isoAnnot = isoAnnot.append(verticalCombined)

    ## add comparison data

    if not quiet:
        print('Step 4 of 6: Calculating DIF ...')

    # log2FC

    ps = foldChangePseudoCount

    isoAnnot["iso_log2_fc"] = np.log2((isoAnnot["iso_value_2"] + ps) / (isoAnnot["iso_value_1"] + ps))

    isoAnnot["gene_log2_fc"] = np.log2((isoAnnot["gene_value_2"] + ps) / (isoAnnot["gene_value_1"] + ps))

    # isoform fraction values

    isoAnnot["IF1"] = np.round(isoAnnot["iso_value_1"] / isoAnnot["gene_value_1"], decimals=4)

    isoAnnot["IF2"] = np.round(isoAnnot["iso_value_2"] / isoAnnot["gene_value_2"], decimals=4)

    isoAnnot["dIF"] = isoAnnot["IF2"] - isoAnnot["IF1"]

    ## filtering for eligible data

    if not quiet:
        print('Step 5 of 6: Filtering for eligible data...')

    ## extract genes with multiple isoforms

    # get unique gene - isoform pairs
    geneIsoOverview = isoAnnot[["isoform_id", "gene_id"]].drop_duplicates()

    # count isoforms per gene
    iso_counter = Counter(geneIsoOverview["gene_id"])

    # select genes with at least 2 isoforms
    genes_of_interest = [x for x in isoAnnot["gene_id"] if iso_counter[x] > 1]

    # reduce data to selected isoforms
    isoAnnot = isoAnnot.loc[isoAnnot["gene_id"].isin(genes_of_interest)]

    ## calculations

    if not quiet:
        print('Step 6 of 6: Adding switch information...')

    # add replicate number

    isoAnnot["nrReplicates_1"] = [conditions["nrReplicates"][conditions["condition"] == x].values[0] for x in
                                    isoAnnot["condition_1"]]

    isoAnnot["nrReplicates_2"] = [conditions["nrReplicates"][conditions["condition"] == x].values[0] for x in
                                    isoAnnot["condition_2"]]

    isoAnnot["switchDirection"] = ["up" if x > 0 else "down" for x in isoAnnot.dIF]

    # write ISOFORM REPORT:
    isoAnnot.to_csv(disease + "ISOFORM_REPORT_" + method + ".tsv", sep="\t")

    return isoAnnot


# switches, switch_pairs = detect_iswitches(method="spada", disease="breast invasive carcinoma")

# example_plot_switches("ABCF3", switch_pairs, switches)

switches, switch_pairs = detect_iswitches(method="spada", disease="breast invasive carcinoma")


#TODO##########################################################################################################################
#TODO##########################################################################################################################

#TODO rename this python file

#TODO outputs:
"""
- ISOFORM_REPORT: 
existing for spada: NO --> do it like for isa
existing for isa: breast_invasive_carcinoma_significantSwitches (früher switches.csv):
    gene_id,isoform_id,iso_value_1,iso_stderr_1,gene_value_1,gene_stderr_1,iso_value_2,iso_stderr_2,gene_value_2,gene_stderr_2,condition_1,condition_2,iso_log2_fc,gene_log2_fc,iso_q_value,gene_q_value,IF1,IF2,dIF,isoform_switch_q_value,gene_switch_q_value,nrReplicates_1,nrReplicates_2,gene_cv1,gene_lower_CI_1,gene_lower_CI_2,gene_var_1,gene_var_2,iso_var_1,iso_var_2,IF_var_1,IF_var_2,dIF_std_err,t_statistics,deg_free,p_value
min cols: 
    isoform_id,gene_id,condition_1,condition_2,IF1,IF2,dIF oder besser IF_case,IF_control ? 
--> for every isoform --> equal number of lines for both methods

- GENE_REPORT.tsv 
existing for spada: geneReport.tsv
    Gene    Symbol  Number-of-switches
not existing for isa, see summary
--> for every gene --> equal number of lines for both methods

- SWITCH_REPORT: 
existing for spada: significantSwitches / switches_spada
    Experiment,GeneId,Symbol,Control_transcript,Case_transcript,CDS_control,CDS_case,CDS_change,5_UTR_change,3_UTR_change,Samples,frac
existing for isa: switch_pairs kidney-clear-cell-carcinoma_isoformSwitchAnalyzePy.tsv: 
    GeneId	Control_transcript	Case_transcript
--> diverging number of lines for each method

- DOMAIN_REPORT:

--> diverging number of lines for each method

"""
