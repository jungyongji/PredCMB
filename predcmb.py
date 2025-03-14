# PredCMB release version 1.0, Dec 19, 2024

from scipy.stats import norm
from scipy.stats import ranksums
import random
import statistics
from statsmodels.stats.multitest import fdrcorrection
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats
import multiprocessing
import argparse
import os

def init_worker(global_producing_only, global_reaction_abundance):
    global isProducingOnly, isReactionAbundanceConsidered
    isProducingOnly = global_producing_only
    isReactionAbundanceConsidered = global_reaction_abundance

def computeMetaboliteZ(pEnzyme, cEnzyme, gfName, gfPValue, gfLog2FoldChange, gfReactionAbundance, mbMaxAbsoluteZ):
    enzymeZ = []
    actualProducingEnzyme = 0
    actualConsumingEnzyme = 0

    for i in range(len(pEnzyme)):
        if pEnzyme[i] in gfName:
            actualProducingEnzyme += 1
            idx = gfName.index(pEnzyme[i])

            if gfLog2FoldChange[idx] > 0:
                newP = gfPValue[idx]/2
            else:
                newP = 1 - gfPValue[idx]/2

            ez = norm.ppf(1 - newP)

            if ez == float('inf'):
                ez = mbMaxAbsoluteZ
            elif ez == -1*float('inf'):
                ez = -1*mbMaxAbsoluteZ

            if isReactionAbundanceConsidered == 'True':
                ez *= gfReactionAbundance[idx]

            enzymeZ.append(ez)

    if isProducingOnly == 'False':
        for i in range(len(cEnzyme)):
            if cEnzyme[i] in gfName:
                actualConsumingEnzyme += 1
                idx = gfName.index(cEnzyme[i])

                if gfLog2FoldChange[idx] > 0:
                    newP = 1 - gfPValue[idx]/2
                else:
                    newP = gfPValue[idx]/2

                ez = norm.ppf(1 - newP)

                if ez == float('inf'):
                    ez = mbMaxAbsoluteZ
                elif ez == -1*float('inf'):
                    ez = -1*mbMaxAbsoluteZ

                if isReactionAbundanceConsidered == 'True':
                    ez *= gfReactionAbundance[idx]

                enzymeZ.append(ez)

    if len(enzymeZ) > 0:
        return sum(enzymeZ)/len(enzymeZ), actualProducingEnzyme, actualConsumingEnzyme
    else:
        return 0, actualProducingEnzyme, actualConsumingEnzyme

def computeMetaboliteZWeighted(pEnzyme, pEnzymeWeight, cEnzyme, cEnzymeWeight, gfName, gfPValue, gfLog2FoldChange, \
                               gfReactionAbundance, mbMaxAbsoluteZ):
    enzymeZ = []
    actualProducingEnzyme = 0
    actualConsumingEnzyme = 0

    for i in range(len(pEnzyme)):
        if pEnzyme[i] in gfName:
            actualProducingEnzyme += 1
            idx = gfName.index(pEnzyme[i])

            if gfLog2FoldChange[idx] > 0:
                newP = gfPValue[idx] / 2
            else:
                newP = 1 - gfPValue[idx] / 2

            ez = norm.ppf(1 - newP)
            ez *= pEnzymeWeight[i]

            if ez == float('inf'):
                ez = mbMaxAbsoluteZ
            elif ez == -1*float('inf'):
                ez = -1*mbMaxAbsoluteZ

            if isReactionAbundanceConsidered == 'True':
                ez *= gfReactionAbundance[idx]

            enzymeZ.append(ez)

    if isProducingOnly == 'False':
        for i in range(len(cEnzyme)):
            if cEnzyme[i] in gfName:
                actualConsumingEnzyme += 1
                idx = gfName.index(cEnzyme[i])

                if gfLog2FoldChange[idx] > 0:
                    newP = 1 - gfPValue[idx] / 2
                else:
                    newP = gfPValue[idx] / 2

                ez = norm.ppf(1 - newP)
                ez *= cEnzymeWeight[i]

                if ez == float('inf'):
                    ez = mbMaxAbsoluteZ
                elif ez == -1*float('inf'):
                    ez = -1*mbMaxAbsoluteZ

                if isReactionAbundanceConsidered == 'True':
                    ez *= gfReactionAbundance[idx]

                enzymeZ.append(ez)

    if len(enzymeZ) > 0:
        return sum(enzymeZ) / len(enzymeZ), actualProducingEnzyme, actualConsumingEnzyme
    else:
        return 0, actualProducingEnzyme, actualConsumingEnzyme

def compute_metabolite_z_weighted(i, mb, pe, pew, ce, cew, fgfn, fgfp, fgfl2fc, fgfra, mbMaxAbsoluteZ):
    metaboliteZ_i, metaboliteProducingEnzymeAmount_i, metaboliteConsumingEnzymeAmount_i = \
        computeMetaboliteZWeighted(pe[i], pew[i], ce[i], cew[i], fgfn, fgfp, fgfl2fc, fgfra, mbMaxAbsoluteZ)
    result = (i, metaboliteZ_i, metaboliteProducingEnzymeAmount_i, metaboliteConsumingEnzymeAmount_i)
    return result

def compute_metabolite_z(i, mb, pe, ce, fgfn, fgfp, fgfl2fc, fgfra, mbMaxAbsoluteZ):
    metaboliteZ_i, metaboliteProducingEnzymeAmount_i, metaboliteConsumingEnzymeAmount_i = \
        computeMetaboliteZ(pe[i], ce[i], fgfn, fgfp, fgfl2fc, fgfra, mbMaxAbsoluteZ)
    result = (i, metaboliteZ_i, metaboliteProducingEnzymeAmount_i, metaboliteConsumingEnzymeAmount_i)
    return result

def evaluateMetaboiteP(mZ, pEnzymeAmount, cEnzymeAmount, gfPValue, gfLog2FoldChange, gfReactionAbundance, permutation, \
                       mbMaxAbsoluteZ):
    if (pEnzymeAmount + cEnzymeAmount) == 0:
        return 1

    if (pEnzymeAmount + cEnzymeAmount) > len(gfPValue):
        print("Too few enzymatic gene families with differentiality. Cannot properly evaluate metaboite p-value.")
        return 1

    permutedMZ = [0]*permutation

    for i in range(permutation):
        idx = random.sample(range(len(gfPValue)), k = (pEnzymeAmount + cEnzymeAmount))
        pIdx = idx[:pEnzymeAmount]
        cIdx = idx[pEnzymeAmount:]
        enzymeZ = []

        for j in range(len(pIdx)):
            if gfLog2FoldChange[pIdx[j]] > 0:
                newP = gfPValue[pIdx[j]]/2
            else:
                newP = 1 - gfPValue[pIdx[j]]/2

            ez = norm.ppf(1 - newP)

            if ez == float('inf'):
                ez = mbMaxAbsoluteZ
            elif ez == -1*float('inf'):
                ez = -1*mbMaxAbsoluteZ

            if isReactionAbundanceConsidered == 'True':
                ez *= gfReactionAbundance[pIdx[j]]

            enzymeZ.append(ez)

        if isProducingOnly == 'False':
            for j in range(len(cIdx)):
                if gfLog2FoldChange[cIdx[j]] > 0:
                    newP = 1 - gfPValue[cIdx[j]]/2
                else:
                    newP = gfPValue[cIdx[j]]/2

                ez = norm.ppf(1 - newP)

                if ez == float('inf'):
                    ez = mbMaxAbsoluteZ
                elif ez == -1*float('inf'):
                    ez = -1*mbMaxAbsoluteZ

                if isReactionAbundanceConsidered == 'True':
                    ez *= gfReactionAbundance[cIdx[j]]

                enzymeZ.append(ez)

        if len(enzymeZ) > 0:
            permutedMZ[i] = sum(enzymeZ)/len(enzymeZ)

    permutedMZMean = statistics.mean(permutedMZ)
    permutedMZStdev = statistics.stdev(permutedMZ)
    mZNullAdjusted = (mZ - permutedMZMean)/permutedMZStdev

    if mZNullAdjusted > 0:
        return 1 - norm.cdf(mZNullAdjusted)
    else:
        return norm.cdf(mZNullAdjusted)

def evaluateMetaboitePWeighted(mZ, pEnzymeAmount, pEnzymeWeight, cEnzymeAmount, cEnzymeWeight, gfPValue, \
                               gfLog2FoldChange, gfReactionAbundance, permutation, mbMaxAbsoluteZ):
    if (pEnzymeAmount + cEnzymeAmount) == 0:
        return 1

    if (pEnzymeAmount + cEnzymeAmount) > len(gfPValue):
        print("Too few enzymatic gene families with differentiality. Cannot properly evaluate metaboite p-value.")
        return 1

    permutedMZ = [0] * permutation

    for i in range(permutation):
        idx = random.sample(range(len(gfPValue)), k=(pEnzymeAmount + cEnzymeAmount))
        pIdx = idx[:pEnzymeAmount]
        cIdx = idx[pEnzymeAmount:]
        pwidx = random.sample(range(len(pEnzymeWeight)), k = pEnzymeAmount)
        cwidx = random.sample(range(len(cEnzymeWeight)), k=cEnzymeAmount)
        enzymeZ = []

        for j in range(len(pIdx)):
            if gfLog2FoldChange[pIdx[j]] > 0:
                newP = gfPValue[pIdx[j]] / 2
            else:
                newP = 1 - gfPValue[pIdx[j]] / 2

            ez = norm.ppf(1 - newP)
            ez *= pEnzymeWeight[pwidx[j]]

            if ez == float('inf'):
                ez = mbMaxAbsoluteZ
            elif ez == -1*float('inf'):
                ez = -1*mbMaxAbsoluteZ

            if isReactionAbundanceConsidered == 'True':
                ez *= gfReactionAbundance[pIdx[j]]

            enzymeZ.append(ez)

        if isProducingOnly == 'False':
            for j in range(len(cIdx)):
                if gfLog2FoldChange[cIdx[j]] > 0:
                    newP = 1 - gfPValue[cIdx[j]] / 2
                else:
                    newP = gfPValue[cIdx[j]] / 2

                ez = norm.ppf(1 - newP)
                ez *= cEnzymeWeight[cwidx[j]]

                if ez == float('inf'):
                    ez = mbMaxAbsoluteZ
                elif ez == -1*float('inf'):
                    ez = -1*mbMaxAbsoluteZ

                if isReactionAbundanceConsidered == 'True':
                    ez *= gfReactionAbundance[cIdx[j]]

                enzymeZ.append(ez)

        if len(enzymeZ) > 0:
            permutedMZ[i] = sum(enzymeZ) / len(enzymeZ)

    permutedMZMean = statistics.mean(permutedMZ)
    permutedMZStdev = statistics.stdev(permutedMZ)
    mZNullAdjusted = (mZ - permutedMZMean) / permutedMZStdev

    if mZNullAdjusted > 0:
        return 1 - norm.cdf(mZNullAdjusted)
    else:
        return norm.cdf(mZNullAdjusted)

def evaluate_metabolite_p_weighted(i, mb, mz, pea, pew, cea, cew, fgfp, fgfl2fc, fgfra, perm, mbMaxAbsoluteZ):
    metaboliteP_i = evaluateMetaboitePWeighted(mz[i], pea[i], pew[i], cea[i], cew[i], fgfp, fgfl2fc, fgfra, perm, mbMaxAbsoluteZ)
    result = (i, metaboliteP_i)
    return result

def evaluate_metabolite_p(i, mb, mz, pea, cea, fgfp, fgfl2fc, fgfra, perm, mbMaxAbsoluteZ):
    metaboliteP_i = evaluateMetaboiteP(mz[i], pea[i], cea[i], fgfp, fgfl2fc, fgfra, perm, mbMaxAbsoluteZ)
    result = (i, metaboliteP_i)
    return result

def process_class(cls, data):
    class_scores = data[data['class'] == cls]['Z-score']
    other_scores = data[data['class'] != cls]['Z-score']
    stat, p_value = ranksums(class_scores, other_scores)
    result = {
        'class': cls,
        'Median of normalized Z': statistics.median(class_scores),
        'p-value': p_value
    }
    return result, p_value

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='PredCMB analysis with optional command line arguments.')
    parser.add_argument('geneFamilyAbundanceFileName', type = str, help = 'Gene family abundance file')
    parser.add_argument('sampleConditionFileName', type = str, help = 'Sample condition file')
    parser.add_argument('referenceCondition', type = str, help = 'Reference condition')
    parser.add_argument('targetCondition', type = str, help = 'Target condition')
    parser.add_argument('outputName', type = str, help = 'Name of result output')
    parser.add_argument('--cpu', type = int, default = 1, help = \
        'Number of CPU cores to use (default: 1)')
    parser.add_argument('--compound_enzyme', type = str, default = 'compound_enzymes_net.txt', help = \
        'A file that contains information of compounds and enzymes that contribute to the abundance of the compounds (default: compound_enzymes_net.txt)')
    parser.add_argument('--pval', type = float, default = 0.05, help = \
        'Threshold for differential abundance p-value of gene family (default: 0.05)')
    parser.add_argument('--log2fc', type = float, default = 1, help = \
        'Threshold for |log2(fold-change ratio)| of differentially abundant gene family (default: 1)')
    parser.add_argument('--cpd_mclass', type = str, default = 'hmdb_class_curation.tsv', help = \
        'A file that maps compounds to metabolite classes (default: hmdb_class_curation.tsv)')
    parser.add_argument('--cpd_name', type = str, default = 'kegg_cpd_name.tsv', help = \
        'A file that maps compounds to metabolite names (default: kegg_cpd_name.tsv)')
    parser.add_argument('--perms', type = int, default = 1000, help = \
        'Number of random permutations in computing p-values of metabolite z-values (default: 1000)')
    parser.add_argument('--maxZ', type = float, default = 7, help = \
        'Limited maximum value of metabolite z-scores (default: 7)')
    parser.add_argument('--enzyme_weighted', type = str, choices = ['True', 'False'], default = 'False', \
                        help = 'Whether the enzymes are weighted based on their stoichiometric coefficient or not (default: False)')
    parser.add_argument('--mclass_size', type = int, default = 5, help = \
        'Number of minimum metabolite per class for class summarization (default: 5)')
    parser.add_argument('--reaction_abundance', type=str, choices=['True', 'False'], default='True', \
                        help='Whether the abundance of each reaction is considered or not (default: True)')
    parser.add_argument('--producing_only', type=str, choices=['True', 'False'], default='False', \
                        help='Whether only the producing enzymes are considered or not (default: False)')
    args = parser.parse_args()
    geneFamilyAbundanceFileName = args.geneFamilyAbundanceFileName
    sampleConditionFileName = args.sampleConditionFileName
    referenceCondition = args.referenceCondition
    targetCondition = args.targetCondition
    outputName = args.outputName
    numberOfCPU = args.cpu
    compoundEnzymeFileName = args.compound_enzyme
    pValueThreshold = args.pval
    absoluteLog2FoldChangeThreshold = args.log2fc
    cpdClassFileName = args.cpd_mclass
    cpdNameFileName = args.cpd_name
    permutationAmount = args.perms
    metaboliteMaxAbsoluteZ = args.maxZ
    isEnzymeWeighted = args.enzyme_weighted
    minMetabolitesPerClass = args.mclass_size
    isReactionAbundanceConsidered = args.reaction_abundance
    isProducingOnly = args.producing_only

    if not os.path.exists(outputName):
        os.makedirs(outputName)

    print('Reading ' + geneFamilyAbundanceFileName)
    df = pd.read_csv(geneFamilyAbundanceFileName, sep = '\t', index_col = 0)
    df_int = df.round(0).astype('int64')
    df_int_T = df_int.T
    print('Reading ' + sampleConditionFileName)
    sc_df = pd.read_csv(sampleConditionFileName, sep = '\t', index_col = 0)
    samples_to_keep = sc_df['condition'].isin([referenceCondition, targetCondition])
    geneFamilyAbundance = df_int_T.loc[samples_to_keep]
    sampleCondition = sc_df.loc[samples_to_keep]
    sampleCondition = sampleCondition.reindex(geneFamilyAbundance.index)

    print('Reference: ' + referenceCondition + '\tTarget: ' + targetCondition + '\t' + str(len(geneFamilyAbundance))\
          + ' of ' + str(len(samples_to_keep)) + ' samples')
    genefamily_to_keep = geneFamilyAbundance.columns[geneFamilyAbundance.sum(axis = 0) >= 10]
    geneFamilyAbundance_filtered = geneFamilyAbundance[genefamily_to_keep]
    print('Gene families: ' + str(geneFamilyAbundance.shape[1]) + ', Gene families after filtering: ' + \
          str(len(genefamily_to_keep)) + ', Removed gene family: ' + str(geneFamilyAbundance.shape[1] - \
                                                                         len(genefamily_to_keep)))
    print('Running PyDESeq2 with ' + str(numberOfCPU) + ' CPU')
    print("Unique conditions in metadata:", sampleCondition['condition'].unique())
    inference = DefaultInference(n_cpus = numberOfCPU)
    dds = DeseqDataSet(counts = geneFamilyAbundance_filtered, metadata = sampleCondition, ref_level = ['condition', \
                                                                                                       referenceCondition], \
                       refit_cooks = True, inference = inference)
    
    dds.deseq2()
    stat_res = DeseqStats(dds, inference = inference)
    stat_res.summary()
    stat_res.results_df.to_csv(outputName + '/DEG.txt', sep = '\t')
    print('PyDESeq2 complete')
    geneFamilyName = stat_res.results_df.index.tolist()
    geneFamilyPValue = stat_res.results_df.padj.tolist()
    geneFamilyLog2FoldChange = stat_res.results_df.log2FoldChange.tolist()
    geneFamilyBaseMean = stat_res.results_df.baseMean.tolist()
    metabolite = []
    producingEnzyme = []
    consumingEnzyme = []
    producingEnzymeWeight = []
    consumingEnzymeWeight = []
    print('Reading ' + compoundEnzymeFileName)
    inFile = open(compoundEnzymeFileName, 'r')
    line = inFile.readline()
    line = inFile.readline()

    while line != '':
        line_str = line.rstrip('\n')
        tk = line_str.split('\t')
        metabolite.append(tk[0])

        if tk[1] == '':
            producingEnzyme.append([])
            producingEnzymeWeight.append([])
        else:
            producingEnzyme.append(tk[1].split(','))
            pew = tk[2].split(',')
            pew_num = [float(x) for x in pew]
            producingEnzymeWeight.append(pew_num)

        if tk[3] == '':
            consumingEnzyme.append([])
            consumingEnzymeWeight.append([])
        else:
            consumingEnzyme.append(tk[3].split(','))
            cew = tk[4].split(',')
            cew_num = [float(x) for x in cew]
            consumingEnzymeWeight.append(cew_num)

        line = inFile.readline()

    print('Reading ' + cpdClassFileName)
    df = pd.read_csv(cpdClassFileName, sep = '\t')
    CPDClassName = dict(zip(df['cpd'], df['class']))
    print('Reading ' + cpdNameFileName)
    df = pd.read_csv(cpdNameFileName, sep = '\t')
    CPDName = dict(zip(df['cpd'], df['name']))
    filteredGeneFamilyName = []
    filteredGeneFamilyPValue = []
    filteredGeneFamilyLog2FoldChange = []
    filteredGeneFamilyBaseMean = []

    for i in range(len(geneFamilyName)):
        if geneFamilyPValue[i] < pValueThreshold and abs(geneFamilyLog2FoldChange[i]) > absoluteLog2FoldChangeThreshold:
            filteredGeneFamilyName.append(geneFamilyName[i])
            filteredGeneFamilyPValue.append(geneFamilyPValue[i])
            filteredGeneFamilyLog2FoldChange.append(geneFamilyLog2FoldChange[i])
            filteredGeneFamilyBaseMean.append(geneFamilyBaseMean[i])

    avr_fgbm = sum(filteredGeneFamilyBaseMean)/len(filteredGeneFamilyBaseMean)
    filteredGeneFamilyReactionAbundance = [x/avr_fgbm for x in filteredGeneFamilyBaseMean]
    print('Differential gene family: ' + str(len(filteredGeneFamilyName)))
    print('Computing metabolite z values with ' + str(numberOfCPU) + ' CPU')
    metaboliteZ = [0]*len(metabolite)
    metaboliteProducingEnzymeAmount = [0]*len(metabolite)
    metaboliteConsumingEnzymeAmount = [0]*len(metabolite)

    if isEnzymeWeighted == 'True':
        args = [(i, metabolite, producingEnzyme, producingEnzymeWeight, consumingEnzyme, consumingEnzymeWeight, \
                 filteredGeneFamilyName, filteredGeneFamilyPValue, filteredGeneFamilyLog2FoldChange, \
                 filteredGeneFamilyReactionAbundance, metaboliteMaxAbsoluteZ) for i in range(len(metabolite))]
    else:
        args = [(i, metabolite, producingEnzyme, consumingEnzyme, filteredGeneFamilyName, filteredGeneFamilyPValue, \
                 filteredGeneFamilyLog2FoldChange, filteredGeneFamilyReactionAbundance, metaboliteMaxAbsoluteZ) for i \
                in range(len(metabolite))]

    with multiprocessing.Pool(numberOfCPU, initializer=init_worker, initargs=(isProducingOnly, isReactionAbundanceConsidered)) as pool:
        if isEnzymeWeighted == 'True':
            results = pool.starmap(compute_metabolite_z_weighted, args)
        else:
            results = pool.starmap(compute_metabolite_z, args)

    for i, metaboliteZ_i, metaboliteProducingEnzymeAmount_i, metaboliteConsumingEnzymeAmount_i in results:
        metaboliteZ[i] = metaboliteZ_i
        metaboliteProducingEnzymeAmount[i] = metaboliteProducingEnzymeAmount_i
        metaboliteConsumingEnzymeAmount[i] = metaboliteConsumingEnzymeAmount_i

    print('Computing metabolite p-value with ' + str(numberOfCPU) + ' CPU')
    metaboliteP = [1]*len(metabolite)

    if isEnzymeWeighted == 'True':
        args = [(i, metabolite, metaboliteZ, metaboliteProducingEnzymeAmount, producingEnzymeWeight, \
                 metaboliteConsumingEnzymeAmount, consumingEnzymeWeight, filteredGeneFamilyPValue, \
                 filteredGeneFamilyLog2FoldChange, filteredGeneFamilyReactionAbundance, permutationAmount, \
                 metaboliteMaxAbsoluteZ) for i in range(len(metabolite))]
    else:
        args = [(i, metabolite, metaboliteZ, metaboliteProducingEnzymeAmount, metaboliteConsumingEnzymeAmount, \
                 filteredGeneFamilyPValue, filteredGeneFamilyLog2FoldChange, filteredGeneFamilyReactionAbundance, \
                 permutationAmount, metaboliteMaxAbsoluteZ) for i in range(len(metabolite))]

    with multiprocessing.Pool(numberOfCPU, initializer=init_worker, initargs=(isProducingOnly, isReactionAbundanceConsidered)) as pool:
        if isEnzymeWeighted == 'True':
            results = pool.starmap(evaluate_metabolite_p_weighted, args)
        else:
            results = pool.starmap(evaluate_metabolite_p, args)


    for i, metaboliteP_i in results:
        metaboliteP[i] = metaboliteP_i

    print('FDR-adjusting metabolite p-values')
    ez_amnt = [a + b for a, b in zip(metaboliteProducingEnzymeAmount, metaboliteConsumingEnzymeAmount)]
    valid_indices = [i for i, value in enumerate(ez_amnt) if value > 0]
    valid_p_values = [metaboliteP[i] for i in valid_indices]
    corrected_p_values = fdrcorrection(valid_p_values)[1]
    metaboliteAdjustedP = metaboliteP.copy()

    for idx, corrected_value in zip(valid_indices, corrected_p_values):
        metaboliteAdjustedP[idx] = corrected_value

    print('Writing metabolite result for ' + outputName)
    metaboliteEnzymeAmount = [0]*len(metabolite)
    outFile = open(outputName + '/metabolite.txt', 'w')
    outFile.write("Metabolite KEGG ID\tMetabolite name\tProducing enzyme\tConsuming enzyme\tTotal enzyme\tMetabolite class name\tMetabolite change Z-score\tp-value\tAdjusted p-value\n")

    for i in range(len(metabolite)):
        outFile.write(metabolite[i] + "\t")
        outFile.write(CPDName.get(metabolite[i], 'N/A') + '\t')
        outFile.write(str(metaboliteProducingEnzymeAmount[i]) + '\t')
        outFile.write(str(metaboliteConsumingEnzymeAmount[i]) + '\t')
        metaboliteEnzymeAmount[i] = metaboliteProducingEnzymeAmount[i] + metaboliteConsumingEnzymeAmount[i]
        outFile.write(str(metaboliteEnzymeAmount[i]) + '\t')
        outFile.write(CPDClassName.get(metabolite[i], 'N/A') + '\t')
        outFile.write(str(metaboliteZ[i]) + "\t")
        outFile.write(str(metaboliteP[i]) + "\t")
        outFile.write(str(metaboliteAdjustedP[i]) + "\n")

    outFile.close()
    metaboliteClassName = []

    for i in range(len(metabolite)):
        metaboliteClassName.append(CPDClassName.get(metabolite[i], 'N/A'))

    print('Summarizing metabolite class statistics with ' + str(numberOfCPU) + ' CPU')
    data2 = pd.DataFrame({'Z-score': metaboliteZ, 'class': metaboliteClassName, 'enzymeAmount': metaboliteEnzymeAmount})
    filteredData2 = data2[data2['enzymeAmount'] > 0]
    unique_classes = filteredData2['class'].unique()
    args = [(cls, filteredData2) for cls in unique_classes if len(filteredData2[filteredData2['class'] == cls]) >= \
            minMetabolitesPerClass]

    with multiprocessing.Pool(numberOfCPU) as pool:
        results_with_pvalues = pool.starmap(process_class, args)

    results, p_values = zip(*results_with_pvalues)
    print('FDR-adjusting metabolite class p-values')
    pvals_corrected = fdrcorrection(p_values)[1]
    print('Writing metabolite class summary for ' + outputName)

    for i in range(len(results)):
        results[i]['Adjusted p-value'] = pvals_corrected[i]

    results_df = pd.DataFrame(results)
    results_df.to_csv(outputName + '/metabolite_class.txt', sep='\t', index=False)