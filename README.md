# PredCMB
**PredCMB** (**Pred**iction of **c**hanges in **m**eta**b**olites from shotgun metagenome) is a tool that predicts changes in individual metabolites from shotgun metagenomic sequencing data utilizing microbial gene-metabolite interaction.

- - -

# Requirments
1. R (version >= 4.0.2)
2. [DESeq2](https://github.com/mikelove/DESeq2) (version == 1.30)
3. [Piano ](https://github.com/varemo/piano) (version == 2.6.0)
4. Python (version >= 3.7.7)
5. [HUMAnN](https://github.com/biobakery/humann) (version == 3.0.0)

- - -

# Data preparation
PredCMB uses **gene families file** and **gene-metabolite interaction file** as input. The gene family file is one of the main output files generated by HUMAnN and represents the abundance of reads mapped to each gene family. Gene-metabolite interaction file is the interaction information constructed by integrating gene family, and the metabolite as results from its enzymatic reaction.  


For normalization of inter-sample depth, re-normalization is performed with copies per million (cpm) units using command provided by HUMAnN. Please refer to the human manual for details. We also screen the gene family representing the total abundance at the community level.  



 1. Re-normalization of genefamilies file to 'cpm' units

     ``$ humann_renorm_table --input /path/genefamilies.tsv --units "cpm" --output /path/genefamilies_cpm.tsv``
 
 2. Selection of the abundance of each gene family in the community

     ``$ humann_split_stratified_table --input /path/genefamilies_cpm.tsv --output /dir/``
 
 3. Conversion of gene family abundance to integer

     ``$ asinteger_genefamilies.py -i /path/dir/genefamilies_cpm_unstratified.tsv -o /path/abun_int_cpm.tsv``  



To assign uniref90 to metabolite, the mapping files of "uniref90 to ec", "uniref90 to ko" provided by the HUMAnN tool, and EC number, orthology, and reaction provided by KEGG were used as bridges. The eQuillator tool was used to determine the direction of the final product of bidirectional KEGG reaction. All of these processes are performed by running the construct_interaction.py file, which uniref_com_kerk_oc.tsv file in the interaction/database directory represents the gene-metabolite interaction file. It also provides the .rds format for easy loading from R.



- - -

# How to run


```
$ Rscript run_prediction.R -i $SAMPLE -m $METADATA -o $OUTPUT -r $CONTROL -gb $INTERACTION -p $PVALUE
```  

**Please refer to the column and row format of the input file in the example file below.**


`$SAMPLE` = Gene family abundance file in .tsv format.

`$METADATA` = Sample groups file in .tsv format.

`$OUTPUT` = The predicted metabolite file in .tsv format.

`$INTERACTION` = Gene family and metabolite interaction file in .tsv format.

`$PVALUE` = Significance level required in Deseq2 analysis.


- - -

# Example

## Example input file ##  


The example_data folder contains three example input files. These files are .tsv format.

**It is recommended that you look inside the input file before using the program.**



### 1. genefamilies_ex ###


The genefamilies_ex.tsv file represents the abundance of gene families for each sample.


```
                    sample1 sample2 sample3 sample4 sample5 sample6 sample7 sample8 sample9 sample10 sample11 sample12
UniRef90_A0A010Z266       1       1       0       4       1       2       0       0       8        0        0        0
UniRef90_A0A014AUH4      11       0      52      19      31      48       3       9      69        0       41        0
UniRef90_A0A015P063       6       0      20      32       3      12       0       2      44        0        3        0
UniRef90_A0A015P476      10       0       9       3       1      16       5       1       0        0        4        0
UniRef90_A0A015PAC1       0       0       0       3       0       0       0       1      27        0       19        0
UniRef90_A0A015PAY5       0       0       2       0       0       0       0       0       0        0        1        0
                    sample13 sample14 sample15 sample16 sample17 sample18 sample19 sample20 sample21 sample22 sample23
UniRef90_A0A010Z266        0        1        0        1        1        1        4       18        0        0        0
UniRef90_A0A014AUH4       11        7       73       44       51       10        0       50        2       71       21
UniRef90_A0A015P063        1        1        0        1        5        0        0        1        1        2        2
UniRef90_A0A015P476        2        0        0        0       35        1        1       18        1        1        1
UniRef90_A0A015PAC1        0        0        0        0      154       15        0        0        1        0        6
UniRef90_A0A015PAY5        0        0        0        0       15        1        0        8        0        0        0
                    sample24 sample25 sample26 sample27 sample28 sample29 sample30 sample31 sample32 sample33 sample34
UniRef90_A0A010Z266        0        0        0        0        0        0        0        0        1        0        0
UniRef90_A0A014AUH4       11       14       39       33        9       24       53        2       34       54       10
UniRef90_A0A015P063        2        3        5        4        0       10        7        1       27        5        0
UniRef90_A0A015P476        1        6        2        6        0        1        2        0        9        9        0
UniRef90_A0A015PAC1        0        3        0        0        0       38       15        1        0        0       11
UniRef90_A0A015PAY5        0        2        0        0        0        0        1        0        3        0        0
                    sample35 sample36 sample37 sample38 sample39 sample40
UniRef90_A0A010Z266        0        0        0        0        0        0
UniRef90_A0A014AUH4       43       21        5        1       16        4
UniRef90_A0A015P063        2       11        4        0        4        1
UniRef90_A0A015P476        0        2        1        0        1        2
UniRef90_A0A015PAC1        7        0        0        0        0        6
UniRef90_A0A015PAY5        0        0        1        0        0        0
```


### 2. metadata_ex ###


The metadata_ex file represents each sample and the group that contains the sample.


```
        diagnosis
sample1   Control
sample2   Control
sample3        CD
sample4        CD
sample5        CD
sample6        CD
```


### 3. interaction_ex ###


The interaction_ex.tsv file represents the gene family and metabolites produced by its enzyme reaction.


```
                uniref    cpd
1 UniRef90_A0A076JL57 C07490
2 UniRef90_A0A174AFY3 C07490
3 UniRef90_A0A174IMF5 C07490
4 UniRef90_A0A174NFL4 C07490
5 UniRef90_A0A174WZP7 C07490
6 UniRef90_A0A174YW45 C07490
```

- - -

## Example run ##


If you execute the command below, you perform a differential analysis of the gene abundance with experiment group based on the control group you specified. And it builds a network based on significant differential gene families and metabolite interactions, and predicts the amount of changes in metabolites based on the amount of changes in enzymatic gene families.


```
$ Rscript run_prediction.R -i genefamilies_ex.tsv -m metadata_ex.tsv -o output.tsv -r "Control" -gb interaction_ex.tsv -p 0.05
```


- - - 

## Example output ##


The output file provides a comprehensive overview of predicted metabolite changes based on microbial gene family abundance. Each entry includes the KEGG compound identity (metabolite_KEGG_ID), the name of the metabolite (metabolite_name), its class (metabolite_class_name), the z-value indicating the magnitude and direction of change (metabolite_z_value), and statistical significance measures (p-value and adjusted p-value).


| metabolite_KEGG_ID | metabolite_name              | metabolite_class_name                       | metabolite_z_value | statistics | metabolite_p-value | metabolite_adjusted_p-value |
|--------------------|------------------------------|---------------------------------------------|--------------------|------------|--------------------|-----------------------------|
| C00005             | NADPH                        | Benzenediols                                | 7.7598             | 4.82E-07   | 0.00015247         | 0.001031                    |
| C00017             | Protein                      | NA                                          | 0.1532             | 0.31934    | 0.45353            | 0.2204                      |
| C00021             | S-Adenosyl-L-homocysteine    | Lactones                                    | -0.58743           | 0.12552    | 0.2204             | 0.0031599                   |
| C00022             | Pyruvate                     | Alpha-keto acids and derivatives            | 6.3454             | 3.00E-05   | 0.0031599          | 0.001031                    |
| C00024             | Acetyl-CoA                   | Carbohydrates and carbohydrate conjugates   | 7.8399             | 6.53E-06   | 0.001031           | 0.001031                    |
| C00025             | L-Glutamate                  | Amino acids, peptides, and analogues        | 7.8941             | 5.66E-06   | 0.001031           | 0.001031                    |
```
| metabolite_KEGG_ID | metabolite_name              | metabolite_class_name                       | metabolite_z_value | statistics | metabolite_p-value | metabolite_adjusted_p-value |
|--------------------|------------------------------|---------------------------------------------|--------------------|------------|--------------------|-----------------------------|
| C00005             | NADPH                        | Benzenediols                                | 7.7598             | 4.82E-07   | 0.00015247         | 0.001031                    |
| C00017             | Protein                      | NA                                          | 0.1532             | 0.31934    | 0.45353            | 0.2204                      |
| C00021             | S-Adenosyl-L-homocysteine    | Lactones                                    | -0.58743           | 0.12552    | 0.2204             | 0.0031599                   |
| C00022             | Pyruvate                     | Alpha-keto acids and derivatives            | 6.3454             | 3.00E-05   | 0.0031599          | 0.001031                    |
| C00024             | Acetyl-CoA                   | Carbohydrates and carbohydrate conjugates   | 7.8399             | 6.53E-06   | 0.001031           | 0.001031                    |
| C00025             | L-Glutamate                  | Amino acids, peptides, and analogues        | 7.8941             | 5.66E-06   | 0.001031           | 0.001031                    |

```



- - -
## Options ##

```
usage: run_prediction.R [-h] [-i INPUT_CPM_FILE] [-g SAMPLE_GROUP_FILE]
                        [-o OUTPUT_FILE] [-r CONTROL_GROUP] [-c NUMBER_CORE]
                        [-gb ENZYME_METABOLITE_FILE]

optional arguments:
  -h, --help            Show this help message and exit
  -i INPUT_CPM_FILE, --input INPUT_CPM_FILE
                        Count data from microbial genefamily [REQUIRED]
  -g SAMPLE_GROUP_FILE, --group SAMPLE_GROUP_FILE
                        The total group of sampels [REQUIRED]
  -o OUTPUT_FILE, --output OUTPUT_FILE
                        Predicted individual metabolites [REQUIRED]
  -r CONTROL_GROUP, --reference CONTROL_GROUP
                        A control group under certein conditions [Default:NA]
  -c NUMBER_CORE, --core NUMBER_CORE
                        The number of R processes [Default:1]
  -gb ENZYME_METABOLITE_FILE, --mgx_mbx ENZYME_METABOLITE_FILE
                        Enzymatic gene-metabolite interaction [Default:NA]
  -p PVALUE, --pvalue PVALUE
                        Threshold of p-value
                   
```
