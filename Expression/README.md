# mRNA-miRNAs-cancer-pathways
This folder contains information to obtain gene and miRNA expression values and their interactions.

##Differential-expression
Inside this folder, you can find the results for a differential expression (DE) analysis, for each of the 15 tumours types. 
These results were obtained by using the DEAnalysis module from miARma-Seq and selecting edgeR as software to compare primary tumour sample against normal tissue samples.
More details are available in the article and inside the miARma-files directory


##miARma-files
Here you can find the ini files needed to process transcriptome raw data samples into differentially expressed genes/miRNAs and miRNA-mRNA target pairs using miARma-Seq. 
For each tumor type, two different ini files are provided, one for small RNASeq samples and the other for RNASeq samples. In both files a readcount process is performed and a DEAnalysis. 
Finally, in the RNASeq ini file a TargetPrediction step is performed to retrieve all possible miRNA-mRNA interactions among DE elements (miRNAs and genes).

##miRGate_miRNA-mRNA-pairs
In this folder, miRNA-mRNA interactions pairs obtained from miRGate for DE elements in each tumor type is stored.

##Scripts
This folder contains scripts to make global analysis combining the results from the individual 15 tumor types. That is: 
Summary_miRGate.pl: To find genes/miRNAs DE in several tumor types
miRNA-mRNA_stats.pl: 
Integrative_pairs_all_tumors.pl: 