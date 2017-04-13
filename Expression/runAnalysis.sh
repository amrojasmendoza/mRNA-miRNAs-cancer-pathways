#!/bin/bash
clear
echo "##########################################################################

Here you can find all data to reproduce the results published in 
https://www.nature.com/articles/srep46101 entitled:
 
\"Novel miRNA-mRNA interactions conserved in essential cancer pathways\"
Eduardo Andrés-León, Ildefonso Cases, Sergio Alonso & Ana M. Rojas

##########################################################################
"
sleep 5

echo "General pairs: "
echo ""
sleep 2
perl src/1.Summary_generalPairs_data.pl -t DE -d RNASeq
echo ""
echo "1) Differentially expressed genes from 15 different tumors and 7 pathways will be studied"
sleep 2
echo "2) Differentially expressed miRNAs from 15 different tumors and 7 pathways will be studied"
perl src/1.Summary_generalPairs_data.pl -t DE -d miRNASeq >> src/log.$$ 2>&1
sleep 2
echo "3) miRNA-mRNA pairs predicted by miRGate from 15 different tumors and DDR genes"
perl src/1.Summary_generalPairs_data.pl -t pairs -net DDR >> src/log.$$ 2>&1
sleep 2
echo "4) miRNA-mRNA pairs predicted by miRGate from 15 different tumors and Cell cycle genes"
perl src/1.Summary_generalPairs_data.pl -t pairs -net "Cell cycle" >> src/log.$$ 2>&1
sleep 2
echo "5) miRNA-mRNA pairs predicted by miRGate from 15 different tumors and DNA replication genes"
perl src/1.Summary_generalPairs_data.pl -t pairs -net "DNA Rep" >> src/log.$$ 2>&1
sleep 2
echo "6) miRNA-mRNA pairs predicted by miRGate from 15 different tumors and Telomeres genes"
perl src/1.Summary_generalPairs_data.pl -t pairs -net Telomeres >> src/log.$$ 2>&1
sleep 2
echo "7) miRNA-mRNA pairs predicted by miRGate from 15 different tumors and Necrosis genes"
perl src/1.Summary_generalPairs_data.pl -t pairs -net Necrosis >> src/log.$$ 2>&1
sleep 2
echo "8) miRNA-mRNA pairs predicted by miRGate from 15 different tumors and Apoptosis genes"
perl src/1.Summary_generalPairs_data.pl -t pairs -net Apoptosis >> src/log.$$ 2>&1
sleep 2
echo "9) miRNA-mRNA pairs predicted by miRGate from 15 different tumors and Senescence genes"
perl src/1.Summary_generalPairs_data.pl -t pairs -net Senescence >> src/log.$$ 2>&1
sleep 2
echo "All temporary results are stored in src/DE_data_TCGA/ and all output is saved in src/log.$$"
sleep 5
echo ""
echo "Obtaining miRNA-mRNA pairs conserved in the majority of tumor types"
perl src/2.Integration_generalPairs_data.pl
echo "Result saved in general.pairs.csv file"

echo ""
echo "Lung pairs: "
echo ""

echo "1) miRNA-mRNA pairs predicted by miRGate from Lung tumors and DDR genes"
perl src/3.Summary_LungPairs_data.pl -t pairs -net DDR >> src/log.$$ 2>&1
sleep 1
echo "2) miRNA-mRNA pairs predicted by miRGate from Lung tumors and Cell cycle genes"
perl src/3.Summary_LungPairs_data.pl -t pairs -net "Cell cycle" >> src/log.$$ 2>&1
sleep 1
echo "3) miRNA-mRNA pairs predicted by miRGate from Lung tumors and DNA Replication genes"
perl src/3.Summary_LungPairs_data.pl -t pairs -net "DNA Rep" >> src/log.$$ 2>&1
sleep 1
echo "4) miRNA-mRNA pairs predicted by miRGate from Lung tumors and Telomere genes"
perl src/3.Summary_LungPairs_data.pl -t pairs -net Telomeres >> src/log.$$ 2>&1
sleep 1
echo "5) miRNA-mRNA pairs predicted by miRGate from Lung tumors and Necrosis genes"
perl src/3.Summary_LungPairs_data.pl -t pairs -net Necrosis >> src/log.$$ 2>&1
sleep 1
echo "6) miRNA-mRNA pairs predicted by miRGate from Lung tumors and Apoptosis genes"
perl src/3.Summary_LungPairs_data.pl -t pairs -net Apoptosis >> src/log.$$ 2>&1
sleep 1
echo "7) miRNA-mRNA pairs predicted by miRGate from Lung tumors and Senescence genes"
perl src/3.Summary_LungPairs_data.pl -t pairs -net Senescence >> src/log.$$ 2>&1
sleep 1

echo ""
echo "Obtaining exclusive miRNA-mRNA pairs in Lung tumors (LUAD and LUSC)"
perl src/4.Integration_LungPairs_data.pl -f src/DE_data_TCGA/TCGA_pairs_DE_pathway_agreement2.xls -t unique -paper
echo "Result saved in lung.pairs.csv file"

echo ""
echo "Please check temporary data and src/DE_data_TCGA and final result files: general.pairs.csv and lung.pairs.csv"
echo ""

