# OptiCon: a general framework for systematic and de novo identification of synergistic regulators as candidates for combination therapy
## OptiCon Framework
Most combination therapies are developed based on targets of existing drugs, which only represent a small portion of the human proteome. We introduce a network controllability-based method, OptiCon, for de novo identification of synergistic regulators as candidates for combination therapy.OptiCon represents a general framework for systemic and de novo identification of synergistic regulators underlying a cellular state transition.

By using gene expression as a constraint in the standard network controllability framework, OptiCon first identifies a set of optimal control nodes (OCNs) in a disease-perturbed gene regulatory network. The identified OCNs exert maximal control over deregulated pathways but minimal control over pathways that are not perturbed by the disease. Next, using a synergy score that combines both genetic mutation and gene functional interaction information, OptiCon identifies a set of synergistic OCNs as key regulators in the disease-perturbed network, which can serve as candidate targets for combination therapy.
## Use OptiCon
### Running environment
operation platform: win10;  
CPU: i7;  
memory: 16G or more;  
storage：10GB (It depends on the amount of data);  
JDK version：12.0.1.
### Create the required folder
1. Create an empty folder, such as a "jar" folder, with a variable name.  
2. In the "jar" folder, create three folders: "CRfile", "inputData" and "output" respectively. Note that the names of these three files must be consistent with the above, or the jar package will fail to run.  
3. Put the "OptiCon_java.jar" into the "jar" folder.
### Prepare input data
Put the following data in the "inputData" folder:  
1. “DScore.txt” contains data in two tab-delimited columns. Each row includes a gene identifier (first column) and its corresponding p-value of differential expression (second column) under two conditions (e.g. diseased vs. healthy). The gene identifier is based on Entrez Gene ID with a prefix “En_”. Duplicated gene identifiers are not allowed.  
2. “GeneExpression.txt” contains data in a matrix format. Rows represent genes and columns represent samples. Each entry in the matrix contains a gene expression value in a specific sample. Duplicated gene identifiers are not allowed.   
3.  “RecurMutant_entrez.txt” contains a list of genes (based on Entrez gene IDs) that are known to harbor recurrent somatic mutations in the cancer type of interest. Duplicated gene identifiers are not allowed.  
4.  If you want to identify synergistic key regulators using a customized directed network, Format your network data into a tab-delimited file “MyGeneNetwork.txt”. Each row represents a directed edge from the node in the first column to the node in the second column. Attention: In identifying synergistic OCNs pairs, we use “CancerCensus_En.txt” (Annotation for cancer (driver) genes were downloaded from the Cancer Gene Census database58, which were confirmed to have recurrent somatic mutations in the specific cancer type using data from the Catalogue of Somatic Mutations in Cancer (COSMIC) database), which names genes used Entrez Gene ID with a prefix “En_”. If you want to use your customized directed network, you should be sure about the genes named the same in all the files.  
### Running OptiCon_java.jar
1. Open the "cmd" command window.  
2. Go to the directory where the “OptiCon_java.jar” package resides.  
3. Enter the following command to run the jar package：java -Xms8g -Xmx8g -jar OptiCon_java.jar [parameter 1] [parameter 2] [parameter 3] [parameter 4] [parameter 5] [parameter 6] [parameter 7] [parameter 8] [parameter 9] [parameter 10] [parameter 11]  
  Parameter 1: The directory of three folders("CRfile", "inputData" and "output")；  
  Parameter 2: The name of gene regulatory network file；  
  Parameter 3: The name of gene deregulation file;  
  Parameter 4: The name of gene expression file;  
  Parameter 5: The name of recurrent somatic mutations file;  
  Parameter 6: The separator of all input file;  
  Parameter 7: The number of SCC;  
  Parameter 8: The value of epsilon used in computing ICV;  
  Parameter 9: The value of lambda used in computing ICV;  
  Parameter 10: The boolean value of true or false, which means using the control range(in “CRfile” folder) that you calculated before or not;  
  Parameter 11: The number of parallel computing thread pools, which is usually set to the number of CPU cores.  
For example: java -Xms8g -Xmx8g -jar OptiCon_java.jar D:\JavaProject\jar MyGeneNetwork.txt DScore.txt GeneExpression.txt RecurMutant_entrez.txt \t 1000 0.18 0.3 true 6
### View the result in “output” folder after running
