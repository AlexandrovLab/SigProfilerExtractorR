
[![Docs](https://img.shields.io/badge/docs-latest-blue.svg)](https://osf.io/t6j7u/wiki/home/) 
[![License](https://img.shields.io/badge/License-BSD\%202--Clause-orange.svg)](https://opensource.org/licenses/BSD-2-Clause)
[![Build Status](https://travis-ci.com/AlexandrovLab/SigProfilerExtractor.svg?branch=master)](https://travis-ci.com/AlexandrovLab/SigProfilerExtractor)

# SigProfilerExtractorR
An R wrapper for running the SigProfilerExtractor framework.

**INTRODUCTION**

The purpose of this document is to provide a guide for using the SigProfilerExtractor framework to extract the De Novo mutational signatures from a set of samples and decompose the De Novo signatures into the COSMIC signatures. An extensive Wiki page detailing the usage of this tool can be found at https://osf.io/t6j7u/wiki/home/. For users that prefer working in a Python environment, the tool is written in Python and can be found and installed from: https://github.com/AlexandrovLab/SigProfilerExtractor

# Table of contents
- [Installation](#installation)
- [Functions](#Functions)
  - [importdata](#importdata)
  - [sigProfilerExtractor](#sigProfilerExtractor)
  - [estimate_solution](#estimate_solution)
  - [decompose](#decompose)
- [Copyright](#copyright)
- [Contact Information](#contact)


## <a name="installation"></a> Installation
**PREREQUISITES**

devtools  (R) 
```
>> install.packages("devtools")
```
reticulate* (R) 
```
>> install.packages("reticulate")  
```

*Reticulate has a known bug of preventing python print statements from flushing to standard out. As a result, some of the typical progress messages are delayed.

**QUICK START GUIDE**

This section will guide you through the minimum steps required to extract mutational signatures from genomes:
1. First, install the python package using pip. The R wrapper still requires the python package:
```
                          pip install sigproextractor
```
2. Open an R session and ensure that your R interpreter recognizes the path to your python installation:
```
$ R
>> library("reticulate")
>> use_python("path_to_your_python")
>> py_config()
python:         /anaconda3/bin/python
libpython:      /anaconda3/lib/libpython3.6m.dylib
pythonhome:     /anaconda3:/anaconda3
version:        3.6.5 |Anaconda, Inc.| (default, Apr 26 2018, 08:42:37)  [GCC 4.2.1 Compatible Clang 4.0.1 (tags/RELEASE_401/final)]
numpy:          /anaconda3/lib/python3.6/site-packages/numpy
numpy_version:  1.16.1
```
If you do not see your python path listed, restart your R session and rerun the above commands in order.

2. Install SigProfilerExtractorR using devtools:
```
>>library("devtools")
>>install_github("AlexandrovLab/SigProfilerExtractorR")
```
3. Load the package in the same R session and install your desired reference genome as follows (available reference genomes are: GRCh37, GRCh38, mm9, and mm10):
```
>> library("SigProfilerExtractorR")
>> install("GRCh37", rsync=FALSE, bash=TRUE)
```

This will install the human 37 assembly as a reference genome. 

**SUPPORTED GENOMES**

Other available reference genomes are GRCh38, mm9 and mm10 (and genomes supported SigProfilerMatrixGenerator. 
Information about supported will be found at https://github.com/AlexandrovLab/SigProfilerMatrixGeneratorR

**Quick Example:**

Signatures can be extracted from vcf files or tab delimited mutational table using the sigprofilerextractor function.
```
>> help(sigprofilerextractor)
```
This will show the details about the sigprofilerextractor funtion.


```
>> library("SigProfilerExtractorR")
>> path_to_example_data <- importdata("matrix")
>> data <- path_to_example_data # here you can provide the path of your own data
>> sigprofilerextractor("matrix", 
                     "example_output", 
                     data, 
                     minimum_signatures=2,
                     maximum_signatures=3,
                     nmf_replicates=5,
                     min_nmf_iterations = 1000,
                     max_nmf_iterations =100000,
                     nmf_test_conv = 1000,
                     nmf_tolerance = 0.00000001)
```

The example file will generated in the working directory. Note that the parameters used in the above example are not optimal to get accurate signatures. Those are used only for a quick example. 

## <a name="functions"></a> Functions
The list of available functions are:
- importdata
- sigprofilerextractor
- estimate_solution
- decompose


### <a name="importdata"></a> importdata
Imports the path of example data.  
importdata(datatype).        
datatype:         	
The type of example data. There are two types: 1."vcf", 2."matrix".     

#### importdata Example

```R
library("SigProfilerExtractorR")
path_to_example_table = importdata("matrix")
data = path_to_example_table 
# This "data" variable can be used as a parameter of the "project" argument of the sigprofilerextractor function.

# To get help on the parameters and outputs of the "importdata" function, please use the following:
help(importdata)
```

### <a name="sigProfilerExtractor"></a> sigProfilerExtractor
    
Extracts mutational signatures from an array of samples.

```R 
sigprofilerextractor(input_type, out_put, input_data, reference_genome="GRCh37", opportunity_genome = "GRCh37", context_type = "default", exome = False, 
                         minimum_signatures=1, maximum_signatures=10, nmf_replicates=100, resample = True, batch_size=1, cpu=-1, gpu=False, 
                         nmf_init="alexandrov-lab-custom", precision= "single", matrix_normalization= "100X", seeds= "none", 
                         min_nmf_iterations= 10000, max_nmf_iterations=1000000, nmf_test_conv= 10000, nmf_tolerance= 1e-15,nnls_add_penalty=0.05,
                         nnls_remove_penalty=0.01, initial_remove_penalty=0.05, de_novo_fit_penalty=0.02,get_all_signature_matrices= False)
```

| Category | Parameter | Variable Type | Parameter Description |
| --------- | --------------------- | -------- |-------- |
| **Input Data** |  |  | |
|  | **input_type** | String | The type of input:<br><ul><li>"vcf": used for vcf format inputs.</li><li>"matrix": used for table format inputs using a tab seperated file.</li></ul> |
|  | **out_put** | String | The name of the output folder. The output folder will be generated in the current working directory.  |   
|  | **input_data** | String | Name of the input folder (in case of "vcf" type input) or the input file (in case of "table"  type input). The project file or folder should be inside the current working directory. For the "vcf" type input, the project has to be a folder which will contain the vcf files in vcf format or text formats. The "text" type projects have to be a file. | 
|  | **reference_genome** | String | The name of the reference genome. The default reference genome is "GRCh37". This parameter is applicable only if the input_type is "vcf". | 
|  | **opportunity_genome** | String | The build or version of the reference signatures for the reference genome. The default opportunity genome is GRCh37. If the input_type is "vcf", the genome_build automatically matches the input reference genome value. | 
|  | **context_type** | String | A string of mutaion context name/names separated by comma (","). The items in the list defines the mutational contexts to be considered to extract the signatures. The default value is "96,DINUC,ID", where "96" is the SBS96 context, "DINUC" is the DINUCLEOTIDE context and ID is INDEL context. | 
|  | **exome** | Boolean | Defines if the exomes will be extracted. The default value is "False".  | 
| **NMF Replicates** |  |  |  | 
|  | **minimum_signatures** | Positive Integer | The minimum number of signatures to be extracted. The default value is 1. | 
|  | **maximum_signatures** | Positive Integer | The maximum number of signatures to be extracted. The default value is 25. | 
|  | **nmf_replicates** | Positive Integer | The number of iteration to be performed to extract each number signature. The default value is 100. | 
|  | **resample** | Boolean | Default is True. If True, add poisson noise to samples by resampling. | 
|  | **seeds** | String | It can be used to get reproducible resamples for the NMF replicates. A path of a tab separated .txt file containing the replicated id and preset seeds in a two columns dataframe can be passed through this parameter. The Seeds.txt file in the results folder from a previous analysis can be used for the seeds parameter in a new analysis. The Default value for this parameter is "none". When "none", the seeds for resampling will be random for different analysis. | 
| **NMF Engines** |  |  |  | 
|  | **matrix_normalization** | String | Method of normalizing the genome matrix before it is analyzed by NMF. Default is value is "gmm". Other options are, "log2", "custom" or "none". | 
|  | **nmf_init** | String | The initialization algorithm for W and H matrix of NMF. Options are 'random', 'nndsvd', 'nndsvda', 'nndsvdar' and 'nndsvd_min'. Default is 'nndsvd_min'. | 
|  | **precision** | String | Values should be single or double. Default is single. | 
|  | **min_nmf_iterations** | Integer | Value defines the minimum number of iterations to be completed before NMF converges. Default is 10000. | 
|  | **max_nmf_iterations** | Integer | Value defines the maximum number of iterations to be completed before NMF converges. Default is 1000000. | 
|  | **nmf_test_conv** | Integer | Value defines the number number of iterations to done between checking next convergence. Default is 10000. | 
|  | **nmf_tolerance** | Float | Value defines the tolerance to achieve to converge. Default is 1e-15. | 
| **Execution** |  |  |  | 
|  | **cpu** | Integer | The number of processors to be used to extract the signatures. The default value is -1 which will use all available processors. | 
|  | **gpu** | Boolean | Defines if the GPU resource will used if available. Default is False. If True, the GPU resources will be used in the computation. | 
|  | **batch_size** | Integer | Will be effective only if the GPU is used. Defines the number of NMF replicates to be performed by each CPU during the parallel processing. Default is 1. | 
| **Solution Estimation Thresholds** |  |  |  | 
|  | **stability** | Float | Default is 0.8. The cutoff thresh-hold of the average stability. Solutions with average stabilities below this thresh-hold will not be considered. | 
|  | **min_stability** | Float | Default is 0.2. The cutoff thresh-hold of the minimum stability. Solutions with minimum stabilities below this thresh-hold will not be considered.  | 
|  | **combined_stability** | Float | Default is 1.0. The cutoff thresh-hold of the combined stability (sum of average and minimum stability). Solutions with combined stabilities below this thresh-hold will not be considered. | 
| **Decomposition** |  |  |  | 
|  | **de_novo_fit_penalty** | Float | Takes any positive float. Default is 0.02. Defines the weak (remove) thresh-hold cutoff to assign denovo signatures to a sample. | 
|  | **nnls_add_penalty** | Float | Takes any positive float. Default is 0.05. Defines the strong (add) thresh-hold cutoff to assign COSMIC signatures to a sample. | 
|  | **nnls_remove_penalty** | Float | Takes any positive float. Default is 0.01. Defines the weak (remove) thresh-hold cutoff to assign COSMIC signatures to a sample. | 
|  | **initial_remove_penalty** | Float | Takes any positive float. Default is 0.05. Defines the initial weak (remove) thresh-hold cutoff to COSMIC assign signatures to a sample. | 
|  | **refit_denovo_signatures** | Boolean | Default is True. If True, then refit the denovo signatures with nnls. | 
|  | **make_decomposition_plots** | Boolean | Defualt is True. If True, Denovo to Cosmic sigantures decompostion plots will be created as a part the results. | 
| **Others** |  |  |  | 
|  | **get_all_signature_matrices** | Boolean | If True, the Ws and Hs from all the NMF iterations are generated in the output. | 
|  | **export_probabilities** | Boolean | Defualt is True. If False, then doesn't create the probability matrix. | 
    
#### sigProfilerExtractor Example
```R    

library("SigProfilerExtractorR")   
# to get input from vcf files.  
path_to_example_folder_containing_vcf_files = importdata("vcf").   
data = path_to_example_folder_containing_vcf_files # you can put the path to your folder containing the vcf samples.  
sigprofilerextractor("vcf", "example_output", data, minimum_signatures=1, maximum_signatures=10)


#Wait untill the excecution is finished. The process may a couple of hours based on the size of the data.
#Check the current working directory for the "example_output" folder.


# to get input from table format (mutation catalog matrix)
path_to_example_table = importdata("matrix")
data = path_to_example_table # you can put the path to your tab delimited file containing the mutational catalog matrix/table
sigprofilerextractor("matrix", "example_output", data, opportunity_genome="GRCh38", minimum_signatures=1,maximum_signatures=10)
```

#### sigProfilerExtractor Output
To learn about the output, please visit https://osf.io/t6j7u/wiki/home/
  

### <a name="estimate_solution"></a> Estimation of the Optimum Solution
Estimate the optimum solution (rank) among different number of solutions (ranks). 

```R
estimate_solution(base_csvfile, 
          All_solution, 
          genomes, 
          output, 
          title,
          stability, 
          min_stability, 
          combined_stability)
```  
    
| Parameter | Variable Type | Parameter Description |
| --------------------- | -------- |-------- |
| **base_csvfile** | String | Default is "All_solutions_stat.csv". Path to a  csv file that contains the statistics of all solutions. |
| **All_solution** | String | Default is "All_Solutions". Path to a folder that contains the results of all solutions. |
| **genomes** | String | Default is Samples.txt. Path to a tab delimilted file that contains the mutation counts for all genomes given to different mutation types. |
| **output** | String | Default is "results". Path to the output folder. |
| **title** | String | Default is "Selection_Plot". This sets the title of the selection_plot.pdf |
| **stability** | Float | Default is 0.8. The cutoff thresh-hold of the average stability. Solutions with average stabilities below this thresh-hold will not be considered. |
| **min_stability** | Float | Default is 0.2. The cutoff thresh-hold of the minimum stability. Solutions with minimum stabilities below this thresh-hold will not be considered. |
| **combined_stability** | Float | Default is 1.0. The cutoff thresh-hold of the combined stability (sum of average and minimum stability). Solutions with combined stabilities below this thresh-hold will not be considered. |
        
#### Estimation of the Optimum Solution Example
```R 
estimate_solution(base_csvfile="All_solutions_stat.csv", 
          All_solution="All_Solutions", 
          genomes="Samples.txt", 
          output="results", 
          title="Selection_Plot",
          stability=0.8, 
          min_stability=0.2, 
          combined_stability=1.25)
```                

#### Estimation of the Optimum Solution Output
The files below will be generated in the output folder:
| File Name | Description |
| ----- | ----- |
| **All_solutions_stat.csv** | A csv file that contains the statistics of all solutions. |
| **selection_plot.pdf** | A plot that depict the Stability and Mean Sample Cosine Distance for different solutions. |


### <a name="decompose"></a> Decompose

Decomposes the De Novo Signatures into COSMIC Signatures and assigns COSMIC signatures into samples

```R 
decompose(signatures, activities, samples,  output, signature_database=None, nnls_add_penalty=0.05, nnls_remove_penalty=0.01, initial_remove_penalty=0.05, de_novo_fit_penalty=0.02, genome_build="GRCh37", refit_denovo_signatures=True, make_decomposition_plots=True, connected_sigs=True, verbose=False)
``` 

| Parameter | Variable Type | Parameter Description |
| --------------------- | -------- |-------- |
| **signatures** | String | Path to a  tab delimited file that contains the signaure table where the rows are mutation types and colunms are signature IDs. |
| **activities** | String | Path to a tab delimilted file that contains the activity table where the rows are sample IDs and colunms are signature IDs. |
| **samples** | String | Path to a tab delimilted file that contains the activity table where the rows are mutation types and colunms are sample IDs. |
| **output** | String | Path to the output folder. |
| **de_novo_fit_penalty** | Float | Takes any positive float. Default is 0.02. Defines the weak (remove) thresh-hold cutoff to be assigned denovo signatures to a sample. Optional parameter. |
| **nnls_add_penalty** | Float | Takes any positive float. Default is 0.01. Defines the weak (remove) thresh-hold cutoff to be assigned COSMIC signatures to a sample. Optional parameter. |
| **nnls_remove_penalty** | Float | Takes any positive float. Default is 0.01. Defines the weak (remove) thresh-hold cutoff to be assigned COSMIC signatures to a sample. Optional parameter. |
| **initial_remove_penalty** | Float | Takes any positive float. Default is 0.05. Defines the initial weak (remove) thresh-hold cutoff to be COSMIC assigned signatures to a sample. Optional parameter. |
| **genome_build** | String | The genome type. Example: "GRCh37", "GRCh38", "mm9", "mm10". The default value is "GRCh37" |
| **verbose** | Boolean | Prints statements. Default value is False.  |
        
#### Decompose Example
```R
signatures = "path/to/De_Novo_Solution_Signatures.txt"
activities="path/to/De_Novo_Solution_Activities.txt"
samples="path/to/Samples.txt"
output="name or path/to/output"
decompose(signatures, activities, samples, output, genome_build="GRCh37", verbose=False)
```   

#### Decompose Output   
Values:
  The files below will be generated in the output folder:
  - comparison_with_global_ID_signatures.csv
  - Decomposed_Solution_Activities.txt
  - Decomposed_Solution_Samples_stats.txt
  - Decomposed_Solution_Signatures.txt
  - decomposition_logfile.txt
  - dendogram.pdf
  - Mutation_Probabilities.txt
  - Signature_assaignment_logfile.txt
  - Signature_plot[MutatutionContext]_plots_Decomposed_Solution.pdf
        
### <a name="plotActivity"></a> Activity Stacked Bar Plot
Generates a stacked bar plot showing activities in individuals


### GPU support

If CUDA out of memory exceptions occur, it will be necessary to reduce the number of CPU processes used (the `cpu` parameter).

#### For more information, help, and examples, please visit: https://osf.io/t6j7u/wiki/home/

## <a name="copyright"></a> Copyright
This software and its documentation are copyright 2018 as a part of the sigProfiler project. The SigProfilerExtractor framework is free software and is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

## <a name="contact"></a> Contact Information
Please address any queries or bug reports to S M Ashiqul Islam (Mishu) at m0islam@ucsd.edu
