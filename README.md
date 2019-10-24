[![Docs](https://img.shields.io/badge/docs-latest-blue.svg)](https://osf.io/s93d5/wiki/home/) [![License](https://img.shields.io/badge/License-BSD\%202--Clause-orange.svg)](https://opensource.org/licenses/BSD-2-Clause) [![Build Status](https://travis-ci.com/AlexandrovLab/SigProfilerMatrixGeneratorR.svg?branch=master)](https://travis-ci.com/AlexandrovLab/SigProfilerMatrixGeneratorR)

# SigProfilerExtractorR
An R wrapper for running the SigProfilerExtractor framework.

**INTRODUCTION**

The purpose of this document is to provide a guide for using the SigProfilerExtractor framework to extract the De Novo mutational signatures from a set of samples and decompose the De Novo signatures into the COSMIC signatures. An extensive Wiki page detailing the usage of this tool can be found at https://osf.io/t6j7u/wiki/home/. 

For users that prefer working in a Python environment, the tool is written in Python and can be found and installed from: https://github.com/AlexandrovLab/SigProfilerExtractor

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

Information about supported will be found at https://github.com/AlexandrovLab/SigProfilerMatrixGeneratorR

**Extraction Signatures**

Signatures can be extracted from vcf files or tab delimited mutational table using the sigprofilerextractor function.
```
>> help(sigprofilerextractor)
```
This will show the details about the sigprofilerextractor funtion.

Example:
```
>> library("SigProfilerExtractorR")
>> path_to_example_data <- importdata("table")
>> data <- path_to_example_data # here you can provide the path of your own data
>> sigprofilerextractor("table", "example_output", data, minsigs=1, maxsigs=3, replicates=10, cpu=-1)
```

The example file will generated in the working directory


**Decomposition of signature**
SigProfilerExtractorR offer a separate function, decomposition, that decomposes a tab delimited file containing a set of De Novo signatures to the COSMIC signatures and attribute them to a given set of sample.
```
help(decomposition)
```
This will show the details about the decomposition function.
## COPYRIGHT
This software and its documentation are copyright 2018 as a part of the sigProfiler project. The SigProfilerExtractor framework is free software and is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

## CONTACT INFORMATION
Please address any queries or bug reports to S M Ashiqul Islam (Mishu) at m0islam.ucsd.edu

