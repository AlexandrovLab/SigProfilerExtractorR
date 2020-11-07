
#' sigprofilerextractor
#' @description Extracts De Novo signature from genomes and Decomposes the De Novo Signatures to Cosmic Signatures
#'
#' @param input_type input_type: A string. Type of input. The type of input should be one of the following:\cr
#' ----- "vcf": used for vcf format inputs.\cr
#' ----- "matrix": used for matrix/table format inputs using a tab seperated file.\cr

#' @param output A string. The name of the output folder. The output folder will be generated in the current working directory.\cr
#' @param inputdata A string. Name/path of the input folder (in case of "vcf" type input) or the input file (in case of "matrix"  type input).\cr
#' @param reference_genome: A string, optional. The name of the reference genome. The default reference genome is "GRCh37". This parameter is applicable only if the input_type is "vcf".\cr
#' @param opportunity_genome: The build or version of the reference signatures for the reference genome. The default opportunity genome is GRCh37. If the input_type is "vcf", the genome_build automatically matches the input reference genome value.\cr
#' @param context_type: A list of strings, optional. The items in the list defines the mutational contexts to be considered to extract the signatures. The default value is "SBS96,DBS78,ID83".\cr
#' @param exome: Boolean, optional. Defines if the exomes will be extracted. The default value is "False".\cr
#' @param minimum_signature: A positive integer, optional. The minimum number of signatures to be extracted. The default value is 1.\cr
#' @param maximum_signatures: A positive integer, optional. The maximum number of signatures to be extracted. The default value is 25.\cr
#' @param nmf_replicates: A positive integer, optional. The number of iteration to be performed to extract each number signature. The default value is 500. \cr
#' @param  resample: Boolean, optional. Default is True. If True, add poisson noise to samples by resampling.\cr
#' @param seeds: Boolean. Default is "random". If random, then the seeds for resampling will be random for different analysis.\cr
#' If not random, then seeds will be obtained from a given path of a .txt file that contains a list of seed. \cr
#' @param matrix_normalization: A string. Method of normalizing the genome matrix before it is analyzed by NMF. Default is "log2". Other options are "gmm", "100X" or "no_normalization". \cr
#' @param nmf_init: A String. The initialization algorithm for W and H matrix of NMF. Options are 'random', 'nndsvd', 'nndsvda', 'nndsvdar' and 'nndsvd_min'\cr
#' Default is 'nndsvd_min'.
#' @param precision: A string. Values should be single or double. Default is single.
#' @param  min_nmf_iterations: An integer. Value defines the minimum number of iterations to be completed before NMF converges. Default is 10000.
#' @param  max_nmf_iterations: An integer. Value defines the maximum number of iterations to be completed before NMF converges. Default is 1000000
#' @param nmf_test_conv: An integer. Value definer the number number of iterations to done between checking next convergence.
#' @param nmf_tolerance: A float. Value defines the tolerance to achieve to converge.
#' @param cpu: An integer, optional. The number of processors to be used to extract the signatures. The default value is -1 which will use all available processors.
#' @param gpu:Boolean, optional. Defines if the GPU resource will used if available. Default is False. If True, the GPU resource will be used in the computation.
#' @param batch_size: An integer. Will be effective only if the GPU is used. Defines the number of NMF replicates to be performed by each CPU during the parallel processing. Default is 1.
#' @param stability: A float. Default is 0.8. The cutoff thresh-hold of the average stability. Solutions with average stabilities below this thresh-hold will not be considered.
#' @param min_stability: A float. Default is 0.2. The cutoff thresh-hold of the minimum stability. Solutions with minimum stabilities below this thresh-hold will not be considered.
#' @param combined_stability: A float. Default is 1.0. The cutoff thresh-hold of the combined stability (sum of average and minimum stability). Solutions with combined stabilities below this thresh-hold will not be considered.
#' @param de_novo_fit_penalty: Float, optional. Takes any positive float. Default is 0.02. Defines the weak (remove) thresh-hold cutoff to be assigned denovo signatures to a sample.
#' @param nnls_add_penalty: Float, optional. Takes any positive float. Default is 0.05. Defines the strong (add) thresh-hold cutoff to be assigned COSMIC signatures to a sample.
#' @param nnls_remove_penalty: Float, optional. Takes any positive float. Default is 0.01. Defines the weak (remove) thresh-hold cutoff to be assigned COSMIC signatures to a sample.
#' @param initial_remove_penalty: Float, optional. Takes any positive float. Default is 0.05. Defines the initial weak (remove) thresh-hold cutoff to be COSMIC assigned signatures to a sample.
#' @param refit_denovo_signatures: Boolean, optional. Default is False. If True, then refit the denovo signatures with nnls.
#' @param make_decomposition_plots: Boolean, optional. Defualt is True. If True, Denovo to Cosmic sigantures decompostion plots will be created as a part the results.
#' @param get_all_signature_matrices: A Boolean. If true, the Ws and Hs from all the NMF iterations are generated in the output.
#' @param export_probabilities: A Boolean. Defualt is True. If False, then doesn't create the probability matrix.
#' @return a folder with results
#' @export sigprofilerextractor
#'
#' @examples
sigprofilerextractor <- function(input_type,
                                 output,
                                 input_data,
                                 reference_genome="GRCh37",
                                 opportunity_genome = "GRCh37",
                                 context_type = "default",
                                 exome = F,
                                 minimum_signatures=1,
                                 maximum_signatures=25,
                                 nmf_replicates=500,
                                 resample = T,
                                 batch_size=1,
                                 cpu=-1,
                                 gpu=F,
                                 nmf_init="random",
                                 precision= "single",
                                 matrix_normalization= "gmm",
                                 seeds= "random",
                                 min_nmf_iterations= 10000,
                                 max_nmf_iterations=1000000,
                                 nmf_test_conv= 10000,
                                 nmf_tolerance= 1e-15,
                                 nnls_add_penalty=0.05,
                                 nnls_remove_penalty=0.01,
                                 de_novo_fit_penalty=0.02,
                                 initial_remove_penalty=0.05,
                                 refit_denovo_signatures=T,
                                 clustering_distance="cosine",
                                 export_probabilities=T,
                                 make_decomposition_plots=T,
                                 stability=0.8,
                                 min_stability=0.2,
                                 combined_stability=1.0,
                                 get_all_signature_matrices= F) {

  sys <- reticulate::import("sys")
  sigpro <- reticulate::import("SigProfilerExtractor.sigpro")


  #minsigs=as.integer(minsigs)
  #maxsigs = as.integer(maxsigs)
  #replicates=as.integer(replicates)
  #exome=F
  #cpu = as.integer(cpu)

  minimum_signatures=as.integer(minimum_signatures)
  maximum_signatures=as.integer(maximum_signatures)
  nmf_replicates=as.integer(nmf_replicates)
  min_nmf_iterations=as.integer(min_nmf_iterations)
  max_nmf_iterations=as.integer(max_nmf_iterations)
  nmf_tolerance=as.numeric(nmf_tolerance)
  nmf_test_conv= as.integer(nmf_test_conv)
  nnls_add_penalty=as.numeric(nnls_add_penalty)
  nnls_remove_penalty=as.numeric(nnls_remove_penalty)
  de_novo_fit_penalty=as.numeric(de_novo_fit_penalty)
  initial_remove_penalty=as.numeric(initial_remove_penalty)
  stability=as.numeric(stability)
  min_stability=as.numeric(min_stability)
  combined_stability=as.numeric(combined_stability)
  batch_size=as.integer(batch_size)
  cpu=as.integer(cpu)


  sigpro$sigProfilerExtractor(input_type,
                              output,
                              input_data,
                              reference_genome=reference_genome,
                              opportunity_genome = opportunity_genome,
                              context_type = context_type,
                              exome = exome,
                              minimum_signatures=minimum_signatures,
                              maximum_signatures=maximum_signatures,
                              nmf_replicates=nmf_replicates,
                              resample = resample,
                              batch_size=batch_size,
                              cpu=cpu,
                              gpu=gpu,
                              nmf_init=nmf_init,
                              precision= precision,
                              matrix_normalization= matrix_normalization,
                              seeds= seeds,
                              min_nmf_iterations= min_nmf_iterations,
                              max_nmf_iterations=max_nmf_iterations,
                              nmf_test_conv= nmf_test_conv,
                              nmf_tolerance= nmf_tolerance,
                              nnls_add_penalty=nnls_add_penalty,
                              nnls_remove_penalty=nnls_remove_penalty,
                              de_novo_fit_penalty=de_novo_fit_penalty,
                              initial_remove_penalty=initial_remove_penalty,
                              refit_denovo_signatures=refit_denovo_signatures,
                              clustering_distance=clustering_distance,
                              export_probabilities=export_probabilities,
                              make_decomposition_plots=make_decomposition_plots,
                              stability=stability,
                              min_stability=min_stability,
                              combined_stability=combined_stability,
                              get_all_signature_matrices= get_all_signature_matrices)
sys$stdout$flush()
}



#' Title
#'
#' @param datatype The type of example data. There are two types: 1."table", 2."vcf"
#'
#' @return The path of the explample data
#' @export importdata
#'
#' @examples
importdata <- function(datatype){
  sigpro <- reticulate::import("SigProfilerExtractor.sigpro")
  project = sigpro$importdata(datatype)
  return(project)
}



#' decomposition
#' @description Decomposes the De Novo Signatures to Cosmic Signatures
#' @param signatures: A string. Path to a  tab delimited file that contains the signaure table where the rows are mutation types and colunms are signature IDs.
#' @param activities: A string. Path to a tab delimilted file that contains the activity table where the rows are sample IDs and colunms are signature IDs.
#' @param samples: A string. Path to a tab delimilted file that contains the activity table where the rows are mutation types and colunms are sample IDs.
#' @param output: A string. Path to the output folder.
#' @param genome_build = A string. The genome type. Example: "GRCh37", "GRCh38", "mm9", "mm10". The default value is "GRCh37"
#' @param verbose = Boolean. Prints statements. Default value is False.
#'
#' @return A folder with decomposition results The files below will be generated in the output folder: \cr
#'    \cr
#'    1.comparison_with_global_ID_signatures.csv \cr
#     2.Decomposed_Solution_Activities.txt \cr
#'    3.Decomposed_Solution_Samples_stats.txt \cr
#'    4.Decomposed_Solution_Signatures.txt \cr
#'    5.decomposition_logfile.txt \cr
#'    6.Mutation_Probabilities.txt \cr
#'    7.Signature_assaignment_logfile.txt \cr
#'    8.Signature_plot[MutatutionContext]_plots_Decomposed_Solution.pdf
#' @export decomposition
#'
#' @examples
decomposition <- function(signatures,
                          activities,
                          samples,
                          output,
                          signature_database=None,
                          nnls_add_penalty=0.05,
                          nnls_remove_penalty=0.01,
                          initial_remove_penalty=0.05,
                          de_novo_fit_penalty=0.02,
                          genome_build="GRCh37",
                          refit_denovo_signatures=T,
                          make_decomposition_plots=T,
                          connected_sigs=T,
                          verbose=F){

  nnls_add_penalty=as.numeric(nnls_add_penalty)
  nnls_remove_penalty=as.numeric(nnls_remove_penalty)
  de_novo_fit_penalty=as.numeric(de_novo_fit_penalty)
  initial_remove_penalty=as.numeric(initial_remove_penalty)

  decomposition <-reticulate::import("SigProfilerExtractor.decomposition")
  result = decomposition$decompose(signatures,
                                   activities,
                                   samples,
                                   output,
                                   signature_database=signature_database,
                                   nnls_add_penalty=nnls_add_penalty,
                                   nnls_remove_penalty=nnls_remove_penalty,
                                   initial_remove_penalty=initial_remove_penalty,
                                   de_novo_fit_penalty=de_novo_fit_penalty,
                                   genome_build="GRCh37",
                                   refit_denovo_signatures=refit_denovo_signatures,
                                   make_decomposition_plots=make_decomposition_plots,
                                   connected_sigs=connected_sigs,
                                   verbose=verbose)
  return(result)
}


#' Model_Selection
#' @description Decomposes the De Novo Signatures to Cosmic Signatures
#' @param signatures: A string. Path to a  tab delimited file that contains the signaure table where the rows are mutation types and colunms are signature IDs.
#' @param activities: A string. Path to a tab delimilted file that contains the activity table where the rows are sample IDs and colunms are signature IDs.
#' @param samples: A string. Path to a tab delimilted file that contains the activity table where the rows are mutation types and colunms are sample IDs.
#' @param output: A string. Path to the output folder.
#' @param genome_build = A string. The genome type. Example: "GRCh37", "GRCh38", "mm9", "mm10". The default value is "GRCh37"
#' @param verbose = Boolean. Prints statements. Default value is False.
#'
#' @return A folder with decomposition results The files below will be generated in the output folder: \cr
#'    \cr
#'    1.comparison_with_global_ID_signatures.csv \cr
#     2.Decomposed_Solution_Activities.txt \cr
#'    3.Decomposed_Solution_Samples_stats.txt \cr
#'    4.Decomposed_Solution_Signatures.txt \cr
#'    5.decomposition_logfile.txt \cr
#'    6.Mutation_Probabilities.txt \cr
#'    7.Signature_assaignment_logfile.txt \cr
#'    8.Signature_plot[MutatutionContext]_plots_Decomposed_Solution.pdf
#' @export decomposition
#'
#' @examples


estimate_solution<- function(base_csvfile="All_solutions_stat.csv",
                  All_solution="All_Solutions",
                  genomes="Samples.txt",
                  output="results",
                  title="Selection_Plot",
                  stability=0.8,
                  min_stability=0.2,
                  combined_stability=1.0,
                  statistics=T,
                  select=None){

  stability=as.numeric(stability)
  min_stability=as.numeric(min_stability)
  combined_stability=as.numeric(combined_stability)

  model_selection <-reticulate::import("SigProfilerExtractor.estimate_best_solution")
  results=model_selection$estimate_solution(base_csvfile=base_csvfile,
                                            All_solution=All_solution,
                                            genomes=genomes,
                                            output=output,
                                            title=title,
                                            stability=stability,
                                            min_stability=min_stability,
                                            combined_stability=combined_stability,
                                            statistics=statistic,
                                            select=select)



  return(result)

}







#' Title
#'
#' @param genome
#' @param custom
#' @param rsync
#' @param bash
#' @param ftp
#'
#' @return
#' @export install
#'
#' @examples
install <- function(genome, custom=F, rsync=F, bash=T, ftp=T){
  os <- reticulate::import("os")
  sys <- reticulate::import("sys")
  genInstall <- reticulate::import("SigProfilerMatrixGenerator.install")
  genInstall$install(genome, custom, rsync,bash,ftp)
  sys$stdout$flush()
}
