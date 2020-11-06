
#' sigprofilerextractor
#' @description Extracts De Novo signature from genomes and Decomposes the De Novo Signatures to Cosmic Signatures
#'
#' @param input_type input_type: A string. Type of input. The type of input should be one of the following:\cr
#' ----- "vcf": used for vcf format inputs.\cr
#' ----- "table": used for table format inputs using a tab seperated file.\cr

#' @param output A string. The name of the output folder. The output folder will be generated in the current working directory.\cr
#' @param inputdata A string. Name of the input folder (in case of "vcf" type input) or the input file (in case of "table"  type input).\cr
#' The project file or folder should be inside the current working directory.\cr
#' For the "vcf" type input,the project has to be a folder which will contain the vcf files in vcf format or text formats.\cr
#' The "text"type projects have to be a file.
#' @param refgen A string, optional. The name of the reference genome. The default reference genome is "GRCh37". This parameter is applicable only if the input_type is "vcf".\cr
#' @param genome_build
#' @param minsigs A positive integer, optional. The minimum number of signatures to be extracted. The default value is 1
#' @param maxsigs A positive integer, optional. The maximum number of signatures to be extracted. The default value is 10
#' @param replicates A positive integer, optional. The number of iteration to be performed to extract each number signature. The default value is 8
#' @param mtype A  strings, optional. The items in the list defines the mutational contexts to be considered to extract the signatures. The default value is "96,DINUCS,ID",\cr
#'  where "96" is the SBS96 context, "DINUC" is the DINULEOTIDE context and ID is INDEL context. Each context should be written in a single string separated by "," without any space.
#' @param init A string. Represents the NMF initialization. Options are "random", "nndsvdar" and "nndsvd"
#' @param exome Boolean, optional. Defines if the exomes will be extracted. The default value is "False".
#' @param cpu An integer, optional. The number of processors to be used to extract the signatures. The default value is -1 which will use all available processors.
#'
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
  nmf_replicates=as.integer(nm_replicates)
  min_nmf_iterations=as.integer(min_nmf_iterations)
  max_nmf_iterations=as.integer(max_nmf_iterations)
  nmf_test_conv= as.integer(nmf_test_conv)
  batch_size=as.integer(batch_size)
  cpu=as.integer(cpu)

  print("v1")
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
#' @param signatures A string. Path to a  tab delimited file that contains the signaure table where the rows are mutation types and colunms are signature IDs.
#' @param activities A string. Path to a tab delimilted file that contains the activity table where the rows are sample IDs and colunms are signature IDs.
#' @param samples  A string. Path to a tab delimilted file that contains the activity table where the rows are mutation types and colunms are sample IDs.
#' @param output A string. Path to the output folder.
#' @param mutation_type A string. The context type. Example: "96", "192", "1536", "6144", "INDEL", "DINUC". The default value is "96".
#' @param genome_build A string. The genome type. Example: "GRCh37", "GRCh38", "mm9", "mm10". The default value is "GRCh37"
#' @param verbose Boolean. Prints statements. Default value is False.
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
decomposition <- function(signatures, activities, samples, output, mutation_type="96", genome_build="GRCh37", verbose=F){
  decompose <-reticulate::import("SigProfilerExtractor.decomposition")
  result = decomposition$decompose(signatures, activities, samples, output, mutation_type="96", genome_build="GRCh37", verbose=F)
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
