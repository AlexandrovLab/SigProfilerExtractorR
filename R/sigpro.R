
#' SigprofilerExtractor Function
#'
#' @param input_type The type of input
#' @param output The name/path of the result folder
#' @param input Input file or folder
#' @param refgen Reference genome
#' @param genome_build Genome opportunity
#' @param minsigs Minimum number of signatures to be extracted
#' @param maxsigs Maximum number of signatures to be extracted
#' @param replicates Number of replication of signature extraction
#' @param mtype Context type
#' @param init NMF initialization
#' @param exome Exome true of False
#' @param cpu Number of CPUs to be used during multiprocessing
#'
#' @return a folder with results
#' @export sigProfilerExtractor
#'
#' @examples
sigProfilerExtractor <- function(input_type, output, input, refgen = 'GRCh37', genome_build="GRCh37",minsigs = 1, maxsigs = 3,  replicates=5, mtype = c('96,DINUCS,ID'), init="random",exome=F, cpu=-1) {

  sys <- reticulate::import("sys")
  sigpro <- reticulate::import("sigproextractor.sigpro")


  minsigs=as.integer(minsigs)
  maxsigs = as.integer(maxsigs)
  replicates=as.integer(replicates)
  exome=FALSE
  cpu = as.integer(cpu)
  sigpro$sigProfilerExtractor(input_type, output, input, refgen = 'GRCh38',genome_build="GRCh38",startProcess=minsigs, endProcess = maxsigs, totalIterations=replicates, exome=exome, init=init)
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
  sigpro <- reticulate::import("sigproextractor.sigpro")
  project = sigpro$importdata("vcf")
  return(project)
}
