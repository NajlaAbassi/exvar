#' Convert a FASTQ file to an indexed BAM
#'
#' This function takes in FASTQ files and performs quality control before
#' aligning to a reference genome. It assumes paired-end samples are of the same
#' file name with an underscore (_) and a number to signify different reads of
#' the same sample. Each sample's outputs will be stored in a separate
#' directory.
#'
#' @param file A list of paths to FASTQ files. If no paths are entered, defaults
#' to all fastq files in dir.
#' @param dir Output directory. By default, it uses the current working directory.
#' @param paired A logical indicating whether the data is single or pair-end.
#' @param threads The number of cores to use in the process.
#' @param molecule A character string indicating either DNA or RNA samples.
#' @param species A character string indicating the ID of the species to use
#' (e.g. "1" for hg19,
#' "2" for hg38...). If NULL, an interactive menu will prompt the user.
#' @return A list of \code{BamFile} objects.
#' @export
processfastq <- function(file = list_files_with_exts(dir = dir, exts = "fastq"),
                         dir = getwd(),
                         paired = FALSE,
                         threads = 4L,
                         molecule = "RNA",
                         species = NULL) {

  if (Sys.info()[['sysname']] != "Linux") {
    stop("This function is only available on Linux.")
  }

  if (is.null(species)) {
    message(
        "These are the species currently supported by exvar: \n",
        "[1] Homo sapiens (hg19) \n",
        "[2] Homo sapiens (hg38) \n",
        "[3] Mus musculus \n",
        "[4] Arabidopsis thaliana \n",
        "[5] Drosophila melanogaster \n",
        "[6] Danio rerio \n",
        "[7] Rattus norvegicus \n",
        "[8] Saccharomyces cerevisiae \n",
        "[9] Caenorhabditis elegans \n"
    )
    species <-
      readline("Type the ID of the species that you would like to use as a reference: ")
  }

  wd <- getwd()

  # get the reference genome

  refgen <- .get_genome_info(species)

  message("Obtaining GSNAP parameters...")
  snapParam <- gmapR::GsnapParam(refgen, unique_only = TRUE,
                                 molecule = molecule, nthreads = threads)

  bams <- c()

  if (isTRUE(paired)) {
    inputFastq <- file
    message("Receiving fastq file...")
    inputFastq <- sort(inputFastq)
    fastqn <- length(inputFastq) / 2
    inputfl <- 1:fastqn
    foldernames <- gsub("\\_.*", "", file_path_sans_ext(basename(inputFastq)))

    for (x in inputfl) {
      dir.create(paste0(dir, "/", foldernames[c(2 * x)]))
    }

    for (x in inputfl) {
      message("Quality checking fastq...")
      fastqPath <- file.path(c(inputFastq[c(2 * x - 1)], inputFastq[c(2 *
                                                                        x)]))
      read1 <- file_path_as_absolute(fastqPath[1])
      read2 <- file_path_as_absolute(fastqPath[2])

      setwd(paste0(dir, "/", foldernames[c(2 * x)]))
      json_report <- Rfastp::rfastp(read1 = read1,
                                    read2 = read2,
                                    outputFastq = paste0(foldernames[2 * x],
                                                         '_quality_checked'),
                                    thread = threads,
                                    maxReadLength = 200L)

      QC <- qcSummary(json_report)

      write.csv(QC, "QC_summary.csv")
      gc()

      ## Alignment using the gsnap function
      ## Iterates over all previously selected files
      ## Output is an indexed bam file
      read1 <- paste0(foldernames[2 * x - 1], '_quality_checked_R1.fastq.gz')
      read2 <- paste0(foldernames[2 * x], '_quality_checked_R2.fastq.gz')

      message("Aligning reads...")
      output <- gmapR::gsnap(read1,
                             read2,
                             params = snapParam,
                             output = paste0(getwd(), "/", foldernames[x]))

      message("Creating bam file...")
      bamfl <- as(output, "BamFile")
      append(bams, bamfl)
      gc()
      setwd(wd)
    }
  } else {
    inputFastq <- file
    message("Receiving fastq file...")
    inputFastq <- sort(inputFastq)
    fastqn <- length(inputFastq)
    inputfl <- 1:fastqn
    foldernames <- gsub("\\_.*", "", file_path_sans_ext(basename(inputFastq)))

    for (x in inputfl) {
      dir.create(paste0(dir, "/", foldernames[c(x)]))
    }

    for (x in inputfl) {
      message("Quality checking fastq...")
      fastqPath <- file.path(inputFastq[c(x)])
      read1 <- file_path_as_absolute(fastqPath)
      setwd(paste0(dir, "/", foldernames[c(x)]))
      json_report <- Rfastp::rfastp(read1 = read1,
                                    outputFastq = paste0(foldernames[x],
                                                         '_quality_checked'),
                                    thread = threads,
                                    maxReadLength = 200L)

      QC <- qcSummary(json_report)

      write.csv(QC, "QC_summary.csv")
      gc()

      ## Alignment using the gsnap function
      ## Iterates over all previously selected files
      ## Output is an indexed bam file
      message("Unzipping fastq.gz")

      R.utils::gunzip(paste0(foldernames[x], '_quality_checked_R1.fastq.gz'))
      read <- paste0(foldernames[x], '_quality_checked_R1.fastq')

      message("Aligning reads...")
      output <- gmapR::gsnap(read,
                             input_b = NULL,
                             params = snapParam,
                             output = paste0(getwd(), "/", foldernames[x]))
      message("Creating bam file...")
      bamfl <- as(output, "BamFile")
      append(bams, bamfl)
      gc()
      setwd(wd)
    }
  }
  setwd(wd)
  return(bams)
}
