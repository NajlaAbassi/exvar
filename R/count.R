#' Count genes
#'
#' This function counts reads of gene regions between sample groups.
#' It assumes that sample BAM files are ordered in a directory structure such as
#' "group/sample/" as processfastq() would order it. It outputs a CSV file
#' showing gene counts. Works similarly to geneExpression(), but outputs count
#' data instead of differential expression data.
#'
#' @param dir The parent directory of the sample groups.
#' @param groups Folder names of the sample groups. The default is all folders in dir.
#' @param outputdir Output directory of CSV file. By default, it uses the current
#' working directory.
#' @param threads Number of cores to use.
#' @param paired Indicates whether the samples are from paired-end reads.
#' @param species A character string indicating the ID of the species to use
#' (e.g. "1" for hg19,
#' "2" for hg38...). If NULL, an interactive menu will prompt the user.
#' @export
#' @return A data frame containing annotated gene counts.

count <- function(dir = getwd(),
                  groups = NULL,
                  outputdir = getwd(),
                  threads = 4L,
                  paired = FALSE,
                  species = NULL) {

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

  bpp = MulticoreParam(threads)

  # sets the reference genome that corresponds to the species chosen by the user
  annots <- .get_annotation_resources(species)
  geneExons <- exonsBy(annots$txdb, by = "gene")

  if (is.null(groups)) {
    folders <- list.dirs(path = dir, full.names = TRUE, recursive = FALSE)
  } else {
    groups <- paste0(dir, "/", groups)
    folders <- normalizePath(groups)
  }
  totaldir <- 1:length(folders)
  bamFilesToCount <- c()
  groupNames <- c()
  samples <- c()

  for (x in totaldir) {
    sampledir <- list.dirs(path = folders[c(x)], full.names = TRUE)
    sampledir <- sampledir[-1]
    totsample <- 1:length(sampledir)

    for (i in totsample) {
      print(basename(dirname(sampledir[c(i)])))
      sampletype <- basename(dirname(sampledir[c(i)]))
      groupNames <- append(groupNames, sampletype)
      print(tail(groupNames, n = 1L))
      bam <- list_files_with_exts(dir = sampledir[c(i)], exts = "bam")
      bamFilesToCount <- append(bamFilesToCount, bam)
      names <- basename(dirname(bam))
      samples <- append(samples, names)
    }

  }

  ##This creates a list of bam files to do gene counting with
  ##The folder where the bam files are found (ie sample folder) will inform the
  ##name of the sample
  ##Samples are grouped according to the condition (named by the condition
  ##folders)

  ##A DESeq object is created from the gene counts
  ##A DESeq object is created from the gene counts
  names(bamFilesToCount) <- samples
  myBams <- BamFileList(bamFilesToCount, yieldSize = 10000, asMates = paired)
  geneCounts <- summarizeOverlaps(geneExons, myBams, ignore.strand = TRUE,
                                  BPPARAM = bpp, singleEnd = isFALSE(paired))
  metaData <- data.frame(Group = groupNames,
                         row.names = colnames(geneCounts))
  countMatrix <- assay(geneCounts)
  countDF <- data.frame(countMatrix)
  annotatedCount <- countDF

  orgdb <- annots$orgdb
  eToSym <- AnnotationDbi::select(get(orgdb, envir = asNamespace(orgdb)),
                                  keys = rownames(countDF),
                                  keytype = "ENTREZID",
                                  columns= c("SYMBOL", "ENSEMBL"))

  annotatedCount <- merge(eToSym,countDF,
                          by.x=1,
                          by.y=0,
                          all.x=FALSE,
                          all.y=TRUE)
  annotatedCount <- annotatedCount[order(as.integer(annotatedCount$ENTREZID)),]

  write.csv(annotatedCount, "Count_Data.csv")
  return(annotatedCount)
}
