.genome_lookup <- list(
  "1" = list(pkg = "BSgenome.Hsapiens.UCSC.hg19", obj = "BSgenome.Hsapiens.UCSC.hg19", dir = "hg19"),
  "2" = list(pkg = "BSgenome.Hsapiens.UCSC.hg38", obj = "BSgenome.Hsapiens.UCSC.hg38", dir = "hg38"),
  "3" = list(pkg = "BSgenome.Mmusculus.UCSC.mm10", obj = "BSgenome.Mmusculus.UCSC.mm10", dir = "mm10"),
  "4" = list(pkg = "BSgenome.Athaliana.TAIR.TAIR9", obj = "BSgenome.Athaliana.TAIR.TAIR9", dir = "TAIR9"),
  "5" = list(pkg = "BSgenome.Dmelanogaster.UCSC.dm6", obj = "BSgenome.Dmelanogaster.UCSC.dm6", dir = "dm6"),
  "6" = list(pkg = "BSgenome.Drerio.UCSC.danRer11", obj = "BSgenome.Drerio.UCSC.danRer11", dir = "danRer11"),
  "7" = list(pkg = "BSgenome.Rnorvegicus.UCSC.rn5", obj = "BSgenome.Rnorvegicus.UCSC.rn5", dir = "rn5"),
  "8" = list(pkg = "BSgenome.Scerevisiae.UCSC.sacCer3", obj = "BSgenome.Scerevisiae.UCSC.sacCer3", dir = "sacCer3"),
  "9" = list(pkg = "BSgenome.Celegans.UCSC.ce11", obj = "BSgenome.Celegans.UCSC.ce11", dir = "ce11")
)

#' List supported genomes in exvar
#'
#' This function returns a data frame with available species IDs and corresponding
#' genome packages
#'
#' @return A data.frame with columns: ID, Package and Directory
#' @export
list_supported_genomes <- function() {
  data.frame(
    ID = names(.genome_lookup),
    Package = sapply(.genome_lookup, function(x) x$pkg),
    Directory = sapply(.genome_lookup, function(x) x$dir),
    stringsAsFactors = FALSE
  )
}


#' @noRd
#' @keywords internal
.get_genome_info <- function(species) {
  info <- .genome_lookup[[species]]
  if (is.null(info)) {
    stop("Invalid species identifier.",
    "Use list_supported_genomes() to view valid options.")
  }

  if (!requireNamespace(info$pkg, quietly = TRUE)) {
    stop(glue::glue(
      "Genome package '{info$pkg}' is not installed.\n",
      "Install it with:\nBiocManager::install('{info$pkg}')"
    ))
  }

  organism <- getExportedValue(info$pkg, info$obj)
  ref_path <- file.path(find.package("exvar"), info$dir)

  if (!dir.exists(ref_path)) {
    message("Reference genome not found. Creating reference. This might take a while...")
    refgen <- gmapR::GmapGenome(organism,
                                directory = find.package("exvar"),
                                create = TRUE)
  } else {
    refgen <- gmapR::GmapGenome(organism, directory = find.package("exvar"))
  }

  return(refgen)
}


#' @noRd
#' @keywords internal
.fetch_counting_params <- function(species) {
  param_lookup <- list(
    "1" = list(
      genome = "BSgenome.Hsapiens.UCSC.hg19",
      txdb = "TxDb.Hsapiens.UCSC.hg19.knownGene",
      orgdb = "org.Hs.eg.db",
      keytype = "ENTREZID"
    ),
    "2" = list(
      genome = "BSgenome.Hsapiens.UCSC.hg38",
      txdb = "TxDb.Hsapiens.UCSC.hg38.knownGene",
      orgdb = "org.Hs.eg.db",
      keytype = "ENTREZID"
    ),
    "3" = list(
      genome = "BSgenome.Mmusculus.UCSC.mm10",
      txdb = "TxDb.Mmusculus.UCSC.mm10.knownGene",
      orgdb = "org.Mm.eg.db",
      keytype = "ENTREZID"
    ),
    "4" = list(
      genome = "BSgenome.Athaliana.TAIR.TAIR9",
      txdb = "TxDb.Athaliana.BioMart.plantsmart28",
      orgdb = "org.At.tair.db",
      keytype = "ENTREZID"
    ),
    "5" = list(
      genome = "BSgenome.Dmelanogaster.UCSC.dm6",
      txdb = "TxDb.Dmelanogaster.UCSC.dm6.ensGene",
      orgdb = "org.Dm.eg.db",
      keytype = "ENTREZID"
    ),
    "6" = list(
      genome = "BSgenome.Drerio.UCSC.danRer11",
      txdb = "TxDb.Drerio.UCSC.danRer11.refGene",
      orgdb = "org.Dr.eg.db",
      keytype = "ENTREZID"
    ),
    "7" = list(
      genome = "BSgenome.Rnorvegicus.UCSC.rn5",
      txdb = "TxDb.Rnorvegicus.UCSC.rn5.refGene",
      orgdb = "org.Rn.eg.db",
      keytype = "ENTREZID"
    ),
    "8" = list(
      genome = "BSgenome.Scerevisiae.UCSC.sacCer3",
      txdb = "TxDb.Scerevisiae.UCSC.sacCer3.sgdGene",
      orgdb = "org.Sc.sgd.db",
      keytype = "ENTREZID"
    ),
    "9" = list(
      genome = "BSgenome.Celegans.UCSC.ce11",
      txdb = "TxDb.Celegans.UCSC.ce11.refGene",
      orgdb = "org.Ce.eg.db",
      keytype = "ENTREZID"
    )
  )

  info <- param_lookup[[species]]
  if (is.null(info)) {
    stop("Invalid species ID. Use list_supported_genomes() to view valid options.")
  }

  pkgs <- c(info$genome, info$txdb, info$orgdb)


  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(glue::glue("Missing required package: {pkg}.
                      Install it using BiocManager::install('{pkg}')"))
      }
  }

  counting_params <- list(
    genome = getExportedValue(info$genome, info$genome),
    txdb = getExportedValue(info$txdb, info$txdb),
    orgdb = info$orgdb,
    keytype = info$keytype
  )

  return(counting_params)
}

#' @noRd
#' @keywords internal
.get_annotation_resources <- function(species) {
  info <- .fetch_counting_params[[species]]
  if (is.null(info)) {
    stop("Invalid species ID.")
  }

  pkgs <- c(info$genome, info$txdb, info$org)
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      stop(glue::glue("Missing required package: {p}.
                      Install it using BiocManager::install('{p}')"))
    }
  }

  list(
    genome = getExportedValue(info$genome, info$genome),
    txdb = getExportedValue(info$txdb, info$txdb),
    orgdb = info$org
  )
}


.fetch_variant_params <- function(species) {
  param_lookup <- list(
    "1" = list(
      genome = "BSgenome.Hsapiens.UCSC.hg19",
      dbsnp = "SNPlocs.Hsapiens.dbSNP144.GRCh37",
      all_snp = "XtraSNPlocs.Hsapiens.dbSNP144.GRCh37",
      bed = "TxDb.Hsapiens.UCSC.hg19.knownGene",
      refgen_name = "hg19"
    ),
    "2" = list(
      genome = "BSgenome.Hsapiens.UCSC.hg38",
      dbsnp = "SNPlocs.Hsapiens.dbSNP144.GRCh38",
      all_snp = "XtraSNPlocs.Hsapiens.dbSNP144.GRCh38",
      bed = "TxDb.Hsapiens.UCSC.hg38.knownGene",
      refgen_name = "hg38"
    ),
    "3" = list(
      genome = "BSgenome.Mmusculus.UCSC.mm10",
      dbsnp = NULL,
      all_snp = NULL,
      bed = NULL,
      refgen_name = "mm10"
    ),
    "4" = list(
      genome = "BSgenome.Athaliana.TAIR.TAIR9",
      dbsnp = NULL,
      all_snp = NULL,
      bed = NULL,
      refgen_name = "TAIR9"
    ),
    "5" = list(
      genome = "BSgenome.Dmelanogaster.UCSC.dm6",
      dbsnp = NULL,
      all_snp = NULL,
      bed = NULL,
      refgen_name = "dm6"
    ),
    "6" = list(
      genome = "BSgenome.Drerio.UCSC.danRer11",
      dbsnp = NULL,
      all_snp = NULL,
      bed = NULL,
      refgen_name = "danRer11"
    ),
    "7" = list(
      genome = "BSgenome.Rnorvegicus.UCSC.rn5",
      dbsnp = NULL,
      all_snp = NULL,
      bed = NULL,
      refgen_name = "rn5"
    ),
    "8" = list(
      genome = "BSgenome.Scerevisiae.UCSC.sacCer3",
      dbsnp = NULL,
      all_snp = NULL,
      bed = NULL,
      refgen_name = "sacCer3"
    ),
    "9" = list(
      genome = "BSgenome.Celegans.UCSC.ce11",
      dbsnp = NULL,
      all_snp = NULL,
      bed = NULL,
      refgen_name = "ce11"
    )
  )

  info <- param_lookup[[species]]
  if (is.null(info)) {
    stop("Invalid species ID. Use list_supported_genomes() to view valid options.")
  }

  # load only non NULL pkgs
  pkgs <- Filter(Negate(is.null), c(info$genome, info$dbsnp, info$all_snp, info$bed))



  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(glue::glue("Missing required package: {pkg}.
                      Install it using BiocManager::install('{pkg}')"))
    }
  }

  variant_params <- list(
    genome = getExportedValue(info$genome, info$genome),
    dbsnp = if (!is.null(info$dbsnp)) getExportedValue(info$dbsnp, info$dbsnp) else NULL,
    all_snp = if (!is.null(info$all_snp)) getExportedValue(info$all_snp, info$all_snp) else NULL,
    bed = if (!is.null(info$bed)) getExportedValue(info$bed, info$bed) else NULL,
    refgen_name = info$refgen_name
  )

  return(variant_params)
}



