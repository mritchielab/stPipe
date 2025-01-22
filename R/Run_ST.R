#####################################################
# Pre-processing sST data
#####################################################

#' @name Run_ST
#' @title Pre-processing function for sequencing-based spatial transcriptomics
#'
#' @description Pre-processing from FASTQ files to gene count matrix and extract spatial location information,
#' processing multiple sample is done via parallel computing. NB: for multiple samples, 'species' and 'technology_version' should be the same.
#' This function processes sequencing-based spatial transcriptomics data using various steps, including BAM to FASTQ conversion, trimming, index building, alignment, and barcode detection.
#' For Slideseq technology, the input should be BAM file and for all other technologies the input should be FASTQ file.
#' @examples
#' \dontrun{
#' Run_ST(config = demo_config_path, show.config = TRUE)
#' }
#' @param config Path to the YAML configuration file.
#' @param show.config Logical value indicating whether to print the configuration. Defaults to TRUE.
#' @return None. Outputs are saved to specified directories.
#' @importFrom Rcpp evalCpp
#' @export

Run_ST <- function(config, show.config = TRUE) {

  config <- yaml::read_yaml(config)
  if (show.config) {
    print(config)
  }

  data_dirs <- strsplit(as.character(config$data_directory), ",\\s*")[[1]]  # Multiple data directories
  out_dirs <- strsplit(as.character(config$output_directory), ",\\s*")[[1]]  # Multiple output directories

  if (length(data_dirs) != length(out_dirs)) {
    stop("The number of data directories and output directories must match!")
  }

  # Parallel processing with progress bar
  pbmcapply::pbmclapply(seq_along(data_dirs), function(i) {
    data_dir <- data_dirs[i]
    out_dir <- out_dirs[i]

    # Create the output directory and set the working directory
    dir.create(out_dir, showWarnings = FALSE)

    technology_version <- as.character(config$technology_version)
    species <- as.character(config$species)
    analysis <- as.character(config$analysis_source)
    index_fa <- as.character(config$index_fa)
    gff3 <- as.character(config$index_gff3)
    scpipe_nthreads <- as.numeric(config$scpipe_nthreads)
    max_reads <- as.numeric(config$max_reads)
    min_count <- as.numeric(config$min_count)
    number_of_locations <- as.numeric(config$number_of_locations)

    # FASTQ read structure design
    read_structure <- list(
      bs1 = config$bs1,
      bl1 = config$bl1,
      bs2 = config$bs2,
      bl2 = config$bl2,
      us = config$us,
      ul = config$ul
    )

    # Python script for BAM to FASTQ conversion
    py_reformat <- "
    import sys
    import pysam
    import os

    def process_bam_to_fastq(input_bam, output_fastq):
        with pysam.AlignmentFile(input_bam, 'rb') as bam_in, pysam.AlignmentFile('temp.bam', 'wb', template=bam_in) as bam_out:
            for read in bam_in:
                barcode = read.get_tag('XC')
                umi = read.get_tag('XM')
                new_read_name = f'{barcode}_{umi}#{read.query_name}'
                read.query_name = new_read_name
                bam_out.write(read)
        with pysam.AlignmentFile('temp.bam', 'rb') as bam_file, open(output_fastq, 'w') as fastq_file:
            for read in bam_file:
                read_name = read.query_name
                sequence = read.query_sequence
                quality = pysam.array_to_qualitystring(read.query_qualities)
                fastq_file.write(f'@{read_name}\\n{sequence}\\n+\\n{quality}\\n')

    input_bam = r.input_bam
    output_fastq = './trimmed.fastq'
    process_bam_to_fastq(input_bam, output_fastq)
    "

    # Process based on the technology version
    if (technology_version == "Slideseq") {

      input_bam <- as.character(list.files(data_dir, pattern = ".*\\.bam$", full.names = TRUE))
      print("Using the following BAM file:")
      print(input_bam)
      reticulate::py_run_string(py_reformat)
      combined_fastq <- file.path

    } else if (technology_version == "Stereoseq") {
      print("Using the following FASTQ files:")

      read_1_fq_path <- list.files(data_dir, pattern = "R1", full.names = TRUE)[1]
      read_2_fq_path <- list.files(data_dir, pattern = "R2", full.names = TRUE)[1]
      print(read_1_fq_path)
      print(read_2_fq_path)

      if (is.null(read_1_fq_path) || is.null(read_2_fq_path)) {
        stop("R1 or R2 FASTQ files not found!")
      }

      output_fq_path <- file.path(out_dir, "demultiplexed.fq.gz")
      h5_mapping_path <- as.character(config$h5_mapping_path)
      bin_size <- as.numeric(config$bin_size)
      n_reads <- max_reads
      coord_bc_start <- config$bs2
      coord_bc_len <- config$bl2
      umi_start <- config$us
      umi_len <- config$ul

      # use Rcpp function RunDemultiplex
      RunDemultiplex(
        read_1_fq_path = read_1_fq_path,
        read_2_fq_path = read_2_fq_path,
        h5_mapping_path = h5_mapping_path,
        output_fq_path = output_fq_path,
        n_reads = n_reads,
        coord_bc_start = coord_bc_start,
        coord_bc_len = coord_bc_len,
        umi_start = umi_start,
        umi_len = umi_len,
        bin_size = bin_size
      )

      combined_fastq_gz <- file.path(out_dir, "demultiplexed.fq.gz")
      combined_fastq <- file.path(out_dir, "demultiplexed.fq")

      if (file.exists(combined_fastq_gz)) {
        system(paste("gunzip -c", shQuote(combined_fastq_gz), ">", shQuote(combined_fastq)))
        message("Decompression completed: .fq.gz converted to .fq")
      } else {
        stop("Error: .fq.gz file not found at the specified location")
      }

    } else if (grepl("Visium", technology_version)) {

      print("Using the following FASTQ files:")
      fq_R2_files <- list.files(data_dir, pattern = "R1", full.names = TRUE)
      print(fq_R2_files)
      fq_R1_files <- list.files(data_dir, pattern = "R2", full.names = TRUE)
      print(fq_R1_files)
      if (length(fq_R1_files) != 1 || length(fq_R2_files) != 1) {
        stop("Expected exactly one R1 and one R2 FASTQ file, but found inconsistencies")
      }
      combined_fastq <- file.path(out_dir, "trimmed.fastq")
      scPipe::sc_trim_barcode(
        outfq = combined_fastq,
        r1 = fq_R1_files[1],
        r2 = fq_R2_files[1],
        read_structure = read_structure
      )
    }

    aligned_bam <- file.path(out_dir, "aligned.bam")
    mapped_bam <- file.path(out_dir, "aligned.mapped.bam")

    # for probe-based technologies
    if (grepl("probe", technology_version)) {

      if (technology_version == "Visium_probe_v1" && species == "human") {
        csv_file <- stPipe::Visium_Human_Transcriptome_Probe_Set_v1_0_GRCh38_2020_A
      } else if (technology_version == "Visium_probe_v2" && species == "human") {
        csv_file <- stPipe::Visium_Human_Transcriptome_Probe_Set_v2_0_GRCh38_2020_A
      } else if (technology_version == "Visium_probe_v1" && species == "mouse") {
        csv_file <- stPipe::Visium_Mouse_Transcriptome_Probe_Set_v1_0_mm10_2020_A
      } else {
        stop("Unsupported combination of technology_version and species for probe-based Visium")
      }

      # Generate FASTA and GFF3 files based on the probe set
      generate_fasta <- function(gene_id, probe_seq) {
        return(paste0(">", gene_id, "\n", probe_seq))
      }
      remove_duplicates <- function(csv_file) {
        unique_data <- csv_file[!duplicated(csv_file[[1]]), ]
        return(unique_data)
      }
      generate_gff3 <- function(csv_file, output_file) {
        unique_rows <- remove_duplicates(csv_file)
        file_conn <- file(output_file, "w")
        writeLines("##gff-version 3", file_conn)
        for (i in 1:nrow(unique_rows)) {
          gene_id <- unique_rows[i, 1]
          gff3_line <- paste(gene_id, ".", "exon", "1", "100000000", ".", "+", ".", paste0("gene_id=", gene_id), sep = "\t")
          writeLines(gff3_line, file_conn)
        }
        close(file_conn)
        cat("GFF3 file generated based on the provided probe set:", output_file, "\n")
      }

      fasta_file <- file.path(out_dir, "probe_set.fa")
      gff3_file <- file.path(out_dir, "probe_set.gff3")

      unique_rows <- remove_duplicates(csv_file)
      fasta_content <- ""
      for (i in 1:nrow(unique_rows)) {
        gene_id <- unique_rows[i, 1]
        probe_seq <- unique_rows[i, 2]
        fasta_content <- paste0(fasta_content, generate_fasta(gene_id, probe_seq), "\n")
      }
      fasta_content <- sub("\n$", "", fasta_content)
      write(fasta_content, file = fasta_file)
      cat("FASTA file generated based on the provided probe set:", fasta_file, "\n")

      generate_gff3(csv_file, gff3_file)

      Rsubread::buildindex(
        out_dir,
        reference = fasta_file
      )
      index_path <- out_dir
      exon_anno <- gff3_file
    } else {
      # For non-probe-based technologies (polyA-based), use the provided GFF3 and FASTA files
      Rsubread::buildindex(
        out_dir,
        reference = index_fa
      )
      index_path <- out_dir
      exon_anno <- gff3
    }

    # Alignment via Rsubread
    Rsubread::align(
      index = index_path,
      readfile1 = combined_fastq,
      output_file = file.path(out_dir, "out.aln.bam"),
      nthreads = scpipe_nthreads
    )

    # Mapping and counting
    barcode_anno <- "sample_index.csv"
    scPipe::sc_exon_mapping(
      file.path(out_dir, "out.aln.bam"),
      file.path(out_dir, "out.map.bam"),
      exon_anno,
      nthreads = scpipe_nthreads
    )
    scPipe::sc_detect_bc(
      infq = combined_fastq,
      outcsv = file.path(out_dir, barcode_anno),
      bc_len = read_structure$bl2,
      max_reads = max_reads,
      min_count = min_count,
      number_of_cells = number_of_locations
    )
    scPipe::sc_count_aligned_bam(
      inbam = file.path(out_dir, "out.aln.bam"),
      outbam = mapped_bam,
      annofn = exon_anno,
      bc_len = read_structure$bl2,
      UMI_len = read_structure$ul,
      outdir = out_dir,
      bc_anno = file.path(out_dir, barcode_anno)
    )

  }, mc.cores = length(data_dirs))  # Number of cores based on data directories
}
