#' Run_ST: Pre-process from FASTQ files to gene count matrix and extract spatial location information
#'
#' This function processes sequencing-based spatial transcriptomics data using various steps, including BAM to FASTQ conversion, trimming, index building, alignment, and barcode detection.
#' For Slideseq technology, the input should be BAM file and for all other technologies the input should be FASTQ file.
#' @param config Path to the YAML configuration file.
#' @param show.config Logical value indicating whether to print the configuration. Defaults to TRUE.
#' @export
#' @examples mousebrain_visium <- Run_ST(config = "~/Desktop/config_stPipe.yml", show.config = TRUE)

Run_ST <- function(config, show.config = TRUE) {

  # Read configuration from the provided YAML file
  config <- yaml::read_yaml(config)
  # Conditionally print config if show_config is TRUE
  if (show.config) {
    print(config)
  }

  # Extract needed configuration parameters
  data_dir <- as.character(config$data_directory)
  out_dir <- as.character(config$output_directory)
  technology_version <- as.character(config$technology_version)
  species <- as.character(config$species)
  analysis <- as.character(config$analysis_source)
  index_fa <- as.character(config$index_fa)
  gff3 <- as.character(config$index_gff3)
  # Numeric parameters from config
  scpipe_nthreads <- as.numeric(config$scpipe_nthreads)
  max_reads <- as.numeric(config$max_reads)
  min_count <- as.numeric(config$min_count)
  number_of_locations <- as.numeric(config$number_of_locations)
  # Set the working directory to the output directory
  dir.create(out_dir, showWarnings = FALSE)
  setwd(out_dir)

  # FASTQ read_structure design from config file
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

  # Process FASTQ file based on provided technology
  if (technology_version == "Slideseq") {
    input_bam <- as.character(list.files(data_dir, pattern = ".*\\.bam$", full.names = TRUE))
    print("Using the following BAM file:")
    print(input_bam)
    reticulate::py_run_string(py_reformat)
    combined_fastq = file.path(out_dir, "trimmed.fastq")
  } else {
    print("Using the following FASTQ files:")
    fq_R1_files <- list.files(data_dir, pattern = "R1", full.names = TRUE)
    print(fq_R1_files)
    fq_R2_files <- list.files(data_dir, pattern = "R2", full.names = TRUE)
    print(fq_R2_files)
    if (length(fq_R1_files) != 1) {
      print(fq_R1_files)
      stop("Expected exactly one R1 FASTQ file, but found ", length(fq_R1_files))
    }
    if (length(fq_R2_files) != 1) {
      print(fq_R2_files)
      stop("Expected exactly one R2 FASTQ file, but found ", length(fq_R2_files))
    }
    # note scPipe has different R1 and R2 setting with current Visium
    fq_R1 <- fq_R2_files[1]
    fq_R2 <- fq_R1_files[1]
    combined_fastq = file.path(out_dir, "trimmed.fastq")
    # trim the FASTQ files
    scPipe::sc_trim_barcode(
      outfq = combined_fastq,
      r1 = fq_R1,
      r2 = fq_R2,
      read_structure = read_structure
    )
  }

  aligned_bam = file.path(out_dir, "aligned.bam")
  mapped_bam = file.path(out_dir, "aligned.mapped.bam")

  if (grepl("probe", technology_version)) {
    if (technology_version == "Visium_probe_v1" && species == "human") {
      csv_file <- stPipe::Visium_Human_Transcriptome_Probe_Set_v1_0_GRCh38_2020_A
    } else if (technology_version == "Visium_probe_v2" && species == "human") {
      csv_file <- stPipe::Visium_Human_Transcriptome_Probe_Set_v2_0_GRCh38_2020_A
    } else if (technology_version == "Visium_probe_v1" && species == "mouse") {
      csv_file <- stPipe::Visium_Mouse_Transcriptome_Probe_Set_v1_0_mm10_2020_A
    } else {
      stop("Unsupported combination of technology_version and species")
    }

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
      cat("GFF3 file generated based on provided probe set:", output_file, "\n")
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
    cat("FASTA file generated based on provided probe set:", fasta_file, "\n")

    generate_gff3(csv_file, gff3_file)

    Rsubread::buildindex(
      out_dir,
      reference = fasta_file
    )
    index_path <- out_dir
    exon_anno <- gff3_file
  } else {
    # if not Visium probe-based will directly use provided gff3 and fa files for Rsubread
    Rsubread::buildindex(
        out_dir,
        reference = ref_index_fa
      )
    index_path <- out_dir
    exon_anno <- ref_gff3
    }
  Rsubread::align(
    index = index_path,
    readfile1 = combined_fastq,
    output_file = file.path(out_dir, "out.aln.bam"),
    nthreads = scpipe_nthreads
  )
  barcode_anno <- "sample_index.csv"
  scPipe::sc_exon_mapping(file.path(out_dir, "out.aln.bam"),
                  file.path(out_dir, "out.map.bam"),
                  exon_anno,
                  nthreads = scpipe_nthreads
  )
  scPipe::sc_detect_bc(
    infq = combined_fastq,
    outcsv = barcode_anno,
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
    bc_anno = barcode_anno
  )
}
