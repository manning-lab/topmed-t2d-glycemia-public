library(stringr)
library(readr)
library(magrittr)
library(dplyr)
library(tidyr)

.check_to_split <- function(desired_columns, to_split) {
  if (all(to_split %in% desired_columns)) {
    return(TRUE)
  } else {
    msg <- paste0(
      "to_split columns missing in desired_columns: ",
      paste(to_split[!to_split %in% desired_columns], collapse = ", ")
    )
    stop(msg)
  }
}

get_fields <- function(source) {
  readfile_con <- gzfile(source, "r")
  header <- suppressWarnings(readLines(readfile_con, n = 1))
  close(readfile_con)
  fields <- str_split(header, "\t", simplify = TRUE)
  return(fields)
}

.check_desired <- function(source_file, desired_columns) {
  all_fields <- get_fields(source_file)
  if (all(desired_columns %in% all_fields)) {
    return(TRUE)
  } else {
    msg <- paste0(
      "Desired columns missing in source file: ",
      paste(desired_columns[!desired_columns %in% all_fields], collapse = ", ")
    )
    stop(msg)
  }
}

.has_header <- function(raw_chunk){
  if (any(str_detect(raw_chunk, "^#?chr\\tpos\\tref\\talt"))){
    1
  } else if (any(str_detect(raw_chunk, "^CHROM\\tPOS\\tREF\\tALT"))){
    2
  } else {
    3
  }
}

.get_header <- function(raw_chunk, head_val){
  if (head_val == 1){
    raw_chunk[str_detect(raw_chunk, "^#?chr\\tpos\\tref\\talt")]
  } else {
    raw_chunk[str_detect(raw_chunk, "^CHROM\\tPOS\\tREF\\tALT")]
  }
}

.is_indel <- function(header){
  any(str_detect(header, "indel_focal_length"))
}

.get_fields_from_chunk <- function(raw_chunk) {
  read_tsv(paste0(raw_chunk, collapse = "\n"),
           col_types = cols(.default = col_character()))
}

.get_list <- function(which_list) {
  if (which_list == "old_names") {
    old_names <- c(
      "#chr",
      "MAP20(+-149bp)",
      "MAP35(+-149bp)",
      "GMS_single-end",
      "GMS_paired-end",
      "H1-hESC_fitCons_score", #nolint
      "H1-hESC_fitCons_rankscore", #nolint
      "H1-hESC_confidence_value", #nolint
      "1000G_strict_masked",
      "1000Gp3_AC",
      "1000Gp3_AF",
      "1000Gp3_AFR_AC",
      "1000Gp3_AFR_AF",
      "1000Gp3_EUR_AC",
      "1000Gp3_EUR_AF",
      "1000Gp3_AMR_AC",
      "1000Gp3_AMR_AF",
      "1000Gp3_EAS_AC",
      "1000Gp3_EAS_AF",
      "1000Gp3_SAS_AC",
      "1000Gp3_SAS_AF",
      "fathmm-MKL_non-coding_score",
      "fathmm-MKL_non-coding_rankscore", #nolint
      "fathmm-MKL_non-coding_group",
      "fathmm-MKL_coding_score",
      "fathmm-MKL_coding_rankscore",
      "fathmm-MKL_coding_pred",
      "fathmm-MKL_coding_group",
      "Eigen-raw",
      "Eigen-phred",
      "Eigen-raw_rankscore",
      "Eigen-PC-raw",
      "Eigen-PC-raw_rankscore",
      "CADDraw",
      "CADDphred",
      "MAP20(+-149bp)_unparsed",
      "MAP35(+-149bp)_unparsed",
      "GMS_single-end_unparsed",
      "GMS_paired-end_unparsed",
      "1000G_strict_masked_unparsed",
      "H1-hESC_fitCons_score_unparsed",
      "H1-hESC_fitCons_rankscore_unparsed",
      "fathmm-MKL_non-coding_score_unparsed",
      "fathmm-MKL_non-coding_rankscore_unparsed",
      "fathmm-MKL_coding_score_unparsed",
      "fathmm-MKL_coding_rankscore_unparsed",
      "Eigen-raw_unparsed",
      "Eigen-phred_unparsed",
      "Eigen-raw_rankscore_unparsed",
      "Eigen-PC-raw_unparsed",
      "Eigen-PC-raw_rankscore_unparsed"
    )
    return(old_names)
  } else if (which_list == "new_names") {
    new_names <- c(
      "chr",
      "MAP20_149bp",
      "MAP35_149bp",
      "GMS_single_end",
      "GMS_paired_end",
      "H1_hESC_fitCons_score", #nolint
      "H1_hESC_fitCons_rankscore",
      #nolint
      "H1_hESC_confidence_value",
      #nolint
      "KGP_strict_masked",
      "KGP3_AC",
      "KGP3_AF",
      "KGP3_AFR_AC",
      "KGP3_AFR_AF",
      "KGP3_EUR_AC",
      "KGP3_EUR_AF",
      "KGP3_AMR_AC",
      "KGP3_AMR_AF",
      "KGP3_EAS_AC",
      "KGP3_EAS_AF",
      "KGP3_SAS_AC",
      "KGP3_SAS_AF",
      "fathmm_MKL_non_coding_score",
      "fathmm_MKL_non_coding_rankscore",
      #nolint
      "fathmm_MKL_non_coding_group",
      "fathmm_MKL_coding_score",
      "fathmm_MKL_coding_rankscore",
      "fathmm_MKL_coding_pred",
      "fathmm_MKL_coding_group",
      "Eigen_raw",
      "Eigen_phred",
      "Eigen_raw_rankscore",
      "Eigen_PC_raw",
      "Eigen_PC_raw_rankscore",
      "CADD_raw",
      "CADD_phred",
      "MAP20_149bp_unparsed",
      "MAP35_149bp_unparsed",
      "GMS_single_end_unparsed",
      "GMS_paired_end_unparsed",
      "KGP_strict_masked_unparsed",
      "H1_hESC_fitCons_score_unparsed",
      "H1_hESC_fitCons_rankscore_unparsed",
      "fathmm_MKL_non_coding_score_unparsed",
      "fathmm_MKL_non_coding_rankscore_unparsed",
      "fathmm_MKL_coding_score_unparsed",
      "fathmm_MKL_coding_rankscore_unparsed",
      "Eigen_raw_unparsed",
      "Eigen_phred_unparsed",
      "Eigen_raw_rankscore_unparsed",
      "Eigen_PC_raw_unparsed",
      "Eigen_PC_raw_rankscore_unparsed"
    )
    return(new_names)
  } else if (which_list == "parseable_fields") {
    parseable_fields <- c(
      "splicing_consensus_ada_score",
      "splicing_consensus_rf_score",
      "MAP20",
      "MAP35",
      "MAP20(+-149bp)",
      "MAP35(+-149bp)",
      "GMS_single-end",
      "GMS_paired-end",
      "1000G_strict_masked",
      "RepeatMasker_masked",
      "phyloP46way_primate",
      "phyloP46way_primate_rankscore",
      "phyloP46way_placental",
      "phyloP46way_placental_rankscore",
      "phyloP100way_vertebrate",
      "phyloP100way_vertebrate_rankscore",
      "phastCons46way_primate",
      "phastCons46way_primate_rankscore",
      "phastCons46way_placental",
      "phastCons46way_placental_rankscore",
      "phastCons100way_vertebrate",
      "phastCons100way_vertebrate_rankscore",
      "GERP_NR",
      "GERP_RS",
      "GERP_RS_rankscore",
      "SiPhy_29way_logOdds",
      "SiPhy_29way_logOdds_rankscore",
      "integrated_fitCons_score",
      "integrated_fitCons_rankscore",
      "integrated_confidence_value",
      "GM12878_fitCons_score",
      "GM12878_fitCons_rankscore",
      "GM12878_confidence_value",
      "H1-hESC_fitCons_score",
      "H1-hESC_fitCons_rankscore",
      "H1-hESC_confidence_value",
      "HUVEC_fitCons_score",
      "HUVEC_fitCons_rankscore",
      "HUVEC_confidence_value",
      "GenoCanyon_score",
      "GenoCanyon_rankscore",
      "DANN_score",
      "DANN_rank_score",
      "fathmm-MKL_non-coding_score",
      "fathmm-MKL_non-coding_rankscore",
      "fathmm-MKL_non-coding_group",
      "fathmm-MKL_coding_score",
      "fathmm-MKL_coding_rankscore",
      "fathmm-MKL_coding_pred",
      "fathmm-MKL_coding_group",
      "Eigen-raw",
      "Eigen-phred",
      "Eigen_coding_or_noncoding",
      "Eigen-raw_rankscore",
      "Eigen-PC-raw",
      "Eigen-PC-raw_rankscore",
      "ENCODE_TFBS_score",
      "ENCODE_Dnase_score",
      "ENCODE_Dnase_cells",
      "EnhancerFinder_general_developmental_enhancer",
      "EnhancerFinder_brain_enhancer",
      "EnhancerFinder_heart_enhancer",
      "EnhancerFinder_limb_enhancer",
      "FANTOM5_enhancer_permissive",
      "FANTOM5_enhancer_robust",
      "FANTOM5_CAGE_peak_permissive",
      "FANTOM5_CAGE_peak_robust"
    )
    return(parseable_fields)
  } else if (which_list == "parse_max") {
    parse_max <- c(
      "splicing_consensus_ada_score",
      "splicing_consensus_rf_score",
      "MAP20",
      "MAP35",
      "MAP20(+-149bp)",
      "MAP35(+-149bp)",
      "GMS_single-end",
      "GMS_paired-end",
      "ENCODE_TFBS_score",
      "ENCODE_Dnase_score",
      "ENCODE_Dnase_cells"
    )
    return(parse_max)
  } else   if (which_list == "parse_triples") {
    parse_triples <- list(
      c(
        "Eigen-raw",
        "Eigen-raw_rankscore",
        "Eigen-phred"
      )
    )
    return(parse_triples)
  } else if (which_list == "parse_pairs") {
    parse_pairs <- list(
      c("phyloP46way_primate",
        "phyloP46way_primate_rankscore"),
      c(
        "phyloP46way_placental",
        "phyloP46way_placental_rankscore"
      ),
      c(
        "phyloP100way_vertebrate",
        "phyloP100way_vertebrate_rankscore"
      ),
      c(
        "phastCons46way_primate",
        "phastCons46way_primate_rankscore"
      ),
      c(
        "phastCons46way_placental",
        "phastCons46way_placental_rankscore"
      ),
      c(
        "phastCons100way_vertebrate",
        "phastCons100way_vertebrate_rankscore"
      ),
      c("SiPhy_29way_logOdds",
        "SiPhy_29way_logOdds_rankscore"),
      c("GenoCanyon_score",
        "GenoCanyon_rankscore"),
      c("DANN_score",
        "DANN_rank_score"),
      c("Eigen-PC-raw",
        "Eigen-PC-raw_rankscore"),
      c("GERP_RS",
        "GERP_RS_rankscore"),
      c(
        "integrated_fitCons_score",
        "integrated_fitCons_rankscore"
      ),
      c(
        "GM12878_fitCons_score",
        "GM12878_fitCons_rankscore"
      ),
      c(
        "H1-hESC_fitCons_score",
        "H1-hESC_fitCons_rankscore"
      ),
      c(
        "HUVEC_fitCons_score",
        "HUVEC_fitCons_rankscore"
      ),
      c(
        "fathmm-MKL_non-coding_score",
        "fathmm-MKL_non-coding_rankscore"
      ),
      c(
        "fathmm-MKL_coding_score",
        "fathmm-MKL_coding_rankscore"
      )
    )
    return(parse_pairs)
  } else if (which_list == "parse_string_yes") {
    parse_string_yes <- c(
      "EnhancerFinder_general_developmental_enhancer",
      "EnhancerFinder_brain_enhancer",
      "EnhancerFinder_heart_enhancer",
      "EnhancerFinder_limb_enhancer",
      "FANTOM5_enhancer_permissive",
      "FANTOM5_enhancer_robust",
      "FANTOM5_CAGE_peak_permissive",
      "FANTOM5_CAGE_peak_robust",
      "RepeatMasker_masked"
    )
    return(parse_string_yes)
  } else if (which_list == "parse_string_no") {
    parse_string_no <- c(
      "1000G_strict_masked"
    )
    return(parse_string_no)
  } else if (which_list == "all_fields") {
    # a list of all possible fields in both SNV and indel annotaiton files -
    # new names, including "_unparsed" string
    all_fields <- c(
      "chr",
      "pos",
      "ref",
      "alt",
      "VEP_ensembl_Consequence",
      "VEP_ensembl_Transcript_ID",
      "VEP_ensembl_Gene_Name",
      "VEP_ensembl_Gene_ID",
      "VEP_ensembl_Protein_ID",
      "VEP_ensembl_CCDS",
      "VEP_ensembl_SWISSPROT",
      "VEP_ensembl_Codon_Change",
      "VEP_ensembl_Distance",
      "VEP_ensembl_Amino_Acid_Change",
      "VEP_ensembl_HGVSc",
      "VEP_ensembl_HGVSp",
      "VEP_ensembl_cDNA_position",
      "VEP_ensembl_CDS_position",
      "VEP_ensembl_Protein_position",
      "VEP_ensembl_Exon_or_Intron_Rank",
      "VEP_ensembl_STRAND",
      "VEP_ensembl_CANONICAL",
      "VEP_ensembl_LoF",
      "VEP_ensembl_LoF_filter",
      "VEP_ensembl_LoF_flags",
      "VEP_ensembl_LoF_info",
      "rs_dbSNP147",
      "splicing_consensus_ada_score",
      "splicing_consensus_rf_score",
      "GWAS_catalog_rs",
      "GWAS_catalog_trait",
      "GWAS_catalog_pubmedid",
      "clinvar_rs",
      "clinvar_clnsig",
      "clinvar_trait",
      "clinvar_golden_stars",
      "GTEx_V6_gene",
      "GTEx_V6_tissue",
      "MAP20",
      "MAP35",
      "MAP20_149bp",
      "MAP35_149bp",
      "GMS_single_end",
      "GMS_paired_end",
      "KGP_strict_masked",
      "RepeatMasker_masked",
      "Ancestral_allele",
      "AltaiNeandertal",
      "Denisova",
      "phyloP46way_primate",
      "phyloP46way_primate_rankscore",
      "phyloP46way_placental",
      "phyloP46way_placental_rankscore",
      "phyloP100way_vertebrate",
      "phyloP100way_vertebrate_rankscore",
      "phastCons46way_primate",
      "phastCons46way_primate_rankscore",
      "phastCons46way_placental",
      "phastCons46way_placental_rankscore",
      "phastCons100way_vertebrate",
      "phastCons100way_vertebrate_rankscore",
      "GERP_NR",
      "GERP_RS",
      "GERP_RS_rankscore",
      "SiPhy_29way_logOdds",
      "SiPhy_29way_logOdds_rankscore",
      "integrated_fitCons_score",
      "integrated_fitCons_rankscore",
      "integrated_confidence_value",
      "GM12878_fitCons_score",
      "GM12878_fitCons_rankscore",
      "GM12878_confidence_value",
      "H1_hESC_fitCons_score",
      "H1_hESC_fitCons_rankscore",
      "H1_hESC_confidence_value",
      "HUVEC_fitCons_score",
      "HUVEC_fitCons_rankscore",
      "HUVEC_confidence_value",
      "GenoCanyon_score",
      "GenoCanyon_rankscore",
      "KGP3_AC",
      "KGP3_AF",
      "KGP3_AFR_AC",
      "KGP3_AFR_AF",
      "KGP3_EUR_AC",
      "KGP3_EUR_AF",
      "KGP3_AMR_AC",
      "KGP3_AMR_AF",
      "KGP3_EAS_AC",
      "KGP3_EAS_AF",
      "KGP3_SAS_AC",
      "KGP3_SAS_AF",
      "RegulomeDB_motif",
      "RegulomeDB_score",
      "Motif_breaking",
      "CADD_raw",
      "CADD_phred",
      "CADD_raw_rankscore",
      "DANN_score",
      "DANN_rank_score",
      "fathmm_MKL_non_coding_score",
      "fathmm_MKL_non_coding_rankscore",
      "fathmm_MKL_non_coding_group",
      "fathmm_MKL_coding_score",
      "fathmm_MKL_coding_rankscore",
      "fathmm_MKL_coding_pred",
      "fathmm_MKL_coding_group",
      "Eigen_raw",
      "Eigen_phred",
      "Eigen_coding_or_noncoding",
      "Eigen_raw_rankscore",
      "Eigen_PC_raw",
      "Eigen_PC_raw_rankscore",
      "ENCODE_TFBS",
      "ENCODE_TFBS_score",
      "ENCODE_TFBS_cells",
      "ENCODE_Dnase_score",
      "ENCODE_Dnase_cells",
      "EnhancerFinder_general_developmental_enhancer",
      "EnhancerFinder_brain_enhancer",
      "EnhancerFinder_heart_enhancer",
      "EnhancerFinder_limb_enhancer",
      "FANTOM5_enhancer_permissive",
      "FANTOM5_enhancer_robust",
      "FANTOM5_enhancer_target",
      "FANTOM5_enhancer_expressed_tissue_cell",
      "FANTOM5_enhancer_differentially_expressed_tissue_cell",
      "FANTOM5_CAGE_peak_permissive",
      "FANTOM5_CAGE_peak_robust",
      "Ensembl_Regulatory_Build_Overviews",
      "Ensembl_Regulatory_Build_TFBS",
      "SIFT4G_AAref", # for SNV file, to pivot
      "SIFT4G_AAalt", # for SNV file, to pivot
      "SIFT4G_AApos", # for SNV file, to pivot
      "SIFT4G_score", # for SNV file, to pivot
      "SIFT4G_pred", # for SNV file, to pivot
      "indel_focal_length",
      "focal_snv_number",
      "splicing_consensus_ada_score",
      "splicing_consensus_rf_score",
      "MAP20_149bp_unparsed",
      "MAP35_149bp_unparsed",
      "GMS_single_end_unparsed",
      "GMS_paired_end_unparsed",
      "KGP_strict_masked_unparsed",
      "fathmm_MKL_non_coding_score_unparsed",
      "fathmm_MKL_non_coding_rankscore_unparsed", #nolint
      "fathmm_MKL_non_coding_group_unparsed",
      "fathmm_MKL_coding_score_unparsed",
      "fathmm_MKL_coding_rankscore_unparsed",
      "fathmm_MKL_coding_pred_unparsed",
      "fathmm_MKL_coding_group_unparsed",
      "Eigen_raw_unparsed",
      "Eigen_phred_unparsed",
      "Eigen_raw_rankscore_unparsed",
      "Eigen_PC_raw_unparsed",
      "Eigen_PC_raw_rankscore_unparsed",
      "splicing_consensus_ada_score_unparsed",
      "splicing_consensus_rf_score_unparsed",
      "MAP20_unparsed",
      "MAP35_unparsed",
      "RepeatMasker_masked_unparsed",
      "phyloP46way_primate_unparsed",
      "phyloP46way_primate_rankscore_unparsed",
      "phyloP46way_placental_unparsed",
      "phyloP46way_placental_rankscore_unparsed",
      "phyloP100way_vertebrate_unparsed",
      "phyloP100way_vertebrate_rankscore_unparsed",
      "phastCons46way_primate_unparsed",
      "phastCons46way_primate_rankscore_unparsed",
      "phastCons46way_placental_unparsed",
      "phastCons46way_placental_rankscore_unparsed",
      "phastCons100way_vertebrate_unparsed",
      "phastCons100way_vertebrate_rankscore_unparsed",
      "GERP_RS_unparsed",
      "GERP_RS_rankscore_unparsed",
      "SiPhy_29way_logOdds_unparsed",
      "SiPhy_29way_logOdds_rankscore_unparsed",
      "integrated_fitCons_score_unparsed",
      "integrated_fitCons_rankscore_unparsed",
      "GM12878_fitCons_score_unparsed",
      "GM12878_fitCons_rankscore_unparsed",
      "H1_hESC_fitCons_score_unparsed",
      "H1_hESC_fitCons_rankscore_unparsed",
      "HUVEC_fitCons_score_unparsed",
      "HUVEC_fitCons_rankscore_unparsed",
      "GenoCanyon_score_unparsed",
      "GenoCanyon_rankscore_unparsed",
      "DANN_score_unparsed",
      "DANN_rank_score_unparsed",
      "ENCODE_TFBS_score_unparsed",
      "ENCODE_Dnase_score_unparsed",
      "ENCODE_Dnase_cells_unparsed",
      "EnhancerFinder_general_developmental_enhancer_unparsed",
      "EnhancerFinder_brain_enhancer_unparsed",
      "EnhancerFinder_heart_enhancer_unparsed",
      "EnhancerFinder_limb_enhancer_unparsed",
      "FANTOM5_enhancer_permissive_unparsed",
      "FANTOM5_enhancer_robust_unparsed",
      "FANTOM5_CAGE_peak_permissive_unparsed",
      "FANTOM5_CAGE_peak_robust_unparsed"
    )
    return(all_fields)
  } else {
    msg <- paste0("which_list must be one of 'old_names', 'new_names', ",
                  "'parseable_fields', 'parse_max', 'parse_pairs', ",
                  "'parse_triples', 'parse_group', 'parse_string_yes', ",
                  "'all_fields', or 'parse_string_no'.")
    stop(msg)
  }
}

.check_names <- function(field_names){
  old_names <- .get_list("old_names")
  any(old_names %in% field_names)
}

.parse_snv_chunk <- function(all_fields,
                             desired_columns,
                             to_split,
                             WGSA_version) {
  
  # pick out the desired columns for further operation
  selected_columns <- all_fields %>%
    select(one_of(desired_columns)) %>% # select fields of interest
    mutate(wgsa_version = WGSA_version) # add wgsa version
  
  # pivot the VEP_* fields
  expanded <- selected_columns %>%
    separate_rows(one_of(to_split), sep = "\\|")
  
  # split the VEP_ensembl_Codon_Change_or_Distance field as follows:
  # if number, put in VEP_ensembl_Distance field
  # if string, put in VEP_ensembl_Codon_Change field
  if ("VEP_ensembl_Codon_Change_or_Distance" %in% desired_columns) {
    expanded <- expanded %>%
      extract(
        col = VEP_ensembl_Codon_Change_or_Distance, #nolint
        into = c("VEP_ensembl_Distance", "VEP_ensembl_Codon_Change"),
        regex = "(\\d*)(\\D*)"
      ) %>%
      mutate_at(vars(one_of(
        c("VEP_ensembl_Distance",
          "VEP_ensembl_Codon_Change")
      )),
      funs(str_replace(
        ., pattern = "^$", replacement = "."
      ))) # fill blanks with "."
  }
  
  # if it's an older version annotation file, rename columns from WGSA fields
  # with weird characters to database column names
  if (.check_names(names(expanded))) {
    names(expanded) <- .fix_names(names(expanded))
  }
  
  expanded <- distinct(expanded)
}

.write_to_file <- function(parsed_lines,
                           destination,
                           desired_columns,
                           header_flag,
                           indel_flag) {
  # desired_columns may have VEP_ensembl_Codon_Change_or_Distance. Fix.
  if ("VEP_ensembl_Codon_Change_or_Distance" %in% desired_columns) {
    desired_columns <- c(
      setdiff(desired_columns,
              "VEP_ensembl_Codon_Change_or_Distance"),
      "VEP_ensembl_Distance",
      "VEP_ensembl_Codon_Change"
    )
  }
  
  # add _unparsed columns to desired_columns for parsed indel chunk
  if (indel_flag) {
    desired_columns <- names(parsed_lines)[names(parsed_lines) %in%
                                             .get_list("all_fields")]
  }
  
  # use select statement to write column headers and make sure of column order.
  if (header_flag) {
    parsed_lines %>%
      select(one_of(c(desired_columns, "wgsa_version"))) %>%
      write_tsv(path = destination, append = FALSE)
  } else {
    parsed_lines %>%
      select(one_of(c(desired_columns, "wgsa_version"))) %>%
      write_tsv(path = destination, append = TRUE)
  }
}

parse_to_file <- function(source_file,
                          destination,
                          desired_columns,
                          to_split,
                          WGSA_version = "WGSA065",
                          chunk_size = 10000,
                          verbose = TRUE) {


  # check that desired_columns and to_split are possible
  if (!.check_to_split(desired_columns, to_split)) {
    stop("all to_split fields must be in desired_columns")
  }


  if (!.check_desired(source_file, desired_columns)) {
    stop("not all desired_columns are in source_file")
  }

  # main loop - read file by chunk, process chunk, write chunk
  readfile_con <- gzfile(source_file, "r")
  index <- 0L
  while (TRUE) {
    # read a raw chunk
    raw_chunk <- suppressWarnings(readLines(readfile_con, n = chunk_size))

    # readLines() returns a zero length result at EOF


    if (length(raw_chunk) == 0) {
      break
    }

    # check for header line and read raw chunk to all_fields tibble
    head_val <- .has_header(raw_chunk)
    if (head_val < 3) {
      found_header <- TRUE
      header_flag <- TRUE
      raw_header <- .get_header(raw_chunk, head_val)
      indel_flag <- .is_indel(raw_header)
      all_fields <- .get_fields_from_chunk(raw_chunk)
    } else {
      # header line should have been read in a previous chunk
      if (!found_header){
        stop("Didn't find header line in source_file!")
      }
      header_flag <- FALSE
      modified_chunk <- c(raw_header, raw_chunk)
      all_fields <- .get_fields_from_chunk(modified_chunk)
    }

    # end iteration if all_fields has 0 observations
    # (to avoid dplyr error arising from empty tibble)
    if (dim(all_fields)[1] == 0) {
      index <- index + 1L
      next
    }

    # parse the all_fields tibble
    if (indel_flag) {
      parsed_lines <- .parse_indel_chunk(all_fields,
                                         desired_columns,
                                         to_split,
                                         WGSA_version)
    } else {
      parsed_lines <- .parse_snv_chunk(all_fields,
                                       desired_columns,
                                       to_split,
                                       WGSA_version)
    }

    # write tibble to tsv file
    .write_to_file(parsed_lines,
                   destination,
                   desired_columns,
                   header_flag,
                   indel_flag)

    # ready for the next chunk!
    index <- index + 1L

    # update progress if desired
    if (verbose) {
      msg <- paste0(
        "Chunks completed: ", index,
        "\n Sourcefile lines processed <= ", chunk_size * index,
        "\n Records in current import: ", dim(parsed_lines)[1]
      )
      message(msg)
    }
  }
  close(readfile_con)
}


args <- commandArgs(trailingOnly=T)
anno.file <- args[1]
out.file <- paste(args[2],".tsv",sep="")
target_columns <- unlist(strsplit(args[3],","))
columns_to_split <- unlist(strsplit(args[4],","))

all_fields <- get_fields(anno.file)
print(all_fields)

parse_to_file(source = anno.file, 
              destination = out.file, 
              desired_columns = target_columns, 
              to_split = columns_to_split, 
              chunk_size = 10000)