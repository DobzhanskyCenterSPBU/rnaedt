#' @keywords internal
is_mismatch <- function(reference, A, C, G, T) {
    switch(
        reference,
        "A" = (C > 0) + (G > 0) + (T > 0) > 0,
        "T" = (A > 0) + (C > 0) + (G > 0) > 0,
        "G" = (A > 0) + (C > 0) + (T > 0) > 0,
        "C" = (A > 0) + (G > 0) + (T > 0) > 0,
        FALSE)
}

#' @keywords internal
is_mismatch_vec <- Vectorize(is_mismatch)

#' Tries to define correct strand of the position via counting reads aligned to both strands
#'
#' @author Irina Shchyukina
#'
#' @keywords internal
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom Rsamtools scanBamFlag ScanBamParam scanBam
get_strand <- function(pos_first, current_chromosome_length, current_chromosome, bam_file, single_end) {
    eps <- 250
    left_range <- max(0, pos_first - eps)
    right_range <- min(current_chromosome_length, pos_first + eps)
    current_range <- GRanges(current_chromosome, IRanges(left_range, right_range))
    if (single_end) {
        bam_flag <- scanBamFlag(
            isSecondaryAlignment = FALSE,
            isNotPassingQualityControls = FALSE,
            isDuplicate = FALSE)
        bam_params <- ScanBamParam(which = current_range, flag = bam_flag, what = c("strand"), mapqFilter = 30)
        first_area_info <- scanBam(bam_file, param = bam_params)[[1]][[1]]
        plus_read_count <- table(first_area_info)['+']
        minus_read_count <- table(first_area_info)['-']
    } else {
        bam_flag <- scanBamFlag(
            isPaired = TRUE,
            isProperPair = TRUE,
            isSecondaryAlignment = FALSE,
            isNotPassingQualityControls = FALSE,
            isDuplicate = FALSE,
            hasUnmappedMate = FALSE,
            isFirstMateRead = TRUE)
        bam_params <- ScanBamParam(which = current_range, flag = bam_flag, what = c("strand"), mapqFilter = 30)
        first_area_info <- scanBam(bam_file, param = bam_params)[[1]][[1]]
        bam_flag <- scanBamFlag(
            isPaired = TRUE,
            isProperPair = TRUE,
            isSecondaryAlignment = FALSE,
            isNotPassingQualityControls = FALSE,
            isDuplicate = FALSE,
            hasUnmappedMate = FALSE,
            isSecondMateRead = TRUE)
        bam_params <- ScanBamParam(which = current_range, flag = bam_flag, what = c("strand"), mapqFilter = 30)
        second_area_info <- scanBam(bam_file, param = bam_params)[[1]][[1]]
        plus_read_count <- table(first_area_info)['+'] + table(second_area_info)['-']
        minus_read_count <- table(first_area_info)['-'] + table(second_area_info)['+']
    }
    if (plus_read_count / minus_read_count > 2) {
        return("+")
    } else if (minus_read_count / plus_read_count > 2) {
        return("-")
    }
    return("*")
}

#' Actually get info from bam for given chromosome
#'
#' @author Irina Shchyukina
#'
#' @keywords internal
#'
#' @importMethodsFrom GenomeInfoDb seqlengths
#' @importFrom magrittr "%>%"
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom Rsamtools pileup scanBamFlag ScanBamParam scanBam getSeq
#' @importFrom lubridate now
#' @importFrom dplyr mutate arrange filter select rename bind_rows bind_cols lead lag row_number
#' @importFrom pryr object_size
#' @importFrom rlang .data
#' @importFrom tidyr spread
fixed_chromosome_positions <- function(
    bam_file,
    current_chromosome,
    bam_index,
    pileup_params,
    fasta_file,
    min_coverage,
    strand_specific,
    single_end)
{
    print(current_chromosome)
    print("Get raw pileups from BAM")
    print(now())
    current_chromosome_length <- seqlengths(bam_file)[current_chromosome]
    current_range <- GRanges(current_chromosome, IRanges(1, current_chromosome_length))
    if (single_end) {
        bam_flag <- scanBamFlag(
            isSecondaryAlignment = FALSE,
            isNotPassingQualityControls = FALSE,
            isDuplicate = FALSE)
        bam_params <- ScanBamParam(which = current_range, flag = bam_flag)
        pileup_table <- pileup(
            file = bam_file,
            index = bam_index,
            pileupParam = pileup_params,
            scanBamParam = bam_params)
        pileup_table$which_label <- NULL
        if (!strand_specific) {
            pileup_table <- mutate(pileup_table, strand = "*")
        }
    } else {
        if (strand_specific) {
            # get raw pileup (strand-specific RNA-Seq)
            bam_flag <- scanBamFlag(
                isPaired = TRUE,
                isProperPair = TRUE,
                isSecondaryAlignment = FALSE,
                isNotPassingQualityControls = FALSE,
                isDuplicate = FALSE,
                hasUnmappedMate = FALSE,
                isFirstMateRead = TRUE)
            bam_params <- ScanBamParam(which = current_range, flag = bam_flag)
            pileup_table_first <- pileup(
                file = bam_file,
                index = bam_index,
                pileupParam = pileup_params,
                scanBamParam = bam_params)
            pileup_table_first$which_label <- NULL
            bam_flag <- scanBamFlag(
                isPaired = TRUE,
                isProperPair = TRUE,
                isSecondaryAlignment = FALSE,
                isNotPassingQualityControls = FALSE,
                isDuplicate = FALSE,
                hasUnmappedMate = FALSE,
                isSecondMateRead = TRUE)
            bam_params <- ScanBamParam(which = current_range, flag = bam_flag)
            pileup_table_second <- pileup(
                file = bam_file,
                index = bam_index,
                pileupParam = pileup_params,
                scanBamParam = bam_params)
            pileup_table_second$which_label <- NULL
            pileup_table_second <- pileup_table_second %>%
                mutate(strand = ifelse(.data$strand == "*", "*", ifelse(.data$strand == "+", "-", "+")))
            print(now())

            pileup_table <- bind_rows(pileup_table_first, pileup_table_second)
            rm(pileup_table_first, pileup_table_second)
        } else {
            # get raw pileup (non-strand-specific RNA-Seq)
            print("Get raw pileups from BAM")
            print(now())
            bam_flag <- scanBamFlag(
                isPaired = TRUE,
                isProperPair = TRUE,
                isSecondaryAlignment = FALSE,
                isNotPassingQualityControls = FALSE,
                isDuplicate = FALSE,
                hasUnmappedMate = FALSE)
            bam_params <- ScanBamParam(which = current_range, flag = bam_flag)
            pileup_table <- pileup(
                file = bam_file,
                index = bam_index,
                pileupParam = pileup_params,
                scanBamParam = bam_params)
            pileup_table$which_label <- NULL
            pileup_table <- mutate(pileup_table, strand = "*")
            print(now())
        }
    }
    print("Combine two pileups")
    print(now())
    pileup_table <- pileup_table %>%
        arrange(.data$pos, .data$strand, .data$nucleotide) %>%
        mutate(is_next_same = (lead(.data$pos) == .data$pos & lead(.data$strand) == .data$strand & lead(.data$nucleotide) == .data$nucleotide)) %>%
        mutate(
            count = ifelse(!is.na(lag(.data$is_next_same)),
            ifelse(lag(.data$is_next_same), .data$count + lag(.data$count), .data$count),
            .data$count)) %>%
        filter(!.data$is_next_same) %>%
        select(-.data$is_next_same)
    print(pryr::object_size(pileup_table))
    if (nrow(pileup_table) == 0) { return(NULL) }
    print(now())

    # get reference nucleotides
    print("Get reference nucleotides")
    print(now())
    positions <- GRanges(pileup_table$seqnames, IRanges(start = pileup_table$pos, end = pileup_table$pos))
    reference_base <- getSeq(fasta_file, positions)
    pileup_table <- pileup_table %>% bind_cols(as.data.frame(reference_base)) %>% rename(reference = .data$x)
    print(pryr::object_size(pileup_table))

    # spread rows by position
    print("Spread rows by position")
    print(now())
    pileup_table <- pileup_table %>% spread(key = .data$nucleotide, value = .data$count, fill = 0)
    print(pryr::object_size(pileup_table))

    # calculate coverage, fraction and filter by coverage
    print("Leave only mismatches, filter by coverage")
    print(now())
    pileup_table <- pileup_table %>%
        mutate(is_mismatch = is_mismatch_vec(.data$reference, .data$A, .data$C, .data$G, .data$T)) %>%
        filter(
            is_mismatch |
            ((lead(.data$pos) == .data$pos) & lead(.data$is_mismatch)) |
            ((lag(.data$pos) == .data$pos) & lag(.data$is_mismatch))) %>%
        mutate(coverage = .data$A + .data$G + .data$T + .data$C) %>%
        filter(.data$coverage > .data$min_coverage)
    if (nrow(pileup_table) == 0) { return(NULL) }
    print(pryr::object_size(pileup_table))
    if (strand_specific) {
        # choose strand
        print("Choose strand to use")
        print("Find paired positions, split data")
        print(now())
        pileup_table <- pileup_table %>% mutate(is_paired = (lead(.data$pos) == .data$pos | lag(.data$pos) == .data$pos))
        if (nrow(pileup_table) > 1) {
            pileup_table[1, "is_paired"] = (pileup_table[1, "pos"] == pileup_table[2, "pos"])
            pileup_table[nrow(pileup_table), "is_paired"] = (pileup_table[nrow(pileup_table), "pos"] == pileup_table[nrow(pileup_table) - 1, "pos"])
        } else {
            pileup_table[1, "is_paired"] = FALSE
        }
        pairs <- pileup_table %>%
            filter(.data$is_paired) %>%
            mutate(is_mess = 0)
        pileup_table <- pileup_table %>% filter(!.data$is_paired)

        # 0 - don't know
        # 1 - correct
        # 2 - delete
        # try coverage ratio
        print("Try coverage ratio")
        print(now())
        pairs <- pairs %>%
            mutate(
                is_mess = ifelse(row_number() %% 2 == 1,
                ifelse(.data$coverage / lead(.data$coverage) > 2,
                1,
                ifelse(lead(.data$coverage) / .data$coverage > 2, 2, 0)),
                0)) %>%
            mutate(
                is_mess = ifelse(row_number() %% 2 == 0,
                ifelse(lag(.data$is_mess) == 1,
                2,
                ifelse(lag(.data$is_mess) == 2, 1, 0)),
                .data$is_mess))
        pileup_table <- bind_rows(pileup_table, filter(pairs, .data$is_mess == 1) %>% select(-.data$is_mess))
        messy_positions <- filter(pairs, .data$is_mess == 2) %>% select(-.data$is_mess, -.data$is_paired)
        pairs <- filter(pairs, .data$is_mess == 0)

        # delete sites with low coverage
        print("Filter by coverage")
        print(now())
        coverage_threshold <- 10
        pairs <- pairs %>%
            mutate(
                is_mess = ifelse(row_number() %% 2 == 1,
                ifelse(.data$coverage < coverage_threshold & lead(.data$coverage) < coverage_threshold, 2, 0),
                ifelse(.data$coverage < coverage_threshold & lag(.data$coverage) < coverage_threshold, 2, 0)))
        messy_positions <- bind_rows(messy_positions, filter(pairs, .data$is_mess == 2) %>% select(-.data$is_paired, -.data$is_mess))
        pairs <- filter(pairs, .data$is_mess == 0)

        # count reads
        print("Count reads around")
        print(now())
        if (nrow(pairs) > 0) {
            for (i in 1:nrow(pairs)) {
                if (i %% 2 == 0) {
                    first_mate_flag <- pairs[i - 1, "is_mess"][[1]]
                    if (first_mate_flag == 0) {
                        pairs[i, "is_mess"] <- 0
                    } else if (first_mate_flag == 1) {
                        pairs[i, "is_mess"] <- 2
                    } else {
                        pairs[i, "is_mess"] <- 1
                    }
                } else {
                    correct_strand <- get_strand(
                        pairs[i, "pos"][[1]],
                        current_chromosome_length,
                        current_chromosome,
                        bam_file,
                        single_end)
                    if (correct_strand == "*") {
                        pairs[i, "is_mess"] <- 0
                    } else if (correct_strand == pairs[i, "strand"][[1]]) {
                        pairs[i, "is_mess"] <- 1
                    } else {
                        pairs[i, "is_mess"] <- 2
                    }
                }
            }
        }
        pileup_table <- bind_rows(pileup_table, filter(pairs, .data$is_mess != 2) %>% select(-.data$is_mess)) %>% select(-.data$is_paired)
        messy_positions <- bind_rows(messy_positions, filter(pairs, .data$is_mess == 2) %>% select(-.data$is_mess, -.data$is_paired)) %>% arrange(.data$pos)
        pairs <- filter(pairs, .data$is_mess == 0) %>% select(-.data$is_mess, -.data$is_paired) %>% arrange(.data$pos)
    }
    pileup_table <- pileup_table %>% filter(.data$is_mismatch) %>% select(-.data$is_mismatch) %>% arrange(.data$pos)

    # output tables
    print("End of processing")
    print(now())

    if (strand_specific) {
        list(pileup_table = pileup_table, messy_positions = messy_positions, pairs = pairs)
    }
    else {
        list(pileup_table = pileup_table)
    }
}

#' @title Reads the BAM file
#'
#' @author Irina Shchyukina
#'
#' @description Separates data from the input BAM file into a number of data frames, each suitable for the further RNA editing analysis
#'
#' @param bam_file_name A \code{character(1)} path to the input BAM file
#' @param reference_fasta A \code{character(1)} path to the reference genome FASTA file
#' @param max_depth An \code{integer(1)} maximum number of nucleotides to be included in pileup
#' @param min_base_quality An \code{integer(1)} minimum QUAL value for each nucleotide in an alignment
#' @param min_mapq An \code{integer(1)} minimum MAPQ value for an alignment to be included in pileup
#' @param min_nucleotide_depth An \code{integer(1)} minimum count of each nucleotide at a given position required for said nucleotide to appear in the result
#' @param min_minor_allele_depth An \code{integer(1)} left undocumented
#' @param min_coverage An \code{integer(1)} minimum coverage for position to be included in final dataset
#' @param use_noncanonical A \code{logical(1)} indicating  whether to use (\code{TRUE}) all chromosomes from BAM including non-canonical ones'
#' @param strand_specific A \code{logical(1)} indicating the presence either of the strand-specific (\code{TRUE}) RNA-seq experiment, or non-strand-specific (\code{FALSE}) one
#' @param single_end A \code{logical(1)} indicating the presence of the single-end (\code{TRUE}) RNA-seq experiment
#'
#' @return On \code{strand_specific == TRUE} the list of three data frames (with names \code{pileup_table}, \code{messy_positions}, and \code{pairs}), otherwise the list of only one data frame (with name \code{pileup_table}) suitable for the further RNA editing analysis
#'
#' @importMethodsFrom GenomeInfoDb seqlengths
#' @importFrom Rsamtools BamFile FaFile PileupParam indexBam
#'
#' @export
bamrd <- function(
    bam_file_name,
    reference_fasta,
    max_depth = 8000,
    min_base_quality = 13,
    min_mapq = 30,
    min_nucleotide_depth = 1,
    min_minor_allele_depth = 0, # get all covered sites
    min_coverage = 0,
    use_noncanonical = TRUE,
    strand_specific = FALSE,
    single_end = TRUE)
{
    # check if bam exists: throw and stop
    if (!file.exists(bam_file_name)) { stop(paste0(bam_file_name, ".bam -- BAM file does not exist")) }

    # check if index exists: create
    bam_index = paste0(bam_file_name, ".bai")
    if (!file.exists(bam_index))

    # check if reference genome exists: throw and stop
    if (!file.exists(reference_fasta)) { stop(paste0(reference_fasta, " -- reference FASTA does not exist")) }

    bam_file <- BamFile(bam_file_name, bam_index)
    fasta_file <- FaFile(file = reference_fasta)
    pileup_params <- PileupParam(
        max_depth = max_depth,
        min_base_quality = min_base_quality,
        min_mapq = min_mapq,
        min_nucleotide_depth = min_nucleotide_depth,
        min_minor_allele_depth = min_minor_allele_depth,
        include_deletions = FALSE)
    chromosome_names <- names(seqlengths(bam_file))
    if (!use_noncanonical) {
        canonical_names <- c(
            "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19",
            "20", "21", "22", "X", "Y", "M", "MT", "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
            "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19",
            "chr20", "chr21", "chr22", "chrX", "chrY", "chrM", "chrMT")
        chromosome_names <- intersect(chromosome_names, canonical_names)
    }
    print("Chromosome names:")
    print(chromosome_names)

    # prepare output data

    if (strand_specific) {
        do.call(
            function(x, y) { list(
                pileup_table = rbind(x$pileup_table, y$pileup_table),
                messy_positions = rbind(x$messy_positions,y$messy_positions),
                ambiguous_file = rbind(x$ambiguous_file,y$ambiguous_file)) },
            lapply(
                chromosome_names,
                function(chr) fixed_chromosome_positions(
                    bam_file, chr, bam_index, pileup_params, fasta_file, min_coverage, strand_specific, single_end)))
    } else {
        do.call(
            function(x, y) { list(pileup_table = rbind(x$pileup_table, y$pileup_table)) },
            lapply(
                chromosome_names,
                function(chr) fixed_chromosome_positions(
                    bam_file, chr, bam_index, pileup_params, fasta_file, min_coverage, strand_specific, single_end)))
    }
}
