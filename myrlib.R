# My R library

#' Count the total number of bases in a Bed file
#' 
#' @param uniqueBedFile A BED file path
#' @return The total number of bps in \code{uniqueBedFile}.
#' @examples
#' TotalBaseCoverage(bedfile)
#' 
#' @export
TotalBaseCoverage <- function(uniqueBedFile) {
    df <- read.delim(uniqueBedFile, header = FALSE)
    #colnames(df[,1:3]) <- c("chr", "start", "end")
    #print(head(df))
    sum(df[,3] - df[,2])
}


#' Depth of coverage by the total bases in a Bed file
#' 
#' @param uniqueBedFile A BED file path
#' @param genomeSize The size of the genome
#' @return The depth of coverage by \code{uniqueBedFile}.
#' @examples
#' GenomeCoverage(bedfile, genomeSize)
#' 
#' @export
GenomeCoverage <- function(uniqueBedFile, genomeSize) {
    TotalBaseCoverage(uniqueBedFile)/genomeSize
}


#' Count the number of fragments in a Bed file
#' Implemented using Unix command wc or native R functions
#' 
#' @param bedFile A BED file path
#' @return The number of fragments in \code{bedFile}.
#' @examples
#' CountFragments(bedfile)
#' 
#' @export
CountFragments <- function(bedFile, useunix = FALSE) {
    fc <- 0
    if (useunix) {
        cmd <- paste0("wc -l ", bedFile)
        #as.integer(unlist(strsplit(trimws(system(cmd, intern = TRUE)), "\\D+")))
        fc <- as.integer(unlist(strsplit(trimws(system(cmd, intern = TRUE)), " "))[1])
    } else {
        fc <- length( readLines(bedFile) )
    }
    fc
}


#' Generate a random DNA sequence
#' 
#' @param len The length of a desired random DNA sequence
#' @param fixed State of seed, default to FALSE
#' @param lowercase Lowercase or not, default to FALSE
#' @param CG The desired CG content, default to 0.5 
#' @return A random DNA sequence with specified parameters
#' @examples
#' RandomDNA(50)
#' 
#' @export
RandomDNA <- function(len, fixed = FALSE, lowercase=FALSE, CG=0.5) {
    alphabet <- c("A", "T", "C", "G")
    if (lowercase) 
        alphabet <- tolower(alphabet)
    if (fixed) 
        set.seed(123)
    AT <- 1-CG
    paste(sample(alphabet, len, replace = TRUE, 
                 prob = c(AT/2, AT/2, CG/2, CG/2)), collapse = "")
}

#' Generate a random DNA sequence
#' 
#' @param len The length of a desired random DNA sequence
#' @param firstLine The length of bases at the first line, default to 0
#' @param perLine The max bases allowed per line, default to 60
#' @param fixed State of seed, default to FALSE
#' @param lowercase Lowercase or not, default to FALSE
#' @param CG The desired CG content, default to 0.5 
#' @return A random DNA sequence with specified parameters.
#' @examples
#' RandomDNA(50)
#' 
#' @export
FormattedRandomDNA <- function(len, firstLine=0, perLine=60, fixed = FALSE, lowercase=FALSE, CG=0.5) {
    stopifnot(firstLine >= 0)
    stopifnot(perLine > 0)
    stopifnot((perLine - firstLine) > 0)
    
    dnaSeq <- RandomDNA(len, fixed, lowercase, CG)
    
    if (firstLine != 0) {
        leadingSpaces <- paste0(rep(" ", perLine - firstLine), collapse = "")
        dnaSeq <- paste0(leadingSpaces, dnaSeq)
        len <- len + perLine - firstLine
    }
    
    if(len %% perLine == 0) {
        substring(dnaSeq, seq(1, len, perLine), seq(perLine, len, perLine))
    } else {
        substring(dnaSeq, seq(1, len, perLine),c(seq(perLine, len, perLine), len))
    }
}

#' Write a shell script that executes a set of commands for each of provided folders
#' 
#' @param commands A vector of shell commands
#' @param data_dirs A vector of data folders
#' @param cmd_file_path The complete file path to be saved
#' @return The path of the created shell script file.
#' @examples
#' write_shell_commands(commands, data_dirs, cmd_file_path)
#' 
#' @export
write_shell_commands<- function(commands, data_dirs, cmd_file_path) {
    
    if (file.exists(cmd_file_path)) {
        stop( paste0("The file already exists: ", cmd_file_path) )
    }
    
    file.create(cmd_file_path)
    file_conn<-file(cmd_file_path, open = "a")
    writeLines("#!/bin/bash", file_conn)
    writeLines("set -eu -o pipefail", file_conn)
    
    # curr_dir <- getwd()
    
    for (d in data_dirs) {
        cd_cmd <- paste("cd", d)
        writeLines("echo =============================================================", 
                   file_conn)
        writeLines(paste("echo", cd_cmd), file_conn)
        writeLines(cd_cmd, file_conn)
        for (cmd in commands) {
            writeLines(paste("echo", cmd), file_conn)
            writeLines(cmd, file_conn)
        }
    }
    writeLines("echo =============================================================", 
               file_conn)
    writeLines('say "done"', file_conn)
    close(file_conn)
    
    # make the file executable
    system(paste0("chmod 755 ", cmd_file_path))
    # reset the working dir back to the saved dir
    # setwd(curr_dir)
    cmd_file_path
}


#' Quantile Normalize a data set to a reference distribution
#' 
#' @param data A vector of real numbers
#' @param ref A vector of data specifying the reference distribution
#' @param method A string that describe how to normalize values between two quantiles
#' @return quantile-normalized data 
#' @examples
#' norm_qutantile(data, ref, method = 'average')
#' 
#' @export
norm_qutantile = function(data, ref, method = 'average') {
    if(! method %in% c('average', 'linear') ) {
        stop("The method should be either average or linear.")
    }
    data.sorted = sort(data, decreasing = TRUE)
    data.rank = rank(-data)
    ref.sorted = sort(ref, decreasing = TRUE)
    data.sorted.q = rank(data.sorted)/length(data.sorted)*100
    ref.sorted.q = rank(ref.sorted)/length(ref.sorted)*100
    
    indexes = vector("list", length(data.sorted.q))
    
    j=1
    is_end_of_ref = FALSE
    for (i in seq_along(data.sorted.q)) {
        while(!is_end_of_ref && data.sorted.q[i] < ref.sorted.q[j]) {
            j = j+1
            if(j>length(ref.sorted.q)) {
                is_end_of_ref = TRUE
            }
        }
        
        if(is_end_of_ref) {
            indexes[[i]] = c(j-1, j-1)
        } else if (data.sorted.q[i] == ref.sorted.q[j]) {
            indexes[[i]] = c(j, j)
        } else{
            indexes[[i]] = c(j-1, j)
        }
    }
    
    norm.sorted.data = sapply(seq_along(data.sorted), function(i) {
        if (indexes[[i]][1] == indexes[[i]][2]) {
            return(ref.sorted[ indexes[[i]][1] ])
        } else {
            # take the average
            if(method=='average') {
                return((ref.sorted[indexes[[i]][1]] + ref.sorted[ indexes[[i]][2] ])/2) 
            } else {
                # linear interpolation
                r = (ref.sorted.q[indexes[[i]][1]] - data.sorted.q[i])/(ref.sorted.q[indexes[[i]][1]] - ref.sorted.q[indexes[[i]][2]])
                return(ref.sorted[indexes[[i]][1]] * (1 - r))
            }
        }
    })
    
    
    norm.data = sapply(data.rank, function(i) {
        norm.sorted.data[i]
    })
    
    return(norm.data)
}


#' Quantile Normalize a data set to a reference distribution not counting zeros
#' Zeros in both data and reference are not used for quantile normalization
#' 
#' @param data A vector of real numbers
#' @param ref A vector of data specifying the reference distribution
#' @param method A string that describe how to normalize values between two quantiles
#' @return quantile-normalized data after removing zeroes from data and reference
#' @examples
#' norm_qutantile_rmzeros(data, ref, method = 'average')
#' 
#' @export
norm_qutantile_rmzeros = function(data, ref, method = 'average') {
    if(! method %in% c('average', 'linear') ) {
        stop("The method should be either average or linear.")
    }
    data2 = data[data!=0]
    ref2 = ref[ref!=0]
    norm.data2 = norm.qutantile(data2, ref2, method=method)
    norm.data = rep(0, length(data))
    norm.data[which(data!=0)] = norm.data2
    return(norm.data)
}


#' Test the significance of two datasets of sequence count data
#' 
#' @param count_1 the first vector of seq count numbers
#' @param count_2 the second vector of seq count numbers
#' @param breaks a vector giving the breakpoints between desired bins
#' @return dataframe about test statistics about the difference between \code{count_1} and \code{count_2}   
#' @examples
#' counts_sig_test(count_1, count_2, breaks)
#' 
#' @export
counts_sig_test = function(count_1, count_2, 
                           breaks = c(-1, 0, 0.25, 0.5, 0.8, 1.2, 1.8, 2.6, 3.6, 7, 9, 10, 15, 25, 1000)) {
    require(dplyr)
    
    percent_std = function(v) {
        s = sd(v)
        if (v[1]==0 & v[2]!=0) {
            return(c(s, s/v[2]))
        } else if (v[2]==0 & v[1]!=0) {
            return(c(s/v[1], s))
        } else if ((v[1]==0 & v[2]==0)) {
            return(c(NA, NA))
        } else {
            return(c(s/v[1], s/v[2]))
        }
    }
    
    df = data.frame(count_1, count_2)
    res = apply(df, 1, percent_std)
    df2 = cbind(df, psd_1 = res[1,], psd_2 = res[2,])
    df2$count_1_group = cut(df2$count_1, breaks = brk, labels = 1:(length(breaks)-1))
    df2$count_2_group = cut(df2$count_2, breaks = brk, labels = 1:(length(breaks)-1))
    
    df2 = df2 %>%
        group_by(count_1_group) %>%
        summarize(mean_psd_1=mean(psd_1, na.rm = TRUE)) %>%
        left_join(df2, ., by="count_1_group")
    
    df2 = df2 %>%
        group_by(count_2_group) %>%
        summarize(mean_psd_2=mean(psd_2, na.rm = TRUE)) %>%
        left_join(df2, ., by="count_2_group") %>%
        mutate(sd_1=ifelse(count_1==0, mean_psd_1, mean_psd_1*count_1)) %>%
        mutate(sd_2=ifelse(count_2==0, mean_psd_2, mean_psd_2*count_2)) %>%
        mutate(Z_score=(count_2-count_1)/sqrt(sd_1*sd_1 + sd_2*sd_2)) %>%
        mutate(logFC=ifelse(count_1==0 & count_2==0, NA, log2(count_2/count_1))) %>%
        mutate(p_value=2*pnorm(-abs(Z_score))) 
    
    # Multiple testing correction BH
    
    # Consider the 0, 0 pairs
    # df2$adj_p_value = p.adjust(df2$p_value, method = "fdr")
    # select(df2, count_1, count_2, logFC, Z_score, p_value, adj_p_value)
    
    # Do not consider the 0, 0 pairs
    df2$adj_p_value = rep(1., nrow(df2))
    df2$adj_p_value[df2$count_1 != 0 | df2$count_2 != 0] = 
        p.adjust(df2$p_value[df2$count_1 != 0 | df2$count_2 != 0], method = "fdr")
    res = select(df2, count_1, count_2, logFC, Z_score, p_value, adj_p_value)
}


#' Differential Enhancer Analysis of paired-end ChIP-Seq data using csaw
#' 
#' @param pebam.files a vector of paired-end bam files with duplicated marked 
#' @param design a design matrix for comparison
#' @param minq the minimum quality valur used to filter out low quality reads 
#' @param window.size the size of the sliding windows
#' @param spacing the distance between two adjacent sliding windows
#' @param bkgrd.width minimum quality valur used to filter out low quality reads
#' @param bkgrd.threshold the filter threshold of fold change over the non-specifc background enrichment
#' @param merge.tol the max distance betwwen two sliding windows allowed for merging
#' @param fdr.threshold the false discovery rate threshold used for BH multi-testing correction
#' @param output.file a file path to save the output 
#' @return None   
#' @examples 
#' diff_enhancer_csaw = function(pebam.files, design, minq=10, window.size=150, spacing=50, bkgrd.width=10000, bkgrd.threshold=5, merge.tol=100, fdr.threshold=0.05, output.file='result.tsv')
#' 
#' @export
diff_enhancer_csaw = function(pebam.files, design, minq=NA, window.size=150, spacing=50, bkgrd.width=10000, bkgrd.threshold=5, merge.tol=100, fdr.threshold=0.05, output.file='result.tsv') {
    require(csaw)
    require(edgeR)
    require(rtracklayer)
    
    # 1. Loading in data from BAM files.
    mm10.backlist = import("~/genome/mm10/mm10.blacklist.bed")
    std.chr = paste0("chr", c(1:19, "X", "Y"))
    max.frag.len = 600
    
    # variable: minq 
    pe.param = readParam(minq = minq, max.frag = max.frag.len, pe="both", restrict = std.chr, discard = mm10.backlist)
    # pe.param = readParam(max.frag = max.frag.len, pe="both", restrict = std.chr, discard = mm10.backlist)
    max.dist = 600
    cr = correlateReads(pebam.files, max.dist, param = reform(pe.param, dedup=TRUE))
    avg.frag.len = maximizeCcf(cr)
    ########## vaiables: width, filter
    data = windowCounts(pebam.files, ext = avg.frag.len, width = window.size, spacing = spacing, filter = 1, param = reform(pe.param, dedup=TRUE))
    
    # -------------histone Marks-----------
    # 2. Filtering out uninteresting windows.
    # By global enrichment
    bin10k = windowCounts(pebam.files, bin=TRUE, width=bkgrd.width, param = reform(pe.param, dedup=TRUE))
    filter.stat = filterWindows(data, bin10k, type = "global")
    
    # threshold of fold chnage over background
    keep = filter.stat$filter > log2(bkgrd.threshold)
    filtered.data = data[keep,]
    
    # 3. Calculating  normalization factors 
    # 3.1 Normalizing IP effiency bias (for histone marks).
    filtered.data = normFactors(filtered.data, se.out=TRUE)
    filtered.data$norm.factors
    # -------------histone Marks-----------
    
    # -------------Txn factors-----------
    #By abundance threshold
    #abundances = aveLogCPM(asDGEList(data))
    #summary(abundances)
    #keep.simple = abundances > -1
    #filtered.data = data[keep.simple,]
    
    # 3. Calculating  normalization factors 
    # 3.2 Normalizing composition bias (for TFs).
    #bin10k = windowCounts(pebam.files, bin=TRUE, width=10000, param = reform(pe.param, dedup=TRUE))
    #filtered.data = normFactor(bin10k, se.out=filtered.data)
    # -------------Txn factors-----------
    
    # glm fit and test
    y = asDGEList(filtered.data)
    y = estimateDisp(y, design)
    fit = glmQLFit(y, design, robust=TRUE)
    results = glmQLFTest(fit, contrast = c(0, 1))
    
    # 5. Correcting for multiple testing.
    ######### variable: tol, <200 or 500-1000; max.width, 2000-10000
    merged = mergeWindows(rowRanges(filtered.data), tol=merge.tol)
    tabcom = combineTests(merged$id, results$table)
    tab.best = getBestTest(merged$id, results$table)
    tabcom$best.logFC = tab.best$logFC
    tabcom$best.start = start(rowRanges(filtered.data))[tab.best$best]
    
    # 6. Output the result
    is.sig = tabcom$FDR <= fdr.threshold
    sig.results = cbind(as.data.frame(merged$region), tabcom)[is.sig,]
    write.table(sig.results, file = output.file, row.names = FALSE, quote = FALSE, sep="\t")
}