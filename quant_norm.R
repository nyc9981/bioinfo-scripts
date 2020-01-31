#!/usr/bin/env Rscript
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])

args = commandArgs(trailingOnly=TRUE)
usage = paste("Usage:", script.name, "Input File", "[Reference File]")
Default.Ref.File = "~/genome/hg19/Norm_Std_H3K27AC_2798.tsv"

if(length(args)==0 || length(args)>2) {
    stop(usage, call.=FALSE)
}

if(length(args) == 1) {
    from.file = args[1]
    std.file = Default.Ref.File
} else {
    from.file = args[1]
    std.file = args[2]
}

if(!file_test("-f", from.file)) stop(paste("File", from.file, "does not exist!"))
if(!file_test("-f", std.file)) stop(paste("File", std.file, "does not exist!"))

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

norm_qutantile_rmzeros = function(data, ref, method = 'average') {
    if(! method %in% c('average', 'linear') ) {
        stop("The method should be either average or linear.")
    }
    data2 = data[data!=0]
    ref2 = ref[ref!=0]
    norm.data2 = norm_qutantile(data2, ref2, method=method)
    norm.data = rep(0, length(data))    
    norm.data[which(data!=0)] = norm.data2
    return(norm.data)
}

std.df = read.table(std.file, sep="\t", header = FALSE)
from.df = read.table(from.file, sep="\t", header = FALSE)
data = rev(from.df)[[1]]
target = std.df[[1]]

res = cbind(from.df, norm_qutantile_rmzeros(data, target, method = 'average'))
write.table(res, file = from.file, sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
