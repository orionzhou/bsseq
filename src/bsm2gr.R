#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser(description = 'Read bismark methylation report and save as GRanges format')
parser$add_argument("bsm", nargs=1, help="bismark methylation report (*.CX_report.txt)")
parser$add_argument("out", nargs=1, help="output file (*.rds)")
args <- parser$parse_args()

fi = args$bsm
fo = args$out

if(file.access(fi) == -1)
    stop(sprintf("file ( %s ) cannot be accessed", fi))

require(DMRcaller)
met = readBismark(fi)
met.sort = sort(met)
saveRDS(met.sort, file=fo)

