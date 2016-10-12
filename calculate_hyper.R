#!/usr/bin/Rscript

# load libraries
#install.packages("optparse", dependencies=TRUE, repos='http://cran.rstudio.com/')
library(optparse)

#> set interface
option_list <- list(

    make_option(c("-f", "--flag"),    action="store_true", default=FALSE, help="Flag [default: %default]"),
    make_option(c("-n", "--numeric"), type="integer",      default=3,     help="Numeric parameter [default: %default]"),
	make_option(c("-i", "--input"),   type="character",    default=FALSE, help="File name"),
	make_option(c("-o", "--output"),  type="character",    default=FALSE, help="File name")
)

parser <- OptionParser(usage="%prog [options]", option_list=option_list)

args <- parse_args(parser, positional_arguments = 0)

opt <- args$options
file <- args$args

#< set interface

# Read table
table <- read.table(opt$input, header=TRUE, sep="\t")

# get columns from table
phylostrata = table$phylostrata
ps_name     = table$phylostrata_name
func_term   = table$Functional.term
quant       = table$quant
sample      = table$sample
hit         = table$hit
total       = table$total

# combine columns from table into dataset
dataset <- cbind.data.frame(phylostrata, ps_name, func_term, quant, sample, hit, total)
odds_sample <-quant / (sample - quant)
odds_rest <- (hit - quant) / (total - hit - sample + quant)
real_log_odds <- log(odds_sample/odds_rest)
dataset <- cbind.data.frame(dataset, odds_sample, odds_rest, real_log_odds)

# hypergeometric calculation
CDFHyper = phyper(dataset$quant, dataset$hit, dataset$total - dataset$hit, dataset$sample)
PDFHyper = dhyper(dataset$quant, dataset$hit, dataset$total - dataset$hit, dataset$sample)
CDFHyperOver = (1 - CDFHyper) + PDFHyper
raw_p_value = pmin(CDFHyper, CDFHyperOver)*2
dataset = cbind.data.frame(dataset, CDFHyper, CDFHyperOver, raw_p_value)

# calculate values for mapping
# sort by raw_p_value to calculate the FDR
dataset$raw_p_value_map = ifelse (dataset$raw_p_value < 0.001, "<0.001", ifelse(dataset$raw_p_value < 0.01, "<0.01", ifelse(dataset$raw_p_value < 0.05, "<0.05", "ns")))
dataset_sorted = dataset[order(dataset$raw_p_value),]
niz = phylostrata
dataset_sorted = cbind.data.frame(dataset_sorted, niz)
dataset_sorted = dataset_sorted[order(dataset_sorted$phylostrata),]
FDR_p_value = dataset_sorted$raw_p_value* max(dataset$phylostrata)/dataset_sorted$niz
dataset_sorted = cbind.data.frame(dataset_sorted, FDR_p_value)
dataset_sorted$for_map_p_value = ifelse (dataset_sorted$FDR_p_value < 0.001, "<0.001", ifelse(dataset_sorted$FDR_p_value < 0.01, "<0.01", ifelse(dataset_sorted$FDR_p_value < 0.05, "<0.05", "ns")))

# Bonferroni correction
Bonferroni_p_value = dataset_sorted$raw_p_value* max(dataset$phylostrata)
dataset_sorted = cbind.data.frame(dataset_sorted, Bonferroni_p_value)
dataset_sorted$for_map_p_value_bon = ifelse (dataset_sorted$Bonferroni_p_value < 0.001, "<0.001", ifelse(dataset_sorted$Bonferroni_p_value < 0.01, "<0.01", ifelse(dataset_sorted$Bonferroni_p_value < 0.05, "<0.05", "ns")))

# write output to file
write.table(dataset_sorted, file = opt$output, quote = FALSE, sep = "\t", row.names=FALSE)


# SYNOPSIS

# ./calculate_hyper.R -i ./t/data/disease.txt -o ./t/data/disease_hyper.txt

# DESCRIPTION

# This script requires input file with this header:
# phylostrata\tphylostrata_name\tFunctional term\tquant\tsample\thit\ttotal
# after header paste your values like this:
# 1\tCellular organisms\tdisease_genes\t970\t8285\t1760\t22845
# ...

# it will then calculate log-odds, hypergeometric test, FDR and Bonferroni correction and write it to tsv file.

# AUTHOR

# <Martin Sebastijan Šestak> (<sestakm@yahoo.com>)

# LICENCE AND COPYRIGHT

# Copyright (c) <2016> <Martin Sebastijan Šestak> (<sestakm@yahoo.com>). All rights reserved.

# This module is free software; you can redistribute it and/or
# modify it under the same terms as Perl itself. See L<perlartistic>.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 


