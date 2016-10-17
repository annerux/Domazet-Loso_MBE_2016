# Domazet-Loso_MBE_2016
Scripts for "No evidence for phylostratigraphic bias impacting inferences on patterns of gene emergence and evolution"

Scripts necessary to reproduce Figures 1C&D, 2, 6 and 7 as well as Tables 1 and S3 were contributed by Anne-Ruxandra Carvunis
These are R scripts that take as input the files "GeneData6.txt" and "dm3_mixed_map.txt". Please modify the address for these files according to where you download them in the #DATA section of the scripts.

### R script "calculate_hyper.R"

was used to calculate log-odds and p-values for Figures 3, 4 and 5. Results of these analyses can be found in Table S2.

#### SYNOPSIS

    # take input file, calculate hypergeometric test and write it all to output file
    ./calculate_hyper.R -i ./t/data/Figure3/ectoderm.txt -o ./t/data/Figure3/ectoderm_hyper.txt

#### DESCRIPTION

This script requires input file in TabSeparated format with this header:

    phylostrata\tphylostrata_name\tFunctional term\tquant\tsample\thit\ttotal

after header paste your values like this:

    1\tCellular organisms\tdisease_genes\t970\t8285\t1760\t22845


it will then calculate log-odds, hypergeometric test, FDR and Bonferroni correction and write it to tsv file.

#### AUTHOR

Copyright (C) <2016> <Martin Sebastijan Šestak> (<sestakm@yahoo.com>). All rights reserved.

#### LICENSE

Copyright (C) Martin Sebastijan Šestak (<sestakm@yahoo.com>)

This program is free software; you can redistribute it and/or
modify it under the same terms as Perl itself.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 

