biobroom: Tidying up computational biology
====================

This package contains methods for converting standard objects constructed by bioinformatics packages, especially those in [BioConductor](http://www.bioconductor.org/), and converting them to [tidy data](http://www.jstatsoft.org/v59/i10). It thus serves as a complement to the [broom package](https://github.com/dgrtwo/broom), and follows the same the tidy/augment/glance division of tidying methods. Tidying data makes it easy to recombine, reshape and visualize bioinformatics analyses.

biobroom implements tidying methods for both S3 and S4 classes. Objects that can be tidied include

* ExpressionSet objects
* per-gene differential expression tests from limma, edgeR, and DESeq2
* [qvalue](http://www.bioconductor.org/packages/release/bioc/html/qvalue.html) multiple hypothesis testing objects

Installation
------------

First install the package's requirements (requires [devtools](https://github.com/hadley/devtools)):

    source("http://bioconductor.org/biocLite.R")
    biocLite(c("Biobase", "limma", "edgeR", "DESeq2", "GenomicRanges"))
    devtools::install_github("dgrtwo/broom")

Then the package can be installed with

    devtools::install_github("dgrtwo/biobroom")

Find out more about the provided methods with:

    library(biobroom)
    ?edgeR_tidiers
    ?DESeq2_tidiers
    ?limma_tidiers
    ?ExpressionSet_tidiers

Note on returned values
------------

All biobroom `tidy` and `augment` methods, since they tend to be large data frames, return a [tbl_df](http://www.inside-r.org/packages/cran/dplyr/docs/tbl_df) by default (this prevents them from printing many rows at once, while still acting like a traditional data.frame). To change this to a data.frame or data.table, you can set the `biobroom.return` option:

    options(biobroom.return = "data.frame")
    options(biobroom.return = "data.table")
