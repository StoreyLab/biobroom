context("limma")

test_that("limma tidier works as expected", {
    dds <- makeExampleDataSet("limma")
    td <- tidy(dds)
    tdp <- tidy(dds, addTargets=TRUE)

    ## Check that both objects have the bare minimum columns with expected
    ## types
    req.cols <- c(gene='character', sample='character', value='numeric')
    for (rc in names(req.cols)) {
        expect_is(td[[rc]], req.cols[[rc]], info=sprintf('td (%s)', rc))
        expect_is(tdp[[rc]], req.cols[[rc]], info=sprintf('tdp (%s)', rc))
    }

    ## tdp has the extra sample data
    pd <- dds$targets
    pd.cols <- setNames(sapply(pd, class), names(pd))
    for (pc in names(pd.cols)) {
        expect_is(tdp[[pc]], pd.cols[[pc]], info=sprintf("tdp pheno (%s)", pc))
    }
})

test_that("voom tidier adds weight column", {
    dds <- makeExampleDataSet("voom")
    td <- tidy(dds)
    tdp <- tidy(dds, addTargets=TRUE)

    ## weights added correctly
    expect_is(td[['weight']], 'numeric')
    expect_is(tdp[['weight']], 'numeric')
    expect_equal(td[['weight']], tdp[['weight']])

    ## vanilla limma tidy has been verified, so ensure that these tidied objects
    ## match base limma tidies objects
    elist <- dds
    elist$weights <- NULL
    elist$design <- NULL

    ld <- tidy(elist)
    ldp <- tidy(elist, addTargets=TRUE)

    expect_equal(transform(td, weight=NULL), ld)
    expect_equal(transform(tdp, weight=NULL), ldp)
})

test_that("voomWithQualityWeights tidier adds weight and sample.weight columns", {
    dds <- makeExampleDataSet("voomWithQualityWeights")
    td <- tidy(dds)
    tdp <- tidy(dds, addTargets=TRUE)

    ## weights and sample.weight added correctly
    expect_is(td[['weight']], 'numeric')
    expect_is(tdp[['weight']], 'numeric')
    expect_equal(td[['weight']], tdp[['weight']])

    expect_is(td[['sample.weight']], 'numeric')
    expect_is(tdp[['sample.weight']], 'numeric')
    expect_equal(td[['sample.weight']], tdp[['sample.weight']])

    ## vanilla limma tidy has been verified, so ensure that these tidied objects
    ## match base limma tidies objects
    elist <- dds
    elist$weights <- NULL
    elist$design <- NULL

    ld <- tidy(elist)
    ldp <- tidy(elist, addTargets=TRUE)

    expect_equal(transform(td, weight=NULL), ld)
    expect_equal(transform(tdp, weight=NULL), ldp)
})
