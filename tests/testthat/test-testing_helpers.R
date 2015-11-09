context("testing helpers")

test_that("makeExampleDataSet produces dataset of the right type", {
    ds.info <- biobroom:::.ds.info
    for (i in 1:nrow(ds.info)) {
        dtype <- ds.info$type[i]
        eclass <- ds.info$class[i]
        obj <- makeExampleDataSet(dtype)
        expect_is(obj, eclass, info=paste("DataSet for type:", dtype))
    }
})
