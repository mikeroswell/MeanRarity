ab <- c(20,8,5,4,2,1)
test_that("output of radplot is stable", {
  vdiffr::expect_doppelganger("a radplot", radplot(ab))
})


