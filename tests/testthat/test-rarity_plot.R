
ab <- c(20,8,5,4,2,1)
test_that("output of rarity_plot is stable", {
  vdiffr::expect_doppelganger("a rarity balance plot", rarity_plot( ab, l = 0))
})
