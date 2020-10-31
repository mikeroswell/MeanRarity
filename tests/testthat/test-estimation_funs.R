test_that("God isn't dumb", {
  ab <-stats::runif(25)
  ell <- -1:1
  expect_equal(sapply(ell, function(l){rarity(ab, l)})
               , sapply(ell, function(l){GUE(ab, ab, l)}))
})
