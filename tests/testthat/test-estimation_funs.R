test_that("God isn't dumb", {
  ab <-stats::runif(25)
  ell <- -1:1
  expect_equal(sapply(ell, function(l){rarity(ab, l)})
               , sapply(ell, function(l){GUE(ab, ab, l)}))
})



# might get this to work with all.equal, which allows a tolerance
# test_that("God handles undersampling", {
#   ab <-stats::runif(25)
#   ell <- -1:1
#   n <- 75
#   expect_equal(sapply(ell, function(l){rarity(ab, l)})
#                , sapply(ell, function(l){mean(replicate(999999, GUE(sample_infinite(ab, n), ab, l)))}))
# })
