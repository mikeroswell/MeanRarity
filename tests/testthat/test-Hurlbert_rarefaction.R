test_that("Hurlbert rarefaction matches vegan", {
  ab <- 1:50
  k <- 27
  expect_equal(hRare(ab, k)
               , as.numeric(vegan::rarefy(ab, k)))
})
