
library('entropart')
test_that("MeanRarity::rarity gives same output as equivalent functions in other
          pacakges", {
  ab <- 1:10
  expect_equal(rarity(ab, l = 0)
               , as.numeric(vegan::renyi(ab, scales = 1, hill = TRUE )))

  expect_equal(rarity(ab, l = 0)
               , as.numeric(entropart::Diversity(ab, 1, Correction = "None")))

  expect_equal(rarity(ab, l = 0)
               , iNEXT::iNEXT(ab
                              , q = 1
                              , size = sum(ab)
                              , se = FALSE)$iNextEst$size_based$qD
               )
})
