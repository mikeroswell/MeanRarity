test_that("MeanRarity::rarity gives same output as equivalent functions in other
          pacakges", {
  expect_equal(rarity(1:10, 0)
               , as.numeric(vegan::renyi(1:10, scales = 1, hill = TRUE )))

  expect_equal(rarity(1:10, 0)
               , as.numeric(entropart::Diversity(1:10, 1, Correction = "None")))

  expect_equal(rarity(1:10, 0)
               , iNEXT::iNEXT(1:10
                              , q = 1
                              , size = 55
                              , se = FALSE)$iNextEst$size_based$qD
               )
})
