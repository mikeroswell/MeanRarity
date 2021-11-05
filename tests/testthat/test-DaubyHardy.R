test_that("ENS_PIE is Hill-Simpson... with finite sample corrections!", {
  ab <- 10:50
  expect_equal(Chao_Hill_abu(ab, -1), hurl(ab, 2))
})
