test_that("simple grid construction from dimensions", {

  # test parameters
  gdim = c(12, 10)
  g = sk(gdim)

  # pass result to sk and get the same thing back
  expect_identical(g, sk(g))

  # supply grid lines as named argument instead and get the same result
  expect_equal(g, sk(gyx=lapply(gdim, function(x) seq(x)-1L)))
})

test_that("single argument (unnamed) can be grid dimensions, with shorthand for square grids", {

  # test parameters
  gdim = c(12, 10)
  g = sk(gdim)

  # single argument (unnamed) can be grid dimensions, with shorthand for square grids
  expect_equal(sk(gdim=c(2,2)), sk(c(2,2)))
  expect_equal(sk(2), sk(gdim=c(2,2)))
})

test_that("using R's default matrix vectorization order", {

  gdim = c(25, 25)
  gyx = as.list(expand.grid(lapply(gdim, seq)))
  eg_vec = as.numeric( gyx[[2]] %% gyx[[1]] )
  eg_mat = matrix(eg_vec, gdim)
  g = sk(eg_mat)

  # this is R's default matrix vectorization order
  expect_equal(eg_vec, as.vector(eg_mat))
  expect_equal(g, sk(gdim=gdim, gval=as.vector(eg_mat)))
})

test_that("idx_grid maps observed to full grid on multi-layer example", {

  # multi-layer example from matrix
  gdim = c(25, 25)
  n_pt = prod(gdim)
  n_layer = 3
  mat_multi = matrix(stats::rnorm(n_pt*n_layer), n_pt, n_layer)
  g_multi = sk(gdim=gdim, gval=mat_multi)

  # repeat with missing data (note all columns must have consistent NA structure)
  mat_multi[sample.int(n_pt, 0.5*n_pt),] = NA

  g_multi_miss = sk(gdim=gdim, gval=mat_multi)

  # only observed data points are stored, idx_grid maps them to the full grid vector
  mat_diff <- g_multi[['gval']] - g_multi_miss[['gval']][g_multi_miss[['idx_grid']],]
  expect_equal(max(abs(mat_diff), na.rm=TRUE), 0)

  # single bracket indexing magic does the mapping automatically
  expect_equal(max(abs(g_multi[] - g_multi_miss[]), na.rm=TRUE), 0)
})


test_that("vals=FALSE drops layer information from multi-layer example", {

  # multi-layer example from matrix
  gdim = c(25, 25)
  n_pt = prod(gdim)
  n_layer = 3
  mat_multi = matrix(stats::rnorm(n_pt*n_layer), n_pt, n_layer)
  g_multi = sk(gdim=gdim, gval=mat_multi)

  expect_true(all(is.na((sk(gdim=gdim, gval=mat_multi, vals=FALSE)))))
})

