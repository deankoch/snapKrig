# sk_corr
test_that("sk_vario_fun constructs the expected variogram from product", {

  # define test distances, grid, and example kernel
  n_test = 100
  d_test = seq(n_test)-1
  g_example = sk(n_test)
  pars = sk_pars(g_example, c('mat', 'gau'))
  pars_x = pars[['x']]
  corr_x_example = sk_corr(pars_x, d=d_test)

  # get the other component correlation, take product
  pars_y = pars[['y']]
  corr_y_example = sk_corr(pars_y, d=d_test)
  corr_example = corr_y_example * corr_x_example

  # variogram
  variogram_example = sk_vario_fun(pars, d=list(y=d_test, x=d_test))
  variogram_compare = 2 * pars$eps + pars$psill * (1 - corr_example)
  expect_equal(variogram_example, variogram_compare)

  # Toeplitz component matrices built entirely from these correlation vectors
  variance_matrix_example = sk_var(g_example, pars, sep=TRUE)
  expect_equal(variance_matrix_example[['y']][,1L], corr_y_example)
  expect_equal(variance_matrix_example[['x']][,1L], corr_x_example)
})

# sk_corr_mat
test_that("sk_corr_mat and sk_var produce the expected variance matrices", {

  # define test distances, grid, and example kernel
  n_test = 10
  g_example = sk(n_test)
  pars = sk_pars(g_example, c('mat', 'gau'))

  # compute the correlation matrices and their kronecker product
  cx = sk_corr_mat(pars[['x']], n=n_test)
  cy = sk_corr_mat(pars[['y']], n=n_test)
  cxy = kronecker(cx, cy)

  # sk_var can return these two matrices in a list
  cxy_list = sk_var(g_example, pars, sep=TRUE)
  expect_equal(cxy_list[['y']], cy)
  expect_equal(cxy_list[['x']], cx)

  # ... or it can compute the full covariance matrix for model pars (default)
  var_matrix = sk_var(g_example, pars, sep=FALSE)
  var_matrix_compare = (pars$psill*cxy) + diag(pars$eps, n_test^2)
  expect_equal(var_matrix, var_matrix_compare)

  # extract a subgrid without computing the whole thing
  cx_sub = sk_corr_mat(pars[['x']], n=n_test, i=2:4, j=2:4)
  expect_equal(cx_sub, cx[2:4, 2:4])

  # gres scales distances. Increasing gres causes correlations to decrease
  cx_long = sk_corr_mat(pars[['x']], n=n_test, gres=2*g_example$gres)
  expect_true(all(cx_long <= cx))
})

# sk_var

# define example grid with NAs and example predictors matrix
gdim = c(12, 13)
n = prod(gdim)
n_obs = floor(n/3)
idx_obs = sort(sample.int(n, n_obs))
g = g_empty = sk(gdim)
g[idx_obs] = stats::rnorm(n_obs)

# example kernel
psill = 0.3
pars = utils::modifyList(sk_pars(g), list(psill=psill))

# get the full covariance matrix with sep=FALSE (default)
V_full = sk_var(g_empty, pars)
V_obs = sk_var(g, pars)

# cholesky factor and eigen-decomposition
V_obs_chol = sk_var(g, pars, fac_method='chol')
V_obs_eigen = sk_var(g, pars, fac_method='eigen')

# random matrix X
nX = 3
X_all = cbind(1, matrix(stats::rnorm(nX * n), ncol=nX))
X_obs = X_all[idx_obs,]
cprod_obs = crossprod(X_obs, chol2inv(chol(V_obs))) %*% X_obs

test_that("sk_var includes variance of observed as sub-matrix", {

  # check that the correct sub-matrix is there
  expect_equal(V_obs, V_full[idx_obs, idx_obs])
})

test_that("sk_var produces 1d correlation matrices consistent with full covariance matrix", {

  # get 1d correlation matrices with sep=TRUE...
  corr_components = sk_var(g_empty, pars, sep=TRUE)

  # ... these are related to the full covariance matrix through psill and eps
  corr_mat = kronecker(corr_components[['x']], corr_components[['y']])
  V_full_compare = pars$psill * corr_mat + diag(pars$eps, n)
  expect_equal(V_full, V_full_compare)
})

test_that("sk_var computes expected quadratic form with random matrix X", {

  cprod_all = crossprod(X_all, chol2inv(chol(V_full))) %*% X_all
  expect_equal(sk_var(g_empty, pars, X=X_all), cprod_all)

  # test products with inverse of quadratic form with X
  mult_test = stats::rnorm(nX + 1)
  cprod_all_inv = chol2inv(chol(cprod_all))
  cprod_all_inv_chol = sk_var(g_empty, pars, X=X_all, scaled=TRUE, fac_method='eigen')
  mult_result = sk_var_mult(mult_test, pars, fac=cprod_all_inv_chol) - cprod_all_inv %*% mult_test
  expect_equal(mult_result, matrix(0, length(mult_result)))

  # repeat with missing data
  expect_equal(sk_var(g, pars, X=X_obs), cprod_obs)

  # verify inverse
  cprod_obs_inv = chol2inv(chol(cprod_obs))
  cprod_obs_inv_chol = sk_var(g, pars, X=X_obs, scaled=TRUE, fac_method='eigen')
  mult_result = sk_var_mult(mult_test, pars, fac=cprod_obs_inv_chol) - cprod_obs_inv %*% mult_test
  expect_equal(mult_result, matrix(0, length(mult_result)))
})

test_that("scaled=TRUE causes sk_var to divide by psill", {

  # `scaled` indicates to divide matrix by psill
  expect_equal(1+pars[['eps']]/pars[['psill']], unique(diag(sk_var(g, pars, scaled=TRUE))))
  expect_equal(sk_var(g, pars), psill * sk_var(g, pars, scaled=TRUE))
  expect_equal(sk_var(g, pars, X=X_obs, scaled=TRUE), cprod_obs/psill)

  # in Cholesky factor this produces a scaling by square root of psill
  expect_equal(V_obs_chol, sqrt(psill) * sk_var(g, pars, fac_method='chol', scaled=TRUE))

  # and in the eigendecomposition, a scaling of the eigenvalues
  vals_scaled = sk_var(g, pars, fac_method='eigen', scaled=TRUE)$values
  expect_equal(sk_var(g, pars, fac_method='eigen')$values, psill*vals_scaled)
})

# sk_var_mult

# define example grid and data
gdim = c(10, 15)
g = sk(gdim)
n = length(g)
g[] = stats::rnorm(n)

# define covariance parameters
pars = utils::modifyList(sk_pars(g, 'gau'), list(psill=2, eps=0.5))

# compute the full covariance matrix
V = sk_var(g, pars, sep=FALSE)

test_that("sk_var_mult calculates the expected product in complete case", {

  # compute the full covariance matrix
  V = sk_var(g, pars, sep=FALSE)
  V_inv = chol2inv(chol(V))
  out_reference = V_inv %*% g[]
  out_reference_quad = t(g[]) %*% out_reference
  expect_equal(sk_var_mult(g, pars), out_reference)
  expect_equal(sk_var_mult(g, pars, quad=TRUE), out_reference_quad)

  # pre-computed factorization on separable components of correlation matrix
  fac_corr = sk_var(utils::modifyList(g, list(gval=NULL)), pars, fac_method='eigen', sep=TRUE)
  expect_equal(sk_var_mult(g, pars, fac=fac_corr), out_reference)
  expect_equal(sk_var_mult(g, pars, fac=fac_corr, quad=TRUE), out_reference_quad)

  # matrix powers
  out_reference = V %*% g[]
  expect_equal(sk_var_mult(g, pars, fac_method='eigen', p=1), out_reference)
  expect_equal(sk_var_mult(g, pars, fac_method='eigen', p=1, quad=TRUE), t(g[]) %*% out_reference)
})

test_that("sk_var_mult calculates the expected product in incomplete case", {

  # randomly sample missing points
  n_sample = floor(n/10)
  idx_sampled = sort(sample.int(n, n_sample))
  g_miss = sk(gdim)
  g_miss[idx_sampled] = g[idx_sampled]
  V = sk_var(g_miss, pars)

  # correctness check (eigen used by default)
  z = matrix(g[idx_sampled], ncol=1)
  V_inv = chol2inv(chol(V))
  out_reference = (V_inv %*% z)
  out_reference_quad = t(z) %*% out_reference
  expect_equal(sk_var_mult(g_miss, pars), out_reference)
  expect_equal(sk_var_mult(g_miss, pars, quad=TRUE), out_reference_quad)

  # check non-default Cholesky method
  expect_equal(sk_var_mult(g_miss, pars, fac_method='chol'), out_reference)
  expect_equal(sk_var_mult(g_miss, pars, quad=TRUE, fac_method='chol'), out_reference_quad)

  # supply data as a vector instead of list by pre-computing factorization
  fac_chol = sk_var(g_miss, pars, scaled=TRUE, fac_method='chol')
  fac_eigen = sk_var(g_miss, pars, scaled=TRUE, fac_method='eigen')
  expect_equal(sk_var_mult(z, pars, fac=fac_chol), out_reference)
  expect_equal(sk_var_mult(g_miss, pars, fac=fac_eigen), out_reference)
  expect_equal(sk_var_mult(z, pars, fac=fac_chol, quad=TRUE), out_reference_quad)
  expect_equal(sk_var_mult(g_miss, pars, fac=fac_eigen, quad=TRUE), out_reference_quad)

  # matrix powers in eigen mode
  out_reference = V %*% z
  expect_equal(sk_var_mult(g_miss, pars, p=1), out_reference)
  expect_equal(sk_var_mult(g_miss, pars, p=1, quad=TRUE), t(z) %*% out_reference)
  expect_equal(sk_var_mult(g_miss, pars, p=2), V %*% out_reference)

  # verify that multiplying g_miss twice by a square root of V is same as multiplying by V
  g_miss_sqrt = g_miss
  g_miss_sqrt[!is.na(g_miss)] = sk_var_mult(g_miss, pars, p=1/2)
  expect_equal(sk_var_mult(g_miss_sqrt, pars, p=1/2), out_reference)
})

# sk_toep_mult
test_that("sk_toep_mult calculates the expected product", {

  # define example matrix from 1D exponential variogram
  n = 10
  y = exp(1-seq(n))
  y_mat = sk_toep_mult(y)
  expect_equal(y_mat, stats::toeplitz(y))

  # multiply by random matrix and compare with default matrix multiply
  z = matrix(stats::rnorm(n^2), n)
  result_default = y_mat %*% z
  expect_equal(result_default, sk_toep_mult(y_mat, z))

  # save memory by passing only the first row of the Toeplitz matrix
  expect_equal(result_default, sk_toep_mult(y, z))

  # sparsify z and repeat
  idx_sparse = sample.int(n^2, n^2 - n)
  z[idx_sparse] = 0
  result_default = y_mat %*% z
  expect_equal(result_default, sk_toep_mult(y, z))

  # right-multiply with another kernel
  x = exp( 2 *( 1-seq(n) ) )
  x_mat = sk_toep_mult(x)
  result_default = result_default %*% x_mat
  expect_equal(result_default, sk_toep_mult(y, z, x))

  # z can also be supplied as vector of nonzero grid values
  idx_obs = which(z != 0)
  gdim = c(y=n, x=n)
  expect_equal(result_default, sk_toep_mult(y, z=z[idx_obs], x, idx_obs, gdim))
})

