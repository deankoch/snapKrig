# randomized grid size

# sk_rescale
test_that("invert sk_rescale (going down then up) and check factor 1 case", {

  # example data
  gdim = seq(1e2) |> sample(2)
  pars = utils::modifyList(sk_pars(gdim), list(eps=1e-2))
  g = sk_sim(gdim, pars)

  # factor 1 rescale does nothing

  # downscale
  sk_rescale(g, down=1) |> expect_equal(g)

  # upscale - temporary: initialize "crs", missing from sk_rescale(up=...)
  sk_rescale(g, up=1) |> utils::modifyList(list(crs=NULL)) |> expect_equal(g)
  sk_rescale(sk_rescale(g, down=c(5,3)), up=c(5,3)) |>
    utils::modifyList(list(crs=NULL)) |>
    expect_equal(g)
})

# sk_mat2vec
test_that("verify sk_mat2vec computes vectorized index wrt expand.grid", {

  # pick a random size grid and random index
  gdim = seq(1e3) |> sample(2)
  ij = c(i=sample(seq(gdim[1]), 1), j=sample(seq(gdim[2]), 1))

  # matrix indices in column-vectorized order
  gyx = expand.grid(i=seq(gdim[1]), j=seq(gdim[2]))
  result = sk_mat2vec(gyx, gdim)
  as.vector(gyx[["i"]] + gdim[1] * (gyx[["j"]] - 1)) |> expect_equal(result)
  gyx |> as.matrix() |> expect_equal(sk_vec2mat(result, gdim))
})

# sk_vec2mat
test_that("verify round trip with sk_vec2mat -> sk_mat2vec", {

  # pick a random size grid and random index
  gdim = seq(1e3) |> sample(2)
  idx = prod(gdim) |> seq() |> sample(1)
  sk_vec2mat(idx, gdim) |> sk_mat2vec(gdim) |> expect_equal(idx)
})

# sk_sub_idx
test_that("verify sk_sub_idx satisfies basic identities wrt expand.grid", {

  # pick a random size grid and a particular grid point
  gdim = seq(1e3) |> sample(2)
  ij_list = c(i=sample(seq(gdim[1]), 1), j=sample(seq(gdim[2]), 1)) |> as.list()

  # a randomly selected sub-grid with top left corner at ij_list
  origin_sg = unlist(ij_list)
  sg_max = gdim - unlist(origin_sg) + 1
  gdim_sg = c(i=sample(seq(sg_max[1]), 1), j=sample(seq(sg_max[2]), 1))

  # sk_sub_idx returns a logical vector indexing the point (or the index itself)
  is_pt = sk_sub_idx(gdim, ij_list)
  idx_sub = sk_sub_idx(gdim, ij_list, idx=TRUE)

  # equivalent call when ij_list is a single point
  sk_mat2vec(ij_list, gdim) |> expect_equal(idx_sub)

  # if i or j are omitted from ij, the function returns the full row or column
  idx_i = expand.grid(i=ij_list[['i']], j=seq(gdim[2])) |> sk_mat2vec(gdim)
  idx_j = expand.grid(i=seq(gdim[1]), j=ij_list[['j']]) |> sk_mat2vec(gdim)
  sk_sub_idx(gdim, ij_list['i']) |> which() |> expect_equal(idx_i)
  sk_sub_idx(gdim, ij_list['j']) |> which() |> expect_equal(idx_j)

  # indices in column-vectorized order
  sk_sub_idx(gdim, ij_list['i']) |> which() |> expect_equal(sk_sub_idx(gdim, ij_list['i'], idx=TRUE))
  sk_sub_idx(gdim, ij_list['j']) |> which() |> expect_equal(sk_sub_idx(gdim, ij_list['j'], idx=TRUE))
  prod(gdim) |> seq() |> expect_equal(sk_sub_idx(gdim, idx=TRUE))

  # bigger sub-grid example
  ij_sg_list = list(i = origin_sg[1] + seq(gdim_sg[1]) - 1, j = origin_sg[2] + seq(gdim_sg[2]) - 1)
  is_sg = sk_sub_idx(gdim, ij=ij_sg_list)
  idx_sg = sk_sub_idx(gdim, ij=ij_sg_list, idx=TRUE)

  # example with j indices supplied in reverse (descending) order
  ij_list_xflip = utils::modifyList(ij_sg_list, list(j=rev(ij_sg_list[['j']])))

  # ordering in ij$i and ij$j doesn't matter `nosort=FALSE` or `idx=FALSE`
  sk_sub_idx(gdim, ij=ij_sg_list, nosort=TRUE) |> expect_equal(is_sg)
  sk_sub_idx(gdim, ij=ij_list_xflip, idx=TRUE) |> expect_equal(which(is_sg))

  # when `nosort=TRUE` and `idx=TRUE` we get the same indices but in a different order
  idx_sg_xflip = sk_sub_idx(gdim, ij=ij_list_xflip, idx=TRUE, nosort=TRUE)
  sort(idx_sg) |> expect_equal(sort(idx_sg_xflip))
})

# sk_sub
test_that("sk_sub indexes subgrids as expected on random input", {

  # make an example grid with at minimum 10 grid lines in each dimension
  gdim = seq(1e2) |> tail(-9) |> sample(2)
  g = sk(gdim)
  g[] = apply(expand.grid(g[['gyx']]), 1, \(z) cos( 2*sum(z^2) ) )

  # randomly selected bottom-right corner for the sub-grid
  imax = seq(gdim[1]) |> tail(-3) |> sample(1)
  jmax = seq(gdim[2]) |> tail(-3) |> sample(1)

  # subset by specifying grid lines to keep
  ij_keep = list(i=seq(1, imax, by=2), j=seq(1, jmax, by=2))
  g_keep = sk_sub(g, ij_keep)

  # get the indices kept and removed
  idx = sk_sub(g, ij_keep, idx=TRUE)

  # equivalent call specifying grid lines to omit
  sk_sub(g, ij_rem=idx[['rem']]) |> expect_equal(g_keep)

  # edge lines to trim
  n_rem = seq(9) |> sample(2) |> lapply(seq)

  # remove data around the edges of the grid
  idx = sk_sub(g, ij_rem=list(i=n_rem[[1]], j=n_rem[[2]]), mirror=TRUE, idx=TRUE)
  idx_y_pts = sk_sub_idx(gdim, idx[['rem']]['i'], idx=TRUE)
  idx_x_pts = sk_sub_idx(gdim, idx[['rem']]['j'], idx=TRUE)
  g[c(idx_y_pts, idx_x_pts)] = NA
  # !! next line produces an error if we omit seq above in n_rem def
  g_sub = sk_sub(g)
  # (non-uniqueness not being handled maybe?)
  # TODO: fix this bug

  # verify interior sub-grid is as large as expected
  idx[['rem']] |>
    sapply(length) |>
    setNames(c("y", "x")) |>
    expect_equal(dim(g) - dim(g_sub))
})

# sk_sub_find
test_that("sk_sub_find finds randomly selected regular sub-grids", {

  # make an example grid with at minimum 10 grid lines in each dimension
  gdim = seq(1e2) |> tail(-9) |> sample(2)
  pars = utils::modifyList(sk_pars(gdim), list(eps=1e-12))

  # generate some random data
  g = sk_sim(gdim, pars)

  # define a super-grid containing the original data and make sure we can find it
  g_big = sk_rescale(g, down=3)
  sk_sub_find(g_big) |> is.null() |> expect_false()

  # define a smaller sub-grid at random
  spacing = sapply(floor(gdim/10), function(x) 1 + sample.int(x, 1))
  gdim_sg = sapply(floor( (gdim - 1) / spacing), function(x) sample.int(x, 1))
  ij_first = sapply(gdim - ( spacing * gdim_sg ), function(x) sample.int(x, 1))

  # find index of sub-grid lines and vectorized index of points
  ij_sg = Map(function(idx, r, n) seq(idx, by=r, length.out=n), idx=ij_first, r=spacing, n=gdim_sg)
  names(ij_sg) = c('i', 'j')
  is_sg = sk_sub_idx(gdim, ij_sg, idx=FALSE)

  # assign values to the sub-grid points
  g_sub = sk(gdim)
  g_sub[is_sg] = g[is_sg]

  # find it!
  sub_result = sk_sub_find(g_sub)
  unname(sub_result[['gdim']]) |> expect_equal(gdim_sg)
  unname(sub_result[['ij']]) |> expect_equal(unname(ij_sg))

  # sub grids with side length 1 have no spacing defined along that dimension
  spacing[gdim_sg==1] = NA

  # check consistency in spacing
  unname(sub_result[['res_scale']]) |> expect_equal(spacing)

  # can also call on the vector and supply gdim separately
  sk_sub_find(!is.na(g_sub), dim(g_sub)) |> expect_equal(sub_result)
})
