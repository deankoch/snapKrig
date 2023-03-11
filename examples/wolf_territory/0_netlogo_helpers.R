# run_netlogo.R
# helper functions to run NetLogo simulations in R
# Dean Koch, Feb 2023

# This defines the two main functions for the simulation experiment:
# run_netlogo executes the simulations (via RNetLogo package), and
# run_analysis takes the resulting territory maps and does the sampling
# study to find bias at different resolutions.

# requires 'terra' and 'RNetLogo' packages, and an installation of NetLogo
# v6.0.4 in the standard location for Windows 10 (check `nl_path` below)

# Note there is a bug preventing 'RNetLogo' from working with newer versions
# of NetLogo. See https://github.com/NetLogo/NetLogo/issues/1782
# For this reason we are using NetLogo v6.0.4, instead of the newer v6.1.1
# recommended by Sells et al. There doesn't seem to be any difference
# in model behaviour or results and in testing, simulations complete
# without error.

# Before running `run_netlogo()` make sure that `proj_dir` contains
# the "Model/Wolf_Territory_Model.nlogo" project file, as well as subfolders
# "Model/Model_outputs" (which can be empty) and "Model/Model_inputs" (with
# covariates). The easiest way to do this is to simply unzip "Model.zip"
# from the Sells et al data supplement to a `proj_dir` of your choice. This
# is done automatically in run_netlogo.R 

# Runs the NetLogo simulator to get a random map of wolf territories (returns nothing)
run_netlogo = function(proj_dir) {
  
  # NetLogo v6.0.4 must be installed separately to this directory
  nl_path = 'C:/Program Files/NetLogo 6.0.4/app'
  nl_name = 'netlogo-6.0.4.jar'
  
  # NetLogo code (called later) reporting the output file name of interest
  asc_report = '(word "Tsize_owners, " word behaviorspace-run-number", "run-ID", "file-id", "where-model",.asc")'
  
  # projection used in grid files
  asc_epsg = 'EPSG:5070'
  
  # workaround for unsupported NetLogo version in project file
  # copy project file and manually revert version number to prevent error
  nl_proj = proj_dir |> file.path('Model/Wolf_Territory_Model.nlogo')
  nl_proj_604 = proj_dir |> file.path('Model/Wolf_Territory_Model_v604.nlogo')
  file.copy(nl_proj, nl_proj_604, overwrite=TRUE)
  gsub('NetLogo 6.1.1', 'NetLogo 6.0.4', readLines(nl_proj_604)) |> writeLines(nl_proj_604)
  
  # create an instance of NetLogo and load the modified project file
  nl_instance = RNetLogo::NLStart(nl_path, gui=FALSE, nl.jarname=nl_name)
  RNetLogo::NLLoadModel(nl_proj_604, nl_instance)
  
  # disable graphics for speed then initialize the model
  RNetLogo::NLCommand('set depict-landscape? false')
  RNetLogo::NLCommand('setup')
  
  # simulations end when this target number of packs established
  npack_target = RNetLogo::NLReport('total-target')
  
  # timer for the whole simulation
  t_start = Sys.time()
  
  # we will evaluate simulation in chunks, checking in periodically to see progress
  iter_chunk = 1e4
  pb = utils::txtProgressBar(0, npack_target, style=3)
  
  # this would loop indefinitely, but an error at the end of the simulation will break it
  cat(paste('simulating', npack_target, 'packs...\n'))
  err_msg = tryCatch({
    while(TRUE) {
      
      # run simulation chunk quietly (if completed, this produces a file not found error)
      RNetLogo::NLDoCommand(iter_chunk, 'go') |> capture.output()
      
      # update user on progress
      npack_current = RNetLogo::NLReport('count packs')
      setTxtProgressBar(pb, npack_current)
      #RNetLogo::NLGetPatches('packs-already-here', as.matrix=TRUE) |> snapKrig::sk() |> plot()
    }
  }, error=\(e) invisible(e), finally=cat('\nsimulation complete.'))
  
  # print time elapsed
  Sys.time() - t_start
  
  # make a unique file name for output raster and replace problematic punctuation
  file_out = paste0('territory_sim_',  as.character(Sys.time()), 'PST.tif')
  file_out = gsub(' ', '_at_', file_out)
  file_out = gsub(':', '', file_out)
  
  # identify the output raster filename and load it as SpatRaster
  output_asc = RNetLogo::NLReport(asc_report)
  output_grid = terra::rast(file.path(proj_dir, 'Model/Model_outputs', output_asc))
  terra::crs(output_grid) = asc_epsg
  names(output_grid) = 'tsize'
  
  # save as tiff
  terra::writeRaster(output_grid, file.path(proj_dir, file_out))
  
  # close NetLogo
  RNetLogo::NLQuit()
}


# Once run_netlogo is finished, you should have a .tif file in your proj_dir
# delineating the territories. You can call run_netlogo multiple times to get
# multiple repetitions of this data. Load all into R with and convert to a 
# logical sk object by doing `g_source = sk(terra::rast(...)) > 0`. Pass this
# g_source to run_analysis to check for bias at resolution scalings up to and
# including n_test.

# We align the test grids so that each test grid cell snugly fits a (square)
# integer number of original cells. This means none of the original cells
# overlap with more than one test cell. Any such overlap would increase the
# amount by which we are overestimating, so this is a best case scenario
# for bias.

# Samples occupancy at a range of resolutions to measure bias in territory estimates 
run_analysis = function(g_source, n_test) {
  
  # g_source: logical multilayer sk grid indicating known territory
  # n_test: the maximum up-scaling factor f. Test grids with f = 1, 2, ..., n_test
  
  # returns three things in a list:
  # df: data frame of bias results with test grid info (layer indicated by sim)
  # g_coarse: nested list of estimated occupancy grids
  # g_fine: list of occupancy stats in single layer
  
  # g_coarse will have as many elements as there are repetitions (layers in g_source). Each
  # element is another list, containing a sequence of logical sk grids (single-layer)
  # indicating occupancy at increasingly coarse resolutions. This appears in the same order
  # as the results in df.
  
  # g_fine is structured the same as g_coarse, but each of its elements are at the same
  # resolution as the source, and contain integer grids with values 1 indicating occupancy
  # at the test resolution, and 2 indicating occupancy at both test and source resolutions. 
  
  # extract grid info 
  gdim = dim(g_source)
  true_area = g_source[] |> apply(2, \(x) sum(x > 0)) |> units::set_units('km2')
  
  ## Set up test grids

  # define sequence of grids at coarser resolution to test 
  n_test = 30L
  results_df = data.frame( scale = rev(seq(n_test)) ) |>
    dplyr::mutate( n_y = floor( gdim['y']/scale ) ) |>
    dplyr::mutate( n_x = floor( gdim['x']/scale ) ) |>
    dplyr::mutate( n = n_x * n_y ) |>
    dplyr::mutate( res_y = scale * g_source[['gres']]['y']) |>
    dplyr::mutate( res_x = scale * g_source[['gres']]['x']) |>
    dplyr::mutate( pixel_area = units::set_units(res_x * res_y / 1e6, 'km2') ) |>
    dplyr::mutate( pixel_sum = NA)
  
  # make a copy of this results data frame for each repetition
  n_rep = dim(output_rast)[3]
  results_list = seq(n_rep) |> lapply(\(x) results_df)
  
  # storage for coarse occupancy grids
  g_coarse_list = seq(n_rep) |> lapply(\(i) vector(mode='list', length=n_test))
  for(r in seq(n_rep)) g_coarse_list[[r]][[n_test]] = g_source
  
  # storage for occupancy counts at fine resolution
  g_fine_list = seq(n_rep) |> lapply(\(i) vector(mode='list', length=n_test))
  occupancy_vec_fine = seq(n_rep) |> lapply(\(i) lapply(seq(n_test), \(j) as.integer(g_source[,i])))
  
  # run the experiment - sample at each test resolution in a loop (takes 2 minutes)
  t_start = Sys.time()
  n_all = results_df[['n']] |> head(-1) |> sum()
  pb = utils::txtProgressBar(0, n_all, style=3)
  iteration = 0L
  for(i_test in seq(n_test-1L)) {

    # enumerate all of the test cells and storage for results
    d = results_df[i_test,]
    gdim_test = as.numeric(d[c('n_y', 'n_x')]) |> setNames(c('y', 'x'))
    n_cell = prod(gdim_test)
    n_occupied = rep(0L, n_rep)
    
    # occupancy grid at coarse resolution
    occupancy_vec = seq(n_rep) |> lapply(\(i) logical(n_cell))
    
    # loop over rows in test grid
    for(i in seq(d[['n_y']])) {
      
      # corresponding rows in fine grid
      i_fetch = seq(d[['scale']]) + (i-1) * d[['scale']]
      
      # loop over columns in test grid
      for(j in seq(d[['n_x']])) {
        
        # corresponding columns in fine grid
        j_fetch = seq(d[['scale']]) + (j-1) * d[['scale']]
        
        # extract the fine grid indices for these rows and columns
        k_idx = sk_sub_idx(gdim, ij=list(i=i_fetch, j=j_fetch), idx=TRUE, nosort=TRUE)
        
        # loop over layers
        for(r in seq(n_rep)) {
          
          # check for nonzero cells (territory, in km2) from each repetition
          is_nz = any(g_source[k_idx, r]) 
          
          # update occupancy vectors at test resolution and original resolution
          occupancy_vec[[r]][sk_mat2vec(c(i,j), gdim_test)] = is_nz
          occupancy_vec_fine[[r]][[i_test]][k_idx] = occupancy_vec_fine[[r]][[i_test]][k_idx] + is_nz
        }
        
        # progress info for user
        iteration = iteration + 1L
        utils::setTxtProgressBar(pb, iteration)
      }
    }
    
    # add the pixel sum, ie # of cells counted as "occupied"
    n_occupied = occupancy_vec |> sapply(sum)
    for(r in seq(n_rep)) results_list[[r]][['pixel_sum']][i_test] = n_occupied[r]
    
    # build occupancy grids at coarse resolution
    for(r in seq(n_rep)) g_coarse_list[[r]][[i_test]] = sk(gdim = gdim_test,
                                                           gval = occupancy_vec[[r]])
    
    # build occupancy counts at fine resolution
    for(r in seq(n_rep)) g_fine_list[[r]][[i_test]] = sk(gdim = gdim,
                                                         gval = occupancy_vec_fine[[r]][[i_test]])
  }
  close(pb)
  t_end = Sys.time()
  print(t_end - t_start)
  
  # add estimates and bias columns
  for(r in seq(n_rep)) {
    
    # fill in the final row (test grid is same as source grid - perfect detection)
    results_list[[r]][['pixel_sum']][n_test] = sum(g_source[,r] > 0) 
    
    # compute estimated area for each test grid and bias
    results_list[[r]] = results_list[[r]] |> 
      dplyr::mutate( est_area = pixel_sum * pixel_area ) |>
      dplyr::mutate(percent_over = (est_area - true_area[r]) / true_area[r] ) |>
      dplyr::mutate(percent_over = percent_over |> units::set_units('%') ) 
  }
  
  # merge results data frames from all repetitions
  results_all_df = do.call(rbind, lapply(seq(n_rep), \(i) {
    
    # adds sim column to indicate layer number for the result
    dplyr::mutate(results_list[[i]], sim=i)
  }))
  
  return(list(df = results_all_df,
              g_coarse = g_coarse_list,
              g_fine = g_fine_list))
}


# Makes a line chart of results from run_analysis
plot_bias = function(bias_df, res_plot, ...) {
  
  # a selection of pixel areas to look at more closely below
  res_survey = units::set_units(600L, 'km2')
  res_plot = res_plot |> units::set_units('km2')
  
  ## plot all statistics in a line chart
  
  # initialize the plot
  plot(percent_over~pixel_area, data = bias_df,
       pch = NA,
       ...)
  
  # indicate Sells et al resolution and test resolutions for raster plots
  abline(v=res_survey, lwd=2, lty='dotted', col=adjustcolor('black', alpha.f=0.3))
  abline(v=res_plot, lwd=1, col=adjustcolor('red', alpha.f=0.3))
  
  # plot point results
  points(percent_over~pixel_area, data = bias_df,
         col = adjustcolor('black', alpha.f=0.3))
  
  # color the points
  n_rep = bias_df[['sim']] |> unique() |> length()
  sim_col = hcl.colors(n_rep, 'Dark 3') |> adjustcolor(alpha.f=0.3)
  sim_col_vec = sim_col[bias_df[['sim']]]
  points(percent_over~pixel_area, data = bias_df,
         pch = 16,
         col = sim_col_vec)
  
  # indicate within-repetition results with lines connecting points
  seq(n_rep) |> lapply(\(i) lines(percent_over~pixel_area,
                                  data = dplyr::filter(bias_df, sim==i),
                                  col = sim_col[i])) |> invisible()
}


# Plot wolf territory maps with three panes showing different resolutions
plot_territory = function(bias_results, res_plot, sim=1, ...) {
  
  n_rep = bias_df[['sim']] |> unique() |> length()
  bias_df = bias_results[['df']]
  
  # find the results corresponding to the resolutions in res_plot
  idx_plot = match(res_plot, dplyr::filter(bias_df, sim==sim)[['pixel_area']]) 
  bias_df[idx_plot,]
  
  # plot all res_plot
  close.screen(all.screens=TRUE)
  n_plot = length(idx_plot)
  for(i in split.screen(c(1, n_plot))) {
    
    g_overlap = bias_results[['g_fine']][[sim]][[idx_plot[i]]]
    g_overlap[] = c('extra', 'truth')[ match(g_overlap[], c(1, 2)) ] |> 
      factor(levels=c('truth', 'extra'))
    
    # omit empty white-space on the right
    gdim = dim(g_overlap)
    j_keep = round(0.7*gdim['x']) |> seq()
    g_overlap = g_overlap |> sk_sub(ij_keep=list(i=seq(gdim['y']), j=j_keep))
    
    screen(i)
    newline = '\n'
    desc_string = as.expression(bquote(.(newline)~cell~area:~.(res_plot[i])~km^2))
    plot(g_overlap,
         zlab = 'territory',
         main = desc_string, 
         axes = FALSE,
         xlab = '',
         ylab = '',
         col_box = NA,
         ...)
  }
}
