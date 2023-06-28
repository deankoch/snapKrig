# demo showing kriging-based interpolation to downscale a raster

# INSTRUCTIONS:
#
# Make sure you have the required packages installed. They are
# listed below in the "library(...)" calls. You will also need
# the "raster" package, which is where the Rlogo.png image comes
# from. eg:
#
# install.packages('raster')
#
# Then do the same for the following: ('shiny', 'terra' etc).
# Then run the whole script in your console, or click the
# "Run App" button in Rstudio

library(shiny)
library(terra)
library(dplyr)
library(shinyWidgets)
library(snapKrig)

# setup

# javascript to check for shape in conditionals below
nm_var = c('exp', 'gau', 'sph', 'gxp', 'mat')
has_shape = sapply(nm_var, \(x) length(sk_kp(x)) > 1)
js_shape = paste0("['", paste(nm_var[has_shape], collapse="', '"), "']")

# open example image
r_path = system.file('ex/logo.tif', package='terra')
g = rast(r_path, lyrs=1) |> scale() |> sk()

# upscale it to make the demo more interesting
g = sk_rescale(g, up=2)

# sample semivariogram from g, set max distance for this example
vg = sk_sample_vg(g, n_bin=50) |> filter(d < 30)

# uncomment to fit a mat x mat kernel to the image
#p_fit = sk_fit(g, 'gau', quiet=TRUE) |> sk_pars_update()
p_fit = c(eps=4e-1, psill=0.8, y.rho=3.4, x.rho=3.4)

# define shiny UI
ui = fluidPage(

  titlePanel('snapKrig demo'),

    sidebarPanel(

      # radio buttons for mode go here
      uiOutput('display_mode'),
      hr(),

      # nugget and sill
      sliderInput('nug', 'nugget', min=1e-3, max=5, value=p_fit['eps'], ticks=F),
      sliderInput('psill', 'parial sill', min=1e-2, max=10, value=p_fit['psill'], ticks=F),
      hr(),

      # neither kernel necessarily has a shape param
      uiOutput('x_kernel'),
      conditionalPanel(condition = paste0(js_shape, '.includes(input.x_var)'),
                       sliderInput('x_shape', 'shape',
                                   min=1e-2, max=25, value=2, ticks=F)),

      sliderInput('x_range', 'range',
                  min=1, max=15, value=p_fit['x.rho'], ticks=F),
      hr(),


      uiOutput('y_kernel'),
      conditionalPanel(condition = paste0(js_shape, '.includes(input.y_var)'),
                       sliderInput('y_shape', 'shape',
                                   min=1e-2, max=25, value=2, ticks=F)),

      sliderInput('y_range', 'range',
                  min=1, max=15, value=p_fit['y.rho'], ticks=F),

      width = 3
  ),

  mainPanel(

    fillPage(plotOutput('plot', width='100%', height='750px')),


    fluidRow(align = 'center',
             conditionalPanel(condition = 'input.mode == "interpolation"',
                              sliderInput('ds',
                                          'downscaling factor',
                                          min=1, max=10, value=1, ticks=F)))
  )
)

server = function(input, output, session) {

  # collect the display mode
  output$display_mode = renderUI({ radioButtons('mode',
                                                label='display',
                                                choices=c('interpolation',
                                                          'model'),
                                                selected='model') })

  # collects the chosen kernel names before displaying sliders in UI
  output$x_kernel = renderUI({ shinyWidgets::pickerInput('x_var',
                                                         label='x variogram',
                                                         choices=nm_var,
                                                         selected='gau',
                                                         width='fit') })

  output$y_kernel = renderUI({ shinyWidgets::pickerInput('y_var',
                                                         label='y variogram',
                                                         choices=nm_var,
                                                         selected='gau',
                                                         width='fit') })

  # draw the covariance plot
  output$plot = renderPlot({

    # required input covariance parameters
    req(input$mode, input$ds,
        input$x_var, input$y_var,
        input$nug, input$psill,
        input$x_range, input$y_range,
        input$x_shape, input$y_shape)

    # shape parameter dropped for certain kernels
    y_shape = NULL
    x_shape = NULL
    if( has_shape[input$x_var] ) x_shape = input$x_shape
    if( has_shape[input$y_var] ) y_shape = input$y_shape

    # define parameter list (NULL values dropped by R)
    p_new = c(input$nug, input$psill, input$y_range, y_shape, input$x_range, x_shape)
    g_pars = sk_pars(g, c(x=input$x_var, y=input$y_var)) |> sk_pars_update(p_new)

    # draw the covariance footprint and semivariogram
    if(input$mode == 'model') {

      # two panes arranged vertically
      split.screen(c(2,1))

      # covariance footprint
      screen(1)
      sk_plot_pars(g_pars, g, cex=1.3)

      # semi-variogram plot
      screen(2)
      sk_plot_semi(vg, g_pars)
      close.screen(all.screens=TRUE)
    }

    # draw the input raster and predict a downscaled version
    if(input$mode == 'interpolation')
    {
      # downscale to requested dimensions then interpolate
      g_down = sk_rescale(g, down=input$ds)
      g_pred = sk_cmean(g_down, g_pars, X=0)

      # set common heatmap scale for both plots
      z_lim = range(g_pred) |> c( range(g) ) |> range()

      # set up titles with resolution
      g_res_txt = paste0('(', paste(dim(g), collapse=' x '), ')')
      pred_res_txt = paste0('(', paste(dim(g_pred), collapse=' x '), ')')
      g_title = paste('rlogo.png red channel', g_res_txt)
      ds_title = paste0('kriging prediction at ', input$ds, 'X resolution ', pred_res_txt)
      if(input$ds == 1) ds_title = 'kriging prediction at original resolution'

      # two panes arranged vertically
      split.screen(c(2,1))

      # plot the source image
      screen(1)
      plot(g, zlim=z_lim, pal='Reds', main=g_title, cex=1.3)

      # plot the interpolated downscaled image
      screen(2)
      plot(g_pred, zlim=z_lim, pal='Reds', main=ds_title, cex=1.3)
      close.screen(all.screens=TRUE)
    }
  })
}


shinyApp(ui=ui, server=server)
