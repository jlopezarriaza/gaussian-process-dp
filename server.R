#server.R
#MultiObjective
source('helpers.R')
# source("GP_Embed.R")
library('ggplot2')
library('akima')
library('MASS')
library('reshape2')
library('shinyjs')
shinyServer(server <- function(input, output) {
  ###########################################################################
  observeEvent(input$submit_loc1, {
    show("submit_loc2")
    show("submit_loc4")
    disable("time.series.length")
    disable("V")
    disable("k_1")
    disable("k_2")
    disable("c1")
    disable("c2")
    disable("r_1")
    disable("r_2")
    disable("num")
  })
  ###########################################################################
  observeEvent(input$submit_loc2, {
    show("submit_loc3")
  })
  ###########################################################################
  observeEvent(input$submit_loc4, {
    enable("submit_loc1")
    enable("time.series.length")
    enable("V")
    enable("k_1")
    enable("k_2")
    enable("c1")
    enable("c2")
    enable("r_1")
    enable("r_2")
    enable("num")
    reset("time.series.length")
    reset("V")
    reset("k_1")
    reset("k_2")
    reset("c1")
    reset("c2")
    reset("r_1")
    reset("r_2")
    reset("num")
    hide("submit_loc2")
    hide("submit_loc3")
    toggle("submit_loc4")
  })
  ###########################################################################
  output$input_ui <- renderUI({
    num <- as.integer(input$num)
    lapply(1:num, function(i) {
      sliderInput(paste0("k_", i),
                  label = withMathJax(paste0("$$k_", i,"$$")),
                  value = 0.5,
                  min = 0,
                  max = 2,
                  step = 0.25)
    })
  })
  ###########################################################################
  output$input_ui2 <- renderUI({
    num <- as.integer(input$num)
    lapply(1:num, function(i) {
      sliderInput(paste0("r_", i),
                  label = paste0("$$r_", i,"$$"), 
                  value = 2.5,
                  min = 1,
                  max = 3,
                  step = 0.25)
    })
  })
  ###########################################################################
  output$input_ui3 <- renderUI({
    num <- as.integer(input$num)
    lapply(1:num, function(i) {
      sliderInput(paste0("c_", i), 
                  label = withMathJax(paste0("$$c_", i,"$$")), 
                  value = 1,
                  min = 0,
                  max = 3,
                  step = 0.25)
    })
  })
  ###########################################################################
  simulate.data <- eventReactive(input$submit_loc1,{
    shinyjs::disable("submit_loc1")
    shinyjs::disable("submit_loc2")
    shinyjs::disable("submit_loc3")
    shinyjs::disable("submit_loc4")
    
    time.series.length = input$time.series.length
    V = input$V
    num <- as.integer(input$num)
    k_ricker = c(sapply(1:num, function(i) {
      input[[paste0("k_", i)]]
    }))
    r_ricker = c(sapply(1:num, function(i) {
      input[[paste0("r_", i)]]
    }))
    competition.parameters = c(input$c1, input$c2)
    
    x = matrix(data = NA, nrow = time.series.length + 1, ncol = num)
    x[1,] = k_ricker * runif(num, min = 1, max = 3)
    for (t in 1:(time.series.length)) {
      x[t + 1,] = generating.models(x[t,], r_ricker, k_ricker, competition.parameters) * 
        exp(rnorm(n = num, mean = 0, sd = V ^ 0.5))
    }
    withProgress(message = 'Generating Data and Fitting GP ...', value = 1, {
      gp.fit = fit.gaussian.process(x)
    })
    shinyjs::enable("submit_loc2")
    shinyjs::enable("submit_loc4")
    list(abundance = x, 
         K = k_ricker, 
         R = r_ricker, 
         num = num,
         gp.fit = gp.fit,
         V = V,
         time.series.length = time.series.length, 
         competition.parameters = competition.parameters)
  })
  ###########################################################################
  dynamic.programming <- eventReactive(input$submit_loc2,{
    shinyjs::disable("submit_loc1")
    shinyjs::disable("submit_loc2")
    shinyjs::disable("submit_loc3")
    shinyjs::disable("submit_loc4")
    
    simulated.data.output = simulate.data()
    
    x = simulated.data.output[['abundance']]
    K = simulated.data.output[['K']]
    R = simulated.data.output[['R']]
    num = simulated.data.output[['num']]
    V = simulated.data.output[['V']]
    
    gp.fit = simulated.data.output$gp.fit
    withProgress(message = 'Finding Optimal Harvest ...', value = 1, {
      dynamic.programming.fit = markov.decision.process(gp.fit)
    })
    gp.fit$K = K
    gp.fit$R = R
    gp.fit$num = num
    gp.fit$V = V
    dynamic.programming.results = append(gp.fit, dynamic.programming.fit)
    shinyjs::enable("submit_loc3")
    shinyjs::enable("submit_loc4")
    
    
    list(abundance = x, 
         K = K,
         R = R, 
         num = simulated.data.output$num, 
         dynamic.programming.results = dynamic.programming.results, 
         competition.parameters = simulated.data.output$competition.parameters)
  })
  
  ###########################################################################
  
  forward.simulations <- eventReactive(input$submit_loc3,{
    shinyjs::disable("submit_loc1")
    shinyjs::disable("submit_loc2")
    shinyjs::disable("submit_loc3")
    shinyjs::disable("submit_loc4")
    withProgress(message = 'Doing Forward Simulations...', value = 1, {
      dynamic.programming.output = dynamic.programming()
    })
    dynamic.programming.results = dynamic.programming.output$dynamic.programming.results
    dynamic.programming.results$competition.parameters = c(input$c1, input$c2)
    Y = forward.iterations(dynamic.programming.results)
    #shinyjs::enable("submit_loc1")
    #shinyjs::enable("submit_loc2")
    #shinyjs::enable("submit_loc3")
    shinyjs::enable("submit_loc4")
    list(dynamic.programming.results = dynamic.programming.results,
         Y = Y)
  })
  
  
  
  ###########################################################################
  
  output$graph1 <- renderPlot({
    simulated.data.output = simulate.data()
    x = simulated.data.output[["abundance"]]
    K = simulated.data.output[['K']]
    num = simulated.data.output[['num']]
    D = as.data.frame(cbind(x, matrix(c(1 : dim(x)[1]))))
    names(D) <- c(1 : dim(x)[2],'Year')
    df_molten = melt(D, id.vars = "Year")
    ggplot(df_molten) +
      geom_line(aes(x = Year,
                    y = value,
                    color = variable),
                size = 2) +
      scale_colour_manual(name  = "",
                          values = c("firebrick", 'royalblue4', "#00FF00"),
                          labels = c('Species 1', 'Species 2', 'Species 3')) +
      ggtitle('Time Series') +
      labs(x = 'Year', y = 'Abundance') +
      theme(axis.text = element_text(size = 16),
            axis.title = element_text(size = 22),
            plot.title = element_text(size = 22),
            legend.text = element_text(size = 20))
  })
  ###########################################################################
  
  output$graph2 <- renderPlot({
    simulated.data.output = simulate.data()
    x = simulated.data.output[["abundance"]]
    num = simulated.data.output[['num']]
    time.series.length = simulated.data.output[['time.series.length']]
    v = list()
    for (i in 1 : num) {
      D = as.data.frame(cbind(x[1:(time.series.length - 1),i],x[2:(time.series.length),i]))
      names(D) <- c('hist', 'obs')
      v[[i]] <- ggplot(D, aes(x = hist, y = obs)) +
        geom_point(color = 'firebrick', size = 5) +
        labs(x = substitute('x'[a][','][t], list(a = i)),
             y = substitute('x'[a][','][t + 1], list(a = i)),
             title = '1-D Map') +
        theme(axis.text = element_text(size = 16),
              axis.title = element_text(size = 22),
              plot.title = element_text(size = 22),
              legend.text = element_text(size = 20))
    }
    if (num == 1)
      multiplot(v[[1]], cols = 1)
    if (num == 2)
      multiplot(v[[1]], v[[2]], cols = 2)
    if (num == 3)
      multiplot(v[[1]], v[[2]], v[[3]], cols = 3)
  })
  ###########################################################################
  
  output$graph3 <- renderPlot({
    dynamic.programming.outputs = dynamic.programming()
    dynamic.programming.results = dynamic.programming.outputs[["dynamic.programming.results"]]
    num = dynamic.programming.outputs[['num']]
    if (num == 1) {
      FitOutput = data.frame(mu = exp(dynamic.programming.results$mu[[1]][,1]))
      FitOutput$x = dynamic.programming.results$state.grid
      
      G = matrix(NA, nrow = length(dynamic.programming.results$state.grid), ncol = 1)
      
      for (dd in 1:length(dynamic.programming.results$state.grid)) {
        G[dd] = generating.models(dynamic.programming.results$state.grid[dd], dynamic.programming.outputs[['R']], dynamic.programming.outputs[['K']])
      }
      RealModel = as.data.frame(cbind(dynamic.programming.results$state.grid,G))
      names(RealModel) = c('x','y')
      
      FitOutput$lower = exp(dynamic.programming.results$mu[[1]][,1] - 2 * sqrt(diag(dynamic.programming.results$sig) + dynamic.programming.results$var[[1]]))
      FitOutput$upper = exp(dynamic.programming.results$mu[[1]][,1] + 2 * sqrt(diag(dynamic.programming.results$sig) + dynamic.programming.results$var[[1]]))
      FitOutput$harvest = dynamic.programming.results$harvest * dynamic.programming.results$state.grid
      DataIn = data.frame(hist = dynamic.programming.results$hist[,1],
                          obs = dynamic.programming.results$obs[,1])
      
      ggplot(FitOutput) +
        geom_line(aes(x = x, y = mu, color = 'Mu')) +
        geom_line(aes(x = x, y = harvest , color = 'Harvest'),
                  size = 2) +
        geom_ribbon(aes(x = x, ymin = lower, ymax = upper, fill = "GP Fit"),
                    alpha = 0.25) +
        geom_point(data = DataIn, aes(x = hist, y = obs, color = 'Data'),
                   size = 5) +
        geom_line(aes(x = RealModel$x, y = RealModel$y, color='True Model')) +
        labs(x = expression(x[1][','][t]), 
             y = expression(x[1][','][t + 1]),
             title = 'Optimization') +
        theme(axis.text = element_text(size = 16),
              axis.title = element_text(size = 22),
              plot.title = element_text(size = 22),
              legend.text = element_text(size = 20)) +
        scale_colour_manual(name = '' ,
                            values = c("Data" = "black",
                                       "Harvest" = '#009E73',
                                       'Mu' = 'red',
                                       'True Model'= 'blue')) +
        scale_fill_manual('', values = '#E69F00') +
        guides(color = guide_legend(override.aes = 
                                      list(shape = c(19, NA, NA, NA),
                                           linetype = c(0, 1, 1, 1)))) 
    } else if (num == 2) {
      v = list()
      for (i in 1:2) {
        D = dynamic.programming.results$interpolated$dynamics[[i]]$z
        D3d <- melt(D)
        names(D3d) <- c("x", "y", "z")
        D3d$x = dynamic.programming.results$interpolated$dynamics[[i]]$x[D3d$x]
        D3d$y = dynamic.programming.results$interpolated$dynamics[[i]]$y[D3d$y]
        v[[i]] <- ggplot(D3d, aes(x, y, z = z)) +
          geom_tile(aes(fill = z)) +
          scale_fill_gradientn(colours = jet.cols) +
          guides(fill = guide_colorbar(barwidth = 0.5, barheight = 10)) +
          labs(x = expression(x[1][','][t]),
               y = expression(x[2][','][t]),
               fill = substitute('x'[a][','][t + 1], list(a = i)),
               title = substitute('Transition Dynamics: x'[a], list(a = i))) +
          theme(axis.text = element_text(size = 16),
                axis.title = element_text(size = 22),
                plot.title = element_text(size = 22),
                legend.text = element_text(size = 20),
                legend.title=element_text(size = 20))
      }
      multiplot(v[[1]], v[[2]], cols = 2)
    }
  })
  ###########################################################################
  
  output$graph5 <- renderPlot({
    dynamic.programming.output = dynamic.programming()
    dynamic.programming.results = dynamic.programming.output[["dynamic.programming.results"]]
    num = dynamic.programming.output[['num']]
    if (num == 1) {
      FitOutput = data.frame(mu = exp(dynamic.programming.results$mu[[1]][,1]))
      FitOutput$x = dynamic.programming.results$state.grid
      
      G = matrix(NA,nrow = length(dynamic.programming.results$state.grid),ncol = 1)
      
      
      for (dd in 1:length(dynamic.programming.results$state.grid)) {
        G[dd] = generating.models(dynamic.programming.results$state.grid[dd],
                       dynamic.programming.output[['R']],
                       dynamic.programming.output[['K']])
      }
      RealModel = as.data.frame(cbind(dynamic.programming.results$state.grid,G))
      names(RealModel) = c('x','y')
      
      FitOutput$lower = exp(dynamic.programming.results$mu[[1]][,1] -
                              2 * sqrt(diag(dynamic.programming.results$sig) + 
                                         dynamic.programming.results$var[[1]]))
      FitOutput$upper = exp(dynamic.programming.results$mu[[1]][,1] +
                              2 * sqrt(diag(dynamic.programming.results$sig) + 
                                         dynamic.programming.results$var[[1]]))
      FitOutput$harvest = dynamic.programming.results$harvest * dynamic.programming.results$state.grid
      DataIn = data.frame(hist = dynamic.programming.results$hist[,1])
      DataIn$obs = dynamic.programming.results$obs[,1]
      
      ggplot(FitOutput) +
        geom_line(aes(x = x,y = mu,color = 'Mu')) +
        geom_line(aes(x = x,y = harvest ,color = 'Harvest'),size = 2) +
        geom_ribbon(aes(
          x = x,ymin = lower,ymax = upper,fill = "GP Fit"
        ),alpha = .21) +
        geom_point(data = DataIn,aes(x = hist, y = obs,color = 'Data'),size =
                     5) +
        labs(x = expression(x[1][','][t]),
             y = expression(x[1][','][t + 1]),title = 'Optimization') +
        theme(
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 22),
          plot.title = element_text(size = 22),
          legend.text = element_text(size = 20),
          legend.title=element_text(size = 20)
        ) +
        scale_colour_manual(
          name = '' ,
          values = c(
            "Data" = "black","Harvest" = '#009E73','Mu' = 'red'
          )
        ) +
        scale_fill_manual('',values = '#E69F00') +
        guides(color = guide_legend(override.aes = list(shape = c(19, NA,NA),
                                                        linetype = c(0,1,1)))) +
        geom_line(aes(x = RealModel$x, y = RealModel$y)) + 
        geom_line(aes(x = RealModel$x,y = RealModel$x))
      
    }else if (num == 2) {
      vv = list()
      for (i in 1:2) {
        D = dynamic.programming.results$interpolated$policy[[i]]$z
        D3d <- melt(D)
        names(D3d) <- c("x", "y", "z")
        D3d$x = dynamic.programming.results$interpolated$policy[[i]]$x[D3d$x]
        D3d$y = dynamic.programming.results$interpolated$policy[[i]]$y[D3d$y]
        vv[[i]] <- ggplot(D3d, aes(x, y, z = z)) +
          geom_tile(aes(fill = z)) +
          #stat_contour(bins=1,aes(x,y,z=z), color="black", size=0.6)+
          scale_fill_gradientn(colours = jet.cols) +
          guides(fill = guide_colorbar(barwidth = 0.5, barheight = 10)) +
          labs(x = expression(x[1][','][t]),
               y = expression(x[2][','][t]),
               fill = 'harvest',
               title = substitute('Optimal harvest: x'[a], list(a = i))
          ) +
          theme(axis.text = element_text(size = 16),
                axis.title = element_text(size = 22),
                plot.title = element_text(size = 22),
                legend.text = element_text(size = 20)
          )
      }
      multiplot(vv[[1]],vv[[2]],cols = 2)
      
    }
    
  })
  ###########################################################################
  
  output$graph4 <- renderPlot({
    forward.simulation.output = forward.simulations()
    X = forward.simulation.output$Y
    dynamic.programming.results = forward.simulation.output$dynamic.programming.results
    vv = list()
    steady = matrix(NA, ncol = dim(X$X)[2])
    if (dim(X$X)[2] == 1)
      steady = dynamic.programming.results$K[1] / dynamic.programming.results$R[1]
    if (dim(X$X)[2] == 2){
      steady[1] = NaN
      steady[2] = NaN
    }
    for (i in 1 : dim(X$X)[2]) {
      M = as.data.frame(cbind(X$X[,i,], seq(1:dim(X$X)[1])))
      num_name = paste("x", 1:dim(X$X)[3], sep = "")
      names(M) = c(num_name, 'Year')
      M_melt = melt(M, 'Year')
      vv[[i]] <- ggplot(M_melt) +
        geom_line(aes(x = Year, y = value, group = variable),
                  color = 'coral2',
                  alpha = 0.5) +
        guides(colour = FALSE) +
        ylim(0, max(X$X) + 0.5*max(X$X)) +
        geom_hline(yintercept = steady[i])+
        ggtitle('Monte Carlo Simulation') +
        labs(x = 'Year', y = 'Abundance')+
        theme(axis.text = element_text(size = 16),
              axis.title = element_text(size = 22),
              plot.title = element_text(size = 22),
              legend.text = element_text(size = 20))
    }
    v = list()
    for (i in 1 : dim(X$X)[2]) {
      M = as.data.frame(cbind(X$H[,i,], seq(1:dim(X$H)[1])))
      num_name = paste0("x", 1 : dim(X$H)[3])
      names(M) = c(num_name, 'Year')
      M_melt = melt(M, 'Year')
      v[[i]] <- ggplot(M_melt) +
        geom_line(aes(x = Year, y = value, group = variable),
                  color = 'firebrick', 
                  alpha = 0.5) +
        ylim(0, max(X$H)) +
        guides(colour = FALSE) +
        ggtitle('Monte Carlo Simulation') +
        labs(x = 'Year', y = 'Harvest') +
        theme(axis.text = element_text(size = 16),
              axis.title = element_text(size = 22),
              plot.title = element_text(size = 22),
              legend.text = element_text(size = 20))
    }
    if (dim(X$X)[2] == 1)
      multiplot(vv[[1]], v[[1]], cols = 2)
    else
      multiplot(vv[[1]], v[[1]], vv[[2]], v[[2]], cols = 2)
    
  })
  
})