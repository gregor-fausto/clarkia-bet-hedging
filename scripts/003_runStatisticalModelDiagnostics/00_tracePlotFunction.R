####
####
# Function to produce trace plots as JPEGs for each parameter, by population
####
####

# function to produce trace plots as JPEGs per population
f = function(x = "parm", model = "modelName", jpeg.quality = 75) {
  
  # get chains for parameter as a list
  chains <- MCMCchains(mcmcSamples, params = x, mcmc.list = TRUE)
  chains.list <- lapply(chains, as.matrix)
  
  n = (dim(chains.list[[1]])[2])/20
  n = ceiling(n)
  counter <- 1
  parm = x
  while (counter <= n) {
    
    # write convergence plot to a directory
    jpeg(filename = paste0(outputDirectory, "convergence-", model, "-", parm, 
                           "-", counter, ".jpeg"), quality = 100)
    # name of parameter
    par.names = colnames(chains[[1]])
    
    # plot dimensions
    par(mfrow = c(4, 5), oma = c(5, 4, 0, 0) + 0.1, mar = c(0, 0, 1, 1) + 0.1)
    
    # run over populations
    for (i in 1:20) {
      
      if((20 * (counter - 1) + i)>dim(chains[[1]])[2]){
        plot(NA,xlim=c(0,1),ylim=c(0,1),axes=FALSE)
        
      } else {
      
      plot(as.vector(chains[[1]][, (20 * (counter - 1) + i)]), type = "n", 
           axes = FALSE, frame = FALSE, xaxt = "n", yaxt = "n", xlab = "", ylab = "")
      
      # each of three chains
      lines(as.vector(chains[[1]][, (20 * (counter - 1) + i)]), col = rgb(0.9882353, 0.5529412, 0.3490196, 0.5))
      lines(as.vector(chains[[2]][, (20 * (counter - 1) + i)]), col = rgb(1,1, 0.75, 0.5))
      lines(as.vector(chains[[3]][, (20 * (counter - 1) + i)]), col = rgb(0.5686275, 0.7490196, 0.8588235, 0.5))
      
      # add parameter name to jpeg
      text(0.5 * length(as.vector(chains[[1]][, (20 * (counter - 1) + i)])), 
           0.9 * max(as.vector(chains[[1]][, (20 * (counter - 1) + i)])), par.names[(20 * (counter - 1) + i)])
      }
    }
    dev.off()
    counter <- counter + 1
  }
  
}