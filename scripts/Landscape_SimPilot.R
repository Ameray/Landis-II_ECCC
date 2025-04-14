################################################################################
## organisation ECCC
### Dominic Cyr /Ameray Abderrahmane
# Clear workspace
rm(list = ls())

# Set working directory
wwd <- "D:/Wena_Project/Testing/Landscape/2025-04-14"
setwd(wwd)
wwd <- getwd()
wwd

# Load simulation information
simInfo <- read.csv("simInfo.csv", colClasses = c(simID = "character"))
simDir <- simInfo$simID
simInfo
# Load required libraries
require(parallel)
require(doSNOW)
# Define the number of simultaneous processes and batch size
maxProcesses <- 2  # Adjust based on system capacity
numSimulations <- length(simDir)
numBatches <- ceiling(numSimulations / maxProcesses)
# Loop through batches
for (batch in 1:numBatches) {
  startIndex <- (batch - 1) * maxProcesses + 1
  endIndex <- min(batch * maxProcesses, numSimulations)
  currentSimulations <- simDir[startIndex:endIndex]
  
  # Create a cluster for the current batch
  cl <- makeCluster(min(maxProcesses, length(currentSimulations)), outfile = "")
  registerDoSNOW(cl)
  
  # Run simulations in parallel
  foreach(i = startIndex:endIndex, .packages = c("parallel")) %dopar% {
    setwd(paste(wwd, simInfo[i, "simID"], sep = "/"))
    
    # Read and update README file
    readmeOld <- readLines("README.txt")
    readmeOld <- readmeOld[1:which(grepl(tail(colnames(simInfo), 1), readmeOld))]
    x <- as.character(shell("landis-ii-extensions.cmd list", intern = TRUE))
    
    # Append details to the README file
    sink(file = "README.txt")
    for (l in seq_along(readmeOld)) {
      cat(readmeOld[l], "\n")
    }
    sink()
    
    sink(file = "README.txt", append = TRUE)
    cat("\n#######################################################################\n")
    cat("########### Installed LANDIS-II extensions\n")
    cat("#######################################################################\n")
    for (l in seq_along(x)) {
      cat(x[l], "\n")
    }
    cat("\n#######################################################################\n")
    cat("########### System Info\n")
    cat(write.table(as.data.frame(Sys.info()), quote = FALSE, col.names = FALSE))
    sink()
    # Run the simulation
    shell("landis-ii-7.cmd scenario.txt", wait = TRUE)
    
  }
  
  # Remove no needed  CSV files for our study, 
  
  #for (sim in currentSimulations) {
   #simPath <- file.path(wwd, sim)
   #csvFiles <- c( "log_Flux.csv", "log_FluxDOM.csv", "log_Pools.csv")
    
  # for (csvFile in csvFiles) {
  #    filePath <- file.path(simPath, csvFile)
  # if (file.exists(filePath)) {
   #   file.remove(filePath)
   #}
   #}
  # Stop the cluster after the batch completes
  stopCluster(cl)
  # Print batch completion message
  cat(sprintf("Batch %d/%d completed and log files deleted.\n", batch, numBatches))
  
  # Perform garbage collection to free up memory
  gc()
  
}
  
