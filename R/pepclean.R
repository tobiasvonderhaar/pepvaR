#' Read in MS results from files or data frames
#'
#' \code{MS.load} reads in peptide information from a dataframe located in the workspace, or from a .csv file located in the working directory.
#'
#' @param MSfilename The name of the dataframe or file containing the MS data.
#' 
#' @param varfilename The name of the file containing the variant peptide data (output of make.var())
#'
#' @param scorecutoff The MASCOT score value below which detected peptides are excluded from the analysis (4.9 by default).
#'
#' @param scorecol The number of the column containg the MASCOT scores.
#'
#' @param sequencecol The number of the column containing the sequences of the detected peptides.
#'
#' @param datarange The numbers of the columns containing the spectral counts for the detected peptides (one column per sample)
#' 
#' @param saveresult Indicates whether to output the result to file, defaults to TRUE
#' 
#' @param outfilename The name of the file that will be generated if saveresult is TRUE. The filename will be constructed based on varfilename if not specified.
#'
#'
#'
#' @return
#'
#'
#'
#' @export

MS.load <- function(MSfilename, varfilename, scorecutoff = 4.9, scorecol = 8, sequencecol = 9, datarange = 11:28, saveresult = TRUE, outfilename='default') {

  if(MSfilename == 'MSdata' | MSfilename %in% ls()) {
    eval(parse(text=paste('MSframe <- ', MSfilename, sep = '')))
  } else {
    if(MSfilename %in% dir()) {
      MSframe <- read.csv(MSfilename)
    } else {
      paste('The specified name does not correspond to a variable in the workspace or a file in the working directory. Names are case sensitive.', sep = ''); return()}
  }

  #Convert values from the Score column from factor to numeric. 
  scores <- as.numeric(levels(MSframe[,scorecol])[MSframe[,scorecol]])
  #Remove the column containing factorised scores and replace with the numeric scores, and only retain scores, sequences and data in MSframe.
  MSframe <- cbind(scores, MSframe[,sequencecol], MSframe[datarange])
  colnames(MSframe)[1:2] <- c('MASCOTscore', 'peptides')
  

  #Remove all rows where score is NA from MSframe and from scores
  MSframe <- MSframe[!is.na(scores),]

  #Remove all rows with score <scoreutoff from MSframe
  MSframe <- MSframe[MSframe$MASCOTscore >= scorecutoff,]

  #The ion search sometimes assigns identical spectra to multiple peptides. Search for these instances. If the highest MASCOT score
  #for each group of peptides is >0.8 above the second highest, retain the peptide with the highest score, if not, reject all peptides
  #as a "true" match cannot be identified.
  
  #make a derivative dataframe which only contains the spectral values
  SConly <- MSframe[,3:ncol(MSframe)]
  uniquerows <- rownames(unique(SConly))
  
  #make a copy of MSframe to receive the cleaned up data
  cleanframe <- MSframe[1,]
  cleanframe[1,1] <- 0
  cleanframe[1,2] <- NA
  cleanframe[1, 3:ncol(cleanframe)] <- 0
  
  #go through each row of MSframe
  for(n in 1:nrow(MSframe)) {
    #if the row is unique just add it to cleanframe
    if(rownames(MSframe[n,]) %in% uniquerows) {
      cleanframe <- rbind(cleanframe, MSframe[n,])
    } else {
      #copy all the duplicated rows into a new data frame
      this.row <- SConly[n,] 
      dupindex <- mat.or.vec(1, nrow(SConly))
      for(x in 1:nrow(SConly)) {if(all(SConly[x,] == this.row)) {dupindex[x] = 1}}
      this.group <- MSframe[dupindex == 1,]
      #add the row with the highest MASCOT score to cleanframe, but only if the highest MASCOT value is >0.7 higher than the next highest one
      this.group[order(-this.group$MASCOTscore),]
      if(this.group$MASCOTscore[1] > this.group$MASCOTscore[2] + 0.7) {cleanframe <- rbind(cleanframe, this.group[1,])}
    }
  }
  #remove row 1 (empty from constructing cleanframe)
  cleanframe <- cleanframe[2:nrow(cleanframe),]
  
  #sum all values for different ions or modification states of the same peptide
  observedpeps <- levels(droplevels(cleanframe$peptides))
  summedpeps <- cleanframe[1,]
  #go through each of the observed peptides
  for(n in 1:length(observedpeps)) {
    #extract all rows for this peptide from cleanframe
    these.peptides <- cleanframe[cleanframe$peptides == observedpeps[n],]
    #Populate summedpeps with the summed spectral intensities for all rows
    summedpeps[n,1:2] <- these.peptides[1,1:2]
    summedpeps[n,3:ncol(summedpeps)] <- colSums(these.peptides[,-(1:2)])
    }
 
  #Add data on nature and site of substitution for each peptide based on varfilename
  variants <- readRDS(varfilename)
  subsframe <- as.data.frame(mat.or.vec(1,3))
  colnames(subsframe) <- c('SubSite','Expected', 'Observed')
  
  for(n in 1:nrow(summedpeps)) {
    subsframe[n,] <- variants[[3]][variants[[3]]$Sequence == summedpeps[n,2], 3:5]
  }
   summedpeps <- cbind(summedpeps[,1:2], subsframe, summedpeps[,3:ncol(summedpeps)])
  
  #renumber row names in sequential order
  row.names(summedpeps) <- 1:nrow(summedpeps)
  
  #If saveresults = TRUE, output results to file
  if(saveresult == TRUE) {
    #prepare a final output list.
  
    #List item 1 is a dataframe for logging all functions that have acted on the list. The first entry is MSload itself. Oher functions in this library add entries when they are called.
    loglist <- data.frame(func = as.character(match.call()[[1]]), date = format.Date(Sys.Date(), "%y%m%d"), time = format(Sys.time(), "%H:%M"), infile = '--', outfile = '--')
  
    #create a list with loglist and summedpeps as the two list items.
    resultslist <- list(loglist, summedpeps)
    names(resultslist) <- c('log', 'peptides')
  
    #if save is TRUE, save the results list as an RDS file. If filename is 'default', construct a based on varfilename, avoiding overwriting existing files by appending underscores
    if(outfilename == 'default') {
      varname <- as.list(match.call())[[2]]
      outfilename <- paste(as.character(format.Date(Sys.Date(), "%y%m%d")), ' ', varname, ' MS.Rds', sep = "")
      while(outfilename %in% dir()) {outfilename <- paste(strsplit(outfilename, '[.]')[[1]][1], '_.', strsplit(outfilename, '[.]')[[1]][2], sep = "")}
      resultslist[[1]]$outfile <- outfilename
      resultslist[[1]]$infile <- varfilename
    }
    
    saveRDS(resultslist, outfilename)
    
  }
  
  
  return(summedpeps)

}


