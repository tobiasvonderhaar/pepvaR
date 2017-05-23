#' Read in MS results from files or data frames
#'
#' \code{MS.read} reads in peptide information from a dataframe located in the workspace, or from a .csv file located in the working directory.
#'
#' @param MSfilename The name of the dataframe or file containing the MS data.
#'
#' @param scorecutoff The MASCOT score value below which detected peptides are excluded from the analysis (4.9 if not specified).
#'
#' @param scorecol The number of the column containg the MASCOT scores.
#'
#' @param sequencecol The number of the column containing the sequences of the detected peptides.
#'
#' @param datarange The numbers of the columns containing the spectral counts for the detected peptides (one column per sample)
#'
#'
#'
#' @return
#'
#'
#'
#' @export

MS.read <- function(MSfilename, scorecutoff = 4.9, scorecol = 8, sequencecol = 9, datarange = 11:28) {

  if(MSfilename == 'MSdata' | MSfilename %in% ls()) {
    eval(parse(text=paste('MSframe <- ', MSfilename, sep = '')))
  } else {
    if(MSfilename %in% dir()) {
      MSframe <- read.csv(MSfilename)
    } else {
      paste('The specified name does not correspond to a variable in the workspace or a file in the working directory. Names are case sensitive.', sep = ''); return()}
  }

  #only retain the columns for scores, peptide sequences, and spectral counts.
  MSframe = MSframe[,c(scorecol, sequencecol, datarange)]

  #Convert values from the Score column from factor to numeric
  scores <- as.numeric(levels(MSframe$Score))[MSframe$Score]

  #Remove all rows where score is NA from MSframe and from scores
  MSframe <- MSframe[!is.na(scores),]
  scores <- na.omit(scores)

  #Remove all rows with score <scoreutoff from MSframe
  MSframe <- MSframe[scores >= scorecutoff,]

  #remove all peptides where multiple identical values have been reported for all samples for multiple peptides (these are artefacts).
  SConly <- MSframe[,3:ncol(MSframe)]
  MSframe <- MSframe[rownames(MSframe) %in% rownames(unique(SConly)),]

  #renumber row names in sequential order
  row.names(MSframe) <- 1:nrow(MSframe)

  return(MSframe)

}


#' Sum spectral counts for different ions of the same peptide.
#'
#' \code{MS.condense} sums spectral counts for differently charged ions for the same peptide.
#'
#' @param MSdf A dataframe containing MS data (One column for Sequence, one or more columns for spectral counts).
#'
#' @param averagerepeats Determines whether repeats should be averaged, default is TRUE.
#'
#' @param repeats The number of repeat measurements for each sample.
#'
#' @param samplenames A vector of names to be used for the samples.
#'
#'
#'
#' @return
#'
#'
#'
#' @export
#'


MS.condense <- function(MSdf, repeats = 3, averagerepeats = TRUE, samplenames = '') {

  #if present remove Score column, which is no longer meaningful for added peptides
  if('Score' %in% colnames(MSdf)) {MSdf <- MSdf[,!(colnames(MSdf) == 'Score')]}



  #sum all values for different ions of the same peptides
  observedpeps <- levels(droplevels(MSdf$Sequence))
  summedpeps <- MSdf[1,]
  for(n in 1:length(observedpeps)) {
    summedpeps[n,1] <- observedpeps[n]
    for(m in 1:ncol(MSdf)) {
      if(is.numeric(MSdf[,m])) {
        summedpeps[n,m] <- sum(MSdf[MSdf$Sequence == observedpeps[n],m])
      }
    }
  }
  rm(observedpeps)
  MSdf <- summedpeps

  #Average data for sample repeats
  if(averagerepeats) {
    #Determine whether the number of repeats makes sense given the number of data columns.
    sampleno <- (ncol(MSdf) - 1) / repeats
    if(round(sampleno) - sampleno != 0) {cat('The number of observatons is not a multiple of the repeat number') } else {
      newframe = as.data.frame(MSdf[,1])
      for(n in seq(2,ncol(MSdf), by =3)) {
        for( m in 1:nrow(MSdf)) {
          newframe[m,(n-1+2)/3+1] <- mean(as.numeric(MSdf[m,n:(n+2)]))
        }
      }
      if(length(samplenames) == 1 && samplenames == '') {samplenames = seq(1:sampleno)}
      colnames(newframe) <- c('Sequence', samplenames)
    }

  }
  return(newframe)

}




#' Add the peptide data to the results list.
#'
#' \code{MS.save} loads a peptide file, appends the MS data as a new list item, and saves the file back to its original filename.
#'
#' @param MSdata The workpace variable containing the clean MS data.
#'
#' @param filename The name of the file to which the cleaned-up MS data will be appended.
#'
#'
#'
#' @return
#'
#'
#'
#' @export


MS.save <- function(MSdata, filename) {
  resultslist = readRDS(filename)
  resultslist[[4]] <- MSdata
  names(resultslist)[4] <- 'MSdata'
  saveRDS(resultslist, filename)

}
