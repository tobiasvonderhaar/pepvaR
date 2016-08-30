#' Read in CSV files containing MS results
#'
#' \code{MS.read} reads in CSV files with lists of detected peptide variants (typically prepared using the Mass Spec-associated software).
#'
#' @param MSfilename The name of the file containing the MS data.
#'
#' @param pepfilename The name of the .Rds file containing the peptide list with which the MS data were generated. MS data are added as a new list item.
#'
#' @param scorecutoff The score value below which detected peptides are excluded from the analysis.
#'
#' @param scorecol The number of the column containg the MASCOT scores.
#'
#' @param sequencecol The number of the column containing the sequences of the detected peptides.
#'
#' @param datarange The numbers of the columns containing the spectral counts for the detected peptides.
#'
#'
#'
#' @return
#'
#'
#'
#' @export

MS.read <- function(MSfilename, scorecutoff = 4.9, scorecol = 9, sequencecol = 10, datarange = 20:37) {

  #conduct quality checks on loadfilename and savefilename: check whether filenames exists in the current working directory
  if(!(MSfilename %in% dir() )) {paste('One or both of the pecified filenames do not exist in the working directory. Note that filenames are case sensitive.', sep = ''); return()}

  MSframe <- read.csv(MSfilename)

  #only retain the columns for score, pepides sequence, and spectral counts.
  MSframe = MSframe[,c(scorecol, sequencecol, datarange)]

  #Convert values from the Score column from facor to numeric
  scores <- as.numeric(levels(MSframe$Score))[MSframe$Score]

  #Remove all rows where score is NA from MSframe
  MSframe <- MSframe[!is.na(scores),]
  scores <- na.omit(scores)

  #Remove all rows with score <5 from MSframe
  MSframe <- MSframe[scores >= scorecutoff,]

  #remove all peptides here for mltiple pepides idenical values have been reported for all samples (these are artefacts).
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
