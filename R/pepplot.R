#' Plot the coverage in terms of wild-type peptides
#' 
#' \code{plot.cover} Plots the coverage of the protein sequence in terms of observed wild-type peptides, as well as locations of observed substitutions.
#' 
#' @param subsdata The name of the file that contains the Mascot peptide sequences and substitution data.
#' 
#' 
#' 
#' @return 
#' 
#' 
#' 
#' @export

plot.cover <- function(subsdata){
 
  wtpeptides <- as.character(subsdata$Sequence[subsdata$is.wt])
  seq.length <- nchar(subsdata$sequence)
  
  #prepare the basic plot window
  plot(c(0,12), c(0,12), type = "n", axes = FALSE, ylab = NA, xlab = NA, asp = 1)
  rect(1,1,11,2, col = 'grey75', border = NA)
  
  #draw boxes for the wild-type coverage
  for(n in 1:length(wtpeptides)) {
   
  #extract the location of this peptide in seq
    location <- gregexpr(wtpeptides[n], subsdata$sequence)
   
  #plot boxes corresponding to each peptide
    plot.start <- location[[1]][1] / seq.length * 10 + 1
    plot.end <- (location[[1]][1] / seq.length + attributes(location[[1]])$match.length / seq.length) * 10 + 1
    rect(plot.start,1,plot.end,2, col= 'black', border = 'white')
    rm(location)
  }
  
  #draw ticks at sites of substitutions
  ticks <- as.numeric(levels(subsdata$MSdata$Site)[subsdata$MSdata$Site])
  ticks <- unique(ticks[ticks != 0])
  rect((ticks) / seq.length *10 + 1, 2.1, (ticks) / seq.length * 10 + 1, 2.6)
  
  #add axes labels
  rect(1,0.5,1,1); text(1,0.4, '0', adj = c(0.5,1))
  rect(11,0.5,11,1); text(11,0.4, as.character(nchar(dats$sequence)), adj = c(0.5,1))
  text(6,0.4, 'Codon Number', adj = c(0.5,1))
  
  #add legend
  rect(2,4,2.5,4.5, col='grey75'); text(2.6, 4.25, 'WT Peptide not detected', adj = c(0,0.5))
  rect(2,5,2.5,5.5, col='black'); text(2.6, 5.25, 'WT Peptide detected', adj = c(0,0.5))
  rect(2.25,6,2.25,6.5); text(2.6,6.25, 'Substitution Observed', adj = c(0,0.5))
}

#' Plot the relative spectral counts for sustitute peptides
#' 
#' \code{plot.counts} Plots cumulative bar graphs for spectral counts of substituted peptides, normalised to the corresponding wild type peptide.
#' 
#' @param filename The name of the file that contains the Mascot peptide sequences and substitution data.
#' 
#' @param excludetop The number of high-spectral-count peptides to exclude, counted from the top
#' 
#' 
#' 
#' @return 
#' 
#' 
#' 
#' @export

plot.counts <- function (filename, excludetop = 0) {
  #load the specified file
  dats <- readRDS(filename)
  
  #determine which substituted peptides correspond to which wild-type petides
  wtpeptidenos <- as.integer(rownames(dats$MSdata[dats$MSdata$is.wt,]))
  groupings <- mat.or.vec(length(dats$MSdata$Sequence),1)
  for(n in 1:length(wtpeptidenos)) {
    for(m in 1:length(dats$MSdata$Sequence)) {
      thiswtpep <- as.character(dats$MSdata$Sequence[wtpeptidenos[n]])
      thisquerypep <- as.character(dats$MSdata$Sequence[m])
      if(nchar(thiswtpep) < nchar(thisquerypep)) {thisquerypep <- substr(thisquerypep,1,nchar(thiswtpep))}
      else if(nchar(thiswtpep) > nchar(thisquerypep)) {thiswtpep <- substr(thiswtpep,1,nchar(thisquerypep))}
      if(m == wtpeptidenos[n]) {groupings[m] <- wtpeptidenos[n]} else if(groupings[m] == 0 && adist(thiswtpep,thisquerypep) < 2) {groupings[m] <- wtpeptidenos[n]}
    }
  }
  groupings <- as.factor(groupings)
  
  #present the available peptides and request input to identify the peptide to be plotted
  #peptidechoice <- "5" #Placeholder!!
  cat(as.character(groupings))
  peptidechoice <- ''
  
  while(!peptidechoice %in% levels(groupings)) {peptidechoice <- as.character(readline('Please enter the number for the wt peptide you wish to plot:'))}
  
  
  #normalise the peptide data in preparation for plotting
  this.subset <- dats$MSdata[groupings == peptidechoice,]
  this.subset.normalised <- this.subset
  cols2process <- colnames(this.subset)[!colnames(this.subset) %in% c('Sequence','is.wt','Site','Substitution')]
  for(n in 1:length(cols2process)){
    for(m in 1:nrow(this.subset)){
      this.subset.normalised[m,cols2process[n]] <- as.numeric(this.subset[m,cols2process[n]])/as.numeric(this.subset[this.subset$is.wt,cols2process[n]])
    }
  }
  this.subset.normalised <- this.subset.normalised[this.subset.normalised$is.wt == FALSE,]
  #dd[with(dd, order(-z, b)), ]
  this.subset.normalised <- this.subset.normalised[with(this.subset.normalised, order(-this.subset.normalised[,5])),]
  rownames(this.subset.normalised) <- seq(1:nrow(this.subset.normalised))
  this.table <- as.matrix(this.subset.normalised[(excludetop + 1):nrow(this.subset.normalised),cols2process])
  barplot(this.table, col = rainbow(nrow(this.subset.normalised))[(excludetop + 1):nrow(this.subset.normalised)], ylab = 'Substituted peptide spectral counts (normalised to unsubstituted)')
  return(this.subset.normalised)
}