#' Plot the coverage in terms of wild-type peptides from a file produced by MS.load()
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

plot.cover <- function(infile){
  
  #check whether there is a file corresponding to infilename in the curren working directory
  
  if(!infilename %in% dir()) {
    cat('There is no file called', infilename, 'in the current working directory')
    return() 
  }
 
  #read in the MS data from the specified file
  subsdata <- readRDS(infile)[[2]]
  #recover the original protein sequence by retrieving it from the variant file used to generate infile
  sequence <- readRDS(readRDS(infile)[[1]]$infile)[[2]]
  
  #recover the wild-type sequences and their average spectral abundancesfrom the MS file (those that have a SubSite value of 0 ie no substitution)
  wtpeptides <- as.character(subsdata$peptides[subsdata$SubSite == 0])
  wtabundances <- log(rowMeans(subsdata[subsdata$SubSite == 0,6:23]))
  wtabundances <- (wtabundances - min(wtabundances))/(max(wtabundances)-min(wtabundances))
  
  #determine the length of the source sequence
  seq.length <- nchar(sequence)
  
  
  #extract the position of the first amino acid of the peptides in the sequence
  pep.start <- mat.or.vec(length(wtpeptides),1)
  for(n in 1:length(wtpeptides)) {
    pep.start[n] <- as.numeric(gregexpr(pattern = wtpeptides[n], sequence))
  }
  #determine the position of the last amino acid of each peptide
  pep.stop <- pep.start + nchar(wtpeptides)
  
  #prepare the basic plot window
  plot(c(0,seq.length * 1.2), c(0,8), type = "n", axes = FALSE, ylab = NA, xlab = NA)
  rect(seq.length * 0.1,1,seq.length * 1.1,2, col = 'grey75', border = NA)
  
  
  #plot boxes corresponding to the peptides, coloured accroding to their abundance values
  colfunc <- colorRamp(c('blue','yellow'))
  colvec <- rgb(colfunc(wtabundances), maxColorValue=255)
  rect(pep.start + seq.length*0.1,1,pep.stop + seq.length*0.1,2, col= colvec, border = 'white')
  
  #draw ticks at sites of substitutions
  ticks <- as.numeric(unique(subsdata$SubSite[!subsdata$SubSite == 0]))
  rect(ticks + (seq.length * 0.1), 2.1, ticks + (seq.length * 0.1), 2.6)
  
  #add axes labels
  rect(seq.length * 0.1,0.5,seq.length * 0.1,1); text(seq.length * 0.1 ,0.4, '0', adj = c(0.5,1))
  rect(seq.length * 1.1,0.5,seq.length * 1.1,1); text(seq.length * 1.1,0.4, as.character(seq.length), adj = c(0.5,1))
  text(seq.length * 0.6, 0.4, 'Codon Number', adj = c(0.5,1))
  
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
