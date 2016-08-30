#' Map sites of observed substitutions from list of detected peptides.
#' 
#' \code{match.subs} takes the list of observed peptides and maps sites of substitutions to the protein sequence.
#' 
#' @param filename The name of the RDS file which holds the petide and MS data.
#' 
#' @param save Indicates whether to save the RDS file with the added data back to disk (default = TRUE).
#' 
#' 
#' 
#' @return 
#' 
#' 
#' 
#' @export

match.subs <- function(filename, save = TRUE){
  
  #load the specified file
  dats <- readRDS(filename)
  
  #check whether MS data exist in the file
  if(!'MSdata' %in% names(dats)) {stop('Cannot recognise MS data in this file')}
  
  #determine for each peptide whether it is wt or not, insert this column into the MSdata dataframe
  is.wt <- mat.or.vec(nrow(dats$MSdata),1)
  is.wt[(lapply(dats$MSdata$Sequence, function(x) grep(x, dats$sequence)) == 1)]= 1
  is.wt <- as.logical(is.wt)
  dats$MSdata <- cbind(dats$MSdata[1], is.wt, dats$MSdata[2:ncol(dats$MSdata)])
  
  #determine for each peptide what the substitution site and type is (site = 0 means no subtitution)
  wtpeptides <- dats$peptides$Sequence[dats$peptides$Type == 'wt']
  seq <- dats$sequence

  substitutions <- mat.or.vec(nrow(dats$MSdata),2)
  
  for(n in 1:length(wtpeptides)) {
    this.wt <- wtpeptides[n]
    for(m in 1:nrow(dats$MSdata)) {
      if(dats$MSdata$is.wt[m]){substitutions[m,1] = 0; substitutions[m,2] = ''} else {
        location <- gregexpr(this.wt, seq)[[1]][1]
        this.detected <- as.character(dats$MSdata$Sequence[m])
        #check whether this.detected is shorter or equal in length to this.wt. If so compare this.wt, if not compare this.wt fused to next peptide
        if(nchar(this.wt) >= nchar(this.detected)) {
          if(as.numeric(adist(this.detected, substr(this.wt, 1, nchar(this.detected)))) == 1) {
            wt.compare <- unlist(strsplit(substr(this.wt, 1, nchar(this.detected)), split=""))
            det.compare <- unlist(strsplit(this.detected, split=""))
            pos <- which(wt.compare != det.compare)
            substitutions[m,1] <- location + pos -1
            substitutions[m,2] <- paste(wt.compare[pos], 'to', det.compare[pos], sep = ' ')
          }
        } else {
          this.extended <- paste(wtpeptides[n], wtpeptides[n+1], sep="")
          location <- gregexpr(this.extended, seq)[[1]][1]
          if(as.numeric(adist(this.extended, this.detected)) == 1 && nchar(this.extended) == nchar(this.detected)) {
            wt.compare <- unlist(strsplit(this.extended, split=""))
            det.compare <- unlist(strsplit(this.detected, split=""))
            pos <- which(wt.compare != det.compare)
            substitutions[m,1] <- location + pos -1
            substitutions[m,2] <- paste(wt.compare[pos], 'to', det.compare[pos], sep = ' ')
          }
        }
      }
    }
  }
  
  colnames(substitutions) <- c('Site', 'Substitution')
  
  dats$MSdata <- cbind(dats$MSdata[1:2], substitutions, dats$MSdata[3:ncol(dats$MSdata)])
  
  #save (data)if selected) and return extended dataframe
  if(save){saveRDS(dats, filename)}
  
  return(dats)
  
}