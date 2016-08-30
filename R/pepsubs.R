#' Make tryptic peptide variants
#'
#' \code{make.var} generates a list of tryptic peptides and all their possible single amino acid variants from an input protein sequence
#'
#' @param sequence Character string containig the protein sequence in one letter format
#'
#' @param cutoff Peptides shorter than the value stored in cutoff are excluded from the list. Default is 4.
#'
#' @param save Boolean defaulting to TRUE. If TRUE, the list structure containing the output file will be written to file.
#'
#' @param filename Output filename to be used if save = TRUE. By default,  this will be the date followed by the name of the analysed variable and '.Rds'.
#'
#'
#'
#' @return a list containing a dataframe with the tryptic peptide sequence variants as the first list element.
#'
#'
#'
#' @export

makevar <- function(sequence, cutoff = 4, export = TRUE, filename = 'default') {

    #convert sequence to uppercase
  sequence <- toupper(sequence)

  #prepare background variables. noncutsites contains all non-cutting aminoacids, cutsites contains the cut sites.
  noncutters <- c('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'L', 'M', 'N', 'P', 'Q', 'S', 'T', 'V', 'W', 'Y')
  cutters <- c('K', 'R')


  #split the sequence into trpytic peptides at all cutsite aminoa acids.

  #To allow the algorithm to run correctly, the beginning and end of the sequence are represented as cutsites at 0 and at length(sequence).
  splitseq <- unlist(strsplit(sequence, ''))
  cutsites <- c(0, which(splitseq %in% cutters), nchar(sequence))

  #if the last amino acid is a K or R this introduces two cutsites at length(sequence), test for this and remove the last site if appropriate
  if(cutsites[length(cutsites)] == cutsites[length(cutsites)-1]) {cutsites <- cutsites[1:(length(cutsites)-1)]}

  wtpeptides <- as.data.frame(mat.or.vec(length(cutsites)-1,3))
  colnames(wtpeptides) <- c('ID', 'Sequence', 'Type')

  for (n in 1:(length(cutsites)-1)) {
    wtpeptides[n,1] <- paste('wt', n, sep = "")
    wtpeptides[n,2] <- substr(sequence, cutsites[n] + 1, cutsites[n+1])
    wtpeptides[n,3] <- 'wt'
  }

  rm(splitseq, n, cutsites)

  # generate the possible single amio acid substitutions for all peptides
  peptidecounter <- 1
  mutpeptides <- wtpeptides[1,]; mutpeptides[1,] <- c('','','')

  #1: process the loop for all wt peptides
  for (n in 1:dim(wtpeptides)[1]) {
    #2: do straightforward 1:1 swaps for all peptides >cutoff, at all position except the terminal position (a change here woudl affect the cut site and this isprocessed separately in step below).
    # Do not process substitutions to K or R which would introduce additional cut sites and would shorten the peptide.
    # Go through m = all letters in the peptide except the last which is K or R.
    for (m in 1:(nchar(wtpeptides[n,2])-1)) {
      #Go through all elements of aminoacids1.
      for (o in 1:length(noncutters)) {
        mutpeptides[peptidecounter,"ID"] <- paste(wtpeptides[n,"ID"], "var", peptidecounter, sep = "")
        mutpeptides[peptidecounter,"Sequence"] <- paste(substr(wtpeptides[n,"Sequence"], 0, m-1), noncutters[o], substr(wtpeptides[n,"Sequence"], m+1, nchar(wtpeptides[n,"Sequence"])), sep="")
        mutpeptides[peptidecounter,"Type"] <- 'variant'
        peptidecounter <- peptidecounter + 1
      }
      # for substitutions to K or R, substitute and then end the peptide as this introduces a new cut site
      for (p in 1:length(cutters)) {
        mutpeptides[peptidecounter,"ID"] <- paste(wtpeptides[n,"ID"], "var", peptidecounter, sep = "")
        mutpeptides[peptidecounter,"Sequence"] <- paste(substr(wtpeptides[n, "Sequence"], 0, m-1), cutters[p], sep="")
        mutpeptides[peptidecounter,"Type"] <- 'variant'
        peptidecounter <- peptidecounter + 1
      }
    }
    #process the last amino acid of each peptide
    #if the last amino acid is not K or R (this is the case for the C-terminal peptide of the protein), plain subsitute this
    m = nchar(wtpeptides[n,"Sequence"])
    if(substr(wtpeptides[n,"Sequence"], nchar(wtpeptides[n,"Sequence"]), nchar(wtpeptides[n,"Sequence"])) %in% noncutters) {
      for (o in 1:length(noncutters)) {
        mutpeptides[peptidecounter,"ID"] <- paste(wtpeptides[n,"ID"], "var", peptidecounter, sep = "")
        mutpeptides[peptidecounter,"Sequence"] <- paste(substr(wtpeptides[n,"Sequence"], 0, m-1), noncutters[o], sep="")
        mutpeptides[peptidecounter,"Type"] <- 'variant'
        peptidecounter <- peptidecounter + 1
      }
      for (p in 1:length(cutters)) {
        mutpeptides[peptidecounter,"ID"] <- paste(wtpeptides[n,"ID"], "var", peptidecounter, sep = "")
        mutpeptides[peptidecounter,"Sequence"] <- paste(substr(wtpeptides[n,"Sequence"], 0, m-1), cutters[p], sep="")
        mutpeptides[peptidecounter,"Type"] <- 'variant'
        peptidecounter <- peptidecounter + 1
      }
    } else {

      #process substitutions of the last amino acid where this is K or R
      #test whether the processed peptide is the last peptide.
      if(n == dim(wtpeptides)[1]) {
        #substitute, do not try to append the next peptide.
        for (o in 1:length(noncutters)) {
          mutpeptides[peptidecounter,"ID"] <- paste(wtpeptides[n,"ID"], "var", peptidecounter, sep = "")
          mutpeptides[peptidecounter,"Sequence"] <- paste(substr(wtpeptides[n,"Sequence"], 0, m-1), noncutters[o], sep="")
          mutpeptides[peptidecounter,"Type"] <- 'variant'
          peptidecounter <- peptidecounter + 1
        }
      }
      #if he processed peptide is not the last peptide
      else {
        #substitute and then append the next peptide as a cut-site has been removed.
        for (o in 1:length(noncutters)) {
          mutpeptides[peptidecounter,"ID"] <- paste(wtpeptides[n,"ID"], "var", peptidecounter, sep = "")
          mutpeptides[peptidecounter,"Sequence"] <- paste(substr(wtpeptides[n,"Sequence"], 0, m-1), noncutters[o], wtpeptides[n + 1,"Sequence"], sep="")
          mutpeptides[peptidecounter,"Type"] <- 'variant'
          peptidecounter <- peptidecounter + 1
        }
      }
        # if the last amino acid is R or K, substitute and end peptide
        for (p in 1:length(cutters)) {
          mutpeptides[peptidecounter,"ID"] <- paste(wtpeptides[n,"ID"], "var", peptidecounter, sep = "")
          mutpeptides[peptidecounter,"Sequence"] <- paste(substr(wtpeptides[n,"Sequence"], 0, m-1), cutters[p], sep="")
          mutpeptides[peptidecounter,"Type"] <- 'variant'
          peptidecounter <- peptidecounter + 1
        }
    }

  }

  #remove all peptides smaller than or equal to the cutoff
  mutpeptides <- mutpeptides[nchar(mutpeptides$Sequence) > cutoff,]
  #remove all wt sequences from mutpeptides
  mutpeptides <- mutpeptides[!mutpeptides$Sequence %in% wtpeptides$Sequence,]
  # rbind wtpeptides and mutpeptides together, renumber rows
  peptides <- rbind(wtpeptides, mutpeptides); row.names(peptides) <- 1:nrow(peptides)

  ##prepare a final output list.

  #List item 1 is a dataframe for logging all functions that have acted on the list. The first entry is trypvar itself.
  loglist <- data.frame(func = as.character(match.call()[[1]]), date = format.Date(Sys.Date(), "%y%m%d"), time = format(Sys.time(), "%H:%M"), infile = '--', outfile = '--')

  #create a list with loglist and peptides as the first two list items. Other items will be populated with other analyses performed with the variant sequences in peptides.
  resultslist <- list(loglist, sequence, peptides)
  names(resultslist) <- c('log', 'sequence', 'peptides')

  #if save is TRUE, save the results list as an RDS file. If filename is 'default', construct a systematic filename, avoiding overwriting existing files by appending asterisks
  if(export) {
    if(filename =='default') {
      varname <- as.list(match.call())[[2]]
      filename = paste(as.character(format.Date(Sys.Date(), "%y%m%d")), ' ', varname, '.Rds', sep = "")
      while(filename %in% dir()) {filename <- paste(strsplit(filename, '[.]')[[1]][1], '_.', strsplit(filename, '[.]')[[1]][2], sep = "")}
      resultslist[[1]]$outfile <- filename
    }

    saveRDS(resultslist, filename)

  }

  return(as.data.frame(resultslist[[3]]))

}


#' Export tryptic peptide variants to fasta file
#'
#' \code{varfasta} writes the petides in a variant list into a fasta file
#'
#' @param sourcevars A list output by the trypvar function.
#'
#' @param filename The filename of the fasta file. If unspecified, a filename will be conructed based on the input variable name.
#'
#'
#'
#' @return None. The function outputs a fasta-formatted file containing the peptide sequences from the specified variable.
#'
#'
#'
#' @export

varfasta <- function(infilename, outfilename = 'default') {

  #check whether filenames exists in the current working directory
  if(!infilename %in% dir()) {paste('\'', infilename, '\' does not exist in the working directory. Note that filenames are case sensitive.', sep = ''); return()}

  #read in the .Rds file
  resultslist <- readRDS(infilename)

  #check if outfilename was specified, if not (ie if filename = default) use the filename of the input file but chnge the extension to .fasta.
  if(outfilename =='default') {
    outfilename = paste(substr(infilename, 0, nchar(varname) - 4), '.fasta', sep = "")
  }


  #generate the .fasta file
  peptides <- as.data.frame(resultslist[[3]])
  sink(outfilename)
  for (n in 1:dim(peptides)[1]) { cat('>', peptides[n,"ID"], '\n', peptides[n,"Sequence"], '\n', sep='') }
  sink()

  #update the loglist.
  this.loglist <- data.frame(log.function = as.character(match.call()[[1]]), log.date = format.Date(Sys.Date(), "%y%m%d"), log.time = format(Sys.time(), "%H:%M"), infile = infilename, outfile = outfilename)
  resultslist[[1]] <- rbind(resultslist[[1]], this.loglist)


  #write the updated .Rds file back to disk
  saveRDS(resultslist, infilename)
}