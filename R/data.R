#' MS dataset of GST-derived peptides containing misincorporation errors.
#'
#' An example MS dataset for the study of amino acdi misincorporation errors generated with
#' tryptic digests of semi-purified GST samples expressed in yeast and E coli, as described
#' in Josse et al (2017).
#'
#' @format A data frame with 603 rows and 28 variables:
#' \describe{
#'   \item{Retention.time..min.}{GC retention time of the peptide in minutes}
#'   \item{Charge}{Measured charge of the peptide ion}
#'   \item{Drift.time..ms.}{Measured drift time in milliseconds}
#'   \item{m.z}{Mass-to-charge ratio of the peptide}
#'   \item{Measured.mass}{Measured mass of the peptide}
#'   \item{Mass.error..u.}{Absolute error of the measured mass}
#'   \item{Mass.error..ppm.}{Realtive error of the measured mass}
#'   \item{Score}{MASCOT quality score as determined by the MS isntrument software}
#'   \item{Sequence}{The sequence of the detected peptide}
#'   \item{Modifications}{Post-translational modifications to the peptide as interpreted by the MS instrument software}
#'   \item{LJ01}{Spectral intensity of peptides for gst-disec sample, replicate 1}
#'   \item{LJ02}{Spectral intensity of peptides for gst-disec sample, replicate 2}
#'   \item{LJ03}{Spectral intensity of peptides for gst-disec sample, replicate 3}
#'   \item{LJ04}{Spectral intensity of peptides for gst-optec sample, replicate 1}
#'   \item{LJ05}{Spectral intensity of peptides for gst-optec sample, replicate 2}
#'   \item{LJ06}{Spectral intensity of peptides for gst-optec sample, replicate 3}
#'   \item{LJ07}{Spectral intensity of peptides for gst-dissc sample, replicate 1}
#'   \item{LJ08}{Spectral intensity of peptides for gst-dissc sample, replicate 2}
#'   \item{LJ09}{Spectral intensity of peptides for gst-dissc sample, replicate 3}
#'   \item{LJ10}{Spectral intensity of peptides for gst-dissc + NAT sample, replicate 1}
#'   \item{LJ11}{Spectral intensity of peptides for gst-dissc + NAT sample, replicate 2}
#'   \item{LJ12}{Spectral intensity of peptides for gst-dissc + NAT sample, replicate 3}
#'   \item{LJ13}{Spectral intensity of peptides for gst-optsc sample, replicate 1}
#'   \item{LJ14}{Spectral intensity of peptides for gst-optsc sample, replicate 2}
#'   \item{LJ15}{Spectral intensity of peptides for gst-optsc sample, replicate 3}
#'   \item{LJ16}{Spectral intensity of peptides for gst-optsc + NAT sample, replicate 1}
#'   \item{LJ17}{Spectral intensity of peptides for gst-optsc + NAT sample, replicate 2}
#'   \item{LJ18}{Spectral intensity of peptides for gst-optsc + NAT sample, replicate 3}
#'   ...
#' }
"MSdata"
