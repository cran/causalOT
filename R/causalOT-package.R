#' An R package to perform causal inference using optimal transport distances.
#'
#' R code to perform causal inference weighting using a variety of methods and optimizers. The code can estimate weights, estimate treatment effects, and also give variance estimates. These methods are described in Dunipace, Eric (2021) <https://arxiv.org/abs/2109.01991>.
#'
#' @author Eric Dunipace
#' @name causalOT
#' @useDynLib causalOT, .registration = TRUE
#' @importFrom Rcpp sourceCpp 
#' @importFrom Rcpp evalCpp
#' @importFrom methods as slotNames
#' @importFrom stats .checkMFClasses .getXlevels as.formula binomial coef contrasts cov delete.response drop.terms formula glm lm median model.frame model.matrix model.response na.pass optim pgamma qchisq qnorm quantile rmultinom runif sd setNames terms var weighted.mean
#' @importFrom utils capture.output data methods setTxtProgressBar txtProgressBar
#' @importFrom R6 R6Class
#' @importFrom methods setOldClass
#' @rdname causalOT-package
#' @keywords internal 
"_PACKAGE"

#' An external control trial of treatments for post-partum hemorrhage
#' 
#' A dataset evaluating treatments for post-partum hemorrhage. The data contain  treatment groups receiving misoprostol vs potential
#' controls from other locations that received only oxytocin. The data is 
#' stored as a numeric matrix.
#'
#' The variables are as follows:
#' \itemize{
#'   \item cum_blood_20m. The outcome variable denoting cumulative blood loss in mL 20 minutes after the diagnosis of post-partum hemorrhage (650 -- 2000).
#'   \item tx. The treatment indicator of whether an individual received misoprostol (1) or oxytocin (0).
#'   \item age. the mother's age in years (15 -- 43).
#'   \item no_educ. whether a woman had no education (1) or some education (0).
#'   \item num_livebirth. the number of previous live births.
#'   \item cur_married. whether a mother is currently married (1 = yes, 0 = no).
#'   \item gest_age. the gestational age of the fetus in weeks (35 -- 43).
#'   \item prev_pphyes. whether the woman has had a previous post-partum hemorrahge.
#'   \item hb_test. the woman's hemoglobin in mg/dL (7 -- 15).
#'   \item induced_laboryes. whether labor was induced (1 = yes, 0 = no).
#'   \item augmented_laboryes. whether labor was augmented (1 = yes, 0 = no).
#'   \item early_cordclampyes. whether the umbilical cord was clamped early (1 = yes, 0 = no).
#'   \item control_cordtractionyes. whether cord traction was controlled (1 = yes, 0 = no).
#'   \item uterine_massageyes. whether a uterine massage was given (1 = yes, 0 = no).
#'   \item placenta. whether placenta was delivered before treatment given (1 = yes, 0 = no).
#'   \item bloodlossattx. amount of blood lost when treatment given (500 mL -- 1800 mL)
#'   \item sitecode. Which site is the individual from? (1 = Cairo, Egypt,  2 = Turkey,        3 = Hocmon, Vietnam,  4 = Cuchi, Vietnam, and 5 Burkina Faso).
#' }
#' 
#' @source Data from the following Harvard Dataverse: 
#' \itemize{
#' \item Winikoff, Beverly, 2019, "Two randomized controlled trials of misoprostol for the treatment of postpartum hemorrhage", https://doi.org/10.7910/DVN/ETHH4N, Harvard Dataverse, V1. 
#' }
#' The data was originally analyzed in 
#' \itemize{
#' \item Blum, J. et al. Treatment of post-partum haemorrhage with sublingual misoprostol versus oxytocin in women receiving prophylactic oxytocin: a double-blind, randomised, non-inferiority trial. The Lancet 375, 217--223 (2010).
#' }
#'
#' @docType data
#' @keywords datasets
#' @name pph
#' @usage data(pph)
#' @format A matrix with 802 rows and 17 variables
NULL