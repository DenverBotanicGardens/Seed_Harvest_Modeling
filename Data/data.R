#' 
#' @param TMx_list A list of annual transition matrices without fecundity (square matrix with N stages).

#' A list of annual transition matrices indicating survival and growth from col to row; fecundity in the first row 
#' (square matrix with N stages).
#' 
#' @format A list of square matrices
#' @source Ellis et al. 2012 "Matrix population models from 20 studies of perennial plant populations. Ecology 93(4) p. 951 
"Mx_all"

#' A list of population size vectors by stage (length N). Stable stage distriubiton or most recent count by stage
#' @format A list of vectors
#' @source Ellis et al. 2012 "Matrix population models from 20 studies of perennial plant populations. Ecology 93(4) p. 951 
"Nx_all"


#'  A list of annual transition matrices without fecundity (square matrix with N stages).
#' 
#' @format A list of square matrices
#' @source Ellis et al. 2012 "Matrix population models from 20 studies of perennial plant populations. Ecology 93(4) p. 951 
"TMx_all"


#' Demographic data (annual census) on Sclerocactus glaucus formatted for an Integral Projection Matrix
#' @format A data frame with 6478 rows and 12 variables:
#' \describe{
#'  \item{Site}{Named sites across the range with permanent macroplots covering dense portions of local population}
#'  \item{Transect}{1 meter wide linear transect within macroplot}
#'  \item{t0}{width, size of marked cactus at time t}
#'  \item{t1}{width, size of marked cactus at time t+1}
#'  \item{sds0}{count of seedlings, individuals less than 0.5 cm wide at time t, t0}
#'  \item{sds1}{count of seedlings, individuals less than 0.5 cm wide at time t+1, t1}
#'  \item{survival}{1 or 0 for survival of the marked individual measured in t0 and t1}
#'  \item{reproyesno}{1 or 0 for reproductive status of marked individual at time t, t0}
#'  \item{year0}{year of t0}
#'  \item{year1}{year of t1}
#'  \item{flowers0}{count of flowers for marked individual t0, estimated from one year of counts}
#'  \item{Region}{ScGl populations cluster into a northern and southern genetic region}
#'  }
#'  @source Denver Botanic Gardens through funding from Colorado Bureau of Land Management
"size_scgl"