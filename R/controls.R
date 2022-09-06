#' Title
#'
#' @param itermax
#' @param K_steps
#' @param n_init
#' @param zero_tol
#' @param zero_tol_1
#' @param zero_tol_2
#' @param trace
#'
#' @return
#' @export
#'
#' @examples
CWRM_control <- function(itermax =100,
                           K_steps =50,
                           n_init =10,
                           zero_tol=1e-16,
                           zero_tol_1=1e-16,
                         zero_tol_2=1e-90,
                         trace=0)
{
  list(itermax =itermax,
         K_steps =K_steps,
       n_init =n_init,
         zero_tol=zero_tol,
       zero_tol_1=zero_tol_1,
       zero_tol_2=zero_tol_2,
       trace=trace)
}
