# Simple wrapper around the original trim.cwm function implemented for

# article{garcia2017robust,
#   title={Robust estimation of mixtures of regressions with random covariates, via trimming and constraints},
#   author={Garc{\'\i}a-Escudero, Luis Angel and Gordaliza, Alfonso and Greselin, Francesca and Ingrassia, Salvatore and Mayo-{\'I}scar, Agust{\'\i}n},
#   journal={Statistics and Computing},
#   volume={27},
#   number={2},
#   pages={377--402},
#   year={2017},
#   publisher={Springer}
# }
# and contained in the CWRM.R file

#' Title
#'
#' @param X
#' @param y
#' @param G
#' @param c_X
#' @param c_y
#' @param alpha
#' @param control_CWRM
#'
#' @return
#' @export
#'
#' @examples
fit_CWRM <- function(X,y,G,c_X,c_y,alpha,control_CWRM=CWRM_control()){

  # Set CWRM controls
  itermax <- control_CWRM$itermax
  K_steps <- control_CWRM$K_steps
  n_init <- control_CWRM$n_init
  zero_tol <- control_CWRM$zero_tol
  zero_tol_1 <- control_CWRM$zero_tol_1
  zero_tol_2 <- control_CWRM$zero_tol_2
  trace <- control_CWRM$trace

  trim.cwm(
    X = X,
    Y = y,
    K = G,
    niter = itermax,
    Ksteps = K_steps,
    maxfactx = c_X,
    maxfacty = c_y,
    zero.tol = zero_tol,
    zero.tol1 = zero_tol_1,
    zero.tol2 = zero_tol_2,
    trace = trace,
    alpha = alpha,
    init = n_init
  )
}
