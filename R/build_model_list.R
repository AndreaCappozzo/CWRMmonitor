#' Title
#'
#' @param data
#' @param alpha_range
#' @param G_range
#' @param c_X_range
#' @param c_y_range
#'
#' @return
#' @export
#'
#' @examples
model_list_builder <-
  function(data,
           alpha_range,
           G_range,
           c_X_range,
           c_y_range) {
    # First variable in data must be the response
    y <- data[, 1]
    X_CWM <- as.matrix(data[, -1, drop = FALSE])

    parameters_grid <-
      expand.grid(
        num_groups = G_range,
        restr_factor_X = c_X_range,
        restr_factor_Y = c_y_range,
        trim_level = alpha_range
      )

    parameters_grid %>%
      dplyr::rowwise() %>%
      dplyr::mutate(model = list(tryCatch(
        trim.cwm(
          X = X_CWM,
          Y = y,
          K = num_groups,
          niter = 100,
          Ksteps = 50,
          maxfactx = restr_factor_X,
          maxfacty = restr_factor_Y,
          zero.tol = 1e-16,
          zero.tol1 = 1e-16,
          zero.tol2 = 1e-90,
          trace = 0,
          alpha = trim_level,
          init = 10
        ),
        error = function(e)
          list(obj = NA)
      ))) %>%
      dplyr::mutate(
        TBIC = -2 * model$obj +
          pnlt_term(restr_factor_X, restr_factor_Y, num_groups, p = p_x) *
          log(N * (1 - trim_level)),
        decomp = ifelse(is.na(model$obj), list(rep(NA,3)),list(var_dec(
          a = model, Y = df_CWM$y, X = as.matrix(df_CWM[,-1])
        ))),
        BSS = decomp[1],
        EWSS = decomp[2],
        RWSS = decomp[3],
        R = EWSS / (EWSS + RWSS)
      ) %>%
      dplyr::select(-decomp) %>%
      dplyr::ungroup()
  }

#' Title
#'
#' @param model_list
#'
#' @return
#' @export
#'
#' @examples
build_list_best_alpha_wise <- function(model_list){

  model_list %>%
    dplyr::group_by(trim_level) %>%
    dplyr::filter(TBIC==min(TBIC,na.rm = TRUE)) %>%
    dplyr::slice(1) %>% # in order to avoid ties
    dplyr::ungroup() %>%
    dplyr::arrange(trim_level)
}
