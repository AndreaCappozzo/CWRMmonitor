#' Title
#'
#' @param data
#' @param CWRM_model
#'
#' @return
#' @export
#'
#' @examples
compute_CWRM_decomp <- function(data,CWRM_model){

  data <- data.matrix(data)
  N <- nrow(data)

  df_decomp <- tibble::tibble(model = list(CWRM_model)) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      decomp = ifelse(is.na(CWRM_model$obj), list(rep(NA, 3)), list(
        var_dec(a = CWRM_model, Y = data[,1], X = as.matrix(data[, -1]))
      )),
      NBSS = decomp[1],
      NEWSS = decomp[2],
      NRWSS = decomp[3],
    ) %>%
    dplyr::select(-decomp) %>%
    dplyr::mutate(
      TSS = sum(NBSS + NEWSS + NRWSS),
      dplyr::across(.cols = -c(TSS,model), .fns = ~ .x / sum(TSS)),
    ) %>%
    dplyr::select(-TSS) %>%
    dplyr::ungroup()

  gg_decomp <- df_decomp %>%
    ggtern::ggtern(ggtern::aes(NBSS,NRWSS,NEWSS)) +
    ggtern::theme_rgbw(
    )+
    ggplot2::geom_point() +
    ggplot2::labs(x="NBSS",y="NRWSS",z="NEWSS")

  print(gg_decomp)
  df_decomp
}
