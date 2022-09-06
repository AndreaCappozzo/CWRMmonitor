#' Title
#'
#' @param data
#' @param CWRM_model
#'
#' @return
#' @export
#'
#' @examples
compute_CWRM_discr_factor <- function(data,CWRM_model){

  data <- data.matrix(data)
  N <- nrow(data)

  G_max <- ifelse(length(CWRM_model$cw)>=3,length(CWRM_model$cw),3) # to avoid unnecessary warnings when defining the palette
  color_palette_set_2 <- RColorBrewer::brewer.pal(n=G_max,name = "Set2")

  cluster_prop_colors <- c("0"="black",color_palette_set_2)
  names(cluster_prop_colors)[-1] <- as.character(1:G_max)

  # Silhouette plots

  # Label invariant
  discr_factor_CWM <- Discrim_Function_CWM(
    y = data[,1],
    X = as.matrix(data[, -1, drop = FALSE]),
    alpha = sum(CWRM_model$assig == 0) / length(CWRM_model$assig),
    assig = CWRM_model$assig,
    tau = CWRM_model$cw,
    beta = t(CWRM_model$b),
    mu = t(CWRM_model$center),
    Sigma = CWRM_model$sigma,
    sigma = CWRM_model$v
  )

  df_4_plot_DF <- tibble::enframe(discr_factor_CWM) %>%
    dplyr::mutate(
      name = paste0(stringr::str_replace(
        string = name,
        pattern = "_",
        replacement = "["
      ), "]"),
      assig = list(factor(CWRM_model$assig)),
      rowid = list(1:N)
    ) %>%
    tidyr::unnest(c(value, assig, rowid)) %>%
    dplyr::arrange(assig, name, value) %>%
    dplyr::mutate(name=dplyr::case_when(name=="DF[CWM]"~"DF",
                          name=="DF[regr]"~"DF[YIX]",
                          name=="DF[X]"~"DF[X]"))


  df_4_plot_DF$rowid <-
    factor(df_4_plot_DF$rowid, levels = unique(df_4_plot_DF$rowid))

  # df_4_plot_DF_new <- df_4_plot_DF
  df_4_plot_DF_new <- df_4_plot_DF %>%
    dplyr::group_by(name) %>%
    dplyr::mutate(
      value = ifelse(value <= stats::median(value), stats::median(value), value),
      DF_is_tresh = ifelse(value <= stats::median(value), "yes", "no")
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(name=forcats::fct_relevel(.f = factor(name),"DF[X]",after = 2))

  gg_silhouette <- df_4_plot_DF_new %>%
    ggplot2::ggplot(ggplot2::aes(
      x = value,
      y = rowid,
      fill = assig,
      col = assig
    )) +
    ggplot2::geom_col(show.legend = FALSE, width = .05) +
    ggplot2::geom_point(data = dplyr::filter(df_4_plot_DF_new, DF_is_tresh == "yes"),
               show.legend = FALSE) +
    ggplot2::facet_grid(~ name, scales = "free", labeller = ggplot2::label_parsed) +
    ggplot2::labs(x = "", y = "") +
    ggplot2::geom_vline(xintercept = log(1/10),col="red",linetype=2, size=.85) +
    ggplot2::scale_color_manual(breaks = names(cluster_prop_colors),
                       values = cluster_prop_colors) +
    ggplot2::scale_fill_manual(breaks = names(cluster_prop_colors),
                      values = cluster_prop_colors) +
    ggplot2::scale_y_discrete(guide = ggplot2::guide_axis(n.dodge = 2))

  print(gg_silhouette)

  df_4_plot_DF %>%
    tidyr::pivot_wider(names_from = name,values_from = value) %>%
    dplyr::mutate(rowid=readr::parse_number(as.character(rowid))) %>%
    dplyr::arrange(rowid)
}
