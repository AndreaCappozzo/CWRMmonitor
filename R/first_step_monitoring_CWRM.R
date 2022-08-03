#' Title
#'
#' @param best_alpha_wise
#' @param relabeler_list
#' @param what
#'
#' @return
#' @export
#'
#' @examples
first_step_monitoring_CWRM <-
  function(best_alpha_wise,
           relabeler_list,
           what = c(
             "everything",
             "groups_proportion",
             "line_plots_parameters",
             "ARI",
             "prop_doubtful_assignment"
           )) {


  G_max <- max(best_alpha_wise$num_groups)
  N <- length(best_alpha_wise$model[[1]]$assig)
  alpha <- best_alpha_wise$trim_level
  color_palette_set_2 <- RColorBrewer::brewer.pal(n=G_max,name = "Set2")

  cluster_prop_colors <- c("0"="black",color_palette_set_2)
  names(cluster_prop_colors)[-1] <- as.character(1:G_max)

  # Cluster proportions including trimming level

  gg_mix <- best_alpha_wise %>%
    tidyr::hoist(.col=model, "csize") %>%
    dplyr::mutate(csize_relabeled= purrr::map(1:nrow(best_alpha_wise), ~csize[[.x]][relabeler_list[[.x]]]/N)) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(csize_relabeled=list(c(trim_level,csize_relabeled))) %>%
    dplyr::ungroup()  %>%
    tidyr::unnest_wider(csize_relabeled, names_sep = ".") %>%
    dplyr::select(trim_level, dplyr::starts_with("csize_relabeled")) %>%
    tidyr::pivot_longer(-trim_level) %>%
    dplyr::mutate(name=factor(as.numeric(stringr::str_extract(name,".$"))-1), title="Groups proportion") %>%
    ggplot2::ggplot(ggplot2::aes(x=trim_level,y=value,fill=name)) +
    ggplot2::geom_col() +
    ggplot2::facet_wrap(~title) +
    ggplot2::scale_x_continuous(breaks = alpha,labels = scales::number_format(accuracy = 0.001)) +
    ggplot2::labs(x=quote(alpha), y="", fill="G") +
    ggplot2::scale_fill_manual(breaks = names(cluster_prop_colors[-1]), values = cluster_prop_colors,
    ) +
    ggplot2::theme(legend.position = "bottom",axis.text.x = ggplot2::element_text(angle = 45,hjust = 1))

  BETA_df <- vector(mode = "list",length = ncol(best_alpha_wise$model[[1]]$b))

  for (m in 1:length(BETA_df)) {
    BETA_df[[m]] <- purrr::map(1:nrow(best_alpha_wise),~best_alpha_wise$model[[.x]]$b[,m]) %>%
      tibble::enframe() %>%
      dplyr::mutate(value_relabeled=purrr::map(1:nrow(best_alpha_wise), ~value[[.x]][relabeler_list[[.x]]])) %>%
      tidyr::unnest_wider(col = value_relabeled, names_sep = paste0("Beta",m-1)) %>%
      dplyr::select(-value)
  }

  sigma_regr_df <- purrr::map(1:nrow(best_alpha_wise),~sqrt(best_alpha_wise$model[[.x]]$v)) %>%
    tibble::enframe() %>%
    dplyr::mutate(value_relabeled=purrr::map(1:nrow(best_alpha_wise), ~value[[.x]][relabeler_list[[.x]]])) %>%
    tidyr::unnest_wider(col = value_relabeled, names_sep = "sigma_regr") %>%
    dplyr::select(-value)

  det_sigma_df <- purrr::map(1:nrow(best_alpha_wise),~apply(best_alpha_wise$model[[.x]]$sigma,3,function(s) det(s)^(1/p_x))) %>%
    tibble::enframe() %>%
    dplyr::mutate(value_relabeled=purrr::map(1:nrow(best_alpha_wise), ~value[[.x]][relabeler_list[[.x]]])) %>%
    tidyr::unnest_wider(col = value_relabeled, names_sep = "det_sigma") %>%
    dplyr::select(-value)

  param_list <- append(BETA_df,list(sigma_regr_df,det_sigma_df))
  param_df <- best_alpha_wise %>%
    dplyr::select(trim_level) %>%
    dplyr::bind_cols(plyr::join_all(param_list,by="name",type="inner")) %>%
    dplyr::select(-name) %>%
    tidyr::pivot_longer(-trim_level) %>%
    dplyr::mutate(name=stringr::str_replace_all(string=name,pattern = "value_relabeledBeta",replacement = "b[g]^")) %>%
    dplyr::mutate(name_minus_last=substring(name,1, nchar(name)-1),param=dplyr::case_when(stringr::str_detect(name,pattern = "b\\[g\\]\\^")~name_minus_last,
                                         stringr::str_detect(name,pattern = "Beta0") ~ "b[g]^0",
                                         stringr::str_detect(name,pattern = "sigma_regr")~"sigma[g]",
                                         stringr::str_detect(name,pattern = "det_sigma")~"abs(Sigma[g])^{1/d}"),
           group=stringr::str_extract(string = name,pattern = ".$"))


  gg_param <- param_df %>%
    dplyr::mutate(param=forcats::fct_relevel(.f = param,"abs(Sigma[g])^{1/d}",after = Inf)) %>%
    ggplot2::ggplot(ggplot2::aes(trim_level,value, group=group, col=group)) +
    ggplot2::geom_point() +
    ggplot2::geom_line(ggplot2::aes(lty=group)) +
    ggplot2::facet_grid(rows = dplyr::vars(param),labeller = ggplot2::label_parsed,scales = "free") +
    ggplot2::scale_x_continuous(breaks = alpha,labels = scales::number_format(accuracy = 0.001)) +
    ggplot2::labs(x=quote(alpha), y="", col="G", lty="G") +
    ggplot2::scale_color_brewer(palette = "Set2") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,hjust = 1),
          strip.text.y = ggplot2::element_text(angle = 0))

  half_tick <- diff(alpha[-length(alpha)])[1]/2

  gg_ARI_05_y <-
    tibble::tibble(
      trim_level = alpha[-length(alpha)]+half_tick,
      ARI = purrr::map_dbl(
        1:(nrow(best_alpha_wise) - 1),
        ~ mclust::adjustedRandIndex(best_alpha_wise$model[[.x]]$assig, best_alpha_wise$model[[.x +1]]$assig)
      ), title="Adjusted Rand Index"
    ) %>%
    ggplot2::ggplot(ggplot2::aes(trim_level, ARI)) +
    ggplot2::geom_line()+
    ggplot2::geom_point()+
    ggplot2::scale_y_continuous(limits = c(0.5,1)) +
    ggplot2::scale_x_continuous(breaks = alpha, limits = range(alpha),labels = scales::number_format(accuracy = 0.001)) +
    ggplot2::facet_wrap(~title) +
    ggplot2::labs(x=quote(alpha), y="")+
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,hjust = 1))

  # Monitoring proportion of doubtful assignment

  gg_doubtful_assignment <- best_alpha_wise %>%
    dplyr::rowwise() %>%
    dplyr::mutate(discr_function = list(
      Discrim_Function_CWM(
        y = df_CWM$y,
        X = as.matrix(df_CWM[, -1, drop = FALSE]),
        alpha = trim_level,
        assig = model$assig,
        tau = model$cw,
        beta = t(model$b),
        mu = t(model$center),
        Sigma = model$sigma,
        sigma = model$v
      )
    )) %>%
    dplyr::mutate(doubtful_assignment_list = list(lapply(discr_function, function(x)
      mean(x > log(
        1 / 10
      )))), title="Proportion of Doubtful Assignments") %>%
    tidyr::unnest_wider(doubtful_assignment_list) %>%
    ggplot2::ggplot(ggplot2::aes(x = trim_level)) +
    ggplot2::geom_line(ggplot2::aes(y = DF_CWM), col="black") +
    ggplot2::geom_point(ggplot2::aes(y = DF_CWM), col="black") +
    ggplot2::scale_x_continuous(breaks = alpha, labels = scales::number_format(accuracy = 0.001)) +
    ggplot2::facet_wrap(~title) +
    ggplot2::labs(x=quote(alpha), y="", color=" ") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,hjust = 1))

  if(what=="everything") {
    list(gg_mix,
    gg_param,
    gg_ARI_05_y,
    gg_doubtful_assignment)
  } else if (what == "groups_proportion") {
    gg_mix
  } else if (what == "line_plots_parameters") {
    gg_param
  } else if (what == "ARI") {
    gg_ARI_05_y
  } else if (what == "prop_doubtful_assignment") {
    gg_doubtful_assignment
  }
}
