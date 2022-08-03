#' Title
#'
#' @param model_list
#' @param selected_alpha
#' @param eta_eta_threshold_ARI
#' @param number_of_optimal_sol
#'
#' @return
#' @export
#'
#' @examples
second_step_monitoring_CWRM <-
  function(model_list,
           selected_alpha,
           eta_threshold_ARI = 0.7,
           max_number_of_optimal_sol=4) {

    model_alpha_fixed <- model_list %>%
      dplyr::filter(trim_level == selected_alpha) %>%
      dplyr::select(1:3, 5:6) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(rank_solution = rank(TBIC, ties.method = "first")) %>%
      tidyr::drop_na()

    output_second_step_LIGHTER <-
      solutions_ranking_identifier_LIGHTER(
        solutions_df = model_alpha_fixed,
        number_of_optimal_sol = max_number_of_optimal_sol,
        epsilon_ARI = eta_threshold_ARI
      )



    gg_second_step <- output_second_step_LIGHTER$df_for_plotting %>%
      dplyr::mutate(
        to_fill = dplyr::case_when((is_stable & is_best_interval) ~ TRUE,
                            is_stable ~ FALSE,
                            TRUE ~ NA
        ),
        n_sol = factor(ifelse(is.na(to_fill), NA, n_sol)),
        num_groups = paste("G = ", num_groups)
      ) %>%
      ggplot2::ggplot(ggplot2::aes(x = restr_factor_X, y = restr_factor_Y)) +
      ggplot2::geom_tile(
        ggplot2::aes(
          fill = n_sol,
          alpha = (to_fill) |
            !is.na(rank_best)
        ),
        colour = "black",
        na.rm = TRUE,
        size = 0.5
      ) +
      ggplot2::geom_text(ggplot2::aes(label = rank_best)) +
      ggplot2::labs(x = quote(c[X]),
           y = quote(c[y]),
           fill = "Optimal solutions") +
      ggplot2::guides(alpha = "none") +
      ggplot2::facet_wrap ( ~ num_groups, nrow = 1) +
      ggplot2::scale_fill_brewer(
        palette = "Dark2",
        na.value = "white",
        breaks = 1:(
          dplyr::n_distinct(output_second_step_LIGHTER$df_for_plotting$n_sol) - 1
        )
      ) +
      ggplot2::scale_alpha_discrete(range = c(0.5, 1)) +
      ggplot2::theme(
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        legend.position = "bottom"
      )
    print(gg_second_step)
    output_second_step_LIGHTER$optimal_solutions
  }
