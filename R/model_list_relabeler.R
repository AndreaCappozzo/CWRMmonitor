#' Title
#'
#' @param model_list
#'
#' @return
#' @export
#'
#' @examples
list_best_alpha_wise_relabeler <- function(best_alpha_wise){

  # Highest trim level for finding deepest obs
  most_robust_model <- best_alpha_wise$model[[nrow(best_alpha_wise)]]

  G_most_robust_modeling <- nrow(most_robust_model$center)

  y_representative_fake_obs <- sapply(1:G_most_robust_modeling,function(g)
    most_robust_model$b[g,] %*% c(1,most_robust_model$center[g,]))
  group_representative_fake_obs <- cbind(y_representative_fake_obs, most_robust_model$center)


  z_CWM <-
    estep_CWM(
      y = group_representative_fake_obs[,1],
      X = group_representative_fake_obs[,-1,drop=FALSE],
      tau = most_robust_model$cw,
      beta = t(most_robust_model$b),
      mu = t(most_robust_model$center),
      Sigma = most_robust_model$sigma,
      sigma = most_robust_model$v
    )

  map_CWM <- mclust::map(z_CWM)

  n_models <- nrow(best_alpha_wise)

  relabeler_list <- vector(mode = "list", length = n_models)
  relabeler_list[[n_models]] <- 1:length(map_CWM)

  most_robust_relabeler <- 1:(dplyr::n_distinct(most_robust_model$assig)-1)

  relabeler_list[[n_models]] <- most_robust_relabeler
  group_representative_fake_obs <- group_representative_fake_obs[most_robust_relabeler,]

  labels_relabeled_df <- purrr::map_dfc(.x = 1:n_models, .f = ~identity(best_alpha_wise$model[[.x]]$assig))

  labels_to_be_changed <- labels_relabeled_df[[n_models]]

  for(k in 1:length(relabeler_list[[n_models]])){
    labels_to_be_changed[labels_relabeled_df[[n_models]]==k] <- which(relabeler_list[[n_models]]==k)
  }

  labels_relabeled_df[[n_models]] <-   labels_to_be_changed

  for(n_mod in (n_models-1):1){

    cluster_relabeler_output <- cluster_relabeler_via_CWM_density(
      X = df_CWM,
      group_representative_fake_obs = group_representative_fake_obs,
      old_partition = labels_relabeled_df[[n_mod+1]],
      new_partition = labels_relabeled_df[[n_mod]],
      old_model = best_alpha_wise$model[[n_mod+1]],
      new_model = best_alpha_wise$model[[n_mod]]
    )
    # Save relabeler
    relabeler_list[[n_mod]] <- cluster_relabeler_output$relabeler

    # Relabel model n_mod
    labels_to_be_changed <- labels_relabeled_df[[n_mod]]
    for(k in 1:length(relabeler_list[[n_mod]])){
      labels_to_be_changed[labels_relabeled_df[[n_mod]]==k] <- which(relabeler_list[[n_mod]]==k)
    }
    labels_relabeled_df[[n_mod]] <-   labels_to_be_changed

    # Update groupwise fake obs
    group_representative_fake_obs <- cluster_relabeler_output$group_representative_fake_obs
  }

  relabeler_list

}
