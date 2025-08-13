#' Posterior Interaction Profile — *Waterfall* style
#'
#' Creates the **waterfall** flavour of a Posterior Interaction Profile (PIP)
#' plot, visualising how a single respondent’s latent position (\eqn{\alpha_p})
#' interacts with every item.
#' The left panel shows the posterior density (and optional HDI) of the chosen
#' respondent’s ability.
#' The right panel (“waterfall”) plots each item’s easiness
#' \eqn{\beta_i} as the starting point of a vertical arrow whose tip marks the
#' personalized easiness
#' \eqn{\delta_{ij} = \beta_i - d_{ij}}, where \eqn{d_{ij}} is the latent
#' distance taken from `distance_mat`. If `gamma` is supplied, distances are
#' scaled before computing deltas, i.e. \eqn{\delta_{ij} = \beta_i - \gamma d_{ij}}.
#' Uncertainty bounds in `distance_low`/`distance_up` are scaled by the same \eqn{\gamma}.
#' Arrows pointing **up** indicate the Arrows pointing **up** indicate the
#' item is *easier* for the focal respondent than for the average person,
#' whereas arrows pointing **down** indicate it is *harder*.
#'
#' @param alpha Numeric vector of length *N*. Posterior means (or draws) of
#'   person ability parameters.
#' @param beta Numeric vector of length *I*. Posterior means of item easiness
#'   parameters.
#' @param distance_mat Numeric matrix *N × I* of latent distances
#'   \eqn{d_{pi}} between persons and items.
#' @param gamma Optional numeric scalar used to multiplicatively rescale
#'   all distances (and `distance_low`/`distance_up`, if provided) before
#'   computing personalized easiness. Defaults to `NULL` (no rescaling).
#' @param alpha_lower,alpha_upper Optional numeric vectors (length *N*) giving
#'   lower/upper bounds (e.g., 95 % HDI) for each person’s \eqn{\alpha_j}.  If
#'   provided, the focal respondent’s interval is shaded in pink.
#' @param distance_low,distance_up Optional matrices the same size as
#'   `distance_mat` providing lower/upper HDI bounds for the distances.
#'   When both are supplied, dotted lines visualize the uncertainty in each
#'   personalised easiness.
#' @param item_group Optional character/factor vector of length *I* assigning
#*   each item to a group.  Used to color the item dots and activate a legend.
#' @param focal_id Integer index (1 ≤ `focal_id` ≤ *N*) of the respondent to
#'   highlight.  Default is `1`.
#' @param density_adjust Positive numeric scalar passed to
#'   `ggplot2::geom_density(adjust = ...)` to control the smoothness of the
#'   left-panel density estimate. Values > 1 increase the bandwidth (smoother
#'   curve); values < 1 decrease it (more detail). Default is `2`.
#' @param y_limits Optional numeric length-2 vector `c(min, max)` that fixes the
#'   y-axis range for both panels.
#'
#' @return A [`patchwork`](https://patchwork.data-imaginist.com) object
#'   containing two `ggplot2` panels.  The plot is also displayed as a side
#'   effect, so the returned object is mainly for further customisation.
#'
#' @seealso [`pip_fountain()`] for the complementary “fountain” layout, and
#'   [`interprofile()`]  for a thin wrapper that chooses between the two styles.
#'
#' @import patchwork
#' @importFrom rlang .data
#'
#' @examples
#' set.seed(42)
#' N <- 6; I <- 10
#' alpha <- rnorm(N)
#' beta  <- rnorm(I, sd = 0.7)
#' dist  <- abs(matrix(rnorm(N * I, sd = 0.8), N, I))  # fake distances
#'
#' # Basic waterfall plot for the first respondent
#' pip_waterfall(alpha, beta, gamma = 1.5, dist, focal_id = 2)
#'
#' # Add grouping and uncertainty bands
#' groups <- rep(c("A", "B"), length.out = I)
#' d_low  <- dist * 0.9; d_up <- dist * 1.1
#' a_l   <- alpha - 0.25; a_u <- alpha + 0.25
#' pip_waterfall(alpha, beta, gamma = 1, dist,
#'               alpha_lower = a_l, alpha_upper = a_u,
#'               distance_low = d_low, distance_up = d_up,
#'               item_group = groups, focal_id = 3)
#'
#' @export
pip_waterfall <- function(alpha, beta,
                          distance_mat,
                          gamma         = NULL,
                          alpha_lower   = NULL, alpha_upper   = NULL,
                          distance_low  = NULL, distance_up   = NULL,
                          item_group    = NULL,
                          focal_id      = 1,
                          density_adjust = 2,
                          y_limits = NULL)
{
        stopifnot(length(alpha)  == nrow(distance_mat),
                  length(beta)   == ncol(distance_mat))

        # --- NEW: optional distance scaling -----------------------------------
        if (!is.null(gamma)) {
                distance_mat <- gamma * distance_mat
                if (!is.null(distance_low)) distance_low <- gamma * distance_low
                if (!is.null(distance_up))  distance_up  <- gamma * distance_up
        }

        N <- length(alpha); I <- length(beta)
        y_max <- max(abs(c(alpha,beta,alpha_lower,alpha_upper,distance_low,distance_up)), na.rm = TRUE)
        y_lims <- if (is.null(y_limits)) c(-y_max, y_max) * 1.2 else y_limits

        ## 1 ─ person table ------------------------------------------------------
        persons <- data.frame(
                id      = 1:N,
                alpha   = -alpha,
                alpha_l = if (is.null(alpha_lower)) NA else -alpha_upper,
                alpha_u = if (is.null(alpha_upper)) NA else -alpha_lower
        )

        ## 2 ─ item colors / groups --------------------------------------------
        if (is.null(item_group)) {
                item_group   <- factor(rep("All", I))           # single dummy level
                fill_map     <- c(All = "grey60")               # one grey swatch
                show_legend  <- FALSE
        } else {
                item_group   <- factor(item_group, levels = unique(item_group))
                n_grp        <- nlevels(item_group)
                palette_cols <- scales::hue_pal()(n_grp)        # distinct hues
                fill_map     <- stats::setNames(palette_cols, levels(item_group))
                show_legend  <- TRUE
        }
        items <- data.frame(
                item_id = factor(1:I),
                beta    = beta,
                group   = item_group
        )

        ## 3 ─ convert distance → personalised easiness δ_ij --------------------
        delta_mat <- sweep(-distance_mat, 2, beta, FUN = "+")      # β_i − d_ij
        if (!is.null(distance_low))
                delta_low <- sweep(-distance_low, 2, beta, FUN = "+")
        if (!is.null(distance_up))
                delta_up  <- sweep(-distance_up,  2, beta, FUN = "+")

        ## 4 ─ long table --------------------------------------------------------
        delta_df <- stats::setNames(as.data.frame(delta_mat), paste0("Item_", 1:I))

        long <- delta_df |>
                dplyr::mutate(id = 1:N) |>
                tidyr::pivot_longer(
                        dplyr::starts_with("Item_"),
                        names_to  = "item_id",
                        values_to = "delta"
                ) |>
                dplyr::mutate(item_id = factor(sub("Item_", "", .data$item_id))) |>
                dplyr::left_join(items,   by = "item_id") |>
                dplyr::left_join(persons, by = "id")
        # attach HDI columns
        if (!is.null(distance_low) && !is.null(distance_up)) {
                long$delta_up <- as.vector(t(delta_low))
                long$delta_low  <- as.vector(t(delta_up))
        }

        sel <- dplyr::filter(long, .data$id == focal_id)[order(as.numeric(items$item_id)), ]

        ## 5 ─ plot limits -------------------------------------------------------
        #y_lims <- range(c(alpha, beta, delta_mat,
        #                  alpha_lower, alpha_upper,
        #                  if (exists("delta_low")) delta_low,
        #                  if (exists("delta_up"))  delta_up), na.rm = TRUE)

        ## 6 ─ left panel --------------------------------------------------------
        left <- ggplot2::ggplot(persons, ggplot2::aes(x = alpha)) +
                {if (!is.na(sel$alpha_l[1]) && !is.na(sel$alpha_u[1]))
                        ggplot2::annotate("rect",
                                          xmin = sel$alpha_l[1], xmax = sel$alpha_u[1],
                                          ymin = -Inf, ymax = Inf,
                                          fill = "pink", alpha = .25)} +
                ggplot2::geom_density(fill = "skyblue", alpha = .5,
                                      adjust = density_adjust) +
                ggplot2::geom_point(data = dplyr::distinct(sel, alpha),
                                    ggplot2::aes(x = alpha, y = 0),
                                    shape = 8, size = 4, color = "darkorange") +
                ggplot2::geom_vline(xintercept = sel$alpha[1], color = "red") +
                ggplot2::coord_flip() +
                ggplot2::scale_x_continuous(limits = y_lims, oob = scales::oob_keep) +
                ggplot2::theme_minimal() + ggplot2::labs(x = "", y = "") +
                ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                               axis.ticks.y = ggplot2::element_blank(),
                               panel.grid   = ggplot2::element_blank())

        ## 7 ─ right panel (Waterfall) ------------------------------------------
        arrow_scale <- ggplot2::scale_color_manual(
                values = c("#2b8cbe", "#e34a33"),
                labels = c("Easier", "Harder"),
                name   = NULL)

        right <- ggplot2::ggplot(items,
                                 ggplot2::aes(x = .data$item_id, y = .data$beta)) +
                ggplot2::geom_point(ggplot2::aes(fill = .data$group),
                                    shape = 21, size = 3, color = "black") +
                ggplot2::scale_fill_manual(values = fill_map,
                                           guide = if (show_legend) "legend" else "none",
                                           name  = if (show_legend) "Item group" else NULL)  +
                ggplot2::geom_segment(data = sel,
                                      ggplot2::aes(xend = .data$item_id, yend = .data$delta,
                                                   color = .data$delta > beta),
                                      arrow = ggplot2::arrow(length = grid::unit(.18, "cm")),
                                      linewidth = .9)

        if (!is.null(distance_low) && !is.null(distance_up)) {
                right <- right +
                        ggplot2::geom_segment(data = sel,
                                              ggplot2::aes(xend = .data$item_id,
                                                           y = delta_low,
                                                           yend = delta_up,
                                                           color = .data$delta > beta),
                                              linetype = "dotted",
                                              linewidth = .25)
        }

        right <- right +
                ggplot2::geom_hline(yintercept = sel$alpha[1], color = "red") +
                {if (!is.na(sel$alpha_l[1]) && !is.na(sel$alpha_u[1]))
                        ggplot2::annotate("rect",
                                          xmin = -Inf, xmax = Inf,
                                          ymin = sel$alpha_l[1], ymax = sel$alpha_u[1],
                                          fill = "pink", alpha = .25)} +
                arrow_scale +
                ggplot2::scale_y_continuous(position = "right",
                                            limits = y_lims,
                                            oob    = scales::oob_keep) +
                ggplot2::labs(x = "Item", y = "Easiness") +
                ggplot2::theme_minimal() +
                ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(),
                               axis.text.x = ggplot2::element_text(size = 10))

        ## 8 ─ combine -----------------------------------------------------------
        (left | right) +
                patchwork::plot_layout(widths = c(.65, 1), guides = "collect") &
                ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0, "cm"),
                               legend.position = "right")
}
