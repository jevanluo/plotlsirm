#' Posterior Interaction Profile — *Fountain* style
#'
#' Generates the **fountain** variant of a Posterior Interaction Profile (PIP)
#' plot.
#' The layout is identical to [`pip_waterfall()`] on the left (posterior
#' density for the focal respondent’s ability) but **inverts** the right-hand
#' panel: each item dot is placed at *- \eqn{\beta_i}{beta[i]}* (the “fountain base”) and an arrow
#' rises to the personalized easiness
#' \deqn{\delta_{ij} = \beta_i - d_{ij}}{delta[ij] = beta[i] - d[ij]}.
#' Distance is taken from `distance_mat`. If `gamma` is supplied, distances are
#' scaled before computing deltas, i.e. \eqn{\delta_{ij} = \beta_i - \gamma d_{ij}}{delta[ij] = beta[i] - gamma * d[ij]}.
#' Uncertainty bounds in `distance_low`/`distance_up` are scaled by the same \eqn{\gamma}{gamma}.
#' Arrows that extend **above** the base indicate the item is *easier* for the
#' respondent than average; arrows that fall short indicate it is *harder*.
#'
#' @param alpha Numeric vector of length *N*. Posterior means (or draws) of
#'   person ability parameters.
#' @param beta Numeric vector of length *I*. Posterior means of item easiness
#'   parameters.
#' @param distance_mat Numeric matrix *N × I* containing the latent distances
#'   \eqn{d_{ij}} between persons and items.
#' @param gamma Optional numeric scalar used to multiplicatively rescale
#'   all distances (and `distance_low`/`distance_up`, if provided) before
#'   computing personalized easiness. Defaults to `NULL` (no rescaling).
#' @param alpha_lower,alpha_upper Optional numeric vectors (length *N*)
#'   providing lower/upper posterior intervals (e.g., 95 % HDI) for each
#'   respondent’s \eqn{\alpha_j}.  The focal respondent’s band is shaded.
#' @param distance_low,distance_up Optional matrices matching `distance_mat`
#'   that give lower/upper HDI bounds for each distance.  When both are
#'   supplied, dotted vertical lines depict the uncertainty in every
#'   personalized easiness.
#' @param item_group Optional character/factor vector of length *I* defining
#'   item groupings.  Enables color coding and a legend.
#' @param focal_id Integer (1 ≤ `focal_id` ≤ *N*) selecting the respondent to
#'   highlight.  Defaults to the first row.
#' @param density_adjust Positive numeric scalar passed to
#'   `ggplot2::geom_density(adjust = ...)` to control the smoothness of the
#'   left-panel density estimate. Values > 1 increase the bandwidth (smoother
#'   curve); values < 1 decrease it (more detail). Default is `2`.
#' @param y_limits Optional numeric length-2 vector `c(min, max)` that fixes the
#'   y-axis range for both panels.
#'
#' @return A `patchwork` object containing the combined left‑ and right‑hand
#'   `ggplot2` panels.  The plot is automatically displayed; the value is
#'   returned invisibly for further tweaking.
#'
#' @seealso [`pip_waterfall()`] for the alternative “waterfall” framing, and
#'   [`interprofile()`] for a wrapper that switches between the two.
#'
#' @import patchwork
#' @importFrom rlang .data
#'
#' @examples
#' # Small simulated example -----------------------------------------
#' set.seed(42)
#' N <- 6; I <- 10
#' alpha <- rnorm(N)
#' beta  <- rnorm(I, sd = 0.7)
#' dist  <- abs(matrix(rnorm(N * I, sd = 0.8), N, I))  # fake distances
#'
#' # Plain fountain plot for respondent 2
#' pip_fountain(alpha, beta, gamma = 1.5, dist, focal_id = 2)
#'
#' # Fountain plot with item groups and uncertainty intervals
#' grp <- rep(c("MCQ", "Essay"), length.out = I)
#' d_lo <- pmax(dist - 0.2, 0);  d_up <- dist + 0.2
#' a_lo <- alpha - 0.3; a_up <- alpha + 0.3
#' pip_fountain(alpha, beta, gamma = 1, dist,
#'              alpha_lower = a_lo, alpha_upper = a_up,
#'              distance_low = d_lo, distance_up = d_up,
#'              item_group = grp, focal_id = 4)
#'
#' @export
pip_fountain <- function(alpha, beta,
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
                alpha   = alpha,
                alpha_l = if (is.null(alpha_lower)) NA else alpha_lower,
                alpha_u = if (is.null(alpha_upper)) NA else alpha_upper
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
        delta_mat <- sweep(distance_mat, 2, -beta, FUN = "+")      # β_i − d_ij
        if (!is.null(distance_low))
                delta_low <- sweep(distance_low, 2, -beta, FUN = "+")
        if (!is.null(distance_up))
                delta_up  <- sweep(distance_up,  2, -beta, FUN = "+")

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
                long$delta_up <- as.vector(t(delta_up))
                long$delta_low  <- as.vector(t(delta_low))
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
                labels = c("Harder", "Easier"),
                name   = NULL)

        # Update item positions for "fountain" layout: bottom at -beta
        items$beta_plot <- -items$beta
        sel$beta_plot <- -sel$beta

        right <- ggplot2::ggplot(items,
                                 ggplot2::aes(x = .data$item_id, y = .data$beta_plot)) +
                ggplot2::geom_point(ggplot2::aes(fill = .data$group),
                                    shape = 21, size = 3, color = "black") +
                ggplot2::scale_fill_manual(values = fill_map,
                                           guide = if (show_legend) "legend" else "none",
                                           name  = if (show_legend) "Item group" else NULL) +

                # Arrows from -beta (base) to delta (height)
                ggplot2::geom_segment(data = sel,
                                      ggplot2::aes(x = .data$item_id,
                                                   y = .data$beta_plot,
                                                   xend = .data$item_id,
                                                   yend = .data$delta,
                                                   color = .data$delta > -beta),
                                      arrow = ggplot2::arrow(length = grid::unit(.18, "cm")),
                                      linewidth = .9) +

                # Uncertainty intervals as vertical dotted lines
                {if (!is.null(distance_low) && !is.null(distance_up))
                        ggplot2::geom_segment(data = sel,
                                              ggplot2::aes(x = .data$item_id,
                                                           y = delta_low,
                                                           xend = .data$item_id,
                                                           yend = delta_up,
                                                           color = .data$delta > -beta),
                                              linetype = "dotted",
                                              linewidth = .25)} +

                # Reference line for person alpha
                ggplot2::geom_hline(yintercept = sel$alpha[1], color = "red") +

                # Optional uncertainty band for alpha
                {if (!is.na(sel$alpha_l[1]) && !is.na(sel$alpha_u[1]))
                        ggplot2::annotate("rect",
                                          xmin = -Inf, xmax = Inf,
                                          ymin = sel$alpha_l[1], ymax = sel$alpha_u[1],
                                          fill = "pink", alpha = .25)} +

                ggplot2::labs(x = "Item", y = "Difficulty") +
                arrow_scale +
                ggplot2::scale_y_continuous(position = "right",
                                            limits = y_lims,
                                            oob    = scales::oob_keep) +
                ggplot2::theme_minimal() +
                ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(),
                               axis.text.x = ggplot2::element_text(size = 10))


        ## 8 ─ combine -----------------------------------------------------------
        (left | right) +
                patchwork::plot_layout(widths = c(.65, 1), guides = "collect") &
                ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0, "cm"),
                               legend.position = "right")
}
