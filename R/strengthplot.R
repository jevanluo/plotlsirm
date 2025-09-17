#' Item-strength profile for a single person
#'
#' For a chosen respondent (`person_index`) this function plots the **strength**
#' (likelihood of endorsement) for every item, defined as
#' \deqn{\exp(-\gamma d_{ij})}{exp(-gamma * d_ij)}, where \eqn{d_{ij}} is the
#' Euclidean distance between the person’s latent position \eqn{z_j}{z[j]} and each
#' item position \eqn{w_i}{w[i]}. When `z` and `w` are supplied as *lists* of
#' matrices (posterior draws), the function summarises the distribution of
#' strengths with medians and a `ci_level` credible interval. Bars can be
#' coloured by an item grouping factor, reordered by decreasing strength, and
#' displayed either vertically or horizontally.
#'
#' When **no `item_group` is provided**, bars are **colour-mapped by a
#' similarity gradient** (low → high) by default. You can disable this behaviour
#' and use a single fill colour instead via `use_gradient = FALSE`.
#'
#' @param z A numeric matrix (*N* × *d*) of person coordinates **or** a *list*
#'   of such matrices representing posterior draws.
#' @param w A numeric matrix (*I* × *d*) of item coordinates **or** a *list* of
#'   such matrices, matching the structure of `z`.
#' @param person_index Integer giving the row of `z` (or each draw in `z`)
#'   corresponding to the focal respondent.
#' @param gamma Positive numeric scalar controlling the decay of strength with
#'   distance; default is `1`.
#' @param item_group Optional character/factor vector of length *I* assigning
#'   each item to a group for colour coding and legend.
#' @param item_names Optional character vector of item labels. If `NULL`
#'   defaults to `"I1"`, `"I2"`, … .
#' @param ci_level Width of the credible interval (between 0 and 1) when
#'   posterior draws are given. Ignored for a single point estimate.
#' @param reorder Logical. Reorder items on the axis by decreasing strength?
#'   Default `FALSE`.
#' @param vertical Logical. `TRUE` (default) draws vertical bars; `FALSE`
#'   flips the axes for a horizontal layout.
#' @param title Optional character string to appear as the plot title.
#'
#' @param use_gradient Logical. When `item_group` is `NULL`, colour bars by a
#'   **strength gradient** (low → high)? Default `TRUE`.
#' @param gradient_low,gradient_high Colours for the gradient when
#'   `use_gradient = TRUE`. Defaults `"#d9f0d3"` (low) to `"#1b7837"` (high).
#' @param show_gradient_legend Logical. Show a legend for the gradient (only
#'   when `item_group` is `NULL` and `use_gradient = TRUE`)? Default `TRUE`.
#' @param single_fill_color Single fill colour used when `use_gradient = FALSE`
#'   and `item_group` is `NULL`. Default `"steelblue"`.
#'
#' @return (Invisibly) a `ggplot` object containing the bar plot. The plot is
#'   also printed.
#'
#' @import ggplot2
#' @importFrom stats median quantile
#' @importFrom rlang .data
#'
#' @examples
#' set.seed(1)
#' z  <- matrix(rnorm(40), ncol = 2)   # 20 persons
#' w  <- matrix(rnorm(30), ncol = 2)   # 15 items
#'
#' ## 1) Point-estimate strengths for person 5 (default gradient, ungrouped)
#' strengthplot(z, w, person_index = 5, gamma = 2,
#'              title = "Strengths for person 5 (gradient)")
#'
#' ## 2) Turn off gradient and use a single colour
#' strengthplot(z, w, person_index = 5, gamma = 2,
#'              use_gradient = FALSE, single_fill_color = "tomato",
#'              title = "Strengths for person 5 (single colour)")
#'
#' ## 3) Posterior example with credible intervals and item groups
#' draws_z <- replicate(50, z + matrix(rnorm(length(z), sd = 0.1),
#'                                     nrow(z), ncol(z)), simplify = FALSE)
#' draws_w <- replicate(50, w + matrix(rnorm(length(w), sd = 0.1),
#'                                     nrow(w), ncol(w)), simplify = FALSE)
#' grp <- rep(c("Core", "Peripheral"), length.out = nrow(w))
#' strengthplot(draws_z, draws_w, person_index = 3,
#'              item_group = grp, ci_level = 0.9, vertical = FALSE,
#'              title = "Posterior strength profile for respondent 3")
#'
#' @export
strengthplot <- function(
                z, w,
                person_index,
                gamma       = 1,
                item_group  = NULL,
                item_names  = NULL,
                ci_level    = 0.95,
                reorder     = FALSE,
                vertical    = TRUE,
                title       = NULL,
                # --- new colouring controls for ungrouped plots ---
                use_gradient         = TRUE,
                gradient_low         = "#d9f0d3",
                gradient_high        = "#1b7837",
                show_gradient_legend = TRUE,
                single_fill_color    = "steelblue"
) {
        # validate gamma
        if (!is.numeric(gamma) || length(gamma) != 1 || !is.finite(gamma) || gamma <= 0)
                stop("`gamma` must be a positive finite scalar.")

        # internal strength fun
        str_fun <- function(d) exp(-gamma * d)
        is_post <- is.list(z) && is.list(w)

        # helper distance
        dist_vec_mat <- function(vec, mat) {
                sqrt(rowSums((mat - matrix(vec, nrow(mat), ncol(mat), byrow = TRUE))^2))
        }

        build_single_df <- function(zm, wm) {
                I <- nrow(wm)
                if (!is.null(item_group) && length(item_group) != I)
                        stop("`item_group` length must equal number of items.")
                nm <- if (is.null(item_names)) paste0("I", seq_len(I)) else item_names
                d  <- dist_vec_mat(zm[person_index, ], wm)
                s  <- str_fun(d)
                ord <- if (reorder) nm[order(-s)] else nm
                data.frame(item = factor(nm, levels = ord),
                           strength = s,
                           group = if (is.null(item_group)) NA_character_ else item_group,
                           stringsAsFactors = FALSE)
        }

        if (!is_post) {
                df <- build_single_df(as.matrix(z), as.matrix(w))

                if (is.null(item_group)) {
                        if (isTRUE(use_gradient)) {
                                p <- ggplot(df, aes(x = .data$item, y = .data$strength, fill = .data$strength)) +
                                        geom_col() +
                                        scale_fill_gradient(
                                                name  = "Strength",
                                                low   = gradient_low,
                                                high  = gradient_high,
                                                guide = if (isTRUE(show_gradient_legend)) "legend" else "none"
                                        )
                        } else {
                                p <- ggplot(df, aes(x = .data$item, y = .data$strength)) +
                                        geom_col(fill = single_fill_color)
                        }
                } else {
                        p <- ggplot(df, aes(x = .data$item, y = .data$strength, fill = .data$group)) +
                                geom_col() +
                                scale_fill_brewer(palette = "Set2", name = "Item group")
                }

        } else {
                # posterior draws
                M <- length(z)
                N <- nrow(as.matrix(z[[1]]))
                I <- nrow(as.matrix(w[[1]]))
                if (person_index > N) stop("`person_index` exceeds number of persons in `z`.")
                if (!is.null(item_group) && length(item_group) != I)
                        stop("`item_group` length must equal number of items.")

                str_mat <- matrix(NA_real_, I, M)
                for (m in seq_len(M)) {
                        zm <- as.matrix(z[[m]])[person_index, ]
                        wm <- as.matrix(w[[m]])
                        str_mat[, m] <- str_fun(dist_vec_mat(zm, wm))
                }

                med <- apply(str_mat, 1, stats::median, na.rm = TRUE)
                lwr <- apply(str_mat, 1, stats::quantile, probs = (1 - ci_level) / 2, na.rm = TRUE)
                upr <- apply(str_mat, 1, stats::quantile, probs = 1 - (1 - ci_level) / 2, na.rm = TRUE)

                nm  <- if (is.null(item_names)) paste0("I", seq_len(I)) else item_names
                ord <- if (reorder) nm[order(-med)] else nm

                df <- data.frame(item  = factor(nm, levels = ord),
                                 median = med,
                                 lower  = lwr,
                                 upper  = upr,
                                 group  = if (is.null(item_group)) NA_character_ else item_group,
                                 stringsAsFactors = FALSE)

                if (is.null(item_group)) {
                        if (isTRUE(use_gradient)) {
                                p <- ggplot(df, aes(x = .data$item, y = .data$median, fill = .data$median)) +
                                        geom_col() +
                                        geom_errorbar(aes(ymin = .data$lower, ymax = .data$upper), width = 0.25) +
                                        scale_fill_gradient(
                                                name  = "Strength",
                                                low   = gradient_low,
                                                high  = gradient_high,
                                                guide = if (isTRUE(show_gradient_legend)) "legend" else "none"
                                        )
                        } else {
                                p <- ggplot(df, aes(x = .data$item, y = .data$median)) +
                                        geom_col(fill = single_fill_color) +
                                        geom_errorbar(aes(ymin = .data$lower, ymax = .data$upper), width = 0.25)
                        }
                } else {
                        p <- ggplot(df, aes(x = .data$item, y = .data$median, fill = .data$group)) +
                                geom_col() +
                                geom_errorbar(aes(ymin = .data$lower, ymax = .data$upper), width = 0.25) +
                                scale_fill_brewer(palette = "Set2", name = "Item group")
                }
        }

        # orientation
        if (!vertical) p <- p + coord_flip()

        # axis & theme
        p <- p +
                labs(x = NULL, y = "Likelihood of endorsement") +
                scale_y_continuous(limits = c(0, 1),
                                   breaks = seq(0, 1, 0.2),
                                   expand = expansion(mult = c(0, 0.02))) +
                theme(
                        panel.background = element_rect(fill = "#F5F5F5"),
                        panel.grid.major = element_line(color = "white"),
                        panel.border     = element_rect(color = "grey60", fill = NA, linewidth = 0.8),
                        axis.text.x      = element_text(angle = if (vertical) 90 else 0, vjust = 0.5, hjust = 1),
                        axis.ticks.y     = element_line(color = "grey50")
                )

        if (!is.null(title)) p <- p + ggtitle(title)

        print(p)
        invisible(p)
}
