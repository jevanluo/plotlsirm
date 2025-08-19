#' Similarity profile for a focal item
#'
#' Plots the similarity between one **focal item** and every other item in latent
#' space, optionally including posterior uncertainty bands when a list of draws
#' is supplied.  Similarity is defined as
#' \deqn{\exp(-\gamma\,d_{ij})}
#' where \eqn{d_{ij}} is the
#' Euclidean distance between items *i* and *j*.  Bars can be color-coded by a
#' grouping factor, reordered by decreasing similarity, displayed horizontally
#' or vertically, and annotated with credible intervals.
#'
#' @param w Numeric **matrix** (*I* × *d*) of item coordinates **or** a *list*
#'   of such matrices (posterior draws).  When a list is given the function
#'   summarises similarity across draws and plots medians with
#'   `ci_level` credible intervals.
#' @param focal_item Index (integer) **or** name (character) of the item whose
#'   similarity profile is to be displayed.
#' @param gamma Positive numeric scalar controlling how quickly similarity
#'   decays with distance.  Default is `1`.
#' @param item_group Optional character/factor vector of length *I* indicating
#'   group membership for each item.  Used for bar colors and legend.
#' @param item_names Optional character vector of item labels (length *I*).
#'   Defaults to `"I1"`, `"I2"`, … if `NULL`.
#' @param ci_level Numeric between 0 and 1 giving the width of the credible
#'   interval when `w` is a posterior list.  Ignored for a single draw.
#' @param reorder Logical.  Reorder items on the axis by decreasing similarity
#'   to the focal item?  Default `FALSE`.
#' @param vertical Logical.  `TRUE` (default) plots vertical bars; `FALSE`
#'   flips the axes for a horizontal layout.
#' @param title Optional character string added as the plot title.
#'
#' @param use_gradient Logical. When `item_group` is `NULL`, color bars by a
#'   **similarity gradient** (low→high). Default `TRUE`.
#' @param gradient_low,gradient_high Colors for the similarity gradient when
#'   `use_gradient = TRUE`. Defaults `"#d9f0d3"` (low) to `"#1b7837"` (high).
#' @param show_gradient_legend Logical. Show legend for the similarity gradient
#'   (only when `item_group` is `NULL` and `use_gradient = TRUE`)? Default `TRUE`.
#' @param single_fill_color Single fill color when `use_gradient = FALSE`
#'   and `item_group` is `NULL`. Default `"steelblue"`.
#'
#' @return (Invisibly) a `ggplot` object.  The plot is also drawn as a side
#'   effect.
#'
#' @import ggplot2
#' @importFrom stats median quantile
#' @importFrom rlang .data
#'
#' @examples
#' set.seed(1)
#' w  <- matrix(rnorm(40), ncol = 2)   # 20 items
#' gp <- sample(c("Math", "Verbal"), nrow(w), replace = TRUE)
#'
#' ## 1) Single estimate, default gradient (ungrouped)
#' itemsimilarity(w, focal_item = 3, gamma = 2,
#'                title = "Similarity to item 3 (gradient)")
#'
#' ## 2) Single estimate, turn off gradient and use one color (ungrouped)
#' itemsimilarity(w, focal_item = 3, gamma = 2,
#'                use_gradient = FALSE, single_fill_color = "tomato",
#'                title = "Similarity to item 3 (single color)")
#'
#' ## 3) Grouped bars (gradient ignored because groups are used)
#' itemsimilarity(w, focal_item = 3, gamma = 2, item_group = gp,
#'                title = "Similarity to item 3 (grouped)")
#'
#' ## 4) Posterior list with credible intervals (ungrouped, gradient)
#' draws <- replicate(100, w + matrix(rnorm(length(w), sd = 0.1),
#'                                    nrow(w), ncol(w)), simplify = FALSE)
#' itemsimilarity(draws, focal_item = "I10", ci_level = 0.9,
#'                vertical = FALSE, show_gradient_legend = FALSE)
#'
#' @export
itemsimilarity <- function(
                w,
                focal_item,
                gamma       = 1,
                item_group  = NULL,
                item_names  = NULL,
                ci_level    = 0.95,
                reorder     = FALSE,
                vertical    = TRUE,
                title       = NULL,
                use_gradient         = TRUE,
                gradient_low         = "#d9f0d3",
                gradient_high        = "#1b7837",
                show_gradient_legend = TRUE,
                single_fill_color    = "steelblue"
) {
        if (!is.numeric(gamma) || length(gamma) != 1 || gamma <= 0)
                stop("`gamma` must be a positive scalar.")

        sim_fun  <- function(d) exp(-gamma * d)
        posterior <- is.list(w)

        dist_rows <- function(mat, i) {
                sqrt(rowSums((mat - matrix(mat[i, ], nrow(mat), ncol(mat), byrow = TRUE))^2))
        }

        ## --- setup --------------------------------------------------------
        I <- if (!posterior) nrow(as.matrix(w)) else nrow(as.matrix(w[[1]]))
        if (is.null(item_names)) item_names <- paste0("I", seq_len(I))
        if (!is.null(item_group) && length(item_group) != I)
                stop("item_group length mismatch")

        focal_index <- if (is.character(focal_item)) match(focal_item, item_names) else focal_item
        if (is.na(focal_index) || focal_index < 1L || focal_index > I)
                stop("focal_item not found / out of bounds")

        ## ================================================================
        ## SINGLE DRAW
        ## ================================================================
        if (!posterior) {
                wmat <- as.matrix(w)
                sim  <- sim_fun(dist_rows(wmat, focal_index))
                sim[focal_index] <- NA

                keep <- -focal_index                    # drop focal itself
                nm   <- item_names[keep]
                grp  <- if (is.null(item_group)) NA_character_ else item_group[keep]
                sim  <- sim[keep]
                ord  <- if (reorder) nm[order(-sim)] else nm

                df <- data.frame(item = factor(nm, levels = ord),
                                 similarity = sim,
                                 group = grp)

                if (is.null(item_group)) {
                        if (isTRUE(use_gradient)) {
                                p <- ggplot2::ggplot(df,
                                                     ggplot2::aes(x = .data$item, y = .data$similarity, fill = .data$similarity)) +
                                        ggplot2::geom_col() +
                                        ggplot2::scale_fill_gradient(
                                                name = "Similarity",
                                                low  = gradient_low,
                                                high = gradient_high,
                                                guide = if (isTRUE(show_gradient_legend)) "legend" else "none"
                                        )
                        } else {
                                p <- ggplot2::ggplot(df,
                                                     ggplot2::aes(x = .data$item, y = .data$similarity)) +
                                        ggplot2::geom_col(fill = single_fill_color)
                        }
                } else {
                        p <- ggplot2::ggplot(df,
                                             ggplot2::aes(x = .data$item, y = .data$similarity, fill = .data$group)) +
                                ggplot2::geom_col() +
                                ggplot2::scale_fill_brewer(palette = "Set3", name = "Item group")
                }

                ## ================================================================
                ## POSTERIOR DRAWS
                ## ================================================================
        } else {
                M <- length(w)
                sim_mat <- matrix(NA_real_, I, M)
                for (m in seq_len(M)) {
                        w_m <- as.matrix(w[[m]])
                        s   <- sim_fun(dist_rows(w_m, focal_index))
                        s[focal_index] <- NA
                        sim_mat[, m] <- s
                }

                med <- apply(sim_mat, 1, stats::median, na.rm = TRUE)
                lwr <- apply(sim_mat, 1, stats::quantile, probs = (1 - ci_level) / 2, na.rm = TRUE)
                upr <- apply(sim_mat, 1, stats::quantile, probs = 1 - (1 - ci_level) / 2, na.rm = TRUE)

                keep <- which(!is.na(med))
                nm   <- item_names[keep]
                grp  <- if (is.null(item_group)) NA_character_ else item_group[keep]
                ord  <- if (reorder) nm[order(-med[keep])] else nm

                df <- data.frame(item  = factor(nm, levels = ord),
                                 median = med[keep],
                                 lower  = lwr[keep],
                                 upper  = upr[keep],
                                 group  = grp)

                if (is.null(item_group)) {
                        if (isTRUE(use_gradient)) {
                                p <- ggplot2::ggplot(df,
                                                     ggplot2::aes(x = .data$item, y = .data$median, fill = .data$median)) +
                                        ggplot2::geom_col() +
                                        ggplot2::geom_errorbar(
                                                ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
                                                width = 0.25
                                        ) +
                                        ggplot2::scale_fill_gradient(
                                                name = "Similarity",
                                                low  = gradient_low,
                                                high = gradient_high,
                                                guide = if (isTRUE(show_gradient_legend)) "legend" else "none"
                                        )
                        } else {
                                p <- ggplot2::ggplot(df,
                                                     ggplot2::aes(x = .data$item, y = .data$median)) +
                                        ggplot2::geom_col(fill = single_fill_color) +
                                        ggplot2::geom_errorbar(
                                                ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
                                                width = 0.25
                                        )
                        }
                } else {
                        p <- ggplot2::ggplot(df,
                                             ggplot2::aes(x = .data$item, y = .data$median, fill = .data$group)) +
                                ggplot2::geom_col() +
                                ggplot2::geom_errorbar(
                                        ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
                                        width = 0.25
                                ) +
                                ggplot2::scale_fill_brewer(palette = "Set3", name = "Item group")
                }
        }

        ## --- orientation -------------------------------------------------
        if (!vertical) {
                p <- p + ggplot2::coord_flip() +
                        ggplot2::labs(x = NULL, y = NULL)  # relabel below
        }

        ## --- axis labels & scales ---------------------------------------
        if (vertical) {
                p <- p +
                        ggplot2::labs(x = NULL, y = "Similarity") +
                        ggplot2::scale_y_continuous(
                                limits = c(0, 1),
                                breaks  = seq(0, 1, 0.2),
                                expand  = ggplot2::expansion(mult = c(0, 0.02))
                        )
        } else {   # coord_flip() moves y → x visually; keep scale on y
                p <- p +
                        ggplot2::labs(y = NULL, x = "Similarity") +
                        ggplot2::scale_y_continuous(
                                limits = c(0, 1),
                                breaks  = seq(0, 1, 0.2),
                                expand  = ggplot2::expansion(mult = c(0, 0.02))
                        )
        }

        ## --- theme & title ----------------------------------------------
        p <- p +
                ggplot2::theme(
                        panel.background = ggplot2::element_rect(fill = "#F5F5F5"),
                        panel.grid.major = ggplot2::element_line(color = "white"),
                        panel.border     = ggplot2::element_rect(color = "grey60", fill = NA, linewidth = 0.8),
                        axis.text.x      = ggplot2::element_text(angle = if (vertical) 90 else 0,
                                                                 vjust = 0.5, hjust = 1),
                        axis.ticks.y     = ggplot2::element_line(color = "grey50")
                )

        if (!is.null(title)) p <- p + ggplot2::ggtitle(title)

        print(p)
        invisible(p)
}
