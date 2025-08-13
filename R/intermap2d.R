#' 2‑D latent‑space interaction map (persons vs. items)
#'
#' Produces a publication‑ready scatter plot of person and item coordinates in a
#' two‑dimensional latent space.  The function is highly configurable: points
#' can be grouped and color‑coded, sized or alpha‑scaled by auxiliary
#' statistics (e.g., \eqn{\theta_p}, \eqn{\beta_i}), labelled, reshaped, and
#' bounded by user‑supplied axis limits.  All graphics are rendered with
#' **ggplot2** and the resulting plot object is returned invisibly after being
#' printed, so you can further modify it if desired.
#'
#' @param z Numeric matrix (*N* × 2). Person coordinates (rows = persons,
#'   cols = dimensions).
#' @param w Numeric matrix (*I* × 2). Item coordinates (rows = items,
#'   cols = dimensions).
#' @param item_group,person_group Optional character/factor vectors defining
#'   group membership for items and persons, respectively.  When supplied,
#'   distinct hues are assigned per group and a legend is shown.
#' @param person_colors,item_colors Optional character vectors of explicit
#'   colors (one per row of `z` or `w`).  Overrides the *_group* mechanism and
#'   suppresses the legend.
#' @param beta,theta Optional numeric vectors of length *I* and *N* used for
#'   mapping item/person statistics to point size or transparency.  Ignored
#'   unless the corresponding `scale_*` argument is `TRUE`.
#' @param scale_z_size,scale_w_size Logical.  Rescale `theta` / `beta` to point
#'   sizes using `z_size_range` / `w_size_range`.
#' @param z_size_range,w_size_range Length‑2 numeric vectors giving the minimum
#'   and maximum point sizes for persons/items when size scaling is enabled.
#' @param alpha_z_scale,alpha_w_scale Logical.  Rescale `theta` / `beta` to
#'   point transparency using `z_alpha_range` / `w_alpha_range`.
#' @param z_alpha_range,w_alpha_range Length‑2 numeric vectors giving the
#'   minimum and maximum alpha (opacity) used when alpha scaling is enabled.
#' @param gamma Numeric scalar.  Optional global stretch factor applied to both
#'   person and item coordinates (useful for visualising different values of a
#'   distance weight \eqn{\gamma}).
#' @param show_ticks Logical.  Draw axis ticks and numeric labels (`TRUE`) or
#'   hide them (`FALSE`, default).
#' @param xlim_range,ylim_range Optional length‑2 numeric vectors that override
#'   the automatically determined symmetric axis limits.
#' @param itemlabels,personlabels Character vectors supplying custom text
#'   labels for items/persons.  Defaults to `"I1"`, `"I2"`, … and
#'   `"P1"`, `"P2"`, … when `NULL`.
#' @param figuretitle Optional character string added as the plot title.
#' @param z_shape,w_shape Point shapes for persons and items (see
#'   `?ggplot2::geom_point` for codes).  Default: 16 (solid circle) for persons
#'   and 17 (solid triangle) for items.
#' @param z_border_width Numeric.  Stroke width for person points (only visible
#'   for filled shapes that have a border aesthetic).
#' @param z_label_size,w_label_size Numeric text sizes for person/item labels
#'   (only used when `show_*_labels` is `TRUE`).
#' @param show_z_labels,show_w_labels Logical flags controlling whether text
#'   labels for persons/items are displayed.
#' @param show_z_shapes,show_w_shapes Logical flags controlling whether the
#'   actual point symbols for persons/items are drawn.  For example, you may
#'   wish to show only item shapes and label persons.
#'
#' @return (Invisibly) a `ggplot` object containing the constructed plot.  The
#'   plot is also sent to the current graphics device as a side effect.
#'
#' @import ggplot2
#' @importFrom scales hue_pal
#' @importFrom rlang .data
#'
#' @examples
#' set.seed(1)
#' z <- matrix(rnorm(40), ncol = 2)   # 20 persons
#' w <- matrix(rnorm(30), ncol = 2)   # 15 items
#'
#' ## Basic map with default colors
#' intermap2d(z, w)
#'
#' ## Grouped items vs. grouped persons, with point size ∝ theta and labels size ∝ beta
#' theta <- rnorm(nrow(z))
#' beta  <- rnorm(nrow(w))
#' intermap2d(
#' z, w,
#' person_group  = rep(c("A", "B"), each = 10),
#' item_group    = rep(letters[1:3], length.out = nrow(w)),
#' theta         = theta, beta = beta,
#' scale_z_size  = TRUE, scale_w_size = TRUE,
#' figuretitle   = "Latent-space interaction map"
#' )
#'
#' @export
intermap2d <- function(
                z, w,
                # ----- latent‑space stretch -----
                gamma           = NULL,
                # ----- color control -----
                person_group    = NULL,
                item_group      = NULL,
                person_colors   = NULL,   # vector length nrow(z)
                item_colors     = NULL,   # vector length nrow(w)
                # ----- size scaling -----
                theta           = NULL, beta = NULL,
                scale_z_size    = FALSE, scale_w_size = FALSE,
                z_size_range    = c(2, 6), w_size_range = c(2, 8),
                # ----- alpha scaling -----
                alpha_z_scale   = FALSE, alpha_w_scale = FALSE,
                z_alpha_range   = c(0.30, 1), w_alpha_range = c(0.30, 1),
                # ----- axes & ticks -----
                show_ticks      = FALSE,
                xlim_range      = NULL, ylim_range = NULL,
                # ----- labels & title -----
                itemlabels      = NULL, personlabels = NULL,
                figuretitle     = NULL,
                # ----- shapes / sizes -----
                z_shape         = 16, w_shape = 17,
                z_border_width  = 0.5,
                z_label_size    = 4,  w_label_size = 4,
                show_z_labels   = FALSE, show_w_labels = TRUE,
                show_z_shapes   = TRUE,  show_w_shapes = FALSE
) {

        ## ---------- data prep ----------
        z_df <- as.data.frame(z); colnames(z_df) <- c("x", "y")
        w_df <- as.data.frame(w); colnames(w_df) <- c("x", "y")

        if (!is.null(gamma)) {
                z_df[, 1:2] <- z_df[, 1:2] * gamma
                w_df[, 1:2] <- w_df[, 1:2] * gamma
        }

        if (is.null(personlabels)) personlabels <- paste0("P", seq_len(nrow(z_df)))
        if (is.null(itemlabels))   itemlabels   <- paste0("I", seq_len(nrow(w_df)))

        ## ---------- size & alpha ----------
        if (scale_z_size) z_df$sz <- rescale_to_range(theta, to = z_size_range)
        if (scale_w_size) w_df$sz <- rescale_to_range(beta,  to = w_size_range)
        z_df$al <- if (alpha_z_scale) rescale_to_range(theta, to = z_alpha_range) else 0.8
        w_df$al <- if (alpha_w_scale) rescale_to_range(beta,  to = w_alpha_range) else 1

        ## ---------- coloring ----------
        default_z_col <- "grey60"   # persons muted
        default_w_col <- "red"      # items highlighted

        person_legend <- FALSE
        item_legend   <- FALSE

        # persons ----
        if (!is.null(person_colors)) {
                stopifnot(length(person_colors) == nrow(z_df))
                z_df$col <- person_colors
        } else if (!is.null(person_group)) {
                z_df$grp <- factor(person_group)
                z_cols <- scales::hue_pal()(nlevels(z_df$grp))
                names(z_cols) <- levels(z_df$grp)
                person_legend <- TRUE
        } else {
                z_df$col <- default_z_col
        }

        # items ----
        if (!is.null(item_colors)) {
                stopifnot(length(item_colors) == nrow(w_df))
                w_df$col <- item_colors
        } else if (!is.null(item_group)) {
                w_df$grp <- factor(item_group)
                w_cols <- scales::hue_pal()(nlevels(w_df$grp))
                names(w_cols) <- levels(w_df$grp)
                item_legend <- TRUE
        } else {
                w_df$col <- default_w_col
        }

        use_legend <- person_legend || item_legend

        ## ---------- base plot ----------
        p <- ggplot2::ggplot() + ggplot2::coord_fixed()

        # limits
        rng <- max(abs(c(z_df$x, z_df$y, w_df$x, w_df$y)), na.rm = TRUE)
        if (is.null(xlim_range)) xlim_range <- c(-rng, rng)
        if (is.null(ylim_range)) ylim_range <- c(-rng, rng)
        p <- p + ggplot2::xlim(xlim_range) + ggplot2::ylim(ylim_range)

        p <- p + ggplot2::theme(
                panel.background=ggplot2::element_rect(fill="#F5F5F5"),
                panel.grid.major=ggplot2::element_line(color="white"),
                panel.border=ggplot2::element_rect(color="grey60", fill=NA, linewidth=0.8),
                #panel.grid       = ggplot2::element_blank(),
                axis.title       = ggplot2::element_blank(),
                axis.ticks       = if (show_ticks) ggplot2::element_line() else ggplot2::element_blank(),
                axis.text        = if (show_ticks) ggplot2::element_text(size = 12) else ggplot2::element_blank(),
                legend.position  = if (use_legend) "right" else "none"
        )

        ## ---------- draw persons ----------
        if (show_z_shapes) {
                if (person_legend) {
                        p <- p + ggplot2::geom_point(
                                data = z_df,
                                ggplot2::aes(.data$x, .data$y, color = .data$grp,
                                             size = if (scale_z_size) .data$sz else NULL,
                                             alpha = .data$al),
                                shape = z_shape, stroke = z_border_width)
                } else {
                        p <- p + ggplot2::geom_point(
                                data = z_df,
                                ggplot2::aes(.data$x, .data$y,
                                             size  = if (scale_z_size) .data$sz else NULL,
                                             alpha = .data$al),
                                shape = z_shape, stroke = z_border_width,
                                color = z_df$col)
                }
        }

        ## ---------- draw items ----------
        if (show_w_shapes) {
                if (item_legend) {
                        p <- p + ggplot2::geom_point(
                                data = w_df,
                                ggplot2::aes(.data$x, .data$y, color = .data$grp,
                                             size = if (scale_w_size) .data$sz else NULL,
                                             alpha = .data$al),
                                shape = w_shape)
                } else {
                        p <- p + ggplot2::geom_point(
                                data = w_df,
                                ggplot2::aes(.data$x, .data$y,
                                             size  = if (scale_w_size) .data$sz else NULL,
                                             alpha = .data$al),
                                shape = w_shape,
                                color = w_df$col)
                }
        }

        ## ---------- legend colors ----------
        if (use_legend) {
                manual_vals <- c(if (person_legend) z_cols else NULL,
                                 if (item_legend)  w_cols else NULL)
                p <- p + ggplot2::scale_color_manual(values = manual_vals, name = "Group")
        } else {
                p <- p + ggplot2::scale_color_identity()
        }

        if (scale_z_size || scale_w_size) p <- p + ggplot2::scale_size_identity()
        if (alpha_z_scale || alpha_w_scale) p <- p + ggplot2::scale_alpha_identity(guide = "none")
        p <- p + ggplot2::guides(alpha = "none")

        ## ---------- labels ----------
        ## ---------- labels ----------
        if (show_z_labels) {
                if (person_legend) {
                        if (scale_z_size) {
                                p <- p + ggplot2::geom_text(
                                        data = z_df,
                                        ggplot2::aes(.data$x, .data$y, label = personlabels,
                                                     color = .data$grp, size = .data$sz),
                                        fontface = "bold")
                        } else {
                                p <- p + ggplot2::geom_text(
                                        data = z_df,
                                        ggplot2::aes(.data$x, .data$y, label = personlabels, color = .data$grp),
                                        size = z_label_size, fontface = "bold")
                        }
                } else {
                        if (scale_z_size) {
                                p <- p + ggplot2::geom_text(
                                        data = z_df,
                                        ggplot2::aes(.data$x, .data$y, label = personlabels, size = .data$sz),
                                        fontface = "bold", color = z_df$col)
                        } else {
                                p <- p + ggplot2::geom_text(
                                        data = z_df,
                                        ggplot2::aes(.data$x, .data$y, label = personlabels),
                                        size = z_label_size, fontface = "bold", color = z_df$col)
                        }
                }
        }

        if (show_w_labels) {
                if (item_legend) {
                        if (scale_w_size) {
                                p <- p + ggplot2::geom_text(
                                        data = w_df,
                                        ggplot2::aes(.data$x, .data$y, label = itemlabels,
                                                     color = .data$grp, size = .data$sz),
                                        fontface = "bold")
                        } else {
                                p <- p + ggplot2::geom_text(
                                        data = w_df,
                                        ggplot2::aes(.data$x, .data$y, label = itemlabels, color = .data$grp),
                                        size = w_label_size, fontface = "bold")
                        }
                } else {
                        if (scale_w_size) {
                                p <- p + ggplot2::geom_text(
                                        data = w_df,
                                        ggplot2::aes(.data$x, .data$y, label = itemlabels, size = .data$sz),
                                        fontface = "bold", color = w_df$col)
                        } else {
                                p <- p + ggplot2::geom_text(
                                        data = w_df,
                                        ggplot2::aes(.data$x, .data$y, label = itemlabels),
                                        size = w_label_size, fontface = "bold", color = w_df$col)
                        }
                }
        }

        if (!is.null(figuretitle)) p <- p + ggplot2::labs(title = figuretitle)

        print(p)
        invisible(p)
}
