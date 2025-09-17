#' 2-D latent-space interaction map (persons vs. items)
#'
#' Persons (\eqn{z}{z}) and items (\eqn{w}{w}) in a 2D latent space with flexible styling.
#' Person parameter \eqn{\alpha_p}{alpha_p}, item parameter \eqn{\beta_i}{beta_i}.
#'
#' @details
#' **Coloring & legends**
#' - A color legend appears when colors are **mapped** via groups or gradients.
#'   For fixed colors/shapes on both layers, the function creates a simple
#'   two-entry legend ("Persons", "Items") using constant mappings; when shapes
#'   are off but labels are on, it builds a labels-only legend with colored
#'   swatches. The legend title is controlled by `legend_title`.
#' - When both persons and items use gradient coloring
#'   (`z_shape_color_gradient = TRUE` and `w_shape_color_gradient = TRUE`),
#'   you can either **share one gradient and legend** by setting
#'   `share_gradient_scale = TRUE` (uses `shape_color_gradient_low/high` and
#'   `legend_title`), or show **separate gradients/legends** for z and w by
#'   keeping `share_gradient_scale = FALSE` (default). In separate mode,
#'   persons use the **color** scale with title `legend_title_z`, and items
#'   use the **fill** scale with title `legend_title_w`. Layer-specific palettes
#'   can be supplied via `z_shape_color_gradient_low/high` and
#'   `w_shape_color_gradient_low/high`; if `NULL`, the global
#'   `shape_color_gradient_low/high` are used as fallbacks.
#' - Opacity (ggplot “alpha”) and size legends are hidden by default. Set
#'   `show_size_legend = TRUE` to show a size legend when size mapping is used.
#' - When `person_colors` / `item_colors` are **vectors** and `z_label_color` /
#'   `w_label_color` are `NULL`, label colors follow those per-observation vectors.
#'
#' **Opacity**
#' - Transparency is called *opacity* to avoid confusion with the statistic
#'   `alpha`. When you supply explicit fixed colors (single color or vector),
#'   default opacity is full (`1`) unless you enable `*_shape_opacity_scale`
#'   or set `*_shape_fixed_opacity`.
#'
#' **Sizes**
#' - Shape sizes are fixed by `z_shape_size` / `w_shape_size` unless
#'   `z_shape_size_scale` / `w_shape_size_scale` are `TRUE`, in which case
#'   sizes are scaled by `alpha` / `beta` into `z_shape_size_range` /
#'   `w_shape_size_range`.
#' - Label text sizes are fixed by `z_label_size` / `w_label_size` unless
#'   `z_label_size_scale` / `w_label_size_scale` are `TRUE`, in which case they
#'   are scaled by `alpha` / `beta` into `z_label_size_range` / `w_label_size_range`.
#'
#' @param z,w Numeric matrices with 2 columns: coordinates for persons/items.
#' @param gamma Optional scalar stretch factor applied to both z and w.
#'
#' @param person_group,item_group Optional factor/character for grouping colors.
#' @param person_colors,item_colors Optional explicit color vectors (length N/I).
#'
#' @param alpha,beta Optional vectors used to scale shape/label sizes, opacity,
#'   and (optionally) gradients for persons/items, respectively.
#'
#' @param z_shape_size_scale,w_shape_size_scale Logical: scale **shape** sizes by alpha/beta.
#' @param z_shape_size_range,w_shape_size_range Length-2 numeric: shape size ranges (when scaling).
#' @param z_shape_size,w_shape_size Numeric: fixed shape sizes when size scaling is OFF.
#'
#' @param z_label_size_scale,w_label_size_scale Logical: scale **label** sizes.
#'   If `NULL`, defaults to the corresponding shape-size flag.
#' @param z_label_size_range,w_label_size_range Length-2 numeric: label size ranges.
#'
#' @param z_shape_opacity_scale,w_shape_opacity_scale Logical: scale **shape opacity** by alpha/beta.
#' @param z_shape_opacity_range,w_shape_opacity_range Length-2 numeric: opacity ranges.
#' @param z_shape_fixed_opacity,w_shape_fixed_opacity Optional constant opacity (0..1) for shapes.
#'
#' @param z_shape_color_gradient,w_shape_color_gradient Logical: color shapes by a gradient
#'   (darker ⇒ higher). Overrides groups/colors for that layer. If both are `TRUE`,
#'   the gradients can be shared (`share_gradient_scale = TRUE`) or separated
#'   (`share_gradient_scale = FALSE`, default).
#' @param z_shape_color_values,w_shape_color_values Optional numeric drivers for gradients
#'   (defaults: `alpha` / `beta` respectively).
#' @param shape_color_gradient_low,shape_color_gradient_high Global colors for the gradient
#'   palette. Used directly when `share_gradient_scale = TRUE`, or as **fallbacks**
#'   for layer-specific palettes when `share_gradient_scale = FALSE`.
#' @param z_shape_color_gradient_low,z_shape_color_gradient_high Optional colors for the
#'   **persons (z)** gradient when using separate scales. If `NULL`, fall back to
#'   `shape_color_gradient_low/high`.
#' @param w_shape_color_gradient_low,w_shape_color_gradient_high Optional colors for the
#'   **items (w)** gradient when using separate scales. If `NULL`, fall back to
#'   `shape_color_gradient_low/high`.
#'
#' @param show_ticks Logical: draw axis ticks/labels.
#' @param xlim_range,ylim_range Optional axis limits (symmetric if `NULL`).
#'
#' @param itemlabels,personlabels Optional labels (defaults: "I1..", "P1..").
#' @param figuretitle Optional plot title.
#'
#' @param z_shape,w_shape ggplot2 shape codes for persons/items (see `?ggplot2::geom_point`).
#' @param z_shape_color,w_shape_color Fixed fallback **shape** colors when not mapping.
#' @param z_border_width Stroke width for z shapes (when applicable).
#'
#' @param z_label_size,w_label_size Fixed **label** sizes when not scaling.
#' @param z_label_color,w_label_color Label colors; if `NULL`, labels follow the
#'   layer’s color (group/gradient) or fixed color/vector as appropriate.
#' @param show_z_labels,show_w_labels Logical: draw labels for z/w.
#' @param show_z_shapes,show_w_shapes Logical: draw shapes for z/w.
#'
#' @param legend_title Character or expression: the legend title (when shown).
#' @param share_gradient_scale Logical. If `TRUE`, persons and items share **one**
#'   gradient/legend (uses `shape_color_gradient_low/high` and `legend_title`).
#'   If `FALSE` (default), persons and items use **separate** gradients/legends:
#'   persons mapped to colour (title `legend_title_z`) and items mapped to fill
#'   (title `legend_title_w`).
#' @param legend_title_z,legend_title_w Titles (character or expressions) for the
#'   separate z and w gradient legends, used only when
#'   `share_gradient_scale = FALSE` and both gradient mappings are enabled.
#' @param show_size_legend Logical: show a size legend (default `FALSE`).
#'
#' @return Invisibly returns a `ggplot` object; also prints the plot.
#' @import ggplot2
#' @importFrom rlang .data
#'
#' @examples
#' ### example data
#' set.seed(1)
#' z <- matrix(rnorm(40), 20, 2)  # persons
#' w <- matrix(rnorm(30), 15, 2)  # items
#' alpha <- rnorm(nrow(z))        # person alpha
#' beta  <- rnorm(nrow(w))        # item beta
#'
#' ### 1) minimal, fixed colors & shapes
#' intermap2d(z, w)
#'
#' ### 2) minimal, fixed shapes for persons and labels for items
#' intermap2d(
#'   z, w,
#'   show_w_shapes = FALSE, show_w_labels = TRUE
#' )
#'
#' ### 3) Grouped colors + sized shapes, formal legend title
#' intermap2d(
#'   z, w,
#'   person_group = rep(c("Cohort A", "Cohort B"), length.out = nrow(z)),
#'   item_group   = rep(c("Domain X", "Domain Y", "Domain Z"), length.out = nrow(w)),
#'   alpha = alpha, beta = beta,
#'   z_shape_size_scale = TRUE, z_shape_size_range = c(2, 6),
#'   w_shape_size_scale = TRUE, w_shape_size_range = c(2, 8),
#'   show_z_shapes = TRUE, show_w_shapes = FALSE,
#'   show_w_labels = TRUE, legend_title = "Cohort / Domain"
#' )
#'
#' ### 4) Gradient for persons only (darker = higher α), labels scaled by alpha
#' intermap2d(
#'   z, w,
#'   alpha = alpha,
#'   z_shape_color_gradient = TRUE,     # z by gradient
#'   w_shape_color          = "red",    # w fixed color
#'   z_label_size_scale     = TRUE,     # label size ∝ α
#'   show_w_shapes = FALSE,
#'   show_w_labels = TRUE,
#'   # shape_color_gradient_low = "grey80", shape_color_gradient_high = "navy",
#'   legend_title = expression(alpha[p])
#' )
#'
#' ### 5) Gradient for both alpha (persons) and beta (items), shared legend
#' intermap2d(
#'   z, w,
#'   alpha = alpha, beta = beta,
#'   z_shape_color_gradient = TRUE, w_shape_color_gradient = TRUE,
#'   shape_color_gradient_low = "grey80", shape_color_gradient_high = "navy",
#'   z_shape_size_scale = TRUE, w_shape_size_scale = TRUE,
#'   show_z_shapes = TRUE, show_w_shapes = TRUE,
#'   show_z_labels = FALSE, show_w_labels = FALSE,
#'   legend_title = "Intensity (alpha persons, beta items)"
#' )
#'
#' ### 6) Explicit per-observation colors (vectors) + label-only
#' z_cols <- ifelse(alpha > 0, "#1f77b4", "#AEC7E8")
#' w_cols <- ifelse(beta  > 0, "#d62728", "#FF9896")
#' intermap2d(
#'   z, w,
#'   person_colors = z_cols, item_colors = w_cols,
#'   show_z_shapes = FALSE, show_w_shapes = FALSE,
#'   show_z_labels = TRUE,  show_w_labels = TRUE
#' )
#'
#' ### 7) Opacity scaling + fixed shape sizes
#' intermap2d(
#'   z, w,
#'   alpha = alpha, beta = beta,
#'   z_shape_opacity_scale = TRUE,  z_shape_opacity_range = c(0.2, 1.0),
#'   w_shape_opacity_scale = TRUE,  w_shape_opacity_range = c(0.4, 1.0),
#'   z_shape_color = "black", w_shape_color = "orange3",
#'   z_shape_size = 3, w_shape_size = 3.5,        # fixed sizes (no size scaling)
#'   show_z_shapes = TRUE, show_w_shapes = TRUE
#' )
#'
#' ### 8) Label-only scaling; shape sizes fixed; custom legend title for groups
#' intermap2d(
#'   z, w,
#'   person_group = rep(c("High", "Low"), length.out = nrow(z)),
#'   alpha = alpha, beta = beta,
#'   z_label_size_scale = TRUE, w_label_size_scale = TRUE,
#'   z_label_size_range = c(3, 7), w_label_size_range = c(3, 7),
#'   show_z_labels = TRUE, show_w_labels = FALSE,
#'   show_z_shapes = FALSE, show_w_shapes = TRUE,
#'   legend_title = "Performance Group"
#' )
#'
#' ### 9) Stretch coordinates + axis ticks + symmetric limits
#' intermap2d(
#'   z, w,
#'   gamma = 1.5,
#'   show_ticks = TRUE,
#'   xlim_range = c(-4, 4), ylim_range = c(-4, 4),
#'   person_group = rep(c("Train","Test"), length.out = nrow(z)),
#'   legend_title = "Set Membership"
#' )
#'
#' @export
intermap2d <- function(
                z, w,
                # stretch
                gamma           = NULL,
                # grouping / explicit colors
                person_group    = NULL,
                item_group      = NULL,
                person_colors   = NULL,
                item_colors     = NULL,
                # stats
                alpha           = NULL, beta = NULL,
                # SHAPE size controls
                z_shape_size_scale = FALSE, w_shape_size_scale = FALSE,
                z_shape_size_range = c(2, 6), w_shape_size_range = c(2, 8),
                z_shape_size       = 2.5,     w_shape_size       = 3.0,
                # LABEL size controls
                z_label_size_scale = NULL, w_label_size_scale = NULL,
                z_label_size_range = c(3, 7), w_label_size_range = c(3, 7),
                z_label_color = "navy", w_label_color = "firebrick",
                z_label_size       = 4,       w_label_size       = 4,
                # SHAPE opacity controls
                z_shape_opacity_scale = FALSE, w_shape_opacity_scale = FALSE,
                z_shape_opacity_range = c(0.30, 1), w_shape_opacity_range = c(0.30, 1),
                z_shape_fixed_opacity = NULL,       w_shape_fixed_opacity = NULL,
                # SHAPE color gradient
                z_shape_color_gradient = FALSE, w_shape_color_gradient = FALSE,
                z_shape_color_values   = NULL,  w_shape_color_values   = NULL,
                shape_color_gradient_low  = "grey80",
                shape_color_gradient_high = "navy",
                z_shape_color_gradient_low  = NULL,
                z_shape_color_gradient_high = NULL,
                w_shape_color_gradient_low  = NULL,
                w_shape_color_gradient_high = NULL,
                # axes
                show_ticks      = FALSE,
                xlim_range      = NULL, ylim_range = NULL,
                # labels & title
                itemlabels      = NULL, personlabels = NULL,
                figuretitle     = NULL,
                # SHAPE codes & colors
                z_shape         = 16, w_shape = 17,               # ggplot shape codes
                z_shape_color   = "navy", w_shape_color = "firebrick",
                z_border_width  = 0.5,
                # draw toggles
                show_z_labels   = FALSE, show_w_labels = FALSE,
                show_z_shapes   = TRUE,  show_w_shapes = TRUE,
                # legends
                legend_title    = "legend",
                show_size_legend = FALSE,
                share_gradient_scale = FALSE,
                legend_title_z  = expression(alpha[p]),
                legend_title_w  = expression(beta[i])
) {
        # helpers
        rescale_to_range <- function(x, to = c(0,1)) {
                if (is.null(x)) return(NULL)
                rng <- range(x, na.rm = TRUE, finite = TRUE)
                if (!is.finite(rng[1]) || !is.finite(rng[2]) || diff(rng) == 0) {
                        return(rep(mean(to), length(x)))
                }
                (x - rng[1]) / diff(rng) * diff(to) + to[1]
        }
        first_non_na <- function(x) x[which(!is.na(x))[1]]

        # data
        z_df <- as.data.frame(z); colnames(z_df) <- c("x","y")
        w_df <- as.data.frame(w); colnames(w_df) <- c("x","y")

        # after creating z_df, w_df (and applying gamma)
        if (!is.null(person_colors)) {
                stopifnot(length(person_colors) == nrow(z_df))
                z_df$lbl_col <- person_colors
        }
        if (!is.null(item_colors)) {
                stopifnot(length(item_colors) == nrow(w_df))
                w_df$lbl_col <- item_colors
        }

        # override defaults ONLY when the user did NOT explicitly pass z_label_color/w_label_color
        use_z_label_vector <- (!is.null(person_colors)) && missing(z_label_color)
        use_w_label_vector <- (!is.null(item_colors))   && missing(w_label_color)


        if (!is.null(gamma)) { z_df[,1:2] <- z_df[,1:2] * gamma; w_df[,1:2] <- w_df[,1:2] * gamma }
        if (is.null(personlabels)) personlabels <- paste0("P", seq_len(nrow(z_df)))
        if (is.null(itemlabels))   itemlabels   <- paste0("I", seq_len(nrow(w_df)))

        # sizes: shapes
        if (z_shape_size_scale) { stopifnot(length(alpha) == nrow(z_df)); z_df$sz_shape <- rescale_to_range(alpha, z_shape_size_range) }
        if (w_shape_size_scale) { stopifnot(length(beta)  == nrow(w_df)); w_df$sz_shape <- rescale_to_range(beta,  w_shape_size_range) }

        # sizes: labels
        if (is.null(z_label_size_scale)) z_label_size_scale <- z_shape_size_scale
        if (is.null(w_label_size_scale)) w_label_size_scale <- w_shape_size_scale
        z_df$sz_label <- if (z_label_size_scale) { stopifnot(length(alpha) == nrow(z_df)); rescale_to_range(alpha, z_label_size_range) } else z_label_size
        w_df$sz_label <- if (w_label_size_scale) { stopifnot(length(beta)  == nrow(w_df)); rescale_to_range(beta,  w_label_size_range) } else w_label_size

        # color modes for shapes
        z_map_color <- FALSE; w_map_color <- FALSE
        if (z_shape_color_gradient) {
                z_df$col_var <- if (!is.null(z_shape_color_values)) z_shape_color_values else alpha
                stopifnot(length(z_df$col_var) == nrow(z_df))
                z_color_mode <- "gradient"; z_map_color <- TRUE
        } else if (!is.null(person_colors)) {
                stopifnot(length(person_colors) == nrow(z_df))
                z_df$col_fix <- person_colors; z_color_mode <- "vector"
        } else if (!is.null(person_group)) {
                z_df$grp <- factor(person_group); z_color_mode <- "group"; z_map_color <- TRUE
        } else {
                z_df$col_fix <- z_shape_color; z_color_mode <- "fixed"
        }
        if (w_shape_color_gradient) {
                w_df$col_var <- if (!is.null(w_shape_color_values)) w_shape_color_values else beta
                stopifnot(length(w_df$col_var) == nrow(w_df))
                w_color_mode <- "gradient"; w_map_color <- TRUE
        } else if (!is.null(item_colors)) {
                stopifnot(length(item_colors) == nrow(w_df))
                w_df$col_fix <- item_colors; w_color_mode <- "vector"
        } else if (!is.null(item_group)) {
                w_df$grp <- factor(item_group); w_color_mode <- "group"; w_map_color <- TRUE
        } else {
                w_df$col_fix <- w_shape_color; w_color_mode <- "fixed"
        }

        # shape opacity (no opacity legend)
        if (z_shape_opacity_scale) {
                z_df$op <- rescale_to_range(alpha, z_shape_opacity_range)
        } else if (!is.null(z_shape_fixed_opacity)) {
                z_df$op <- rep(z_shape_fixed_opacity, nrow(z_df))
        } else if (z_color_mode %in% c("fixed","vector")) {
                z_df$op <- 1
        } else {
                z_df$op <- 0.8
        }
        if (w_shape_opacity_scale) {
                w_df$op <- rescale_to_range(beta, w_shape_opacity_range)
        } else if (!is.null(w_shape_fixed_opacity)) {
                w_df$op <- rep(w_shape_fixed_opacity, nrow(w_df))
        } else if (w_color_mode %in% c("fixed","vector")) {
                w_df$op <- 1
        } else {
                w_df$op <- 1.0
        }

        # base plot
        p <- ggplot2::ggplot() + ggplot2::coord_fixed()
        rng <- max(abs(c(z_df$x, z_df$y, w_df$x, w_df$y)), na.rm = TRUE)
        if (is.null(xlim_range)) xlim_range <- c(-rng, rng)
        if (is.null(ylim_range)) ylim_range <- c(-rng, rng)
        p <- p + ggplot2::xlim(xlim_range) + ggplot2::ylim(ylim_range)
        p <- p + ggplot2::theme(
                panel.background = ggplot2::element_rect(fill = "#F5F5F5"),
                panel.grid.major = ggplot2::element_line(color = "white"),
                panel.border     = ggplot2::element_rect(color = "grey60", fill = NA, linewidth = 0.8),
                axis.title       = ggplot2::element_blank(),
                axis.ticks       = if (show_ticks) ggplot2::element_line() else ggplot2::element_blank(),
                axis.text        = if (show_ticks) ggplot2::element_text(size = 12) else ggplot2::element_blank(),
                legend.position  = "right"
        )

        # helper: both layers are single fixed shape colors & shapes (simple legend)
        both_fixed_single <- (
                z_color_mode == "fixed" && w_color_mode == "fixed" &&
                        isTRUE(show_z_shapes) && isTRUE(show_w_shapes) &&
                        is.character(z_shape_color) && length(z_shape_color) == 1 &&
                        is.character(w_shape_color) && length(w_shape_color) == 1
        )

        # draw shapes: persons (z)
        if (show_z_shapes) {
                if (z_map_color) {
                        if (z_shape_size_scale) {
                                p <- p + ggplot2::geom_point(
                                        data = z_df,
                                        ggplot2::aes(.data$x, .data$y,
                                                     color = if (z_color_mode == "group") .data$grp else .data$col_var,
                                                     size  = .data$sz_shape, alpha = .data$op
                                        ),
                                        shape = z_shape, stroke = z_border_width
                                )
                        } else {
                                p <- p + ggplot2::geom_point(
                                        data = z_df,
                                        ggplot2::aes(.data$x, .data$y,
                                                     color = if (z_color_mode == "group") .data$grp else .data$col_var,
                                                     alpha = .data$op
                                        ),
                                        shape = z_shape, size = z_shape_size, stroke = z_border_width
                                )
                        }
                } else if (both_fixed_single) {
                        if (z_shape_size_scale) {
                                p <- p + ggplot2::geom_point(
                                        data = z_df,
                                        ggplot2::aes(.data$x, .data$y, color = "Persons", shape = "Persons",
                                                     size = .data$sz_shape, alpha = .data$op),
                                        stroke = z_border_width
                                )
                        } else {
                                p <- p + ggplot2::geom_point(
                                        data = z_df,
                                        ggplot2::aes(.data$x, .data$y, color = "Persons", shape = "Persons", alpha = .data$op),
                                        size = z_shape_size, stroke = z_border_width
                                )
                        }
                } else {
                        if (z_shape_size_scale) {
                                p <- p + ggplot2::geom_point(
                                        data = z_df,
                                        ggplot2::aes(.data$x, .data$y, size = .data$sz_shape, alpha = .data$op),
                                        shape = z_shape, stroke = z_border_width, color = z_df$col_fix
                                )
                        } else {
                                p <- p + ggplot2::geom_point(
                                        data = z_df,
                                        ggplot2::aes(.data$x, .data$y, alpha = .data$op),
                                        shape = z_shape, size = z_shape_size, stroke = z_border_width, color = z_df$col_fix
                                )
                        }
                }
        }

        # draw shapes: items (w)
        if (show_w_shapes) {
                if (w_map_color) {
                        if (w_shape_size_scale) {
                                # NEW: use separate fill gradient for items when both layers use gradients
                                if (w_color_mode == "gradient" && z_shape_color_gradient && w_shape_color_gradient && !share_gradient_scale) {
                                        p <- p + ggplot2::geom_point(
                                                data = w_df,
                                                ggplot2::aes(.data$x, .data$y,
                                                             fill  = .data$col_var,
                                                             size  = .data$sz_shape,
                                                             alpha = .data$op),
                                                shape = 21, color = "black"
                                        )
                                } else {
                                        p <- p + ggplot2::geom_point(
                                                data = w_df,
                                                ggplot2::aes(.data$x, .data$y,
                                                             color = if (w_color_mode == "group") .data$grp else .data$col_var,
                                                             size  = .data$sz_shape, alpha = .data$op),
                                                shape = w_shape
                                        )
                                }
                        } else {
                                if (w_color_mode == "gradient" && z_shape_color_gradient && w_shape_color_gradient && !share_gradient_scale) {
                                        p <- p + ggplot2::geom_point(
                                                data = w_df,
                                                ggplot2::aes(.data$x, .data$y,
                                                             fill = .data$col_var,
                                                             alpha = .data$op),
                                                shape = 21, size = w_shape_size, color = "black"
                                        )
                                } else {
                                        p <- p + ggplot2::geom_point(
                                                data = w_df,
                                                ggplot2::aes(.data$x, .data$y,
                                                             color = if (w_color_mode == "group") .data$grp else .data$col_var,
                                                             alpha = .data$op),
                                                shape = w_shape, size = w_shape_size
                                        )
                                }
                        }
                } else if (both_fixed_single) {
                        if (w_shape_size_scale) {
                                p <- p + ggplot2::geom_point(
                                        data = w_df,
                                        ggplot2::aes(.data$x, .data$y, color = "Items", shape = "Items",
                                                     size = .data$sz_shape, alpha = .data$op)
                                )
                        } else {
                                p <- p + ggplot2::geom_point(
                                        data = w_df,
                                        ggplot2::aes(.data$x, .data$y, color = "Items", shape = "Items", alpha = .data$op),
                                        size = w_shape_size
                                )
                        }
                } else {
                        if (w_shape_size_scale) {
                                p <- p + ggplot2::geom_point(
                                        data = w_df,
                                        ggplot2::aes(.data$x, .data$y, size = .data$sz_shape, alpha = .data$op),
                                        shape = w_shape, color = w_df$col_fix
                                )
                        } else {
                                p <- p + ggplot2::geom_point(
                                        data = w_df,
                                        ggplot2::aes(.data$x, .data$y, alpha = .data$op),
                                        shape = w_shape, size = w_shape_size, color = w_df$col_fix
                                )
                        }
                }
        }

        # labels: persons (z)
        if (show_z_labels) {
                if (z_map_color) {
                        p <- p + ggplot2::geom_text(
                                data = z_df,
                                ggplot2::aes(.data$x, .data$y, label = personlabels,
                                             color = if (z_color_mode == "group") .data$grp else .data$col_var,
                                             size  = if (z_label_size_scale) .data$sz_label else NULL),
                                fontface = "bold", show.legend = TRUE
                        )
                } else if (use_z_label_vector) {
                        p <- p + ggplot2::geom_text(
                                data = z_df,
                                ggplot2::aes(.data$x, .data$y, label = personlabels,
                                             color = .data$lbl_col,
                                             size  = if (z_label_size_scale) .data$sz_label else NULL),
                                fontface = "bold", show.legend = FALSE
                        )
                } else {
                        p <- p + ggplot2::geom_text(
                                data = z_df,
                                ggplot2::aes(.data$x, .data$y, label = personlabels,
                                             size = if (z_label_size_scale) .data$sz_label else NULL),
                                fontface = "bold",
                                color = z_label_color,
                                show.legend = FALSE
                        )
                }
        }

        # labels: items (w)
        if (show_w_labels) {
                if (w_map_color) {
                        # NEW: when using separate gradients, keep item labels fixed (no second color scale)
                        if (w_color_mode == "gradient" && z_shape_color_gradient && w_shape_color_gradient && !share_gradient_scale) {
                                p <- p + ggplot2::geom_text(
                                        data = w_df,
                                        ggplot2::aes(.data$x, .data$y, label = itemlabels,
                                                     size  = if (w_label_size_scale) .data$sz_label else NULL),
                                        color = w_label_color, fontface = "bold", show.legend = FALSE
                                )
                        } else {
                                p <- p + ggplot2::geom_text(
                                        data = w_df,
                                        ggplot2::aes(.data$x, .data$y, label = itemlabels,
                                                     color = if (w_color_mode == "group") .data$grp else .data$col_var,
                                                     size  = if (w_label_size_scale) .data$sz_label else NULL),
                                        fontface = "bold", show.legend = TRUE
                                )
                        }
                } else if (use_w_label_vector) {
                        p <- p + ggplot2::geom_text(
                                data = w_df,
                                ggplot2::aes(.data$x, .data$y, label = itemlabels,
                                             color = .data$lbl_col,
                                             size  = if (w_label_size_scale) .data$sz_label else NULL),
                                fontface = "bold", show.legend = FALSE
                        )
                } else {
                        p <- p + ggplot2::geom_text(
                                data = w_df,
                                ggplot2::aes(.data$x, .data$y, label = itemlabels,
                                             size = if (w_label_size_scale) .data$sz_label else NULL),
                                fontface = "bold",
                                color = w_label_color,
                                show.legend = FALSE
                        )
                }
        }

        # scales & guides
        p <- p + ggplot2::scale_alpha_identity(guide = "none")
        if (isTRUE(show_size_legend) &&
            (isTRUE(z_shape_size_scale) || isTRUE(w_shape_size_scale) ||
             isTRUE(z_label_size_scale) || isTRUE(w_label_size_scale))) {
                p <- p + ggplot2::scale_size_continuous(name = legend_title)
        } else {
                p <- p + ggplot2::scale_size_identity(guide = "none")
        }

        any_gradient <- z_shape_color_gradient || w_shape_color_gradient
        any_grouping <- (!any_gradient) && (z_map_color || w_map_color)

        if (both_fixed_single) {
                p <- p +
                        ggplot2::scale_color_manual(
                                name   = legend_title,
                                values = c(Persons = z_shape_color, Items = w_shape_color)
                        ) +
                        ggplot2::scale_shape_manual(
                                name   = legend_title,
                                values = c(Persons = z_shape, Items = w_shape)
                        )
        } else if (any_gradient) {
                # NEW: if both layers use gradients and we don't share the scale, split into colour (z) and fill (w)
                if (z_shape_color_gradient && w_shape_color_gradient && !share_gradient_scale) {
                        # choose layer-specific palettes with fallback to global palette
                        z_low <- if (!is.null(z_shape_color_gradient_low)) z_shape_color_gradient_low else shape_color_gradient_low
                        z_high <- if (!is.null(z_shape_color_gradient_high)) z_shape_color_gradient_high else shape_color_gradient_high
                        w_low <- if (!is.null(w_shape_color_gradient_low)) w_shape_color_gradient_low else shape_color_gradient_low
                        w_high <- if (!is.null(w_shape_color_gradient_high)) w_shape_color_gradient_high else shape_color_gradient_high
                        p <- p +
                                ggplot2::scale_color_gradient(
                                        name = legend_title_z,
                                        low  = z_low,
                                        high = z_high
                                ) +
                                ggplot2::scale_fill_gradient(
                                        name = legend_title_w,
                                        low  = w_low,
                                        high = w_high
                                )
                } else {
                        # shared single palette (keeps current behaviour)
                        p <- p + ggplot2::scale_color_gradient(
                                name = legend_title,
                                low  = shape_color_gradient_low,
                                high = shape_color_gradient_high
                        )
                }
        } else if (any_grouping) {
                p <- p + ggplot2::scale_color_discrete(name = legend_title)
        } else {
                # possible labels-only legend
                p <- p + ggplot2::scale_color_identity(guide = "none")
                fill_vals <- c()
                if (show_z_labels && !show_z_shapes) {
                        z_leg_col <- if (!is.null(z_label_color)) z_label_color else if (!is.null(person_colors)) first_non_na(person_colors) else z_shape_color
                        fill_vals["Persons"] <- z_leg_col
                }
                if (show_w_labels && !show_w_shapes) {
                        w_leg_col <- if (!is.null(w_label_color)) w_label_color else if (!is.null(item_colors)) first_non_na(item_colors) else w_shape_color
                        fill_vals["Items"] <- w_leg_col
                }
                if (length(fill_vals)) {
                        p <- p +
                                ggplot2::geom_point(
                                        data = data.frame(x = Inf, y = Inf, who = names(fill_vals)),
                                        ggplot2::aes(.data$x, .data$y, fill = .data$who),
                                        shape = 22, size = 4, alpha = 0, inherit.aes = FALSE, show.legend = TRUE
                                ) +
                                ggplot2::scale_fill_manual(
                                        name   = legend_title,
                                        values = fill_vals
                                ) +
                                ggplot2::guides(
                                        fill   = ggplot2::guide_legend(override.aes = list(alpha = 1, size = 5, shape = 22)),
                                        size   = "none",
                                        alpha  = "none",
                                        colour = "none"
                                )
                }
        }

        if (!is.null(figuretitle)) p <- p + ggplot2::labs(title = figuretitle)
        print(p); invisible(p)
}
