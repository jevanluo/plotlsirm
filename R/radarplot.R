#' Radar plot of branch‐specific abilities (single or multiple subjects)
#'
#' Draws a radar / spider chart in which each axis (“branch”) represents a
#' domain‑specific ability and the radial extent marks the attained score.
#' *Single‑subject mode* shades the overall‑ability circle and fills the polygon
#' formed by the branch scores.
#' *Multi‑subject mode* overlays several polygons on a common background,
#' optionally coloring, labelling, and annotating each subject.
#'
#' @param data Numeric **vector** (one subject) or **matrix/data‑frame**
#'   (*m* × *n*) of scores, where *n* =`length(labels)` is the number of
#'   branches.  Rows correspond to subjects when `data` is two‑dimensional.
#' @param labels Character vector of length *n* giving the label for each
#'   branch / axis.
#' @param max_radius Numeric.  Maximum drawing radius in plot units after the
#'   scores in `ability_range` have been linearly rescaled.  Default `100`.
#' @param branch_max Numeric.  Reserved for future branch‑specific scaling.
#'   Currently ignored.
#' @param overallAbility Numeric scalar *or* length‑*m* vector giving the
#'   overall ability for each subject.  Used to set the radius of the shaded
#'   background circle for that subject.  If `NA` the maximum of
#'   `ability_range` is used.
#' @param ability_range Numeric length‑2 vector `[min, max]` defining the scale
#'   of the input scores.  These limits are mapped to the interval
#'   `[0, max_radius]`.
#' @param abilityCutoffs Numeric.  Reserved for future color gradations; not
#'   used in the current version.
#' @param bgColors Character vector of colors for the shaded background
#'   circles.  Only the first element is used at present.
#' @param markerInd Numeric vector (0 = hollow, 1 = solid) of length *n*
#'   indicating the point style for each branch.  When `data` is a matrix this
#'   can be supplied as an *m* × *n* matrix so each subject has its own marker
#'   pattern.
#' @param point_cex Numeric point size for branch markers.  Default `4`.
#' @param subjectLabels Optional character vector of length *m* naming each
#'   subject in multi‑subject mode.  Row names of `data` are used when
#'   available; otherwise `"Subject 1"`, `"Subject 2"`, … are generated.
#' @param sampleColors Character vector of length *m* giving the polygon/point
#'   r for each subject.  Defaults to a distinct hue palette if `NULL`.
#' @param showOverallAbility Logical.  If `TRUE`, prints each subject’s overall
#'   ability beneath their label.
#' @param title Optional plot title.
#' @param plot_margin A `ggplot2::margin()` object controlling the outer
#'   whitespace around the figure.  Default adds 20 pt on every side.
#' @param label_angle_offset Numeric scalar or length‑*n* vector specifying (in
#'   degrees) how much to rotate each branch label relative to its default
#'   tangential orientation.  Useful for fine‑tuning readability.
#'
#' @return A `ggplot` object representing the radar chart (also printed as a
#'   side effect).
#'
#' @import ggplot2
#' @importFrom scales hue_pal
#'
#' @examples
#' ## Single subject -------------------------------------------------
#' #### the distance from a person to all items/item clusters
#' dist_z_w <- c(item1 = 1.6, item2 = 0.8, item3 = 1.9, item4 = 2.5, item5 = 0.4)
#' #### transform distance to strength
#' strength <- exp(-dist_z_w)
#' #### plot the radar
#' radarplot(strength, labels = names(strength),
#'           overallAbility = 1.8, showOverallAbility = TRUE,
#'           title = "Student A profile")
#'
#' ## Multiple subjects ---------------------------------------------
#' set.seed(1)
#' #### strength for 3 persons on 5 items
#' dat <- matrix(rnorm(15,3,1), nrow = 3,
#'               dimnames = list(NULL, c("item1", "item2", "item3", "item4", "item5")))
#' radarplot(dat, labels = colnames(dat),
#'           overallAbility = c(-1.8,0.5,2.5),
#'           subjectLabels  = c("Alice", "Bob", "Cara"),
#'           sampleColors   = c("#1b9e77", "#d95f02", "#7570b3"),
#'           showOverallAbility = TRUE,
#'           title = "Class-level comparison")
#'
#' @export

radarplot <- function(data = NULL,         # Numeric vector (one subject) or matrix/data frame (multiple subjects)
                             labels,              # Branch (axis) labels (character vector)
                             max_radius = 100,    # Maximum drawing radius; the ability_range is rescaled to [0, max_radius]
                             branch_max = 10,     # (Not used in this version; branch scaling is relative to overall ability)
                             overallAbility = NA, # Overall ability value. For multiple subjects supply either a single value or a vector (one per subject)
                             ability_range = c(-3, 3),  # The numeric range that defines the radar background grid
                             abilityCutoffs = c(-1.68, 1.68), # (Not used in this version for coloring)
                             bgColors = c("red", "yellow", "green"),  # (Not used here for multiple subjects)
                             markerInd = NULL,  # Indicator vector (0: hollow markers, 1: solid markers) for branch markers
                             point_cex = 4,       # Marker (point) size
                             subjectLabels = NULL,  # Optional subject labels (for multiple subjects)
                             sampleColors = NULL,   # Optional vector of colors (one per subject) for the radar polygons/points
                             showOverallAbility = FALSE,  # If TRUE, annotate each subject with its overall ability value
                             title = NULL,        # Plot title (optional)
                             plot_margin = margin(t = 20, r = 20, b = 20, l = 20),  # Custom plot margins
                             label_angle_offset = 0) {  # New parameter: additional angle(s) (in degrees) to add to each branch label

        n <- length(labels)  # number of branches

        # Check data: if data is a vector, then single subject; if a matrix/data frame, then multiple subjects.
        if (is.null(data))
                stop("No data provided!")

        if (is.null(dim(data))) {
                # Single subject: ensure length matches number of labels.
                if (length(data) != n)
                        stop("Data vector must have the same length as labels.")
                m <- 1
        } else {
                m <- nrow(data)
                if (ncol(data) != n)
                        stop("For multiple subjects, the number of columns in data must equal the number of labels.")
        }

        if (is.null(markerInd)) {
                markerInd <- if (m == 1L) {
                        rep(1, n)
                } else {
                        matrix(1, nrow = m, ncol = n)
                }
        }

        ### A. Build the common radar background using ability_range
        ability_min <- ability_range[1]
        ability_max_val <- ability_range[2]

        # Helper to rescale a value from ability_range to the drawing radius [0, max_radius].
        rescale <- function(value) {
                (value - ability_min) / (ability_max_val - ability_min) * max_radius
        }

        # Compute 6 grid values (including endpoints)
        grid_vals <- seq(ability_min, ability_max_val, length.out = 6)
        grid_radii <- sapply(grid_vals, rescale)

        circle_data <- data.frame()
        for (r in grid_radii) {
                theta <- seq(0, 2 * base::pi, length.out = 100)
                circle_data <- rbind(circle_data,
                                     data.frame(x = r * cos(theta),
                                                y = r * sin(theta),
                                                group = r))
        }

        # Define the radial axes (all subjects share these)
        angles <- seq(0, 2 * base::pi, length.out = n + 1)[1:n]
        axes_df <- data.frame(axis = labels,
                              angle = angles,
                              x = max(grid_radii) * cos(angles),
                              y = max(grid_radii) * sin(angles),
                              stringsAsFactors = FALSE)

        # Allow label_angle_offset to be a vector (one per branch) or a single value.
        if (length(label_angle_offset) == 1) {
                label_angle_offset <- rep(label_angle_offset, n)
        } else if (length(label_angle_offset) != n) {
                stop("label_angle_offset must be either a single value or a vector with length equal to the number of labels.")
        }

        # Compute the tangent angle (in degrees) for each branch and add the offset.
        # The tangent angle is computed by converting the branch angle (radians) to degrees and adding 90°.
        axes_df$label_angle <- axes_df$angle * 180 / base::pi + 90 + label_angle_offset
        # Adjust angles to keep text orientation sensible.
        axes_df$label_angle <- ifelse(axes_df$label_angle > 180, axes_df$label_angle - 180, axes_df$label_angle)

        # Build the base radar background.
        p <- ggplot() +
                # Background grid circles:
                geom_path(data = circle_data, aes(x = .data$x, y = .data$y, group = .data$group),
                          color = "gray", linetype = "dashed") +
                # Radial axes lines:
                geom_segment(data = axes_df, aes(x = 0, y = 0, xend = .data$x, yend = .data$y),
                             color = "gray", linetype = "dashed") +
                # Axis labels (placed slightly outside the outermost grid) with adjusted angles:
                geom_text(data = axes_df,
                          aes(x = 1.1 * max(grid_radii) * cos(.data$angle),
                              y = 1.1 * max(grid_radii) * sin(.data$angle),
                              label = .data$axis,
                              angle = .data$label_angle),
                          size = 4) +
                coord_fixed() +
                theme_void() +
                theme(plot.margin = plot_margin)

        if (!is.null(title))
                p <- p + ggtitle(title)

        ### B. Single Subject Case
        if (m == 1) {
                overall_val <- if (!is.na(overallAbility)) overallAbility else ability_max_val
                R <- rescale(overall_val)

                theta_bg <- seq(0, 2 * base::pi, length.out = 100)
                bg_polygon <- data.frame(x = R * cos(theta_bg),
                                         y = R * sin(theta_bg))

                col <- if (!is.null(sampleColors) && length(sampleColors) >= 1) sampleColors[1] else "blue"

                max_val <- max(data, na.rm = TRUE)
                scaling_factor <- ifelse(max_val > 0, R / max_val, 1)
                data_scaled <- data * scaling_factor

                data_polygon <- data.frame(axis = labels,
                                           value = data_scaled,
                                           marker = markerInd,
                                           angle = angles,
                                           x = data_scaled * cos(angles),
                                           y = data_scaled * sin(angles),
                                           stringsAsFactors = FALSE)
                data_polygon <- rbind(data_polygon, data_polygon[1, ])

                p <- p +
                        geom_polygon(data = bg_polygon, aes(x = .data$x, y = .data$y),
                                     fill = col, alpha = 0.3) +
                        geom_polygon(data = data_polygon, aes(x = .data$x, y = .data$y),
                                     fill = col, color = col, alpha = 0.7, linewidth = 1) +
                        geom_point(data = data_polygon[1:n, ],
                                   aes(x = .data$x, y = .data$y),
                                   shape = ifelse(markerInd == 0, 1, 16),
                                   size = point_cex, color = col)

                if (showOverallAbility) {
                        overall_df <- data.frame(x = 0,
                                                 y = - R - 5,
                                                 label = paste0("alpha: ", overall_val))
                        p <- p + geom_text(data = overall_df, aes(x = .data$x, y = .data$y, label = .data$label),
                                           color = "black", size = 5, fontface = "italic")
                }

                return(p)
        }

        ### C. Multiple Subjects Case (Overlay All Subjects on One Radar)
        if (is.null(sampleColors)) {
                sampleColors <- scales::hue_pal()(m)
        } else {
                if (length(sampleColors) < m)
                        stop("Not enough sampleColors provided for the number of subjects.")
        }

        if (is.null(subjectLabels)) {
                if (!is.null(rownames(data))) {
                        subjectLabels <- rownames(data)
                } else {
                        subjectLabels <- paste("Subject", 1:m)
                }
        }

        polygon_centroid <- function(x, y) {
                if (abs(x[1] - x[length(x)]) < 1e-6 && abs(y[1] - y[length(y)]) < 1e-6) {
                        x <- x[-length(x)]
                        y <- y[-length(y)]
                }
                npts <- length(x)
                A <- 0
                Cx <- 0
                Cy <- 0
                for (i in 1:npts) {
                        j <- ifelse(i == npts, 1, i + 1)
                        cross <- x[i] * y[j] - x[j] * y[i]
                        A <- A + cross
                        Cx <- Cx + (x[i] + x[j]) * cross
                        Cy <- Cy + (y[i] + y[j]) * cross
                }
                A <- A / 2
                Cx <- Cx / (6 * A)
                Cy <- Cy / (6 * A)
                return(c(Cx, Cy))
        }

        overall_labels_df <- data.frame()

        for (i in 1:m) {
                subj_data <- as.numeric(data[i, ])
                markerInd_i <- as.numeric(markerInd[i, ])
                if (length(overallAbility) == m) {
                        overall_i <- overallAbility[i]
                } else if (length(overallAbility) == 1) {
                        overall_i <- overallAbility
                } else {
                        overall_i <- ability_max_val
                }
                if (is.na(overall_i))
                        overall_i <- ability_max_val
                R_i <- rescale(overall_i)

                theta_bg <- seq(0, 2 * base::pi, length.out = 100)
                subj_bg <- data.frame(x = R_i * cos(theta_bg),
                                      y = R_i * sin(theta_bg))

                max_val <- max(subj_data, na.rm = TRUE)
                scaling_factor <- ifelse(max_val > 0, R_i / max_val, 1)
                subj_scaled <- subj_data * scaling_factor

                subj_poly <- data.frame(axis = labels,
                                        value = subj_scaled,
                                        marker = markerInd_i,
                                        angle = angles,
                                        x = subj_scaled * cos(angles),
                                        y = subj_scaled * sin(angles),
                                        stringsAsFactors = FALSE)
                subj_poly <- rbind(subj_poly, subj_poly[1, ])

                p <- p +
                        geom_polygon(data = subj_bg, aes(x = .data$x, y = .data$y),
                                     fill = sampleColors[i], alpha = 0.2) +
                        geom_polygon(data = subj_poly, aes(x = .data$x, y = .data$y),
                                     fill = NA, color = sampleColors[i], alpha = 0.6, linewidth = 1) +
                        geom_point(data = subj_poly[1:n, ],
                                   aes(x = .data$x, y = .data$y),
                                   shape = ifelse(markerInd_i == 0, 1, 16),
                                   size = point_cex, color = sampleColors[i])

                label_df <- data.frame(x = 4 * i,
                                       y = R_i + 1,
                                       label = subjectLabels[i])
                p <- p + geom_text(data = label_df, aes(x = .data$x, y = .data$y, label = .data$label),
                                   color = sampleColors[i], size = 5, fontface = "bold")

                if (showOverallAbility) {
                        oa_df <- data.frame(x = 4 * i,
                                            y = (R_i - 4),
                                            label = paste0("alpha: ", overall_i))
                        overall_labels_df <- rbind(overall_labels_df, oa_df)
                }
        }

        if (showOverallAbility && nrow(overall_labels_df) > 0) {
                p <- p + geom_text(data = overall_labels_df, aes(x = .data$x, y = .data$y, label = .data$label),
                                   color = "black", size = 4, fontface = "italic")
        }

        return(p)
}
