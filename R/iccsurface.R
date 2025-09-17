#' Latent‑Space Item Characteristic *Surface*
#'
#' Evaluates the LSIRM probability
#' \deqn{P(Y_{pi}=1)=\operatorname{logit}^{-1}\!\bigl(\alpha+\beta-\gamma\,d\bigr)}
#' on a rectangular grid of *ability* (\eqn{\alpha}) and *person–item
#' distance* (\eqn{d}) values and, by default, renders the resulting
#' surface interactively with **plotly**.
#'
#' ## Colour options
#' * **Uniform** – a single‐colour surface (`colour_mode = "uniform"`,
#'   default).  The colour is set by `surface_col`.
#' * **Gradient** – a continuous *plotly* palette (`colour_mode = "gradient"`);
#'   the palette is chosen via `palette`, and `dark_high = TRUE` reverses the
#'   scale so *higher* probabilities appear darker.
#'
#' ## Wire‑frame
#' Setting `show_grid = TRUE` overlays a black wire‑frame every
#' `grid_step` rows/columns to emphasise the surface curvature.
#'
#' The scene’s aspect is fixed to a cube, and bold zero‑lines are drawn on
#' the \eqn{\alpha}‑ and *distance* axes so their origins align visually.
#'
#' @param beta            Numeric scalar \eqn{\beta_i} shifting the surface
#'   along the probability axis (item “easiness”).
#' @param alpha_lim,n_alpha Numeric. Range `c(min, max)` and grid size for the
#'   ability axis.  Default is \eqn{\alpha\in[-4,4]} with 60 points.
#' @param dist_lim,n_dist  Numeric. Range `c(min, max)` and grid size for the
#'   distance axis.  Default is \eqn{d\in[0,4]} with 60 points.
#' @param gamma           Positive scalar controlling how strongly the
#'   probability decays with distance.
#' @param colour_mode     `"uniform"` (default) or `"gradient"`.
#' @param surface_col     Single colour used when `colour_mode = "uniform"`.
#' @param palette         Character name of a *plotly* continuous palette
#'   (e.g. `"Viridis"`, `"Hot"`, `"Blues"`); only used when
#'   `colour_mode = "gradient"`.
#' @param dark_high       Logical.  If `TRUE` (default) reverses the palette so
#'   high probabilities map to darker shades.
#' @param surface_opacity Numeric in (0, 1].  By default the function chooses
#'   `1` for uniform and `0.9` for gradient surfaces so the wire‑frame remains
#'   visible.
#' @param show_grid       Logical.  Overlay a wire‑frame?  Default `TRUE`.
#' @param grid_step       Positive integer: draw every `grid_step`‑th row/column
#'   when `show_grid = TRUE`.
#' @param plot            Logical.  If `TRUE` (default) return an interactive
#'   **plotly** surface; if `FALSE` return the raw numeric grid.
#'
#' @return
#' * **`plot = TRUE`** – a *plotly* **htmlwidget** (prints automatically).
#' * **`plot = FALSE`** – a list with components `alpha`, `distance`, and
#'   `prob` (an `n_alpha × n_dist` matrix of probabilities).
#'
#' @section Dependencies:
#' Rendering the surface requires the **plotly** package
#' (`install.packages("plotly")`).  No external packages are needed when
#' `plot = FALSE`.
#'
#' @examples
#' ## Numeric output only
#' surf <- iccsurface(beta = 0, n_alpha = 21, n_dist = 21, plot = FALSE)
#' str(surf)
#'
#' \dontrun{
#' ## Interactive surfaces
#' ## 1. Uniform single‑colour
#' iccsurface(beta = -0.5)
#'
#' ## 2. Gradient "Hot" palette, darker = high P
#' iccsurface(beta = 0.3,
#'            colour_mode = "gradient",
#'            palette     = "Hot")
#' }
#'
#' @export
iccsurface <- function(beta,
                       alpha_lim  = c(-4, 4),   n_alpha = 60,
                       dist_lim   = c( 0, 4),   n_dist  = 60,
                       gamma      = 1,
                       colour_mode = c("uniform", "gradient"),  # default = "uniform"
                       surface_col = "steelblue",
                       palette     = "Viridis",  # any Plotly palette name
                       dark_high   = TRUE,       # reverse palette so high P = darker
                       surface_opacity = NULL,   # auto‑set later
                       show_grid   = TRUE,
                       grid_step   = 5,
                       plot        = TRUE) {

        `%||%` <- function(x, y) if (is.null(x)) y else x   # helper

        ## ---------------- sanity checks & grid ------------------------------
        stopifnot(length(beta) == 1,
                  all(diff(alpha_lim) > 0),
                  all(diff(dist_lim)  > 0),
                  grid_step >= 1)

        logistic  <- function(x) 1 / (1 + exp(-x))
        alpha_seq <- seq(alpha_lim[1], alpha_lim[2], length.out = n_alpha)
        dist_seq  <- seq(dist_lim [1], dist_lim [2], length.out = n_dist)

        prob_mat <- outer(alpha_seq, dist_seq,
                          function(th, d) logistic(th + beta - gamma * d))

        if (!plot) {
                return(list(alpha    = alpha_seq,
                            distance = dist_seq,
                            prob     = prob_mat))
        }

        if (!requireNamespace("plotly", quietly = TRUE)) {
                stop("Package 'plotly' is required for plotting; please install it ",
                     "or set plot = FALSE.", call. = FALSE)
        }

        ## ---------------- colour handling -----------------------------------
        colour_mode   <- match.arg(colour_mode)
        palette_flag  <- if (colour_mode == "gradient") palette else NULL
        reverse_flag  <- if (colour_mode == "gradient" && dark_high) TRUE else FALSE
        showscale     <- colour_mode == "gradient"
        surface_opacity <- surface_opacity %||%
                if (colour_mode == "gradient") 0.9 else 1

        ## ---------------- base surface --------------------------------------
        p <- plotly::plot_ly(
                x            = dist_seq,
                y            = alpha_seq,
                z            = prob_mat,
                type         = "surface",
                colorscale   = if (colour_mode == "uniform")
                        list(list(0, surface_col), list(1, surface_col))
                else
                        palette_flag,          # palette name string
                reversescale = reverse_flag,
                showscale    = showscale,
                opacity      = surface_opacity
        )

        ## ---------------- optional wire‑frame -------------------------------
        if (show_grid) {
                row_idx <- seq(1, n_alpha, by = grid_step)
                col_idx <- seq(1, n_dist,  by = grid_step)

                # row lines (θ fixed, distance varies)
                for (r in row_idx) {
                        p <- plotly::add_trace(
                                p,
                                x    = dist_seq,
                                y    = rep(alpha_seq[r], n_dist),
                                z    = prob_mat[r, ],
                                type = "scatter3d",
                                mode = "lines",
                                line = list(color = "black", width = 1),
                                inherit = FALSE, showlegend = FALSE
                        )
                }
                # column lines (distance fixed, θ varies)
                for (c in col_idx) {
                        p <- plotly::add_trace(
                                p,
                                x    = rep(dist_seq[c], n_alpha),
                                y    = alpha_seq,
                                z    = prob_mat[, c],
                                type = "scatter3d",
                                mode = "lines",
                                line = list(color = "black", width = 1),
                                inherit = FALSE, showlegend = FALSE
                        )
                }
        }

        ## ---------------- axes, aspect, zero‑lines --------------------------
        p <- plotly::layout(
                p,
                scene = list(
                        aspectmode = "cube",
                        xaxis = list(title = "d (distance)",
                                     zeroline = TRUE, zerolinecolor = "black", tick0 = 0),
                        yaxis = list(title = "alpha (ability)",
                                     zeroline = TRUE, zerolinecolor = "black", tick0 = 0),
                        zaxis = list(title = "Pr(Y = 1)", zeroline = FALSE)
                )
        )

        p
}
