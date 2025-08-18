#' Latent-Space Item Characteristic Curve (ICC)
#'
#' Plots the LSIRM ICC for **one item** on a user-defined grid of ability
#' values (`alpha_grid`). The function supports two ways of supplying inputs:
#'
#' * **Point-estimate inputs (default)** — leave `posterior = NULL` and
#'   supply deterministic `beta`, `gamma`, `w_pos`, and (if needed) `z_pos`.
#'   A single curve per requested group is drawn (no ribbon).
#'
#' * **Draws-based inputs** — supply a `posterior` list with draws of
#'   `beta`, `gamma`, `w`, and optionally `z`. The plotted curve is the
#'   **posterior-predictive mean probability** at each θ (i.e., average over
#'   draws). Optionally add a credible ribbon via `credibleRibbon = TRUE`
#'   with width `cred_level`.
#'
#' The probability model is
#'
#' \deqn{P(Y_{ij}=1 \mid \theta_j,d_{ij})
#'       = \operatorname{logit}^{-1}\!\bigl(\theta_j + \beta_i - \gamma\,d_{ij}\bigr),}
#'
#' where \eqn{d_{ij}=\lVert z_j-w_i\rVert}. Choice of the reference position
#' (`reference = "item"`, `"origin"`, or `"person-global"`) determines how
#' \eqn{d_{ij}} is computed for the *baseline* (grey) curve.
#'
#' @section Curve types:
#' * **Reference curve** – distance is computed from the chosen reference
#'   position to the item for every posterior draw (or once with point-estimate
#'   inputs). Shown unless `compare = FALSE`.
#' * **Person curve(s)** – distance is computed from the latent position(s) of
#'   respondent(s) listed in `person_id`. Requires `z` (posterior draws or
#'   point estimates).
#'
#' @param item_id   Scalar index of the item to plot.
#' @param posterior Optional list of draws with components
#'   \describe{
#'     \item{`beta`}{\code{M × I} matrix of item intercepts.}
#'     \item{`gamma`}{Length-\code{M} vector of distance weights.}
#'     \item{`w`}{Either an \code{M × I × D} array or a length-\code{M} list
#'       of \code{I × D} matrices of item coordinates.}
#'     \item{`z`}{Optional array/list of person coordinates (same format as
#'       \code{w}). Only needed if you request \code{person_id} or
#'       \code{reference = "person-global"}.}
#'   }
#' @param beta,gamma Numeric point estimates used **only** when
#'   \code{posterior = NULL}.  \code{beta} may be scalar or length-\code{I}.
#' @param w_pos,z_pos Matrices of point estimates for item (\code{I × D}) and
#'   person (\code{N × D}) coordinates, used only when \code{posterior = NULL}.
#' @param alpha_grid Numeric vector of ability values (default
#'   \code{seq(-4, 4, length.out = 201)}).
#' @param person_id  \code{NULL} (no person curves) or integer vector of
#'   respondent indices to overlay.
#' @param compare    Logical. If \code{TRUE} (default) the reference curve is
#'   drawn in addition to any person curves; if \code{FALSE} only person curves
#'   appear.
#' @param credibleRibbon Logical. Draw the credible ribbon for draws-based
#'   inputs? Ignored (forced \code{FALSE}) when \code{posterior = NULL}.
#'   Default \code{FALSE} to keep plots uncluttered.
#' @param cred_level Width of the credible ribbon (e.g., \code{0.95}).
#' @param reference  One of \code{"item"}, \code{"origin"}, or
#'   \code{"person-global"}; see Details.
#' @param ref_col    Colour for the reference curve.
#' @param person_cols Optional vector of colours for person curves; recycled or
#'   auto-generated as needed.
#'
#' @return (Invisibly) a \pkg{ggplot2} object; the plot is displayed as a
#'   side-effect.
#'
#' @examples
#' ## ---- reproducible demonstration ------------------------------------
#' set.seed(1)
#' I <- 6; N <- 40; D <- 2; M <- 300           # toy dimensions
#'
#' ## 1. Point-estimate inputs (default) ---------------------------------
#' beta_hat  <- 0.3
#' gamma_hat <- 1.2
#' w_hat     <- matrix(rnorm(I * D), I, D)
#' z_hat     <- matrix(rnorm(N * D), N, D)
#'
#' # population curve + one person (no ribbon in point-estimate usage)
#' lsirmicc(item_id  = 4,
#'          beta     = beta_hat,
#'          gamma    = gamma_hat,
#'          w_pos    = w_hat,
#'          z_pos    = z_hat,
#'          person_id = 7)
#'
#' ## 2. Draws-based inputs (posterior list) ------------------------------
#' w_base <- matrix(0, I, D); w_base[, 1] <- seq(-1.2, 1.2, length.out = I)
#' z_base <- matrix(0, N, D); z_base[, 1] <- rep(c(-0.6, 0.6), length.out = N)
#'
#' posterior <- list(
#'   beta  = matrix(rnorm(M * I, 0, 0.25), M, I),
#'   gamma = rgamma(M, shape = 300, rate = 300),
#'   w     = array(rep(w_base, each = M), c(M, I, D)) +
#'           array(rnorm(M * I * D, sd = 0.12), c(M, I, D)),
#'   z     = array(rep(z_base, each = M), c(M, N, D)) +
#'           array(rnorm(M * N * D, sd = 0.12), c(M, N, D))
#' )
#'
#' # posterior-predictive mean curve with ribbon and two people
#' lsirmicc(item_id   = 2,
#'          posterior = posterior,
#'          person_id = c(22, 31),
#'          credibleRibbon = TRUE,
#'          cred_level = 0.95,
#'          person_cols = c("red", "blue"))
#'
#' @export
lsirmicc <- function(item_id,
                     posterior   = NULL,
                     beta        = NULL,
                     gamma       = NULL,
                     w_pos       = NULL,
                     z_pos       = NULL,
                     alpha_grid  = seq(-4, 4, length.out = 201),
                     person_id   = NULL,
                     compare     = TRUE,
                     credibleRibbon = FALSE,
                     cred_level  = 0.95,
                     reference   = c("item", "person-global", "origin"),
                     ref_col     = "grey40",
                     person_cols = NULL)
{
        # ------------------------------------------------------------------------
        det_mode <- is.null(posterior)     # TRUE → using point-estimate inputs
        # ------------------------------------------------------------------------

        if (!det_mode) {
                ## ── input shape detection & helpers for draws-based inputs ───────
                w_is_array <- is.array(posterior$w)
                z_is_array <- !is.null(posterior$z) && is.array(posterior$z)

                M <- if (w_is_array) dim(posterior$w)[1] else length(posterior$w)
                D <- if (w_is_array) dim(posterior$w)[3] else ncol(posterior$w[[1]])

                get_wm <- function(m) {
                        if (w_is_array) posterior$w[m, item_id, ] else posterior$w[[m]][item_id, ]
                }
                get_bm <- function(m) posterior$beta[m, item_id]
                get_gm <- function(m) posterior$gamma[m]

                # ── person-specific helpers (draws-based) ─────────────────────────
                if (!is.null(person_id)) {
                        if (is.null(posterior$z))
                                stop("person_id supplied but posterior$z is NULL.")

                        N <- if (z_is_array) dim(posterior$z)[2] else nrow(posterior$z[[1]])
                        if (any(person_id < 1 | person_id > N))
                                stop("person_id outside 1...N.")
                }

                get_zm <- function(m, pid) {
                        if (z_is_array) posterior$z[m, pid, ] else posterior$z[[m]][pid, ]
                }
        } else {
                ## ── point-estimate branch (deterministic inputs) ─────────────────
                if (any(sapply(list(beta, gamma, w_pos), is.null)))
                        stop("Provide beta, gamma, and w_pos when posterior = NULL.")

                if (length(gamma) != 1)
                        stop("gamma must be a single numeric value.")

                I  <- if (length(beta) == 1) nrow(w_pos) else length(beta)
                D  <- ncol(w_pos)
                M  <- 1                                   # single “draw”

                get_wm <- function(m) w_pos[item_id, ]
                get_bm <- function(m) if (length(beta) == 1) beta else beta[item_id]
                get_gm <- function(m) gamma

                z_is_array <- FALSE                       # for later helpers
                if (!is.null(person_id) || reference == "person-global") {
                        if (is.null(z_pos))
                                stop("z_pos must be provided for person-global reference or person_id curves.")
                        get_zm <- function(m, pid) z_pos[pid, ]
                }
        }
        logistic <- function(x) 1 / (1 + exp(-x))
        summarise_probs <- function(mat, lab) {
                alpha <- (1 - cred_level) / 2

                # helper that degrades gracefully to the mean if length == 1
                q_safe <- function(x, prob) {
                        if (length(x) <= 1) mean(x) else stats::quantile(x, prob, names = FALSE)
                }

                as.data.frame(mat) |>
                        stats::setNames(seq_along(alpha_grid)) |>
                        dplyr::mutate(draw = dplyr::row_number()) |>
                        tidyr::pivot_longer(
                                - .data$draw,
                                names_to        = "t_idx",
                                values_to       = "p",
                                names_transform = list(t_idx = as.integer)
                        ) |>
                        dplyr::mutate(
                                theta = alpha_grid[.data$t_idx],
                                group = lab
                        ) |>
                        dplyr::group_by(.data$theta, .data$group) |>
                        dplyr::summarise(
                                med = mean(.data$p, na.rm = TRUE),
                                lo  = q_safe(.data$p, alpha),
                                hi  = q_safe(.data$p, 1 - alpha),
                                .groups = "drop"
                        )
        }

        # ---------------------------------------------------------------- reference choice
        reference <- match.arg(reference)

        ## helper: distance from the chosen reference to the item in draw m ---------------
        dist_ref <- switch(reference,

                           "item" = function(m) 0,                      # anchor ICC at item position

                           "origin" = {
                                   z_ref <- rep(0, D)
                                   function(m) sqrt(sum((z_ref - get_wm(m))^2))
                           },

                           "person-global" = {                           # reference to grand mean of persons
                                   if (det_mode && is.null(z_pos))
                                           stop("reference = 'person-global' requires z_pos")
                                   if (!det_mode && is.null(posterior$z))
                                           stop("reference = 'person-global' requires posterior$z")

                                   z_ref <- if (det_mode) {
                                           colMeans(z_pos)                           # grand mean persons
                                   } else if (z_is_array) {
                                           apply(posterior$z, 3, mean)
                                   } else {
                                           Reduce(`+`, lapply(posterior$z, colMeans)) / M
                                   }
                                   function(m) sqrt(sum((z_ref - get_wm(m))^2))
                           }
        )

        ## ------------------------------------------------------ probability fillers
        fill_probs <- function(dist_fun) {
                out <- matrix(NA_real_, M, length(alpha_grid))
                for (m in seq_len(M)) {
                        d_m     <- dist_fun(m)
                        shift_m <- get_bm(m) - get_gm(m) * d_m
                        out[m, ] <- logistic(outer(alpha_grid, shift_m, "+"))
                }
                out
        }

        ## reference probs ---------------------------------------------------------
        probs_ref <- fill_probs(dist_ref)
        df_ref    <- summarise_probs(probs_ref,
                                     paste0("Reference (", reference, ")"))

        ## person probs ------------------------------------------------------------
        person_dfs <- list()
        if (!is.null(person_id)) {
                if (is.null(person_cols))
                        person_cols <- scales::hue_pal()(length(person_id))
                person_cols <- rep(person_cols, length.out = length(person_id))

                person_dfs <- lapply(seq_along(person_id), function(k) {
                        pid <- person_id[k]
                        probs_k <- fill_probs(function(m)
                                sqrt(sum((get_zm(m, pid) - get_wm(m))^2)))
                        summarise_probs(probs_k, paste0("Person ", pid))
                })
                names(person_dfs) <- as.character(person_id)
        }

        ## --------------------------------------------------------------- ggplot ---
        p <- ggplot2::ggplot()

        # person layers ----------------------------------------------------------
        if (!is.null(person_id)) {
                for (k in seq_along(person_id)) {
                        dfk <- person_dfs[[k]]
                        if (credibleRibbon && M > 1) {
                                p <- p + ggplot2::geom_ribbon(
                                        data = dfk,
                                        ggplot2::aes(x = .data$theta, ymin = .data$lo, ymax = .data$hi, fill = .data$group),
                                        alpha = 0.18
                                )
                        }
                        p <- p + ggplot2::geom_line(
                                data = dfk,
                                ggplot2::aes(x = .data$theta, y = .data$med, colour = .data$group),
                                linewidth = 1
                        )
                }
        }

        # reference layer --------------------------------------------------------
        if (is.null(person_id) || compare) {
                if (credibleRibbon && M > 1) {
                        p <- p + ggplot2::geom_ribbon(
                                data = df_ref,
                                ggplot2::aes(x = .data$theta, ymin = .data$lo, ymax = .data$hi, fill = .data$group),
                                alpha = 0.25
                        )
                }
                p <- p + ggplot2::geom_line(
                        data = df_ref,
                        ggplot2::aes(x = .data$theta, y = .data$med, colour = .data$group),
                        linewidth = 0.9
                )
        }

        # ───────────────────────── scale colours / fills ─────────────────────────
        ref_lab <- paste0("Reference (", reference, ")")

        # 1) always list the person curves; add the reference only if it will be plotted
        legend_breaks <- c(paste0("Person ", person_id),
                           if (compare || is.null(person_id)) ref_lab)

        # 2) make (or recycle) enough distinct colours for the people
        if (length(person_id)) {
                if (is.null(person_cols))
                        person_cols <- scales::hue_pal()(length(person_id))
                person_cols <- rep(person_cols, length.out = length(person_id))
        }

        # 3) build the colour map in EXACT legend order
        col_map <- stats::setNames(c(person_cols, if (ref_lab %in% legend_breaks) ref_col),
                                   legend_breaks)

        # enforce factor levels for every plotted data frame
        df_ref$group <- factor(df_ref$group, levels = legend_breaks)
        for (k in seq_along(person_dfs)) {
                person_dfs[[k]]$group <- factor(person_dfs[[k]]$group,
                                                levels = legend_breaks)
        }

        # manual scales with limits guarantee legend order
        p <- p +
                ggplot2::scale_colour_manual(values = col_map,
                                             limits = legend_breaks, name = NULL) +
                ggplot2::scale_fill_manual(values   = col_map,
                                           limits = legend_breaks, name = NULL) +
                ggplot2::scale_y_continuous(limits = c(0, 1)) +
                ggplot2::labs(
                        title = paste("LSIRM ICC - Item", item_id),
                        x     = "\u03b1",
                        y     = "Pr(Response=1)"
                ) +
                ggplot2::theme_minimal(base_size = 12) +
                ggplot2::theme(legend.position = c(0.02, 0.98),
                               legend.justification = c(0, 1))

        print(p)
        invisible(p)
}
