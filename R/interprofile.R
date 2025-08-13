#' Draw a Posterior Interaction Profile in either style
#'
#' Convenience wrapper that calls **[`pip_fountain()`]** (default) or
#' **[`pip_waterfall()`]** depending on the `style` argument.  All additional
#' arguments are forwarded unchanged to the selected function, so you can pass
#' `alpha`, `beta`, `distance_mat`, HDI bounds, grouping factors, and so on in
#' exactly the same way as you would when calling the underlying plotting
#' functions directly.
#'
#' @param style Character string choosing the layout.  Accepts `"fountain"`
#'   (default) or `"waterfall"`.  Matching is case‑insensitive and only the
#'   first few letters are required (e.g., `"wat"`).
#' @param ... Further arguments passed on to either `pip_fountain()` or
#'   `pip_waterfall()`.  See those functions for a full description of valid
#'   parameters.
#'
#' @return Whatever the chosen PIP function returns: a `patchwork` object that
#'   combines the two `ggplot2` panels, invisibly returned after being printed.
#'
#' @seealso
#' * [`pip_fountain()`] — “base at −β, arrow up” style
#' * [`pip_waterfall()`] — “base at β, arrow down” style
#'
#' @examples
#' # Small simulated example -----------------------------------------
#' set.seed(42)
#' N <- 6; I <- 10
#' alpha <- rnorm(N)
#' beta  <- rnorm(I, sd = 0.7)
#' dist  <- abs(matrix(rnorm(N * I, sd = 0.8), N, I))  # fake distances
#'
#' # pip_profile() defaults to the fountain view
#' interprofile(alpha = alpha,
#'             beta  = beta,
#'             distance_mat = dist,
#'             focal_id = 2)
#'
#' # Switch to waterfall with the same data
#' interprofile("waterfall",
#'             alpha = alpha,
#'             beta  = beta,
#'             distance_mat = dist,
#'             item_group = rep(LETTERS[1:2], length.out = length(beta)),
#'             y_limits=c(-3,2))
#'
#' @export
interprofile <- function(style = c("fountain", "waterfall"), ...) {
        style <- tolower(style[1])          # allow partial/upper‑case matches
        if (startsWith(style, "water")) {
                pip_waterfall(...)
        } else {
                pip_fountain(...)
        }
}
