#' Rescale a numeric vector to a new range
#'
#' Linearly transforms the values in `x` so they fall within a specified
#' interval. Useful, for example, when mapping latent‐space coordinates to
#' aesthetic ranges (point sizes, color scales, etc.) in a plot.
#'
#' @param x Numeric vector. The data to be rescaled.
#' @param to Numeric vector of length 2 giving the lower and upper limits of
#'   the target range. Defaults to `c(0, 1)`.
#' @param na.rm Logical. Should missing values be ignored when computing the
#'   source range? Defaults to `TRUE`. Any `NA`s in `x` are returned unchanged.
#'
#' @return A numeric vector the same length as `x`, with values rescaled to lie
#'   within `to[1]` and `to[2]`. The function preserves the positions of
#'   `NA`s.
#'
#' @details
#' If all non‑missing values in `x` are identical, the function returns the
#' midpoint of the target range (`mean(to)`) for those elements to avoid
#' division by zero.
#'
#' @examples
#' set.seed(123)
#' x <- rnorm(5)
#'
#' # Default 0–1 range
#' rescale_to_range(x)
#'
#' # Rescale to −1–1
#' rescale_to_range(x, to = c(-1, 1))
#'
#' # Preserve NAs but ignore them when determining the range
#' x_with_na <- c(x, NA, 10)
#' rescale_to_range(x_with_na, to = c(0, 100))
#'
#' @export

rescale_to_range <- function(x, to = c(0, 1), na.rm = TRUE) {
        # use to scale point size to a specified range in a latent space
        rng <- range(x, na.rm = na.rm)
        (x - rng[1]) / diff(rng) * diff(to) + to[1]
}
