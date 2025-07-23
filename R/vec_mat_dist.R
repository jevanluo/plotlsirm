#' Euclidean distance from a single vector to each row of a matrix
#'
#' Calculates the Euclidean distance between a reference vector
#' `v` and every row of a matrix `mat`.  This is a thin wrapper around
#' `rowSums()` and avoids an explicit loop, so it is fast even for large
#' matrices.
#'
#' @param v Numeric vector of length *d*. The reference point in *d*‑dimensional
#'   space.
#' @param mat Numeric matrix with *n* rows and *d* columns. Each row is treated
#'   as a point whose distance from `v` is to be computed.
#'   The number of columns in `mat` must match `length(v)`.
#'
#' @return A numeric vector of length *n* containing the Euclidean distance
#'   between `v` and each corresponding row of `mat`.
#'
#' @details
#' Internally the function replicates `v` into an *n × d* matrix, subtracts it
#' from `mat`, squares the element‑wise differences, sums across columns, and
#' finally takes the square root—i.e.
#' \deqn{d_i = \sqrt{\sum_{k=1}^d (m_{ik} - v_k)^2}}
#' for each row *i*.
#' Because the computation is fully vectorised it is considerably faster than a
#' simple `apply()` or a for‑loop implementation.
#'
#' @examples
#' # Two‑dimensional example
#' v   <- c(0, 0)
#' mat <- matrix(c(1, 0,
#'                 0, 2,
#'                 3, 4),
#'               ncol = 2, byrow = TRUE)
#'
#' vec_mat_dist(v, mat)
#' #> [1] 1 2 5
#'
#' @export
vec_mat_dist <- function(v, mat) {
        sqrt(rowSums((mat - matrix(v, nrow(mat), ncol(mat), byrow = TRUE))^2))
}
