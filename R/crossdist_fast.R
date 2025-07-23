#' Fast pairwise (cross) Euclidean distances
#'
#' Computes the Euclidean distance between every row of a “person” matrix
#' (`z`, size *N* × *d*) and every row of an “item” matrix (`w`, size *I* × *d*).
#' An optional `item_labels` argument lets you collapse items into groups first,
#' replacing each group with its centroid before distances are calculated—
#' handy when you only need person–group distances.
#'
#' @param z Numeric matrix of shape *N* × *d*. Each row is a point in *d*‑dimensional
#'   latent space representing a person.
#' @param w Numeric matrix of shape *I × d*. Each row is a point in *d*‑dimensional
#'   latent space representing an item.
#' @param item_labels Optional character or factor vector of length *I* giving a
#'   group label for each item row in `w`.
#'   * If supplied, the function first replaces the items in each group with
#'   their centroid (mean position) so the output contains distances to those
#'   centroids rather than to individual items.
#'   * If `NULL` (default), distances are computed to every item.
#'
#' @return A numeric distance matrix.
#' * When `item_labels` is `NULL`, the result is *N × I*: distance from every
#'   person to every item.
#' * When `item_labels` is provided, the result is *N × G*, where *G* is the
#'   number of distinct groups.
#'
#' @details
#' The computation exploits the identity
#' \deqn{\|z_j-w_i\|^2 = \|z_j\|^2 + \|w_i\|^2 - 2 z_j^\top w_i,}{%
#' ‖z_j−w_i‖² = ‖z_j‖² + ‖w_i‖² − 2 z_jᵀ w_i}
#' allowing all pairwise squared distances to be obtained with a single
#' matrix multiplication. Negative rounding errors are clipped at zero before
#' taking the square root.
#'
#' @examples
#' set.seed(42)
#' z <- matrix(rnorm(15), nrow = 5)   # 5 persons in 3‑D
#' w <- matrix(rnorm(30), nrow = 10)  # 10 items   in 3‑D
#'
#' # Person–item distances
#' d_full <- crossdist_fast(z, w)
#'
#' # Person–group distances (items grouped into two sets)
#' grp <- rep(c("A", "B"), each = 5)
#' d_group <- crossdist_fast(z, w, item_labels = grp)
#'
#' @export
crossdist_fast <- function(z, w, item_labels = NULL) {
        # z: N x d matrix (person positions)
        # w: I x d matrix (item positions)
        # item_labels: optional vector of length I specifying item group labels

        if (!is.null(item_labels)) {
                # Ensure item_labels length matches number of items
                stopifnot(length(item_labels) == nrow(w))

                # Compute group centroids by label
                label_levels <- unique(item_labels)
                group_centroids <- do.call(rbind, lapply(label_levels, function(lbl) {
                        colMeans(w[item_labels == lbl, , drop = FALSE])
                }))
                rownames(group_centroids) <- label_levels

                # Use centroids as new "w"
                w <- group_centroids
        }

        # Precompute squared norms
        z_sq <- rowSums(z^2)        # length N
        w_sq <- rowSums(w^2)        # length I or number of groups

        # Compute distance matrix
        dist_sq <- outer(z_sq, w_sq, "+") - 2 * (z %*% t(w))
        dist_mat <- sqrt(pmax(dist_sq, 0))  # ensure non-negative

        return(dist_mat)
}
