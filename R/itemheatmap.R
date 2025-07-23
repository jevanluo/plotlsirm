#' Heat‑map of item–item similarity in latent space
#'
#' Draws a lower‑triangle heat‑map (including the main diagonal) of the
#' similarity between item positions in a latent space.  Similarity is defined
#' as \deqn{\exp(-\gamma\,d_{ij})} where \eqn{d_{ij}} is the Euclidean distance
#' between items *i* and *j*, and \eqn{\gamma>0} is a scale parameter
#' controlling how quickly similarity decays with distance.  The function can
#' optionally reorder items via hierarchical clustering so that similar items
#' are placed next to one another, making block‑structure easier to see.
#'
#' @param w Numeric matrix or data frame with one row per item and two (or more)
#'   columns giving the latent coordinates of each item.
#' @param gamma Positive numeric scalar.  Controls the steepness of the
#'   similarity decay; larger values make similarity drop off more quickly.
#'   Default is `1`.
#' @param item_names Optional character vector of item labels.  Must have the
#'   same length as `nrow(w)`.  Defaults to `"I1"`, `"I2"`, … if `NULL`.
#' @param reorder Logical.  If `TRUE` (default is `FALSE`) the heat‑map is reordered using
#'   hierarchical clustering of the distance matrix so that similar items are
#'   grouped along the diagonal.
#' @param digits Integer.  Number of decimal places used when printing similarity
#'   values inside the cells.  Default is `2`.
#' @param title Optional character string for the plot title.
#'
#' @return (Invisibly) a `ggplot` object containing the heat‑map. The plot is
#'   also displayed as a side effect.
#'
#' @import ggplot2
#' @importFrom stats dist hclust as.dist
#' @importFrom rlang .data
#'
#' @examples
#' set.seed(123)
#' w <- matrix(rnorm(40), ncol = 2)   # 20 items in 2‑D latent space
#'
#' # Default heat‑map
#' itemheatmap(w)
#'
#' # Stronger decay (gamma = 3) and custom item names
#' itemheatmap(w, gamma = 3, item_names = paste("Item", 1:nrow(w)))
#'
#' # Turn off re‑ordering
#' itemheatmap(w, reorder = FALSE, title = "Fixed item ordering")
#'
#' @export
itemheatmap <- function(
                w,
                gamma       = 1,
                item_names  = NULL,
                reorder     = FALSE,
                digits      = 2,
                title       = NULL
) {
        w_mat <- as.matrix(w)
        I     <- nrow(w_mat)
        if (is.null(item_names)) item_names <- paste0("I", seq_len(I))

        # similarity matrix ------------------------------------------------
        dist_mat <- as.matrix(dist(w_mat))
        sim_mat  <- exp(-gamma * dist_mat)
        diag(sim_mat) <- 1

        # optional ordering -------------------------------------------------
        ord <- if (reorder) hclust(as.dist(dist_mat))$order else seq_len(I)
        sim_ord <- sim_mat[ord, ord]
        lab_ord <- item_names[ord]

        # keep cells where column >= row (upper triangle incl. diag) -------
        keep <- upper.tri(sim_ord, diag = TRUE)
        df <- data.frame(
                item_row = rep(lab_ord, each = I)[keep],   # row index varies slowest
                item_col = rep(lab_ord,       I)[keep],   # column index varies fastest
                sim      = sim_ord[keep]
        )
        # axis factor levels: X left→right I1..I24, Y bottom→top I1..I24
        df$item_row <- factor(df$item_row, levels = lab_ord)
        df$item_col <- factor(df$item_col, levels = lab_ord)

        # plot -------------------------------------------------------------
        p <- ggplot(df, aes(.data$item_col, .data$item_row, fill = .data$sim)) +
                geom_tile(color = "white", linewidth = 0.4) +
                geom_text(aes(label = sprintf(paste0("%.", digits, "f"), .data$sim)), size = 3) +
                scale_fill_gradientn(colors = c("#fbe9ff", "#d07fff", "#9d1dff"),
                                     limits = c(0, 1), name = "Similarity") +
                coord_fixed() +
                theme_minimal(base_size = 11) +
                theme(panel.grid = element_blank(), axis.title = element_blank(),
                      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
        if (!is.null(title)) p <- p + ggtitle(title)
        print(p)
        invisible(p)
}
