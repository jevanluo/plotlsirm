#' Item‑strength profile for a single person
#'
#' For a chosen respondent (`person_index`) this function plots the **strength**
#' (probability of endorsement) for every item, defined as
#' \deqn{\exp(-\gamma\,d_{ij})}{exp(-γ·d_ij)}, where \eqn{d_{ij}} is the
#' Euclidean distance between the person’s latent position \eqn{z_j} and each
#' item position \eqn{w_i}.  When `z` and `w` are supplied as *lists* of
#' matrices (posterior draws), the function summarises the distribution of
#' strengths with medians and a `ci_level` credible interval.  Bars can be
#' coloured by an item grouping factor, reordered by decreasing strength, and
#' displayed either vertically or horizontally.
#'
#' @param z A numeric matrix (*N* × *d*) of person coordinates **or** a *list*
#'   of such matrices representing posterior draws.
#' @param w A numeric matrix (*I* × *d*) of item coordinates **or** a *list* of
#'   such matrices, matching the structure of `z`.
#' @param person_index Integer giving the row of `z` (or each draw in `z`)
#'   corresponding to the focal respondent.
#' @param gamma Positive numeric scalar controlling the decay of strength with
#'   distance; default is `1`.
#' @param item_group Optional character/factor vector of length *I* assigning
#'   each item to a group for colour coding and legend.
#' @param item_names Optional character vector of item labels.  If `NULL`
#'   defaults to `"I1"`, `"I2"`, … .
#' @param ci_level Width of the credible interval (between 0 and 1) when
#'   posterior draws are given.  Ignored for a single point estimate.
#' @param reorder Logical.  Reorder items on the axis by decreasing strength?
#'   Default `FALSE`.
#' @param vertical Logical.  `TRUE` (default) draws vertical bars; `FALSE`
#'   flips the axes for a horizontal layout.
#' @param title Optional character string to appear as the plot title.
#'
#' @return (Invisibly) a `ggplot` object containing the bar plot.  The plot is
#'   also printed.
#'
#' @import ggplot2
#' @importFrom stats median quantile
#' @importFrom rlang .data
#'
#' @examples
#' set.seed(1)
#' z  <- matrix(rnorm(40), ncol = 2)   # 20 persons
#' w  <- matrix(rnorm(30), ncol = 2)   # 15 items
#'
#' ## Point‑estimate strengths for person 5
#' strengthplot(z, w, person_index = 5, gamma = 2)
#'
#' ## Posterior example with credible intervals and item groups
#' draws_z <- replicate(50, z + matrix(rnorm(length(z), sd = 0.1),
#'                                     nrow(z), ncol(z)), simplify = FALSE)
#' draws_w <- replicate(50, w + matrix(rnorm(length(w), sd = 0.1),
#'                                     nrow(w), ncol(w)), simplify = FALSE)
#' grp <- rep(c("Core", "Peripheral"), length.out = nrow(w))
#' strengthplot(draws_z, draws_w, person_index = 3,
#'              item_group = grp, ci_level = 0.9, vertical = FALSE,
#'              title = "Posterior strength profile for respondent 3")
#'
#' @export
strengthplot <- function(
                z, w,
                person_index,
                gamma       = 1,
                item_group  = NULL,
                item_names  = NULL,
                ci_level    = 0.95,
                reorder     = FALSE,
                vertical    = TRUE,
                title       = NULL
) {
        # internal strength fun
        str_fun <- function(d) exp(-gamma * d)
        is_post <- is.list(z) && is.list(w)

        # helper distance
        dist_vec_mat <- function(vec, mat) sqrt(rowSums((mat - matrix(vec, nrow(mat), ncol(mat), byrow = TRUE))^2))

        build_single_df <- function(zm, wm) {
                I <- nrow(wm)
                if (is.null(item_names)) nm <- paste0("I", seq_len(I)) else nm <- item_names
                d <- dist_vec_mat(zm[person_index, ], wm)
                s <- str_fun(d)
                ord <- if (reorder) nm[order(-s)] else nm
                data.frame(item = factor(nm, levels = ord), strength = s,
                           group = if (is.null(item_group)) NA_character_ else item_group)
        }

        if (!is_post) {
                df <- build_single_df(as.matrix(z), as.matrix(w))
                p <- ggplot(df, aes(x = .data$item, y = .data$strength,
                                    fill = if (is.null(item_group)) .data$strength else .data$group)) +
                        geom_col()
                if (is.null(item_group)) {
                        p <- p + scale_fill_gradientn(colours = c("yellow", "orange", "red"), guide = "none")
                } else {
                        p <- p + scale_fill_brewer(palette = "Set2", name = "Item group")
                }
        } else {
                M <- length(z)
                N <- nrow(as.matrix(z[[1]]))
                I <- nrow(as.matrix(w[[1]]))
                if (person_index > N) stop("person_index exceeds N")
                if (!is.null(item_group) && length(item_group)!=I) stop("item_group length mismatch")
                str_mat <- matrix(NA_real_, I, M)
                for (m in seq_len(M)) {
                        str_mat[, m] <- str_fun(dist_vec_mat(as.matrix(z[[m]])[person_index, ], as.matrix(w[[m]])))
                }
                med <- apply(str_mat,1,stats::median)
                lwr <- apply(str_mat,1,stats::quantile,probs=(1-ci_level)/2);
                upr <- apply(str_mat,1,stats::quantile,probs=1-(1-ci_level)/2)
                nm <- if (is.null(item_names)) paste0("I", seq_len(I)) else item_names
                ord <- if (reorder) nm[order(-med)] else nm
                df <- data.frame(item=factor(nm,levels=ord), median=med, lower=lwr, upper=upr,
                                 group = if(is.null(item_group)) NA_character_ else item_group)

                p <- ggplot(df, aes(x=.data$item, y=.data$median,
                                    fill = if (is.null(item_group)) .data$median else .data$group)) +
                        geom_col() +
                        geom_errorbar(aes(ymin=.data$lower, ymax=.data$upper), width=0.25)
                if (is.null(item_group)) {
                        p <- p + scale_fill_gradientn(colours=c("yellow","orange","red"), guide="none")
                } else {
                        p <- p + scale_fill_brewer(palette="Set2", name="Item group")
                }
        }

        # coord flip for horizontal
        if (!vertical) p <- p + coord_flip()

        # axis & theme
        p <- p + labs(x=NULL, y="Likelihood of endorsement") +
                scale_y_continuous(limits=c(0,1), breaks=seq(0,1,0.2), expand=expansion(mult=c(0,0.02))) +
                theme(panel.background=element_rect(fill="#F5F5F5"),
                      panel.grid.major=element_line(color="white"),
                      panel.border=element_rect(color="grey60", fill=NA, linewidth=0.8),
                      axis.text.x=element_text(angle= if(vertical) 90 else 0, vjust=0.5,hjust=1),
                      axis.ticks.y=element_line(color="grey50"))
        if(!is.null(title)) p <- p + ggtitle(title)
        print(p)
        invisible(p)
}
