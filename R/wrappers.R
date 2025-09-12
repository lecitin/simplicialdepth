input_check <- function(X, x=NULL) {
    if (!is.matrix(X) || !is.numeric(X))
        stop("X must be a numeric matrix")

    n <- nrow(X)
    d <- ncol(X)

    if (!is.null(x)) {
        if ((!is.vector(x) && !is.matrix(x)) || !is.numeric(x))
            stop("x must be a numeric matrix or a numeric vector.")
        if (is.matrix(x)) {
            if (ncol(x) != d) {
                stop("x must have the same number of columns as X (i.e. ncol(x) = ncol(X)).")
            }
        } else {
            if (length(x) != d) {
                stop("x must have length equal to the number of columns of X (i.e. length(x) = ncol(X)).")
            }
        }
    }
}

# Constructor for depth_result
depth_result <- function(depth, max_depth, max_point, max_index) {
  structure(
    list(
      depth     = depth,
      max_depth = max_depth,
      max_point = max_point,
      max_index = max_index
    ),
    class = "depth_result"
  )
}

# Print method
print.depth_result <- function(x, ...) {
  cat("Depth result\n")
  cat("  Number of points:", length(x$depth), "\n")
  cat("  Maximum depth   :", x$max_depth, "\n")
  cat("  At point        :", paste(x$max_point, collapse = ", "), "\n")
  cat("  With index      :", paste(x$max_index, collapse = ", "), "\n")
  invisible(x)
}

spherical_asd_on_one_point <- function(X, i) {
    new_X <- X[-i,]
    x <- X[i,]
    return(spherical_asd(new_X, x))
}

#' Angular Simplicial Depth
#'
#' Computes the angular simplicial depth for 2D or 3D data, either for all
#' points in a sample or for user-specified query points.
#'
#' @param X A numeric matrix of size \eqn{n \times d} containing the sample points, where each row is an observation and each column is a coordinate. 
#'   All points are first projected onto a unit circle/sphere and then containting arcs/spherical triangles are counted.
#'   The dimension \eqn{d} must be either 2 or 3.
#' @param x Optional. A query point or set of points:
#'   \itemize{
#'     \item If `NULL` (default):
#'       \itemize{
#'         \item In 2D, computes the angular simplicial depth for all points 
#'         in `X`, along with arc depths and median information.
#'         \item In 3D, computes the depth of each row of `X` with respect to 
#'         the remaining points.
#'       }
#'     \item If a numeric vector of length \eqn{d}, computes the depth of this 
#'       single point with respect to `X`.
#'     \item If a numeric matrix with \eqn{d} columns, computes the depth of 
#'       each row of `x` with respect to `X`.
#'   }
#'
#' @details
#' In 2D, the function normalizes arc and endpoint depths by the appropriate
#' binomial coefficients. It also identifies the point (the "median") with
#' maximum angular simplicial depth.  
#'
#' In 3D, depths are normalized by \eqn{choose(n, 3)} or \eqn{choose(n-1, 3)} 
#' depending on the case.
#'
#' @return
#' A numeric value, numeric vector, or a list depending on the input:
#' \itemize{
#'   \item If `x` is `NULL` and \eqn{d=2}, returns a list with components:
#'     \describe{
#'       \item{arcs}{Data frame with arc depths and endpoint depths.}
#'       \item{median_index}{Index of the point with maximum depth.}
#'       \item{median}{The corresponding arc endpoint.}
#'       \item{max_depth}{Maximum endpoint depth.}
#'     }
#'     In this case the whole run takes \eqn{O(n\log(n))}.
#'   \item Otherwise, returns a numeric value or vector of depths.
#' }
#'
#' @examples
#' # Simulate some 2D data
#' set.seed(123)
#' X <- matrix(rnorm(20), ncol = 2)
#'
#' # Depths for each point in X
#' angularsimplicialdepth(X)
#'
#' # Depth for a single query point
#' angularsimplicialdepth(X, c(0,0))
#'
#' # Depths for multiple query points
#' Y <- matrix(rnorm(6), ncol = 2)
#' angularsimplicialdepth(X, Y)
#'
#' @export
angularsimplicialdepth <- function(X, x = NULL) {
    input_check(X, x)

    n <- nrow(X)
    d <- ncol(X)
    norms <- sqrt(rowSums(X^2))
    X <- X / norms

    if (d == 2) {
        if (is.null(x)) {
            if (n < 3)
                stop("For x=NULL, X must have at least 3 rows.")
            # do something (e.g. call your C++ function for global depth)
            result <- circular_asd_all_arcs(X)
            result$depth = result$depth / choose(n, 2)
            result$end_point_depth = result$end_point_depth / choose(n-1, 2)
            max_depth_index <- which.max(result$end_point_depth)
            # median_row <- result[max_depth_index, ]
            median <- c(result$right_point_x[max_depth_index], result$right_point_y[max_depth_index])
            max_depth <- result$end_point_depth[max_depth_index]
            # return(list(arcs=result,
            #             median_index=max_depth_index,
            #             median=median,
            #             max_depth=max_depth))
            return(depth_result(depth=result$end_point_depth, max_depth=max_depth, max_point=median, max_index=max_depth_index))
        } else {
            if (is.matrix(x)) {
                result <- apply(x, 1, function(row) circular_asd(X, row)) / choose(n, 2)
                max_depth_index <- which.max(result)
                return(depth_result(depth=result,
                                    max_depth=result[max_depth_index],
                                    max_point=x[max_depth_index, ],
                                    max_index=max_depth_index))
            }
            # The only option left is that x is a numeric vector of proper length
            result <- circular_asd(X, x) / choose(n, 2)
            return(depth_result(depth=c(result),
                                max_depth=result,
                                max_point=x,
                                max_index=1))
        }
    }
    if (d == 3) {
        if (is.null(x))
            if (n < 4)
                stop("For x=NULL, X must have at least 3 rows.")
            # stop("For spherical (3D) angular simplicial depth, x cannot be NULL.")
            result <- vapply(seq_len(nrow(X)), function(i) spherical_asd_on_one_point(X, i), numeric(1)) / choose(n-1, 3)
            max_depth_index <- which.max(result)
            return(depth_result(depth=result,
                                max_depth=result[max_depth_index],
                                max_point=X[max_depth_index, ],
                                max_index=max_depth_index))

        if (is.matrix(x)) {
            result <- apply(x, 1, function(row) spherical_asd(X, row)) / choose(n, 3)
            max_depth_index <- which.max(result)
            return(depth_result(depth=result,
                                max_depth=result[max_depth_index],
                                max_point=x[max_depth_index, ],
                                max_index=max_depth_index))
        }
        result <- spherical_asd(X, x) / choose(n, 3)
        return(depth_result(depth=c(result),
                            max_depth=result,
                            max_point=x,
                            max_index=1))
    }
    stop("Not implemented: only 2D or 3D simplicial depth is supported")
}

sd2d_on_one_point <- function(X, i) {
    new_X <- X[-i,]
    x <- X[i,]
    return(simplicial_depth_2d(new_X, x))
}

sdk_on_one_point <- function(X, i, k) {
    new_X <- X[-i,]
    x <- X[i,]
    return(SDk_parallel(new_X, x, k))
}

#' Simplicial Depth
#'
#' Computes the simplicial depth for 2D or 3D data, either for all points in a 
#' sample or for user-specified query points.
#'
#' @param X A numeric matrix of size \eqn{n \times d} containing the sample 
#'   points, where each row is an observation and each column is a coordinate. 
#'   The dimension \eqn{d} must be either 2 or 3.
#' @param x Optional. A query point or set of points:
#'   \itemize{
#'     \item If `NULL` (default):
#'       \itemize{
#'         \item In 2D, computes the simplicial depth for each point in `X` 
#'         with respect to the rest of the sample.
#'         \item In 3D, computes the simplicial depth for each point in `X`, 
#'         provided there are at least 5 sample points. If `n \le 4`, returns 
#'         a zero vector of length `n`.
#'       }
#'     \item If a numeric vector of length \eqn{d}, computes the depth of this 
#'       single point with respect to `X`.
#'     \item If a numeric matrix with \eqn{d} columns, computes the depth of 
#'       each row of `x` with respect to `X`.
#'   }
#'
#' @details
#' In both 2D and 3D, depths are normalized by the corresponding binomial 
#' coefficients:
#' \itemize{
#'   \item In 2D: division by \eqn{choose(n-1, 3)} or \eqn{choose(n, 3)} 
#'   depending on whether depths are computed for sample points or query points.
#'   The whole run takes \eqn{O(n^2\log(n))} time.
#'   \item In 3D: division by \eqn{choose(n-1, 4)} or \eqn{choose(n, 4)}. 
#'   In this case the whole run takes \eqn{O(n^3\log(n))}.
#' }
#'
#' @return A numeric vector of depths (one per query point or per row of `X`), 
#' depending on the input.  
#'
#' @examples
#' # Simulate some 2D data
#' set.seed(123)
#' X <- matrix(rnorm(20), ncol = 2)
#'
#' # Depths for each point in X
#' simplicialdepth(X)
#'
#' # Depth for a single query point
#' simplicialdepth(X, c(0,0))
#'
#' # Depths for multiple query points
#' Y <- matrix(rnorm(6), ncol = 2)
#' simplicialdepth(X, Y)
#'
#' # Simulate some 3D data
#' X3 <- matrix(rnorm(30), ncol = 3)
#' simplicialdepth(X3)  # depths for each point in X3
#'
#' @export
simplicialdepth <- function(X, x=NULL) {
    input_check(X, x)

    n <- nrow(X)
    d <- ncol(X)

    if (d == 2) {
        if (is.null(x)) {
            if (n < 4)
                stop("For x=NULL, X must have at least 4 rows.")
            result <- vapply(seq_len(nrow(X)), function(i) sd2d_on_one_point(X, i), numeric(1)) / choose(n-1, 3)
            max_depth_index <- which.max(result)
            return(depth_result(depth=result,
                                max_depth=result[max_depth_index],
                                max_point=X[max_depth_index, ],
                                max_index=max_depth_index))
        }
        if (is.matrix(x)) {
            result <- apply(x, 1, function(row) simplicial_depth_2d(X, row)) / choose(n, 3)
            max_depth_index <- which.max(result)
            return(depth_result(depth=result,
                                max_depth=result[max_depth_index],
                                max_point=x[max_depth_index, ],
                                max_index=max_depth_index))
        }
        # The only remaining option is that x is a numeric vector
        result <- simplicial_depth_2d(X,x) / choose(n, 3)
        return(depth_result(depth=c(result),
                            max_depth=result,
                            max_point=x,
                            max_index=1))
    }
    if (d == 3) {
        if (is.null(x)) {
            if (n < 5)
                stop("For x=NULL, X must have at least 5 rows.")
            result <- vapply(seq_len(nrow(X)), function(i) sdk_on_one_point(X, i, 4), numeric(1)) / choose(n-1, 4)
            max_depth_index <- which.max(result)
            return(depth_result(depth=result,
                                max_depth=result[max_depth_index],
                                max_point=X[max_depth_index, ],
                                max_index=max_depth_index))
        }
        if (is.matrix(x)) {
            result <- apply(x, 1, function(row) SDk_parallel(X, row, 4)) / choose(n, 4)
            max_depth_index <- which.max(result)
            return(depth_result(depth=result,
                                max_depth=result[max_depth_index],
                                max_point=x[max_depth_index, ],
                                max_index=max_depth_index))
        }
        result <- SDk_parallel(X, x, 4) / choose(n, 4)
        return(depth_result(depth=c(result),
                            max_depth=result,
                            max_point=x,
                            max_index=1))
    }
}

#' k-Hull Depth
#'
#' Computes the \eqn{k}-hull depth for 2D or 3D data, 
#' either for all points in a sample or for user-specified query points.
#'
#' @param X A numeric matrix of size \eqn{n \times d} containing the sample 
#'   points, where each row is an observation and each column is a coordinate.
#' @param x Optional. A query point or set of points:
#'   \itemize{
#'     \item If `NULL` (default): computes the \eqn{k}-hull depth for each row 
#'       of `X` with respect to the rest of the sample.
#'     \item If a numeric vector of length \eqn{d}: computes the depth of this 
#'       single point with respect to `X`.
#'     \item If a numeric matrix with \eqn{d} columns: computes the depth of 
#'       each row of `x` with respect to `X`.
#'   }
#' @param k Integer, the hull size parameter. Must satisfy \eqn{k > d} and 
#'   \eqn{k \ge n} to yield nontrivial depths. If \eqn{k \le d} or 
#'   \eqn{k < n}, the function returns zero depths.
#'
#' @details
#' The \eqn{k}-hull depth is a generalization of simplicial depth, replacing 
#' simplices by convex hulls with \eqn{k} vertices. Depths are normalized by the 
#' corresponding binomial coefficients:
#'
#' @return A numeric vector of depths, one per query point (or one per row of 
#' `X` if `x = NULL`). If `k \le d` or `k < n`, depths are identically zero.
#'
#' @examples
#' # Simulate some 2D data
#' set.seed(123)
#' X <- matrix(rnorm(20), ncol = 2)
#'
#' # Depths for each point in X with k = 4
#' khulldepth(X, k = 4)
#'
#' # Depth for a single query point
#' khulldepth(X, c(0,0), k = 4)
#'
#' # Depths for multiple query points
#' Y <- matrix(rnorm(6), ncol = 2)
#' khulldepth(X, Y, k = 4)
#'
#' @export
khulldepth <- function(X, x=NULL, k) {
    input_check(X, x)

    n <- nrow(X)
    d <- ncol(X)

    if (k <= d || n < k) {
        if (is.null(x)) {
            if (n < k+1)
                return(depth_result(depth=numeric(n), max_depth=0, max_point=X[1,], max_index=1))
        } else { 
            if (is.matrix(x))
                return(depth_result(depth=numeric(nrow(x)), max_depth=0, max_point=X[1,], max_index=1))
            # The only option left is that x is a vector
            return(depth_result(depth=numeric(1), max_depth=0, max_point=X[1,], max_index=1))
        }
    }

    if (is.null(x)) {
        if (n<k+1)
            stop("X must have at least k+1 rows.")
        result <- vapply(seq_len(nrow(X)), function(i) sdk_on_one_point(X, i, k), numeric(1)) / choose(n-1, k)
        max_depth_index <- which.max(result)
        return(depth_result(depth=result,
                            max_depth=result[max_depth_index],
                            max_point=X[max_depth_index, ],
                            max_index=max_depth_index))
    }
    if (is.matrix(x)) {
        result <- apply(x, 1, function(row) SDk_parallel(X, row, k)) / choose(n, k)
        max_depth_index <- which.max(result)
        return(depth_result(depth=result,
                            max_depth=result[max_depth_index],
                            max_point=x[max_depth_index, ],
                            max_index=max_depth_index))
    }
    # The only option left is that x is a numeric vector of length equal to d
    result <- SDk_parallel(X, x, k) / choose(n, k)
    return(depth_result(depth=c(result),
                        max_depth=result,
                        max_point=x,
                        max_index=1))
}
