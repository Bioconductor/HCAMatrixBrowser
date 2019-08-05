#' @export
.HCAMatrix <- setClass("HCAMatrix", contains = "Service")

#' API Entry function for the Human Cell Atlas Matrix service
#'
#' This function allows the use of the HCA Matrix API
#'
#' @return An object of class 'cBioPortal'
#'
#' @importFrom AnVIL Service
#'
#' @export
HCAMatrix <- function() {
    .HCAMatrix(
        Service(
            service = "HCAMatrix",
            host = "matrix.dev.data.humancellatlas.org",
            config = httr::config(ssl_verifypeer = 0L, ssl_verifyhost = 0L,
                http_version = 0L),
            package = "HCAMatrixBrowser",
            schemes = "https"
        )
    )
}

#' @export
available_filters <- function(hca) {
    unlist(
        httr::content(
            hca$matrix.lambdas.api.v1.core.get_filters()
        )
    )
}

#' @export
filter_detail <- function(hca, filter = "genes_detected") {
    httr::content(
        hca$matrix.lambdas.api.v1.core.get_filter_detail(filter_name = filter)
    )
}

