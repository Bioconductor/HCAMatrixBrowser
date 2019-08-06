#' @export
.HCAMatrix <- setClass("HCAMatrix", contains = "Service",
    representation(filter = "list"))

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

setGeneric("filter", function(x, ...) { standardGeneric("filter") })

setGeneric("filters", function(x, value) { standardGeneric("filters") })

setGeneric("filters<-", function(x, value) { standardGeneric("filters<-") })

#' @export
setMethod("filter", "HCAMatrix", function(x, expr) {
    filt <- make_filter(enquo(expr), available_filters(x))
    filters(x) <- filt
    x
})

#' @export
setMethod("filters", "HCAMatrix", function(x) {
    x@filter
})

#' @export
setReplaceMethod("filters", c("HCAMatrix", "list"), function(x, value) {
    current <- filters(x)
    if (length(current))
        slot(x, "filter") <- list( op = "and", content = list(current, value) )
    else
        slot(x, "filter") <- value
    x
})

make_filter <- function(expr, available_filters) {
    afilters <- setNames(available_filters, available_filters)
    filterlist <- c(.oplookup, afilters)
    rlang::eval_tidy(expr, data = filterlist)
}
