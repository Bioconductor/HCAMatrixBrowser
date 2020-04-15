#' @name HCAMatrix-class
#' @title A class for representing the HCAMatrix API
#'
#' @description The `HCAMatrix` class is a representation of the `HCAMatrix`
#' API protocol via OAS version 2.0. The original version OAS 3 was converted
#' using the APIMatic converter (\url{apimatic.io}).
#'
#' @importFrom methods new
#'
#' @seealso \link{HCAMatrix}, \link[AnVIL]{Service}
#'
#' @examples
#'
#' HCAMatrix()
#'
#' @export
.HCAMatrix <- setClass("HCAMatrix", contains = "Service",
    representation(filter = "list"))

#' @rdname HCAMatrix
#'
#' @aliases HCAMatrix
#'
#' @title API Entry function for the Human Cell Atlas Matrix service
#'
#' @description This function allows the use of the HCA Matrix API
#'
#' @param api An HCAMatrix API object
#'
#' @param filter_name character(1) The name of the filter to get more detail on
#'
#' @param format_name character(1) The format for which to obtain more detail on
#'
#' @param feature_name character(1) The feature for which to obtain more detail
#'     on
#'
#' @return An object of class 'HCAMatrix'
#'
#' @importFrom AnVIL Service
#'
#' @examples
#' hca <- HCAMatrix()
#'
#' @export
HCAMatrix <- function() {
    .HCAMatrix(
        Service(
            service = "HCAMatrix",
            host = "matrix.dev.data.humancellatlas.org",
            config = httr::config(ssl_verifypeer = 0L, ssl_verifyhost = 0L,
                http_version = 0L),
            authenticate = FALSE,
            api_url =
                system.file("service/swagger.yaml",
                    package = "HCAMatrixBrowser", mustWork = TRUE),
            package = "HCAMatrixBrowser",
            schemes = "https"
        )
    )
}

#' @name HCAMatrix
#'
#' @section Filters:
#'     * available_filter - Get a list of filters within the API to filter with
#'     * filter_detail - Obtain more detail on a particular filter name
#'
#' @export
available_filters <- function(api) {
    unlist(
        httr::content(
            api$matrix.lambdas.api.v1.core.get_filters()
        )
    )
}

#' @name HCAMatrix
#'
#' @export
filter_detail <- function(api, filter_name) {
    stopifnot(is.character(filter_name), length(filter_name) == 1L,
        !is.na(filter_name), !is.logical(filter_name))

    httr::content(
        api$matrix.lambdas.api.v1.core.get_filter_detail(
            filter_name = filter_name
        )
    )
}

#' @name HCAMatrix
#'
#' @section Formats:
#'     * available_formats - Get a list of matrix format outputs
#'     * format_detail - Obtain more detail on a particular matrix file format
#'
#' @export
available_formats <- function(api) {
    unlist(
        httr::content(
            api$matrix.lambdas.api.v0.core.get_formats()
        )
    )
}

#' @name HCAMatrix
#'
#' @export
format_detail <- function(api, format_name) {
    stopifnot(is.character(format_name), length(format_name) == 1L,
        !is.na(format_name), !is.logical(format_name))

    temphtml <- tempfile(fileext = ".html")
    writeLines(
        httr::content(
                api$matrix.lambdas.api.v1.core.get_format_detail(
                format_name = format_name
            )
        ), con = file(temphtml)
    )
    utils::browseURL(temphtml)
}

#' @name HCAMatrix
#'
#' @section Features:
#'     * available_features - Get a list of feature outputs, either genes or
#'         transcripts
#'     * feature_detail - Obtain more information on a matrix feature type
#'
#' @export
available_features <- function(api) {
    unlist(
        httr::content(
            api$matrix.lambdas.api.v1.core.get_features()
        )
    )
}

#' @name HCAMatrix
#'
#' @export
feature_detail <- function(api, feature_name) {
    httr::content(
        api$matrix.lambdas.api.v1.core.get_feature_detail(
            feature_name = feature_name
        )
    )
}

#' @name filtering
#'
#' @aliases filter filters
#'
#' @title Manipulating HCAMatrix filters
#'
#' @return A \code{\link{HCAMatrix}} object with the filter
#' field replaced by the specified filter expression
#'
#' @examples
#' # make an HCAMatrix object to start
#' hca <- HCAMatrix()
#'
#' head(available_filters(hca))
#'
#' hca1 <- filter(hca, genes_detected >= 500)
#' filters(hca1)
#'
#' @section filter:
#'     The \code{filter} is a convenient setter for the filter
#'     element in \code{\link{HCAMatrix}} objects.
#' @section filters:
#'     The \code{filters} (plural) function is a safe accessor for the filters
#'     already present in the `HCAMAtrix` API object. The filter can also be
#'     set using the `filters<-` function setter (advanced use).
#'
#' @param x the object on which to set the filter list
#' member
#'
#' @param expr a filter expression in the form of
#' the right hand side of a formula, where bare names
#' (without quotes) are allowed if they are available
#' fields associated with the HCAMAtrix object, \code{x}
#'
#' @param value A list of structured filters (internal use)
#'
#' @note Filtering documentation provided by the `GenomicDataCommons` package
#'
#' @export
setGeneric("filter", function(x, expr) { standardGeneric("filter") })

#' @rdname filtering
#'
#' @export
setGeneric("filters", function(x) { standardGeneric("filters") })

#' @rdname filtering
#'
#' @export
setGeneric("filters<-", function(x, value) { standardGeneric("filters<-") })

#' @rdname filtering
#'
#' @export
setMethod("filter", c("HCAMatrix", "ANY"), function(x, expr) {
    filt <- try({
        if (rlang::is_formula(expr))
            make_filter(expr, available_filters(x))
        }, silent = TRUE)
    if (inherits(filt, "try-error"))
        filt <- make_filter(rlang::enquo(expr), available_filters(x))
    filters(x) <- filt
    x
})

#' @rdname filtering
#'
#' @exportMethod filters
setMethod("filters", "HCAMatrix", function(x) {
    x@filter
})

#' @rdname filtering
#'
#' @export
setReplaceMethod("filters", c("HCAMatrix", "ANY"), function(x, value) {
    current <- filters(x)
    if (length(current))
        BiocGenerics:::replaceSlots(
            x,
            filter = list(op = "and", value = list(current, value))
        )
    else
        BiocGenerics:::replaceSlots(x, filter = value)
})

make_filter <- function(expr, available_filters) {
    afilters <- stats::setNames(available_filters, available_filters)
    filterlist <- c(.oplookup, afilters)

    if (rlang::is_formula(expr))
        rlang::eval_tidy(rlang::f_rhs(expr), data = filterlist,
            env = rlang::f_env(expr))
    else
        rlang::eval_tidy(expr, data = filterlist)
}
