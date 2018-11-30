#' @importFrom BiocFileCache bfcquery bfcrpath bfcadd bfcdownload
NULL

.matrix_query_url <- "https://matrix.data.humancellatlas.org/v0/matrix"

.initialize_download <- function(fq_ids, query_url = .matrix_query_url,
    format = "loom")
{
    if (missing(fq_ids))
        stop("<internal> Provide valid fq_ids vector")
    body <- list(bundle_fqids = fq_ids, format = format)
    header <- list(`Content-Type` = "application/json",
        Accept = "application/json")
    response <- httr::POST(query_url, header, body = body, encode = "json",
        httr::verbose())
    httr::content(response)[["request_id"]]
}

#' @name HCAMatrix
#' @title Obtain expression matrix data from the Human Cell Atlas API service
#'
#' @description Using a vector of data bundle identifiers (`bundle_fqids`), users
#' can request the associated matrix of expression values. The query submitted
#' by `getHCAMatrixID` can take some time to be completed. Once the query is
#' completed, users can use the `request_id` to load the dataset as a
#' `LoomExperiment` object using the `loadHCAMatrix` function.
#'
#' @details The matrix_query_url value points to
#'     \url{https://matrix.data.humancellatlas.org/v0/matrix}
#'
#' @param bundle_fqids A character vector of bundle identifiers
#' @param matrix_query_url A single character vector of the API endpoint for
#' downloading matrix data (defaults to HCA matrix endpoint)
#' @param verbose logical (default TRUE) whether to output stepwise messages
#'
#' @return getHCAMatrixID: A request identifier as a single string
#'
#' @examples
#'
#' bundle_fqids <-
#'     c("980b814a-b493-4611-8157-c0193590e7d9.2018-11-12T131442.994059Z",
#'     "7243c745-606d-4827-9fea-65a925d5ab98.2018-11-07T002256.653887Z")
#'
#' req_id <- getHCAMatrixID(bundle_fqids)
#'
#' loadHCAMatrix(req_id)
#'
#' @export
getHCAMatrixID <- function(bundle_fqids, matrix_query_url = .matrix_query_url,
    verbose = TRUE)
{
    req_id <- .initialize_download(bundle_fqids, query_url = matrix_query_url)
    if (verbose)
        message("Query request ID: ", req_id)
    req_address <- file.path(matrix_query_url, req_id)

    while (
        identical(httr::content(httr::GET(req_address))[["status"]],
            "In Progress")
    )
        Sys.sleep(10)
    get_response <- httr::GET(req_address)
    matrix_url <- httr::content(get_response)[["matrix_location"]]

    bfc <- .get_cache()
    rid <- bfcquery(bfc, req_id, "rname")$rid
    if (!length(rid)) {
        rid <- names(bfcadd(bfc, req_id, matrix_url, download = FALSE))
    }
    if (!.cache_exists(bfc, req_id)) {
        if (verbose)
            message("Downloading matrix data for request ID: ", req_id)
            bfcdownload(bfc, rid, ask = FALSE)
    } else
        message("Matrix data already in cache: ", req_id)

    req_id
}

#' @rdname HCAMatrix
#'
#' @param request_id A single string obtained from the getHCAMatrixID function
#' @return loadHCAMatrix: A LoomExperiment object
#'
#'
#' @export
loadHCAMatrix <- function(request_id) {
    bfc <- .get_cache()
    rid <- bfcquery(bfc, request_id, "rname")$rid
    if (!length(rid))
        stop("'request_id' not found in cache")
    mat_loc <- bfcrpath(bfc, rids = rid)
    LoomExperiment::import(mat_loc)
}
