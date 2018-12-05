#' @importFrom BiocFileCache bfcquery bfcrpath bfcadd bfcdownload
NULL

.matrix_query_url <- "https://matrix.data.humancellatlas.org/v0/matrix"

.rname_digest <- function(fq_idvec) {
    digest::digest(sort(fq_idvec), algo = "md5")
}

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
#' @description Using a vector of data bundle identifiers (`bundle_fqids`),
#' users can request the associated matrix of expression values. The query
#' submitted by `loadHCAMatrix` may take some time to be completed. Once the
#' query is completed, a `LoomExperiment` object is loaded.
#'
#' @details The matrix_query_url value points to
#'     \url{https://matrix.data.humancellatlas.org/v0/matrix}
#'
#' @param bundle_fqids A character vector of bundle identifiers
#' @param matrix_query_url A single character vector of the API endpoint for
#' downloading matrix data (defaults to HCA matrix endpoint)
#' @param verbose logical (default FALSE) whether to output stepwise messages
#'
#' @return A LoomExperiment object
#'
#' @examples
#'
#' bundle_fqids <-
#'     c("980b814a-b493-4611-8157-c0193590e7d9.2018-11-12T131442.994059Z",
#'     "7243c745-606d-4827-9fea-65a925d5ab98.2018-11-07T002256.653887Z")
#'
#' loadHCAMatrix(bundle_fqids)
#'
#' @export
loadHCAMatrix <- function(bundle_fqids, matrix_query_url = .matrix_query_url,
    verbose = FALSE)
{
    rname_digest <- .rname_digest(bundle_fqids)
    bfc <- .get_cache()
    rid <- bfcquery(bfc, rname_digest, "rname")$rid
    if (!length(rid)) {
        req_id <-
            .initialize_download(bundle_fqids, query_url = matrix_query_url)
        if (verbose)
            message("Matrix query request_id: ", req_id)

        req_address <- file.path(matrix_query_url, req_id)

        pb <- progress::progress_bar$new(
            format = "  (:spin) waiting for query completion:dots :elapsedfull",
            total = NA, clear = FALSE)

        while (
            identical(httr::content(httr::GET(req_address))[["status"]],
                      "In Progress")
        ) {
            pb$tick(tokens = list(dots = "    "))
            Sys.sleep(1/8)
            pb$tick(tokens = list(dots = ".   "))
            Sys.sleep(1/8)
            pb$tick(tokens = list(dots = "..  "))
            Sys.sleep(1/8)
            pb$tick(tokens = list(dots = "... "))
            Sys.sleep(1/8)
            pb$tick(tokens = list(dots = "...."))
            Sys.sleep(1/8)
        }

        get_response <- httr::GET(req_address)
        matrix_url <- httr::content(get_response)[["matrix_location"]]
        rid <- names(bfcadd(bfc, rname_digest, matrix_url))
    }

    if (verbose)
        message("Matrix data in cache with 'rname': ", rname_digest)

    mat_loc <- bfcrpath(bfc, rids = rid)

    LoomExperiment::import(mat_loc)
}
