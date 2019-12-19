#' @importFrom BiocFileCache bfcquery bfcrpath bfcadd bfcdownload
NULL

.matrix_query_url <- "https://matrix.data.humancellatlas.org/v0/matrix"

.rname_digest <- function(fq_idvec) {
    digest::digest(sort(fq_idvec), algo = "md5")
}

.initialize_download <- function(fq_ids, query_url = .matrix_query_url,
    format = "loom", verbose)
{
    if (missing(fq_ids))
        stop("<internal> Provide valid fq_ids vector")
    body <- list(bundle_fqids = fq_ids, format = format)
    header <- list(`Content-Type` = "application/json",
        Accept = "application/json")
    response <- httr::POST(query_url, header, body = body, encode = "json",
        if (verbose) httr::verbose())
    httr::content(response)[["request_id"]]
}

.api_download <- function(api, fq_ids) {
    bundles <- list(bundle_fqids = fq_ids)
    response <- api$matrix.lambdas.api.v0.core.post_matrix(bundles)
    httr::content(response)[["request_id"]]
}

.bind_content <- function(x) {
    dplyr::bind_rows(
        httr::content(x)
    )
}

.dotter <- function(ndots, maxlength) {
    paste0(
        paste0(rep(".", times = ndots), collapse = ""),
        paste0(rep(" ", times = maxlength-ndots), collapse = ""),
        collapse = ""
    )
}

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
#' @param api An API object of class `HCAMatrix` from the `HCAMatrix`
#'     function
#'
#' @param bundle_fqids A character vector of bundle identifiers
#'
#' @param verbose logical (default FALSE) whether to output stepwise messages
#'
#' @param names.col character (default "CellID") The column name in the
#'     `colData`` metadata to use as column names of the
#'      \linkS4class{LoomExperiment} object
#'
#' @return A \linkS4class{LoomExperiment} object
#'
#' @import LoomExperiment
#'
#' @examples
#'
#' hca <- HCAMatrix()
#'
#' bundle_fqids <-
#'     c("ffd3bc7b-8f3b-4f97-aa2a-78f9bac93775.2019-05-14T122736.345000Z",
#'     "f69b288c-fabc-4ac8-b50c-7abcae3731bc.2019-05-14T120110.781000Z",
#'     "f8ba80a9-71b1-4c15-bcfc-c05a50660898.2019-05-14T122536.545000Z",
#'     "fd202a54-7085-406d-a92a-aad6dd2d3ef0.2019-05-14T121656.910000Z",
#'     "fffe55c1-18ed-401b-aa9a-6f64d0b93fec.2019-05-17T233932.932000Z")
#'
#' loadHCAMatrix(hca, bundle_fqids)
#'
#' @export loadHCAMatrix
loadHCAMatrix <-
    function(api, bundle_fqids, verbose = FALSE, names.col = "CellID")
{
    stopifnot(is.character(names.col), length(names.col) == 1L,
        !is.na(names.col), !is.logical(names.col))

    rname_digest <- .rname_digest(bundle_fqids)
    bfc <- .get_cache()
    rid <- bfcquery(bfc, rname_digest, "rname")$rid
    if (!length(rid)) {
        req_id <- .api_download(api, bundle_fqids)
        if (verbose)
            message("Matrix query request_id: ", req_id)
        getResponse <- function(req) {
            httr::content(
                api$matrix.lambdas.api.v1.core.get_matrix(request_id = req)
            )
        }
        pb <- progress::progress_bar$new(
            format = "  (:spin) waiting for query completion:dots :elapsedfull",
            total = NA, clear = FALSE)

        while (
            identical(getResponse(req_id)[["status"]], "In Progress")
        ) {
            for (ndot in 0:10) {
                pb$tick(tokens = list(dots = .dotter(ndot, 10)))
                Sys.sleep(2/8)
            }
            cat("\n")
        }

        response_obj <- getResponse(req_id)
        if (identical(response_obj[["status"]], "Failed"))
            stop(.msg(response_obj[["message"]]))
        matrix_url <- response_obj[["matrix_url"]]
        rid <- names(bfcadd(bfc, rname_digest, matrix_url))
    }

    if (verbose)
        message("Matrix data in cache with 'rname': ", rname_digest)

    mat_loc <- bfcrpath(bfc, rids = rid)

    lex <- LoomExperiment::import(mat_loc)
    idcol <- lex[[names.col]]
    if (length(idcol))
        colnames(lex) <- idcol
    lex
}
