#' @importFrom BiocFileCache bfcquery bfcrpath bfcadd bfcdownload
NULL

.matrix_query_url <- "https://matrix.data.humancellatlas.org/v0/matrix"

.rname_digest <- function(fq_idvec, fmt) {
    inlist <- list(bundles = sort(fq_idvec), format = fmt)
    digest::digest(inlist, algo = "md5")
}

.filter_digest <- function(api, fmt, feat) {
    inlist <- list(filter = filters(api), format = fmt, feature = feat)
    digest::digest(inlist, algo = "md5")
}

.api_download <- function(api, v1query, fq_ids, fmt, feat) {
    if (v1query)
        args <- list(filter = filters(api), format = fmt, feature = feat)
    else
        args <- list(bundle_fqids = fq_ids, format = fmt)

    endpoint <- paste0("matrix.lambdas.api.",
        if (v1query) "v1" else "v0", ".core.post_matrix")
    httr::content(
        do.call(
            .invoke_fun, c(api = api, name = endpoint, args)
        )
    )[["request_id"]]
}

.bind_content <- function(x) {
    dplyr::bind_rows(
        httr::content(x)
    )
}

.invoke_fun <- function(api, name, ...) {
    if (!is(api, "HCAMatrix"))
        stop("Provide a 'HCAMatrix' class API object")
    ops <- names(AnVIL::operations(api))
    if (!name %in% ops)
        stop("<internal> operation name not found in API")

    do.call(`$`, list(api, name))(...)
}


.getResponse <- function(api, req, v1) {
    endpoint <- paste0("matrix.lambdas.api.",
        if (v1) "v1" else "v0", ".core.get_matrix")
    httr::content(.invoke_fun(api, endpoint, request_id = req))
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
#' @param bundle_fqids character() v0 Bundle identifiers
#'
#' @param verbose logical (default FALSE) whether to output stepwise messages
#'
#' @param names.col character (default "CellID") The column name in the
#'     `colData`` metadata to use as column names of the
#'     \linkS4class{LoomExperiment} object when `format = "loom"`
#'
#' @param format character(1) Data return format, one of:
#'     c("loom", "mtx", "csv"); (default: "loom")
#'
#' @param feature character(1) Provide either cell by "gene" or "transcript"
#'     matrices (default: "gene")
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
#' @export
loadHCAMatrix <-
    function(api, bundle_fqids, verbose = FALSE, names.col = "CellID",
        format = c("loom", "mtx", "csv"), feature = c("gene", "transcript")
    )
{
    stopifnot(is.character(names.col), length(names.col) == 1L,
        !is.na(names.col), !is.logical(names.col))
    v1q <- missing(bundle_fqids)
    format <- match.arg(format)
    feature <- match.arg(feature)

    rname_digest <-
        if (v1q)
            .filter_digest(api, format, feature)
        else
            .rname_digest(bundle_fqids, format)
    bfc <- .get_cache()
    rid <- bfcquery(bfc, rname_digest, "rname")$rid
    if (!length(rid)) {
        req_id <- .api_download(api, v1q, bundle_fqids, format, feature)
        if (verbose)
            message("Matrix query request_id: ", req_id)
        pb <- progress::progress_bar$new(
            format = "  (:spin) waiting for query completion:dots :elapsedfull",
            total = NA, clear = FALSE)

        while (
            identical(.getResponse(api, req_id, v1q)[["status"]], "In Progress")
        ) {
            for (ndot in 0:10) {
                pb$tick(tokens = list(dots = .dotter(ndot, 10)))
                Sys.sleep(2/8)
            }
            cat("\n")
        }

        response_obj <- .getResponse(api, req_id, v1q)
        if (identical(response_obj[["status"]], "Failed"))
            stop(.msg(response_obj[["message"]]))
        matrix_url <- response_obj[["matrix_url"]]
        rid <- names(bfcadd(bfc, rname_digest, matrix_url))
    }

    if (verbose)
        message("Matrix data in cache with 'rname': ", rname_digest)

    mat_loc <- bfcrpath(bfc, rids = rid)

    if (identical(format, "loom")) {
        .checkPkgsAvail("LoomExperiment")
        lex <- LoomExperiment::import(mat_loc)
        idcol <- lex[[names.col]]
        if (length(idcol))
            colnames(lex) <- idcol
        lex
    } else if (identical(format, "mtx")) {
        .checkPkgsAvail("HCAmtxzip")
        HCAmtxzip::import_mtxzip(mat_loc)
    } else {
        .checkPkgsAvail("readr")
        readr::read_csv(mat_loc)
    }
}
