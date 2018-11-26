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

getHCAMatrixID <- function(bundle_fqids, matrix_query_url = .matrix_query_url,
    verbose = TRUE, force = FALSE)
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
    if (!.cache_exists(bfc, req_id) || force) {
        if (verbose)
            message("Downloading matrix data for request ID: ", req_id)
            bfcdownload(bfc, rid, ask = FALSE)
    } else
        message("Matrix data already in cache: ", req_id)

    req_id
}

loadHCAMatrix <- function(request_id) {
    bfc <- .get_cache()
    rid <- bfcquery(bfc, request_id, "rname")$rid
    if (!length(rid))
        stop("'request_id' not found in cache")
    mat_loc <- bfcrpath(bfc, rids = rid)
    LoomExperiment::import(mat_loc)
}
