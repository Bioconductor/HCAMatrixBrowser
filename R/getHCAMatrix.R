.matrix_query_url <- "https://matrix.data.humancellatlas.org/v0/matrix"

.download_matrix <- function(fq_ids, format = "loom") {
    body <- list(bundle_fqids = fq_ids, format = format)
    header <- .build_header(FALSE)
    response <- httr::POST(.matrix_query_url, header, body = body, encode = "json", httr::verbose())
    req_id <- httr::content(response)$request_id
}
