.MTX_ARCHIVE_FILES <- c(
    "features.tsv.gz", "matrix.mtx.gz", "cells.tsv.gz", "barcodes.tsv.gz"
)

.read_tsv <-
    function(path, ..., stringsAsFactors = FALSE, sep = "\t", header = FALSE)
{
    read.delim(
        path, ...,
        stringsAsFactors = stringsAsFactors,
        sep = sep, header = header
    )
}

#' @importFrom Matrix sparseMatrix
.read_mtx <-
    function(path, verbose = FALSE)
{
    headers <- readLines(path, 2L)
    dims <- as.integer(strsplit(headers[2], " ")[[1]][1:2])
    !verbose || .message("dim: ", dims[1], " ", dims[2])
    v <- scan(
        path, list(integer(), integer(), numeric()), skip = 2,
        quiet = !verbose
    )
    sparseMatrix(v[[1]], v[[2]], x = v[[3]], dims = dims)
}

#' Import Human Cell Atlas '.mtx.zip' archives
#'
#' @param path character(1) the path to the remote (`http://` or
#'     `https://`) or local `.zip` archive or, when `.data` is
#'     present, the column name (default `"path"`) in which the path
#'     is found.
#'
#' @param ... additional arguments, not supported.
#'
#' @param exdir character(1) directory in which to extract .zip archive.
#'
#' @param overwrite logical(1) overwrite existing files in `exdir`?
#'
#' @param verbose logical(1) report progress using `message()`.
#'
#' @return `SingleCellExperiment()` representing the expression
#'     data. Rows represent features and columns cells; counts are
#'     represented in a sparse matrix.
#'
#' @importFrom utils read.delim unzip
#' @importFrom SingleCellExperiment SingleCellExperiment
#'
#' @examples
#' \dontrun{
#'     import_mtxzip("path/to/mtx/archive.zip")
#' }
#'
#' @keywords internal
import_mtxzip <-
    function(path, ...,
        exdir = tempfile(), overwrite = FALSE, verbose = FALSE)
{
    stopifnot(
        length(list(...)) == 0L,
        .is_scalar_character(exdir),
        .is_scalar_logical(overwrite),
        overwrite || !file.exists(exdir),
        .is_scalar_logical(verbose)
    )

    if (missing(path))
        stop("Provide a valid path from the HCA Matrix service.")

    !verbose || .message("unzip")
    unzip(path, exdir = exdir, overwrite = overwrite, junkpaths = TRUE)
    path <- exdir

    ar <- file.path(exdir, .MTX_ARCHIVE_FILES)
    names(ar) <- sub(".(tsv|mtx).gz", "", .MTX_ARCHIVE_FILES)

    ## rowData / rowRanges
    !verbose || .message("rowData")
    features <- .read_tsv(ar[["features"]], row.names = 1)

    ## colData
    !verbose || .message("colData")
    cells <- .read_tsv(ar[["cells"]], header = TRUE, row.names = "cellkey")
    if (!"barcode" %in% names(cells)) {
        barcodes <- .read_tsv(
            ar[["barcodes"]], blank.lines.skip = FALSE, col.names = "barcode"
        )
        cells <- cbind(cells, barcodes)
    }

    ## assays
    !verbose || .message("assays")
    counts <- .read_mtx(ar[["matrix"]], verbose)

    ## return value
    !verbose || .message("SingleCellExperiment")
    SingleCellExperiment(
        assays = list(counts = counts),
        colData = cells,
        rowData = features
    )
}
