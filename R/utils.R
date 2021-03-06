.getAnswer <- function(msg, allowed)
{
    if (interactive()) {
        repeat {
            cat(msg)
            answer <- readLines(n = 1)
            if (answer %in% allowed)
                break
        }
        tolower(answer)
    } else {
        "n"
    }
}

.get_cache <- function() {
    cache <- getOption("HCAMatrixBrowser_cache",
        setCache(directory = tools::R_user_dir("HCAMatrixBrowser", "cache"))
    )
    BiocFileCache::BiocFileCache(cache)
}

.cache_exists <- function(bfc, rname) {
    file.exists(bfcrpath(bfc, rname, exact = TRUE))
}

setCache <-
    function(directory = tools::R_user_dir("HCAMatrixBrowser", "cache"),
        verbose = TRUE, ask = interactive())
{
    stopifnot(is.character(directory), identical(length(directory), 1L),
        !is.na(directory))

    if (!dir.exists(directory)) {
        if (ask) {
            qtxt <- sprintf(
                "Create HCAMatrixBrowser cache at \n    %s? [y/n]: ",
                directory
            )
            answer <- .getAnswer(qtxt, allowed = c("y", "Y", "n", "N"))
            if ("n" == answer)
                stop("'HCAMatrixBrowser_cache' directory not created. Use 'setCache'")
        }
        dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    }
    options("HCAMatrixBrowser_cache" = directory)

    if (verbose)
        message("HCAMatrixBrowser cache directory set to:\n    ", directory)
    invisible(directory)
}

.msg <-
    function(fmt, ..., width=getOption("width"))
{
    txt <- strwrap(sprintf(fmt, ...), width=width, exdent=2)
    paste(txt, collapse="\n")
}

.checkPkgsAvail <- function(pkgnames) {
    vapply(pkgnames, function(pkgname) {
        if (!requireNamespace(pkgname, quietly = TRUE)) {
            func <- as.character(sys.call(1L)[[1L]])
            func <- func[!(func %in% c("::", "HCAMatrixBrowser"))]
            stop("Install the '", pkgname, "' package to use '", func, "'",
                call. = FALSE)
        } else {
            TRUE
        }
    }, logical(1L))
}

.is_scalar_character <-
    function(x, allow.na = FALSE, allow.zchar = FALSE)
{
    is.character(x) && length(x) == 1L && (allow.na || !is.na(x)) &&
        (allow.zchar || nzchar(x))
}

.is_scalar_integer <-
    function(x, allow.na = FALSE)
{
    is.integer(x) && length(x) == 1L && (allow.na || !is.na(x))
}

.is_scalar_logical <-
    function(x, allow.na = FALSE)
{
    is.logical(x) && length(x) == 1L && (allow.na || !is.na(x))
}

.message <-
    function(...)
{
    ## FIXME: use futile.logger?
    message(...)
    TRUE
}

