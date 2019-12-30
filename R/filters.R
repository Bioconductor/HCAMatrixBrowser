# code adapted from GenomicDataCommons
.comparison_filter <- function(sep)
{
    force(sep)
    function(e1, e2) {
        force(e1)
        force(e2)
        list(
            op = jsonlite::unbox(sep),
            field = e1,
            value = e2
        )
    }
}

.logical_filter <- function(sep)
{
    force(sep)
    function(e1, e2) {
        force(e1)
        force(e2)
        list(
            op = jsonlite::unbox(sep),
            value = rbind.data.frame(e1, e2, stringsAsFactors = FALSE)
        )
    }
}

.logicfilters <- c("and", "or", "not")

# create table of operations
.names_table <- data.frame(
    op = c("=", "!=", ">", "<", ">=", "<=", "in", "and", "or", "not"),
    name = c("==", "!=", ">", "<", ">=", "<=", "%in%", "&", "|", "!"),
    stringsAsFactors = FALSE
)

.setOps <- function() {
    x <- stats::setNames(.names_table[["op"]], .names_table[["name"]])
    lapply(x, function(operation) {
        if (operation %in% .logicfilters)
            .logical_filter(operation)
        else
            .comparison_filter(operation)
    })
}

.oplookup <- .setOps()
