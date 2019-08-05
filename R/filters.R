# code from GenomicDataCommons
.binary_op <- function(sep)
{
    force(sep)
    function(e1, e2) {
        force(e1)
        force(e2)
        list(
            op = jsonlite::unbox(sep),
            content = list(value = e2, field = e1)
        )
    }
}

.combine_op <- function(sep)
{
    force(sep)
    function(e1, e2) {
        force(e1)
        list(
            op = jsonlite::unbox(sep),
            content =
                if (!identical(sep, "not")) {
                    force(e2)
                    list(e1, e2)
                } else {
                    list(e1)
                }
        )
    }
}

.combined_ops <- c("or", "and", "not")

# create table of operations
.names_table <- data.frame(
    op = c("=", "!=", "<", ">", "and", "or", "<=", ">=", "in", "not"),
    name = c("==", "!=", "<", ">", "&", "|", "<=", ">=", "%in%", "!"),
    stringsAsFactors = FALSE
)

.setOps <- function() {
    x <- setNames(.names_table[["op"]], .names_table[["name"]])
    lapply(x, function(operation) {
        if (operation %in% .combined_ops)
            .combine_op(operation)
        else
            .binary_op(operation)
    })
}

.oplookup <- .setOps()
