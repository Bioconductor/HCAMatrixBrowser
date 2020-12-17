
# HCAMatrixBrowser <a href='https://bioconductor.org/packages/HCAMatrixBrowser/'><img src='https://raw.githubusercontent.com/Bioconductor/BiocStickers/master/HCAMatrixBrowser/HCAMatrixBrowser.png' align="right" height="139" /></a>

## Overview

`HCAMatrixBrowser` gives users access to the matrix download service for
the [Human Cell Atlas](https://data.humancellatlas.org/) project.

It provides a main function for downloading data matrices from either a
set of `bundle_fqids` or a data filters.

  - `loadHCAMatrix` allows either a vector of `bundle_fqids` or an added
    filter operation to the HCA API object

## API object

``` r
library(HCAMatrixBrowser)
#> Loading required package: AnVIL
#> Loading required package: dplyr
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
#> 
#> Attaching package: 'HCAMatrixBrowser'
#> The following object is masked from 'package:dplyr':
#> 
#>     filter
#> The following object is masked from 'package:stats':
#> 
#>     filter

hca <- HCAMatrix()
```

# Quick Start

## with `bundle_fqids`

``` r
my_fqids <-
    c("ffd3bc7b-8f3b-4f97-aa2a-78f9bac93775.2019-05-14T122736.345000Z",
    "f69b288c-fabc-4ac8-b50c-7abcae3731bc.2019-05-14T120110.781000Z")

loadHCAMatrix(hca, bundle_fqids = my_fqids)
#> class: LoomExperiment 
#> dim: 58347 2 
#> metadata(0):
#> assays(1): matrix
#> rownames: NULL
#> rowData names(9): Accession Gene ... genus_species isgene
#> colnames(2): 3c2180aa-0aa4-411f-98dc-73ef87b447ed
#>   1cfe9423-21d1-4281-9f9d-3aaa07b8e1e8
#> colData names(38): CellID barcode ...
#>   specimen_from_organism.provenance.document_id total_umis
#> rowGraphs(0): NULL
#> colGraphs(0): NULL
```

## with filter operations

``` r
hca1 <- filter(hca, bundle_uuid == "ffd3bc7b-8f3b-4f97-aa2a-78f9bac93775")
## filter print out
filters(hca1)
#> $op
#> [x] "="
#> 
#> $field
#> [x] "bundle_uuid"
#> 
#> $value
#> [x] "ffd3bc7b-8f3b-4f97-aa2a-78f9bac93775"

loadHCAMatrix(hca1)
#> class: LoomExperiment 
#> dim: 58347 1 
#> metadata(0):
#> assays(1): matrix
#> rownames: NULL
#> rowData names(9): Accession Gene ... genus_species isgene
#> colnames(1): 3c2180aa-0aa4-411f-98dc-73ef87b447ed
#> colData names(44): CellID analysis_protocol.protocol_core.protocol_id
#>   ... specimen_from_organism.provenance.document_id total_umis
#> rowGraphs(0): NULL
#> colGraphs(0): NULL
```
