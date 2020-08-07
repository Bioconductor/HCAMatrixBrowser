## Changes in version 1.0.0

### New features

* `HCAMatrixBrowser` finally on Bioconductor!
* `HCAMatrixBrowser` uses OpenAPI Specification version 2 and `rapiclient`
to provide R API representations.
* `loadHCAMatrix` provides users with matrix data given a set of 'bundle_fqids'
* Filtering on the main API object is supported see `HCAMatrix` for details
* Representations in all formats is supported this includes (.csv, .mtx, and
.loom)
* Caching implemented using `BiocFileCache`
* MTX format support provided by Martin @mtmorgan
* Support for LOOM files provided by `LoomExperiment`
* CSV format files given as a `tibble` list

### Bug fixes and minor improvements

* Updated vignettes to include API changes
* Allow for singleton `bundle_fqid` queries for v0 endpoint (@dvantwisk, #2)
