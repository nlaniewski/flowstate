# A `flowstate` containing spectrally-associated bead/cellular events.

This function is entirely dependent – as of now – on the following:

- Spectroflo software raw reference controls – tested using SpectroFlo
  3.3.0.

- Spectroflo software naming convention:

  - Directory – './Raw/Reference Group/...'

  - File names – 'MARKER FLUOR (Beads\|Cells).fcs' (literal space
    characters separating)

  - Quality reference controls – well resolved, dominant scatter
    populations

## Usage

``` r
reference.group.spectral.events(
  raw.reference.group,
  transform = FALSE,
  quantiles = NULL,
  population.markers = NULL,
  top.N.percent = 1,
  height.cut.scatter = 0.25,
  n.detector.events = 100,
  detector.override = NULL,
  n.spectral.events = 200,
  n.spectral.events.unstained = 1000,
  n.spectral.events.override = NULL,
  plot = FALSE,
  plot.dir = "$PROJ",
  return.spectral.events = TRUE,
  verbose = FALSE
)
```

## Arguments

- raw.reference.group:

  [file.path](https://rdrr.io/r/base/file.path.html) – one of:

  - a reference group directory containing raw single-color/universal
    negative control .fcs files; all will be processed.

  - a vector of files within that directory; used to process a subset of
    controls.

- transform:

  Logical – default `FALSE`; if `TRUE`, parameters (detectors) having a
  keyword-value pair of `TYPE/Raw_Fluorescence` will be transformed
  using asinh(x/5000).

- quantiles:

  Named [list](https://rdrr.io/r/base/list.html) – default `NULL`; if
  defined, the list must conform to the following:

  - `list('scatter' = probs,'fluor' = probs)`, where
    [probs](https://rdrr.io/r/stats/quantile.html) defines a numeric
    vector of length 2.

    - e.g., `list('scatter' = c(0,0.99),'fluor' = c(0,0.999))`; as
      defined, 'scatter' (forward and side) and 'fluor' (fluorescence)
      associated events will be excluded from processing based on their
      respective `probs`.

- population.markers:

  Named character vector – default `NULL`; `population.markers` must be
  defined for the function to be successful (for cellular controls). The
  named character vector should take the following form:

  - `c(population.name1 = 'sample.name1',population.name2 = 'sample.name2',...)`,
    where `population.name(s)` are lineage cell types and
    `'sample.name(s)'` are lineage markers used to stain those
    respective cell types.

    - e.g.,
      `c(lymphocytes = 'CD3 BV510 (Cells)', monocytes = 'CD14 SB550 (Cells)'`.

  The defined `population.markers` will be used as 'cellular anchors':
  top expressing events (peak detector) in the named sample will be used
  to 'anchor' scatter distributions for proper assignment of the named
  population.

- top.N.percent:

  Numeric – default `1`; defines the number of top expressing events
  (peak detector) used in conjunction with `population.markers` for
  detecting and assigning the named populations.

- height.cut.scatter:

  Numeric – default `0.25`; defines the value at which scatter peak
  heights will be cut for fine tuning selected populations (as defined
  through `population.markers`).

- n.detector.events:

  Numeric – default `100`; the number of maximally expressing events
  (sorted vector) used to auto-detect reference control-specific peak
  detectors.

- detector.override:

  Named character vector – default `NULL`; if defined, the supplied
  detector name (value) will override the auto-detected peak detector on
  a reference control-specific (name) basis. See example.

- n.spectral.events:

  Numeric – default `200`; the number of maximally expressing events
  (sorted vector) used to define 'spectral events'.

- n.spectral.events.unstained:

  Numeric – default `1000`; the number of maximally expressing events
  (sorted vector) used to define '(universal) negative' events.

- n.spectral.events.override:

  Named numeric vector – default `NULL`; if defined, the supplied
  numeric (value) will override `n.spectral.events` on a sample-specific
  (name) basis. See example.

- plot:

  Logical – default `FALSE`; plots diagnostic/QC output for evaluating
  function performance.

- plot.dir:

  Character vector – default `'$PROJ'`; the value associated with the
  keyword '\$PROJ' (experiment name) will be used to construct a plot
  output directory.

- return.spectral.events:

  Logical – default `TRUE`; the returned `flowstate` will be subset to
  include only 'spectral events' – top expressing events used to
  calculate medians to define spectra. If `FALSE`, the returned
  `flowstate` will contain both 'spectral events' and non-'spectral
  events' – useful for diagnostic/QC purposes.

- verbose:

  Logical – default `FALSE`; print function-associated messages to the
  console.

## Value

A `flowstate` containing spectrally-associated bead/cellular events.

## Examples

``` r
if (FALSE) { # \dontrun{
##single-color reference controls
raw.ref.files <- system.file("extdata", package = "flowstate") |>
list.files(full.names = TRUE, pattern = "Beads")

##simulate a SpectroFlo directory structure
raw.ref.dir <- file.path(tempdir(),"Raw","Reference Group")
if(!dir.exists(raw.ref.dir)) dir.create(raw.ref.dir,recursive = TRUE)
file.copy(
  from = raw.ref.files,
  to = file.path(raw.ref.dir,basename(raw.ref.files))
) |> invisible()

##generate spectral events
ref <- reference.group.spectral.events(
  raw.reference.group.directory = raw.ref.dir,
  cluster.types = 'beads',
  plot.select = TRUE,
  verbose = TRUE
)

##a flowstate containing spectral events
class(ref)
ref[['data']]

##peak detectors
ref$data[,.SD,.SDcols = is.factor] |> unique()

##number of spectral events as defined using defaults
ref$data[,.N,by = sample.id]

##saved plot output
list.files(tempdir()) |> grep(pattern = "select.*.pdf",value = TRUE)

##detector.override
##for the sake of example
ref <- reference.group.spectral.events(
  raw.reference.group.directory = raw.ref.dir,
  cluster.types = 'beads',
  detector.override = c("Unstained (Beads)" = "UV7")
)
ref$data[,.SD,.SDcols = is.factor] |> unique()

##n.spectral.events.override
##for the sake of example
ref <- reference.group.spectral.events(
  raw.reference.group.directory = raw.ref.dir,
  cluster.types = 'beads',
  n.spectral.events.override = c("TNFa PE (Beads)" = 500)
)

##number of spectral events as defined using n.spectral.events.override
ref$data[,.N,by = sample.id]
} # }
```
