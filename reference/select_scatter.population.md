# Add scatter-specific populations to `[['data']]`

Effective processing of single-color controls first involves defining
dominant scatter-specific populations that are free of debris/aggregates
and contain max-expressing events. To efficiently select these
populations in cellular controls, a defined set of cluster of
differentiation (CD)/lineage surface markers are used to restrict
peak-finding/cutting methods to the scatter profile of specific immune
cell types. Once these profiles/populations are defined, they are
applied to all of `[['data']]`.

## Usage

``` r
select_scatter.population(
  flowstate,
  population.marker,
  top.N.percent = 10,
  scatter.peak.cut = 0.25,
  plot = F,
  plot.output.dir = tempdir()
)
```

## Arguments

- flowstate:

  A flowstate object as returned from
  [read.flowstate](https://nlaniewski.github.io/flowstate/reference/read.flowstate.md).

- population.marker:

  Named character vector – default `NULL`; `population.marker` must be
  defined for the function to be successful (for cellular controls). The
  named character vector should take the following form:

  - `c(population.name1 = 'sample.name1',population.name2 = 'sample.name2',...)`,
    where `population.name(s)` are CD/lineage cell types and
    `'sample.name(s)'` are CD/lineage markers used to stain those
    respective cell types.

    - e.g.,
      `c(lymphocytes = 'CD3 BV510 (Cells)', monocytes = 'CD14 SB550 (Cells)'`.

- top.N.percent:

  Numeric – default `10`; defines the percentage of top expressing
  events (peak detector) used in conjunction with `population.marker`
  for detecting and assigning the named populations.

- scatter.peak.cut:

  Numeric – default `0.25`; defines the value at which scatter peak
  heights will be cut for fine tuning the scatter profile of selected
  populations (as defined through `population.marker`).

- plot:

  Logical – default `FALSE`; plots diagnostic/QC output for evaluating
  function performance.

- plot.output.dir:

  A character vector of a full path name – default
  [tempdir()](https://rdrr.io/r/base/tempfile.html); any/all generated
  plots will be saved to this directory.

## Value

UPDATES BY REFERENCE and invisibly returns `flowstate`; a new column
(factored) named `population` is added to `[['data']]`.
