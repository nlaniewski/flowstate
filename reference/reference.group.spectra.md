# Title

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
reference.group.spectra(
  raw.reference.group,
  name.fix = NULL,
  population.marker = NULL,
  top.N.percent = 10,
  scatter.peak.cut = 0.25,
  top.n = 200,
  top.n.override = NULL,
  output.dirs = "derived",
  plot = F
)
```

## Arguments

- raw.reference.group:

  [file.path](https://rdrr.io/r/base/file.path.html) – one of:

  - a reference group directory containing raw single-color/universal
    negative control .fcs files; all will be processed.

  - a vector of files within that directory; used to process a subset of
    controls.

- name.fix:

  Named character vector – default `NULL`; if defined, vector names
  should match to `sample.id`(s) and their respective values name
  replacements.

  - e.g., `c('cd8alpha BV786 (Cells)' = 'CD8a BV786')`.

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

- top.n:

  Numeric – default `200`; the number of maximally expressing events
  (sorted vector) used to define 'spectral events' – events used to
  define detector medians.

- top.n.override:

  Named numeric vector – default `NULL`; if defined, the supplied
  numeric (value) will override `top.n` on a sample-specific (name)
  basis. See example.

- output.dirs:

  Character vector/an output directory – default `derived`; the default
  will derive output directories based on the location of
  `raw.reference.group` file(s).

- plot:

  Logical – default `FALSE`; plots diagnostic/QC output for evaluating
  function performance.

## Value

A `flowstate` containing spectrally-associated bead/cellular events.
