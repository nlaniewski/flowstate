# Generate normalized medians/spectra from reference group 'spectral events'

The 'spectral events' returned from
[reference.group.spectral.events](https://nlaniewski.github.io/flowstate/reference/reference.group.spectral.events.md)
are used to generate normalized `[0,1]` medians/spectra for use in
spectral unmixing.

## Usage

``` r
reference.group.medians(
  flowstate.object.reference,
  name.fix = NULL,
  syntactically.valid = FALSE
)
```

## Arguments

- flowstate.object.reference:

  A `flowstate` object as returned from
  [reference.group.spectral.events](https://nlaniewski.github.io/flowstate/reference/reference.group.spectral.events.md).

- name.fix:

  Named character vector – default `NULL`; if defined, vector names
  should match to `sample.id`(s) and their respective values name
  replacements.

  - e.g., `c('Cd8 BV786 (Cells)' = 'CD8 BV786')`.

- syntactically.valid:

  Logical – default `FALSE`; if `TRUE`, spaces, dashes, and dots are
  removed from strings.

## Value

A
[data.table](https://rdatatable.gitlab.io/data.table/reference/data.table.html)
containing normalized reference control medians.
