# Unmix a `flowstate`

A raw/overdetermined `flowstate` is unmixed using a corresponding
`ref.medians`. Data is unmixed using [Ordinary Least
Squares](https://rdrr.io/r/stats/lsfit.html).

## Usage

``` r
flowstate.unmix(flowstate.object.raw, ref.medians, hash.historic = NULL)
```

## Arguments

- flowstate.object.raw:

  A `flowstate` – the return of
  [read.flowstate](https://nlaniewski.github.io/flowstate/reference/read.flowstate.md);
  the `flowstate` must be raw/overdetermined.

- ref.medians:

  A
  [data.table](https://rdatatable.gitlab.io/data.table/reference/data.table.html)
  containing normalized `[0,1]` reference control medians – the return
  of
  [reference.group.medians](https://nlaniewski.github.io/flowstate/reference/reference.group.medians.md).

- hash.historic:

  Character string – default `NULL`; for internal/reproducibility
  purposes, a hash string of the 'unmixing matrix' can be defined.

## Value

A `flowstate` containing unmixed data.
