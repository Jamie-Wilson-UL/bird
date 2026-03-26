## Test environments

* Local macOS (aarch64-apple-darwin20), R 4.5.2
* GitHub Actions:
  * macOS-latest (release)
  * windows-latest (release)
  * ubuntu-latest (devel, release, oldrel-1)

## R CMD check results

Local `devtools::check(cran = TRUE, manual = FALSE)`:

* 0 ERROR
* 1 WARNING
* 2 NOTE

Environment/tooling-related items:

* `checking top-level files ... WARNING`:
  * "A complete check needs the 'checkbashisms' script."
  * This is due to the local system not having `checkbashisms` installed.
* `checking for future file timestamps ... NOTE`:
  * "unable to verify current time"
  * Local environment issue.

## Fortran compiled-code NOTE

`R CMD check` reports Fortran I/O symbols (`__gfortran_st_open/read/write/...`) from
vendored legacy DPpackage Fortran routines used by the nonparametric LDDP engine.

This is expected for the inherited upstream code path. In this package:

* file I/O is redirected to a temporary directory created at runtime (`tempfile("dppackage_")`);
* the working directory is restored after execution;
* temporary files/directories are removed via `on.exit(..., add = TRUE)`.

The temporary-directory handling is implemented in the R wrapper around the Fortran call
(`R/LDDPsurvival.R`, inside `LDDPsurvival.default`).

No persistent user files are written by default during standard usage.
