# Releasing sfcqr

The package repository is https://github.com/dingketang/sfcqr. Run these
steps from a source checkout before publishing a versioned release.

## Check the source

Confirm the intended version in `DESCRIPTION` and summarize changes in
`NEWS.md`. Install the development dependencies in R:

```r
install.packages(c("quantreg", "survival", "testthat", "remotes"))
```

From the repository root:

```sh
Rscript -e 'testthat::test_local(".")'
R CMD INSTALL .
Rscript examples/quickstart.R
R CMD build .
R CMD check --no-manual sfcqr_0.1.2.tar.gz
```

Use the version from `DESCRIPTION` in the archive name. Inspect
`sfcqr.Rcheck/00check.log`, then commit and push the reviewed source changes.
Wait for the `R-CMD-check` workflow to pass on Linux, macOS, and Windows at
that same commit. Record the R and dependency versions used for the analysis.

## Publish a checked version

Check that the intended remote tag is available before creating it:

```sh
git ls-remote --tags origin refs/tags/v0.1.2
```

If it does not exist and the checked package version is 0.1.2, create an
annotated tag from the checked commit and push it:

```sh
git tag -a v0.1.2 -m 'sfcqr 0.1.2'
git push origin v0.1.2
```

Create the GitHub release from that tag and attach the checked R source
archive. Verify installation in a fresh R session:

```r
install.packages("remotes")
remotes::install_github("dingketang/sfcqr", ref = "v0.1.2")
library(sfcqr)
```

After the tag is available and installation succeeds, change the README
installation ref from `main` to `v0.1.2`. Record the exact tagged commit in
the manuscript's reproduction archive and link the study scripts there.

## Subsequent releases

Update `DESCRIPTION` and `NEWS.md`, rerun the example and package checks,
and publish a new tag after all platform checks pass. Keep existing release
tags fixed so previously cited versions remain installable.
