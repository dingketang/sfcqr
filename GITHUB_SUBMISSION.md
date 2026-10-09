# Publishing sfcqr on GitHub

This guide publishes the package source as a GitHub repository and a versioned release. The delivered files are local; a remote repository has not been created or uploaded. Submission to CRAN is a separate process with its own policies and checks.

## 1. Open the package source directory

Extract the delivered GitHub source archive into an independent directory. Open the package directory containing `DESCRIPTION`, `R`, `man`, `tests`, and `.github`; it should be the repository root. Keep the hidden `.github` directory and `.gitignore` file when copying the source.

```sh
cd /path/to/sfcqr
```

Run the following Git commands from this directory. Upload this package source rather than the manuscript project or an installed R library directory.

## 2. Complete the release metadata and run the checks

Confirm the release metadata before publishing:

- `DESCRIPTION` identifies Dingke Tang (`dtang@uottawa.ca`) as the author and maintainer; `LICENSE` and `LICENSE.md` identify Dingke Tang as the 2026 copyright holder. The delivered source uses the MIT license.
- Replace `USERNAME/sfcqr` in the README with your actual repository path.

Once your GitHub username is known, add these fields to `DESCRIPTION`:

```text
URL: https://github.com/USERNAME/sfcqr
BugReports: https://github.com/USERNAME/sfcqr/issues
```

Install the runtime and test dependencies in R:

```r
install.packages(c("quantreg", "survival", "testthat"))
```

From the independent copy, install the package, run the README example script, and build and check the source package:

```sh
Rscript -e 'testthat::test_local(".")'
R CMD INSTALL .
Rscript examples/quickstart.R
cd ..
R CMD build sfcqr
R CMD check --no-manual sfcqr_0.1.2.tar.gz
cd sfcqr
```

The quickstart script uses the manuscript's main setting described in the README. Review `../sfcqr.Rcheck/00check.log` and resolve any errors or warnings before pushing. Rerun the checks after changing the release metadata.

## 3. Create an empty GitHub repository

Sign in to GitHub and create a repository named `sfcqr`, with your preferred public or private visibility. Leave the options to generate a README, license, or `.gitignore` unchecked, because these files are already included locally. These steps follow [GitHub's guide to adding locally hosted code](https://docs.github.com/en/migrations/importing-source-code/using-the-command-line-to-import-source-code/adding-locally-hosted-code-to-github).

## 4. Initialize Git and make the first commit

Run the following commands in the package source directory:

```sh
git init -b main
git status
git add .
git diff --cached --stat
git commit -m 'Initial release of sfcqr 0.1.2'
```

Confirm that the staged file list contains only package files. If Git reports that your commit identity is missing, set it for this repository and retry the commit:

```sh
git config user.name 'YOUR NAME'
git config user.email 'YOUR EMAIL'
git commit -m 'Initial release of sfcqr 0.1.2'
```

## 5. Connect the remote repository and push

For an HTTPS remote, replace `USERNAME` with your GitHub username:

```sh
git remote add origin https://github.com/USERNAME/sfcqr.git
git remote -v
git push -u origin main
```

Authenticate using a personal access token, Git Credential Manager, or GitHub CLI. GitHub account passwords are not accepted for Git operations. If you use GitHub CLI, run `gh auth login` first. If SSH authentication is already configured, use `git@github.com:USERNAME/sfcqr.git` as the remote URL. See [GitHub's authentication documentation](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/about-authentication-to-github).

## 6. Verify GitHub Actions and create a release

After pushing, open the repository's **Actions** tab and inspect the `R-CMD-check` workflow. The included workflow installs dependencies and checks the package with the current R release on Linux, macOS, and Windows. When all three platform checks pass, tag the release:

```sh
git tag v0.1.2
git push origin v0.1.2
```

In GitHub's **Releases** section, create a release from `v0.1.2`. Add a short description of the package and attach the `sfcqr_0.1.2.tar.gz` file produced by `R CMD build`. This archive is an R source package; users need R and the package dependencies to install it.

Other users can then install the tagged version in R:

```r
install.packages("remotes")
remotes::install_github("USERNAME/sfcqr", ref = "v0.1.2")
library(sfcqr)
```

## Subsequent updates

After editing the code, documentation, or tests, update the version in `DESCRIPTION` and the release notes in `NEWS.md`. Run the README example and `R CMD check` again before committing and pushing:

```sh
git add .
git commit -m 'Describe the change'
git push
```

Use a new tag for each published version, such as `v0.1.3` for the next release. The GitHub Actions configuration follows the [official r-lib/actions R package check examples](https://github.com/r-lib/actions/tree/v2/examples).
