# Model and parameter recovery for choice models

A worked example of model and parameter recovery for choice models, on intertemporal choice. Six models of money-earlier-or-later decisions, three based on delay discounting and three on attribute heuristics, are compared by cross-validation on the public data of Ericson et al. (2015), in the adjusted form of Wulff and van den Bos (2018). The same selection procedure is then run on data simulated from every model, which gives the confusion matrix (how often each generating model is recovered) and its inversion (how far a selection should be believed), and each generating model is refit to its own synthetic participants to see whether its parameters come back. On 25 choices per participant the procedure recovers the generating model for 56% of exponential-discounting agents and 84% of trade-off agents, but for only 10% to 26% of the other four, and when it selects the empirical winner, exponential discounting, the chance that this model generated the data is 25%, barely above a uniform guess. The report is `report.qmd`; it reads only frozen results and renders in well under a minute. The rendered report is published at <https://01a08ba1-9c90-e3b1-f8f2-23f1133a22c4.share.connect.posit.cloud>.

## Provenance

- Data: Marzilli Ericson, K. M., White, J. M., Laibson, D., & Cohen, J. D. (2015). Money earlier or later? Simple heuristics explain intertemporal choices better than delay discounting does. *Psychological Science*, 26(6), 826–833. [doi:10.1177/0956797615572232](https://doi.org/10.1177/0956797615572232). The export in `data/raw/choices.csv` is the `choices.csv` of OSF project [9uxve](https://osf.io/9uxve/), licensed CC BY 4.0.
- Model definitions: Wulff, D. U., & van den Bos, W. (2018). Modeling choices in delay discounting. *Psychological Science*, 29(11), 1890–1894. [doi:10.1177/0956797616664342](https://doi.org/10.1177/0956797616664342).

## Layout

| Path | Contents |
|---|---|
| `report.qmd` | The report. It reads only the frozen summary and the raw export, so it renders in well under a minute; the rendered file is not tracked. |
| `R/models.R` | The six models: choice probability, box constraints, start boxes, parameter names. |
| `R/fit.R` | Sample rules, fitting, cross-validation with common splits, simulation, recovery, and the keyed recovery metrics. |
| `R/prepare_data.R` | Builds `data/prepared/choices.rds` from the raw export. |
| `R/run_analysis.R` | Computes every result in resumable stages and writes `output/summary.rds`. |
| `tests/test-recovery.R` | Machinery tests: model sanity, keyed joins, recovery on synthetic truth. |
| `data/raw/` | The untouched OSF export. |
| `data/prepared/` | The analysis sample, one row per answered trial. |
| `output/` | `summary.rds`, the frozen artifact the report reads; `raw/` holds per-repetition stage output and is not tracked. |

## Reproduce

Requirements: R 4.5 with dplyr, tidyr, purrr, readr, forcats, ggplot2, patchwork, furrr, future, withr, digest, knitr, rmarkdown, and testthat, plus Quarto for the render. The package versions used are listed in the report's session information.

1. Prepare the sample from the raw export.

   ```bash
   Rscript R/prepare_data.R
   ```

2. Run the machinery tests.

   ```bash
   Rscript -e "testthat::test_dir('tests')"
   ```

3. Compute the results. This takes several hours on seven cores at the default 100 cross-validation repetitions; the run resumes stage by stage if interrupted, and `CMR_N_SUBJECTS=20 CMR_N_REPS=3` gives a one-minute smoke run.

   ```bash
   Rscript R/run_analysis.R
   ```

4. Render the report from the frozen summary.

   ```bash
   quarto render report.qmd
   ```

## License

MIT, see [LICENSE](LICENSE). The data carry their own license, CC BY 4.0, from the source above.
