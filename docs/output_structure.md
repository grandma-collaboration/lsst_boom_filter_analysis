# Output Structure

The pipeline writes one self-contained run directory for each execution.

```text
boom_filter_outputs/
  run_<creation-time>_range_<start>_to_<end>/
    cands.json
    objs.json
    cands_sum.csv
    objs_sum.csv
    plots/
      lc_<n>_object_<objectId>.png
      delta_metrics.csv
    pop_plots/
      delta_distribution.png
      highlighted_fraction_by_band.png
      highlighted_fraction_global.png
      lsst_points_by_band.png
      lsst_points_global.png
    query_summary.md
    pipeline_manifest.json
```

## Files

`cands.json`: full candidate records returned by the BOOM candidate query.

`objs.json`: full object records returned by the BOOM object query.

`cands_sum.csv`: compact table of candidate IDs, object IDs, detection time, magnitude, and SNR.

`objs_sum.csv`: compact table of object-level photometry counts, time baselines, coordinates, magnitudes, SNRs, and per-band counts.

`plots/`: light-curve PNGs. Each plot shows LSST detections, highlighted points that passed the filter, LSST forced photometry, and ZTF crossmatch detections when enabled and available.

`plots/delta_metrics.csv`: per-object relative time range of interest measurements, highlighted-point fractions, and status labels used to build population-level plots.

`pop_plots/`: population-level PNGs for the run.

`pop_plots/delta_distribution.png`: histogram of usable multi-point relative time range of interest values. Single-point cases are excluded from the histogram, and the annotation reports `#used / (#used + #singlepoint)`.

`pop_plots/highlighted_fraction_by_band.png`: six-panel distribution of highlighted LSST points divided by total LSST detection points in each band.

`pop_plots/highlighted_fraction_global.png`: distribution of highlighted LSST points divided by total LSST detection points across all bands.

`pop_plots/lsst_points_by_band.png`: six-panel distribution of the number of LSST detection points in each band.

`pop_plots/lsst_points_global.png`: distribution of the total number of LSST detection points across all bands.

`query_summary.md`: markdown summary of the number of candidates and objects retrieved, plus a table of objects with at least N filter-passing points.

`pipeline_manifest.json`: machine-readable record of the command options and artifact paths for the run.
