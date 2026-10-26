# Output Structure

The pipeline writes one self-contained run directory for each execution.

```text
boom_filter_outputs/
  run_<creation-time>_range_<start>_to_<end>/
    cands.json
    objs.json
    cands_sum.csv
    objs_sum.csv
    band_rates.json
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

When an excluded time range is set, every object with any point (LSST detection or forced photometry, in any band) inside it is dropped from the run, along with its candidates. All the files below then only describe the objects that were kept.

## Files

`cands.json`: full candidate records returned by the BOOM candidate query.

`objs.json`: full object records returned by the BOOM object query.

`cands_sum.csv`: compact table of candidate IDs, object IDs, detection time, magnitude, and SNR.

`objs_sum.csv`: compact table of object-level photometry counts, time baselines, coordinates, magnitudes, SNRs, and per-band counts.

`band_rates.json`: for each object and band, BOOM's rising and fading rates (mag/day) as of the object's most recent LSST detection in that band, whether or not that detection passed the filter. Taken from `properties.photstats` of the corresponding `LSST_alerts` document. The most recent detection is the latest one BOOM has at query time, which can be later than the end of the queried range.

`plots/`: light-curve PNGs, one panel per band. Each plot shows LSST detections, highlighted points that passed the filter, LSST forced photometry, and ZTF crossmatch detections when enabled and available. The legend of each panel gives that band's rising/fading rate from `band_rates.json`, and the excluded time range, if any, is shaded red where it overlaps the plotted time range. The title gives the object's T0, its relative time range of interest, its number of filter-passing LSST points, and its RA and Dec.

`plots/delta_metrics.csv`: per-object relative time range of interest measurements, highlighted-point fractions, and status labels used to build population-level plots.

`pop_plots/`: population-level PNGs for the run.

`pop_plots/delta_distribution.png`: histogram of the relative time range of interest of objects with at least two filter-passing points and a baseline under 20 days. The annotations report the number of objects used out of all objects, and how many were left out for a baseline over 20 days, for a single filter-passing point, or for both.

`pop_plots/highlighted_fraction_by_band.png`: six-panel distribution of highlighted LSST points divided by total LSST detection points in each band.

`pop_plots/highlighted_fraction_global.png`: distribution of highlighted LSST points divided by total LSST detection points across all bands.

`pop_plots/lsst_points_by_band.png`: six-panel distribution of the number of LSST detection points in each band.

`pop_plots/lsst_points_global.png`: distribution of the total number of LSST detection points across all bands.

`query_summary.md`: markdown summary of the number of candidates and objects retrieved (after dropping objects with points in the excluded range), plus a table of objects with at least N filter-passing points, and the objects in the first bin of the relative time range of interest (baseline under 20 days and relative time range of interest under 0.05).

`pipeline_manifest.json`: machine-readable record of the parameters and artifact paths for the run. When `COMPUTE_FILTER_PASS_RATES` is enabled, `filter_result.filter_pass_rates` holds the number of filter-passing alerts and objects, the total number of LSST alerts and objects in the queried range, the resulting pass fractions, and the count queries used for the totals.
