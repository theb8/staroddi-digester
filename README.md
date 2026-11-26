# staroddi-digester

R functions to process continuous (daily) data from a Star-Oddi oceanographic data logger
mounted on a fishing rig (dredge or trawl).

The logger records depth, temperature, salinity, etc. all day long. This project looks at the
depth + timestamp series and segments it into **sets** (tows): the gear goes down to the bottom,
fishes, and comes back up to the ship. For each set it then computes summary values
(start/end time, time on bottom, median temperature, depth and salinity).

> Status: experimental / work in progress.

## Project layout

```
R/func.R     all the functions
R/run.R      example script: read one daily file and run the full pipeline
data/        logger files (not tracked by git), e.g. data/iml202619/2026-06-12.CSV
calib/       reference data for depth validation (not tracked by git)
```

## Input data

Star-Oddi CSV export, one file per day:

- separator `;`, encoding `latin1`
- decimal comma (`12,5`), converted with `charcomma2numdecimal()`
- `Date.Time` formatted as `dd.mm.yyyy HH:MM:SS`
- columns used: `Date.Time`, `Depth.m.`, `Temp..C.`, `Salinity.psu.`
  (also accepted by `average_x()`: `Conduct.mS.cm.`, `Velocity.m.sec.`)

## Usage

Open `staroddi-digester.Rproj` (or set the working directory to the project root), put a file in
`data/`, adjust `fname` in `R/run.R`, then:

```r
source("R/run.R")
```

This produces:

| object         | content |
|----------------|---------|
| `main_table`   | the full file, with `tow_id` (set number, `NA` between sets) and `depth_corr` (depth relative to the surface) |
| `summary_df`   | one row per set: `datetime_start`, `datetime_end`, `dur`, `bottom_start`, `bottom_end`, `bottom_dur`, median `Temp..C.`, `depth_corr`, `Salinity.psu.` |
| `details_df`   | list with one data frame per set |
| `drift_report` | surface level found between each pair of sets (`local_zeros`, `local_spread`, `n_skipped`) |

## Pipeline

1. **Read & parse** (`run.R`): numeric columns, `POSIXct` datetime.
2. **Clean** – `clean_logger_data()`: drops rows with a missing depth or datetime (warning + row
   numbers in `attr(x, "dropped_rows")`). Later steps refuse `NA` values.
3. **Zero the depth** – `find_surface()`: the surface level is the median of the shallowest points
   (within 3 m of the 5th percentile). `depth_corr = depth - surface`.
4. **Find the sets** – `find_tow()`, hysteresis approach:
   - a candidate set is a run of points deeper than `exit` (1.5 m) that reaches below `enter` (7 m);
   - candidates less than `min_gap` (2 min) apart are merged;
   - sets shorter than `min_dur` (2 min) are dropped;
   - a file with no set returns an all-`NA` `tow_id` with a warning.
5. **Drift check** – `drift_check()`: compares the surface level of each break between sets with
   the global surface level, and warns if one is off by more than 0.8 m (sensor drift). The 30 s
   after and 4 min before each set are ignored; breaks too short to use are skipped (`NA`).
6. **Summarise** – `file_summary()`, `file_details()`, `average_x()`: per-set values are medians
   over the **on-bottom** part of the set, selected by `bottom_mask()` (points within `band` = 2 m
   of the set's bottom depth, its 90th percentile).

## Main parameters

| function       | parameter | default | meaning |
|----------------|-----------|---------|---------|
| `find_tow`     | `enter`   | 7 m     | a set must go deeper than this |
| `find_tow`     | `exit`    | 1.5 m   | a set ends when the gear is shallower than this |
| `find_tow`     | `min_gap` | 2 min   | sets closer than this are merged |
| `find_tow`     | `min_dur` | 2 min   | shorter sets are dropped |
| `bottom_mask`  | `band`    | 2 m     | depth band below the bottom depth counted as "on bottom" |
| `drift_check`  | `DRIFT_TH`| 0.8 m   | max allowed difference of a local surface level (hard coded) |

## Known issues / to do

- **Time zone**: timestamps are parsed in the computer's local time zone, while the summary table
  uses `America/Toronto`. The logger clock has no daylight saving time, so parsing should use a
  fixed zone (UTC or a fixed offset) matching how the logger was programmed.
- **`difftime` units**: the gap between sets in `find_tow()` and `dur` in `file_summary()` are
  computed without explicit units, so R may pick seconds/minutes/hours.
- **Depth validation**: comparing the logger depth with the reference `true_depth` from the ANDES
  sensor data (`calib/`) is not written yet (`depth_validation()`), and the end of `run.R`
  references an undefined `fname_depthref`.

## Author

atheberge
