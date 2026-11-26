# ======================================= -
# func.R
# functions to digest Star-Oddi logger data (continuous daily files)
# mounted on a fishing rig (dredge or trawl) :
# segment the depth/time series into tows (sets) and summarise each tow
#
# typical order (see R/run.R) :
#   charcomma2numdecimal() -> clean_logger_data() -> find_surface() -> find_tow()
#   -> drift_check() -> file_summary() / file_details()
# ======================================= -


# ----------------------- -
# charcomma2numdecimal()
# convert numbers written with a decimal comma ("12,5") to numeric (12.5)
# only the first comma is replaced; unreadable values become NA (with a warning from as.numeric)
# will not fail if x is already numeric
#
# x : character (or numeric) vector
# returns : numeric vector, same length as x
# ----------------------- -
charcomma2numdecimal <- function(x){
  x_dot <- sub(pattern = ",", replacement = ".", x = x)
  x_num <- as.numeric(x_dot)
}


# ----------------------- -
# clean_logger_data()
# drop rows with a missing depth or datetime (unreadable values become NA when parsing)
# must be called after parsing (numeric depth, POSIXct datetime) and before find_surface()/find_tow()
#
# input_table : data frame of the file (Date.Time as POSIXct, Depth.m. as numeric)
# returns : the same data frame without the bad rows
#           the original row numbers of the dropped rows are in attr(x, "dropped_rows")
# ----------------------- -
clean_logger_data <- function(input_table){

  bad <- is.na(input_table$Depth.m.) | is.na(input_table$Date.Time)
  dropped_rows <- which(bad)

  if(length(dropped_rows) > 0){
    warning(sprintf("%d rows dropped (missing depth/datetime), see attr(x, \"dropped_rows\")",
                    length(dropped_rows)))
  }

  clean_table <- input_table[!bad, ]
  attr(clean_table, "dropped_rows") <- dropped_rows
  return(clean_table)
}


# ----------------------- -
# find_surface()
# estimate the depth reading of the sensor at the surface (its "zero")
# = median of the shallowest points : those within 3 m of the 5th percentile of depth
# used to correct the depth : depth_corr = depth - find_surface(depth)
#
# depth : raw depth values of the whole file (m), no NA
# returns : single numeric value, the surface level (m)
# ----------------------- -
find_surface <- function(depth){

  if(anyNA(depth)){
    stop("depth contains NA, run clean_logger_data() first")
  }

  # find surface
  surface_anchor <- quantile(depth, probs = 0.05) #5th percentile to
  surface_level <- 3
  surface_mask <- depth < (surface_anchor + surface_level)
  surface_points <- depth[surface_mask]
  zero <- median(surface_points)
  return(zero)
}

# --------------------- -
# find_tow()
# find the tow events (sets) of the instrument
# "hysteresis" approach to avoid ambiguous depth near the surface :
#   1. a candidate tow is a run of points deeper than `exit` that reaches deeper than `enter`
#   2. candidates separated by less than `min_gap` are merged into one tow
#   3. tows shorter than `min_dur` are dropped
# a file without any tow returns an all-NA vector (with a warning)
#
# depth : corrected depth (m, relative to the surface), no NA
# datetime : POSIXct timestamps, same length as depth, no NA
# enter (m) : depth a tow must reach
# exit (m) : depth threshold to end a tow
# min_gap (min) : min time gap between two tows; closer tows are merged
# min_dur (min) : minimum duration of a tow
# returns : integer vector, same length as depth : tow number (1, 2, ...) or NA between tows
# --------------------- -
find_tow <- function(depth, datetime, 
                     enter=7, exit=1.5,
                     min_gap=2, min_dur=2){
  
  n <- length(depth)
  
  if(length(datetime) != n){
    stop("depth and datetime must be the same length")
  }

  # NA would silently split/truncate tows
  if(anyNA(depth) || anyNA(datetime)){
    stop("depth or datetime contains NA, run clean_logger_data() first")
  }
  
  # identifity points just below the surface of the water (default is 1.5 meters)
  low <- depth >= exit
  
  # group points together when they belong to a single event/tow
  # finmd their respective start and end indices
  r <- rle(low)                           # r$values = TRUE/FALSE, r$lengths = run sizes
  ends   <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1
  
  # Max depth of each tow/event 
  run_id  <- rep(seq_along(r$lengths), r$lengths)
  run_max <- tapply(depth, run_id, max)
  
  # keep the event only if it went deeper than `enter` (7 meters default)
  keep <- r$values & run_max > enter
  tow_indices <- data.frame(start = starts[keep], end = ends[keep])

  # no tow in this file (e.g. bad weather, day at port) : return an all-NA map
  tow_id <- rep(NA_integer_, n)
  if(nrow(tow_indices) == 0){
    warning("no tow found in this file")
    return(tow_id)
  }

  # --- -
  # check 1 : time gap checkup
  # tow events were identified 
  # now let's take a look at the time difference identify false tow events
  # --- -
  
  # compute time gap for all tow events
  gap <- as.numeric(
    difftime(datetime[tow_indices$start[-1]], datetime[tow_indices$end[-nrow(tow_indices)]], units = "min")
    )
  
  # create groups based on time difference
  time_group <- cumsum(c(TRUE, gap >= min_gap))
  
  # within each time-group, find the min/max indice from the first identification 
  tow_indices  <- data.frame(start = tapply(tow_indices$start, time_group, min),
                              end   = tapply(tow_indices$end,   time_group, max))
  
  # --- -
  # check 2 : short tows
  # now let's drop short OR shallow tow events
  # --- -
   
  # Drop: too short or too shallow
  dur  <- as.numeric(
    difftime(datetime[tow_indices$end], datetime[tow_indices$start], units = "mins")
    )
  
  tow_indices   <- tow_indices[dur >= min_dur,]
  
  n_tows <- dim(tow_indices)[1]

  if(n_tows == 0){
    warning("no tow found in this file (all candidate tows shorter than min_dur)")
    return(tow_id)
  }

  for(i in seq_len(n_tows)){
    
    tow_id[tow_indices$start[i]:tow_indices$end[i]] <- i
    
  }
  
  # return a map of the tow id 
  # same length as depth/datetime
  return(tow_id)
  
}


# ----------------------- -
# drift_check()
# Try to find if the depth is drifting
# looking at all "surface events" individually and compare the local zero found
# Will inform the user if a drift has been found in the file
#  
# Note : If a signicant drift is found that would make things more complex to identify when the rig is actually fishing
#
# Note : a "surface event" is the moment between two fishing tow
# when the rig is stand by (at the very surface of the water or on the deck).
#
# the 30 s after and 4 min before each tow are ignored (gear handled at the surface);
# surface events too short to use (< 10 points left) are skipped (NA) with a warning
#
# depth_raw : raw depth values (m, NOT corrected), no NA
# datetime : POSIXct timestamps, same length as depth_raw
# tow_id : output of find_tow(), same length as depth_raw
# returns : list with
#   indices : data frame (start, end) of the row indices of each surface event
#   local_zeros : median depth of each surface event (NA if skipped)
#   local_spread : max - min depth of each surface event (NA if skipped)
#   n_skipped : number of surface events skipped
# ----------------------- -
drift_check <- function(depth_raw, datetime, tow_id){
  
  # hard coded thresholds
  TIME_TH_LO <-  30 # will discard 30 seconds after each previous tow (unit sec)
  TIME_TH_UP <- 4*60 # will discard 4 minutes just before the next tow; since the gear is being moved around at the surface (unit sec)
  DRIFT_TH <- 0.8 # hard threshold to determine if sensor is drifting (unit meters)
  
  #find points between tows
  surface_mask <- is.na(tow_id)
  
  # regroup in event and get the indice
  r <- rle(surface_mask)
  ends   <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1
  keep <- r$values == T
  surface_indices <- data.frame(start = starts[keep], end = ends[keep])
  
  # surface events found
  n_groups <- dim(surface_indices)[1]
  
  # also compute the global zero to compare values
  global_zero <- find_surface(depth_raw)
  
  #preallocate result
  file_zeros <- rep(NA_real_, n_groups)
  file_spread  <- rep(NA_real_, n_groups)

  n <- length(depth_raw)
  MIN_PTS <- 10 # min number of points for a surface event to be used

  for(i in seq_len(n_groups)){

    indice_lo <- surface_indices$start[i]
    indice_up <- surface_indices$end[i]

    start_time <- datetime[indice_lo]
    end_time <- datetime[indice_up]

    # only trim the sides that border a tow
    # (the first/last surface event of the file has no tow before/after it)
    start_time_th <- if(indice_lo > 1) start_time + TIME_TH_LO else start_time
    end_time_th <- if(indice_up < n) end_time - TIME_TH_UP else end_time

    #mask using timestamp (empty if the event is shorter than the trims)
    keep_time <- (datetime >= start_time_th) & (datetime <= end_time_th)
    keep_time[is.na(keep_time)] <- FALSE
    n_kept <- sum(keep_time) #number of points for this "surface event"

    # too short : skip this surface event (left as NA)
    if(n_kept < MIN_PTS){
      next
    }

    # the depth values that are kept
    depth_window <- depth_raw[keep_time]

    # store computed values (zero found for this surface event)
    file_zeros[i] <- median(depth_window, na.rm = TRUE)
    file_spread[i] <- diff(range(depth_window, na.rm = TRUE))

  }

  # report skipped surface events
  skipped <- is.na(file_zeros)
  n_skipped <- sum(skipped)
  if(n_groups == 0 || n_skipped == n_groups){
    warning("no usable surface event found\n  drift could not be checked")
  } else if(n_skipped > 0){
    warning(sprintf("%d of %d surface events too short (< %d points after trimming) and skipped in drift check",
                    n_skipped, n_groups, MIN_PTS))
  }

  # check if a surface event is outside our depth threshold (skipped events are ignored)
  zero_ok <- abs(file_zeros - global_zero) <= DRIFT_TH

  if(!all(zero_ok, na.rm = TRUE)){
    warning_msg <- sprintf("some surface events are outside threshold of %s meters\n  please verify if the instrument is drifting", DRIFT_TH)
    warning(warning_msg)
  }

  list2return <- list(indices=surface_indices,
                      local_zeros=file_zeros,
                      local_spread=file_spread,
                      n_skipped=n_skipped)
  return(list2return)
  
}

# ----------------------- -
# average_x()
# compute the average (median) value of variable x for one tow,
# using only the on-bottom points of the tow (see bottom_mask())
#
# tow : tow number (sequential value starting at 1)
# x : name of the target column, one of : "Date.Time", "Temp..C.", "Depth.m.", "Salinity.psu.",
#     "Conduct.mS.cm.", "Velocity.m.sec.", "depth_corr"
# main_table : data frame of the file (formatted with date and proper numeric),
#              with the tow_id and depth_corr columns
# band : depth band (m) below the bottom depth used to select on-bottom points (see bottom_mask)
# returns : median of x over the on-bottom points of the tow
# ----------------------- -
average_x <- function(tow, x, main_table, band=2){
  
  accepted_vars <- c("Date.Time",
                     "Temp..C.",
                     "Depth.m.",
                     "Salinity.psu.",
                     "Conduct.mS.cm.",
                     "Velocity.m.sec.",
                     "depth_corr"
                     )
  
  # check arguments and flag error
  cols <- colnames(main_table)
  if(!(x %in% accepted_vars)){
    error_msg <- sprintf("variable name is not valid")
    stop(error_msg)
  }
  unique_tows <- unique(main_table$tow_id)
  unique_tows <- unique_tows[!is.na(unique_tows)]
  if(!any(tow == unique_tows)){
    error_msg <- sprintf("selected tow is not available for this file")
    stop(error_msg)
  }
  
  # columns indice for x
  col_indice <- which(cols ==x)
  
  # create a subset for the current tow
  mask <- main_table$tow_id == tow
  mask[is.na(mask)] <- F #set NA to False
  subset_table <- main_table[mask,]

  # keep only the on-bottom part of the tow
  # the goal is to drop values when the instrument is going up/down or moving around at subsurface
  keep <- bottom_mask(subset_table$depth_corr, band = band)

  x_values <- subset_table[keep, col_indice]
  avg <- median(x_values, na.rm = TRUE)
  return(avg)

}


# ----------------------- -
# bottom_mask()
# select the on-bottom points of a single tow, based on depth
# (instead of trimming a fixed time at the start/end of the tow)
#
# depth : depth values of ONE tow (corrected depth, m)
# band : points deeper than (bottom depth - band) are kept (m)
# min_pts : if fewer points than this are kept, fall back to the whole tow (with a warning)
# returns : logical vector, same length as depth
# ----------------------- -
bottom_mask <- function(depth, band=2, min_pts=5){

  # robust "bottom depth" : 90th percentile rather than max (avoids spikes)
  bottom <- quantile(depth, probs = 0.9, na.rm = TRUE)
  mask <- depth >= (bottom - band)
  mask[is.na(mask)] <- FALSE

  if(sum(mask) < min_pts){
    warning("few on-bottom points for this tow; using the whole tow")
    mask <- !is.na(depth)
  }

  return(mask)
}




# ------------------ -
# file_summary()
# build the summary table of a file : one row per tow
#
# main_table : data frame of the file (formatted with date and proper numeric),
#              with the tow_id and depth_corr columns
# returns : data frame, one row per tow, with columns
#   set : tow number
#   datetime_start, datetime_end, dur : start/end time and duration of the whole tow
#   bottom_start, bottom_end, bottom_dur : start/end time and duration on the bottom (min)
#   Temp..C., depth_corr, Salinity.psu. : median value on the bottom (see average_x())
#   (0 rows if the file has no tow)
# ------------------ -
file_summary <- function(main_table){

  na_mask <- is.na(main_table$tow_id)
  tow_unique <- unique(main_table$tow_id[!na_mask])
  n_tow <- length(tow_unique)
  
  vars2output <- c("Temp..C.", "depth_corr", "Salinity.psu.")
  n_var <- length(vars2output)
  
  summary_df <- preallocate_summary_table(main_table$tow_id, vars2output)
  
  for(i in seq_len(n_tow)){
    
    #subset tow
    mask <- main_table$tow_id == i
    mask[is.na(mask)] <- F #set NA to False
    subset_table <- main_table[mask,]
    n_subset <- dim(subset_table)[1]
    
    #find start and end of tow 
    start_time <- subset_table$Date.Time[1]
    end_time <- subset_table$Date.Time[n_subset]
    
    summary_df$datetime_start[i] <- start_time
    summary_df$datetime_end[i] <- end_time
    summary_df$dur[i] <- end_time - start_time

    # on-bottom part of the tow
    on_bottom <- bottom_mask(subset_table$depth_corr)
    bottom_times <- subset_table$Date.Time[on_bottom]
    bottom_start <- min(bottom_times)
    bottom_end <- max(bottom_times)
    summary_df$bottom_start[i] <- bottom_start
    summary_df$bottom_end[i] <- bottom_end
    summary_df$bottom_dur[i] <- as.numeric(difftime(bottom_end, bottom_start, units = "mins"))

    for(j in seq_len(n_var)){
      
      col_indice <- which(vars2output[j] == colnames(summary_df))
      
      cavg <- average_x(tow = i, x = vars2output[j], main_table)
      
      summary_df[i,col_indice] <- cavg # (j+1 to skip col 1)
      
    }
  }
  return(summary_df)
}


# --------------- -
# file_details()
# split the file into one data frame per tow (all points of the tow, all columns)
#
# main_table : data frame of the file, with the tow_id column
# returns : list of data frames, element i = tow i (empty list if the file has no tow)
# --------------- -
file_details <- function(main_table){
  
  na_mask <- is.na(main_table$tow_id)
  tow_unique <- unique(main_table$tow_id[!na_mask])
  n_tow <- length(tow_unique)
  
  vars2output <- c("Temp..C.", "depth_corr", "Salinity.psu.")
  n_var <- length(vars2output)
  
  #preallocate a list to store all the data frames
  detail_list <- vector(mode='list', length = n_tow)
  
  for(i in seq_len(n_tow)){
    
    #subset tow
    mask <- main_table$tow_id == i
    mask[is.na(mask)] <- F #set NA to False
    subset_table <- main_table[mask,]
    
    # store values of this tow
    detail_list[[i]] <- subset_table 
    
  }

  return(detail_list)
  
}


# --------------- -
# preallocate_summary_table()
# create the empty (NA-filled) summary table used by file_summary()
# times are in the hard-coded "America/Toronto" time zone
#
# tow_id : output of find_tow() (tow number or NA for each point)
# vars2output : names of the variable columns to add (numeric, NA)
# returns : data frame, one row per tow : set, datetime_start, datetime_end, dur,
#           bottom_start, bottom_end, bottom_dur, then one column per vars2output
# --------------- -
preallocate_summary_table <- function(tow_id, vars2output){
  
  na_mask <- is.na(tow_id)
  tow_unique <- unique(tow_id[!na_mask])
  n_tow <- length(tow_unique)
  n_var <- length(vars2output)
  
  #hardcoded time zone
  TZ <- "America/Toronto"
  
  # preallocate a generic posix vector
  dummy_posix <- as.POSIXct(
    rep(NA_real_, n_tow),
    origin = "1970-01-01",
    tz = TZ
  )
  
  # preallocate duration object
  dummy_difftime <- as.difftime(
    rep(NA_real_, n_tow),
    units = "mins")
  
  # preallocate output dataframe
  summary_df <- data.frame(set=tow_unique,
                           datetime_start=dummy_posix,
                           datetime_end=dummy_posix,
                           dur=dummy_difftime,
                           bottom_start=dummy_posix,
                           bottom_end=dummy_posix,
                           bottom_dur=dummy_difftime)
  
  # add variables as cols
  for(i in seq_len(n_var)){
    summary_df[[vars2output[i]]] <- rep(NA_real_, n_tow)
  }  
  
  return(summary_df)
  
}


# ----------------------- -
# get_true_depth()
# keep only the "true_depth" rows of the ANDES sensor data
# (reference depth, to validate the logger depth)
#
# sens_data : ANDES sensor data (data frame with a `type` column)
# returns : data frame, the "true_depth" rows of sens_data
# ----------------------- -
get_true_depth <- function(sens_data){
  mask <- sens_data$type == "true_depth"
  depth_data <- sens_data[mask,]
}


# ----------------------- -
# depth_validation()
# (!) not written yet
# goal : compare the logger depth with the reference depth from ANDES (get_true_depth())
# ----------------------- -
depth_validation <- function(){
  
}








