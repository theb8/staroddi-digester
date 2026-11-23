#
#
#
#
#

# x is a char vector with a single comma
# will not fail if x is numeric
charcomma2numdecimal <- function(x){
  x_dot <- sub(pattern = ",", replacement = ".", x = x)
  x_num <- as.numeric(x_dot)
}


# will return a depth value
find_surface <- function(depth){
  # find surface
  surface_anchor <- quantile(depth, probs = 0.05) #5th percentile to
  surface_level <- 3
  surface_mask <- depth < (surface_anchor + surface_level)
  surface_points <- depth[surface_mask]
  zero <- median(surface_points)
  return(zero)
}

# --------------------- -
# find_tow
# will find tow events of the instrument
# "hysteresis" approach to avoid ambiguous depth
#
# ENTER (m) : depth threshold to start a tow
# EXIT (m) : depth threshold to end a tow
# MIN_GAP : min timegap to merge two tows (min)
# MIN_DUR : minimum duration (min)
# returns :
# ---------------------
find_tow <- function(depth, datetime, 
                     enter=5, exit=1.5,
                     min_gap=2, min_dur=2){
  
  n <- length(depth)
  
  if(length(datetime) != n){
    stop("depth and datetime must be the same length")
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
  
  # keep the event only if it towd more than `enter`; which is 5 meters default.
  keep <- r$values & run_max > enter
  tow_indices <- data.frame(start = starts[keep], end = ends[keep])
 
  # --- -
  # check 1 : time gap checkup
  # tow events were identified 
  # now let's take a look at the time difference identify false tow events
  # --- -
  
  # compute time gap for all tow events
  gap <- as.numeric(difftime(datetime[tow_indices$start[-1]], datetime[tow_indices$end[-nrow(tow_indices)]]))
  
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
  dur  <- as.numeric(difftime(datetime[tow_indices$end], datetime[tow_indices$start], units = "mins"))
  tow_indices   <- tow_indices[dur >= min_dur,]
  
  n_tows <- dim(tow_indices)[1]
  
  tow_id <- rep(NA_integer_, n)
                
  for(i in 1:n_tows){
    
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
# depth_raw : 
# datetime :
# time_id : 
# returns : 
# ----------------------- -
drift_check <- function(depth_raw, datetime, tow_id){
  
  # hard coded thresholds
  TIME_TH_LO <-  30 # will discard 30 seconds after each previous tow (unit sec)
  TIME_TH_UP <- 4*60 # will discard 5 minutes just before the next tow; since the gear is being moved around at the surface (unit sec)
  DRIFT_TH <- 0.5 # hard threshold to determine if sensor is drifting (unit meters)
  
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
  file_zeros <- rep(NA, n_groups)
  file_spread  <- rep(NA, n_groups)
  
    
  for(i in 1:n_groups){
    
    indice_lo <- surface_indices$start[i]
    indice_up <- surface_indices$end[i]
    
    #trim 3 minutes on each side
    start_time <- datetime[indice_lo]
    end_time <- datetime[indice_up]
    
    start_time_th <- start_time+TIME_TH_LO
    end_time_th <- end_time - TIME_TH_UP
    
    # flag
    if(end_time_th<start_time_th){
      stop("end time cannot be before start time.\n a tow must have beend misidentified.")
    }
    
    #mask using timestamp
    keep_time <- (datetime >= start_time_th) & (datetime <= end_time_th)
    n_kept <- sum(keep_time) #number of points for this "surface event"
    
    # another flag
    if(n_kept < 10){
      warning("number of points is too low for this surface event\n drift check might not work properly")
    }
    
    # the depth values that are kept
    depth_window <- depth_raw[keep_time]
    
    # zero found for this surface event
    local_zero <- median(depth_window)
    
    # store computed values
    file_zeros[i] <- local_zero
    file_spread[i] <- max(depth_window)-min(depth_window)
    
  }
  
  # check if a surface event is outside our depth threshold
  zero_ok <- (file_zeros <= (global_zero + DRIFT_TH)) & (file_zeros >= (global_zero - DRIFT_TH))
  
  if(!all(zero_ok)){
    warning_msg <- sprintf("some surface events are outside threshold of %s meters\n  please verify if the instrument is drifting", DRIFT_TH)
    warning(warning_msg)
  }
  
  list2return <- list(indices=surface_indices,
                      local_zeros=file_zeros,
                      local_spread=file_spread)
  return(list2return)
  
}
