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
                     enter=7, exit=1.5,
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

# ----------------------- -
# average_x()
# compute average (median) value of variable x
#
# x :
# main_table : data frame of the file (formatted with date and proper numeric)
# tow : tow number (sequential value starting at 1 )
# x : target variable (i.e. date, temp, depth, sal, conduct, vel)
# time_th_lo :  (sec)
# time_th_up : (sec)
# returns : 
# ----------------------- -
average_x <- function(tow, x, main_table, time_th_lo=30, time_th_up=240){
  
  # check arguments and flag error
  cols <- colnames(main_table)
  if(!any(x == cols[2:7])){
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
  n<- dim(subset_table)[1] # dimension of the 
  
  # cut the very start and end of the tow 
  # the goal is to drop values when the instrument is being handled (going up/down or moving around at subsurface)
  datetime<-subset_table$Date.Time
  start_time <- datetime[1]
  end_time <- datetime[n]
  start_time_th <- start_time + time_th_lo
  end_time_th <- end_time - time_th_up
  
  # flag
  if(end_time_th<start_time_th){
    stop("end time cannot be before start time.\n a tow must have beend misidentified.")
  }
  

  keep_time <- (datetime >= start_time_th) & (datetime <= end_time_th)
  
  x_values <- subset_table[keep_time, col_indice]
  avg <- median(x_values)
  
}




# ------------------ -
#
#
# ------------------ -
file_summary <- function(main_table){

  na_mask <- is.na(main_table$tow_id)
  tow_unique <- unique(main_table$tow_id[!na_mask])
  n_tow <- length(tow_unique)
  
  vars2output <- c("Temp..C.", "Depth.m.", "Salinity.psu.")
  n_var <- length(vars2output)
  
  summary_df <- preallocate_summary_table(main_table$tow_id, vars2output)
  
  for(i in 1:n_tow){
    
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
    
    for(j in 1:n_var){
      
      col_indice <- which(vars2output[j] == colnames(summary_df))
      
      cavg <- average_x(tow = i, x = vars2output[j], main_table)
      
      summary_df[i,col_indice] <- cavg # (j+1 to skip col 1)
      
    }
  }
  return(summary_df)
}


# --------------- -
#
#
#
# --------------- -
file_details <- function(main_table){
  
  na_mask <- is.na(main_table$tow_id)
  tow_unique <- unique(main_table$tow_id[!na_mask])
  n_tow <- length(tow_unique)
  
  vars2output <- c("Temp..C.", "Depth.m.", "Salinity.psu.")
  n_var <- length(vars2output)
  
  #preallocate a list to store all the data frames
  detail_list <- vector(mode='list', length = n_tow)
  
  for(i in 1:n_tow){
    
    #subset tow
    mask <- main_table$tow_id == i
    mask[is.na(mask)] <- F #set NA to False
    subset_table <- main_table[mask,]
    
    # store values of this tow
    detail_list[[i]] <- subset_table 
    
  }

  return(detail_list)
  
}


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
                           dur=dummy_difftime)
  
  # add variables as cols
  for(i in 1:n_var){
    summary_df[[vars2output[i]]] <- NA
  }  
  
  return(summary_df)
  
}














