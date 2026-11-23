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
# will find "tow events" of the instrument
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


drift_check <- function(){
  
  
  
}
