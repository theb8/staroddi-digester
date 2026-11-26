# experimental script to digest star odi data
#
# (!) add description when it's stable
#
# author : atheberge

# what to do :
# - open file & read data
# - identify which point belong to which set
# - give absolute and relative set ID
# - travailler en relatif pour commencer
#
# - à la fin un module pour calibrer ou comparer la prof

source("R/func.R")

enc<-"latin1"
fname <- file.path("data/iml202619/2026-06-12.CSV")

input_table <- read.csv(file = fname, header = T, sep = ";", fileEncoding = enc)

# make cols numeric instead of char
col2makenum <- c(3:7)
for(i in col2makenum){
  x <- input_table[,i]
  input_table[,i] <- charcomma2numdecimal(x)
}

datetime <- as.POSIXct(input_table$Date.Time, format="%d.%m.%Y %H:%M:%S")
input_table$Date.Time <- datetime

depth_raw <- input_table$Depth.m.
offset <- find_surface(depth_raw)
depth <- depth_raw - offset

tow_id <- find_tow(depth, datetime)

drift_report <- drift_check(depth_raw, datetime, tow_id)

# a more complete data frame table
main_table <- input_table
main_table$tow_id <- tow_id
main_table$depth_corr <- depth

summary_df <- file_summary(main_table)
details_df <- file_details(main_table)

# WIP ---- -
#fname_sensors <- file.path("calib/IML2026019F_sensor_data.csv")
#data_sensors <- read.csv(fname_depthref, encoding = enc)
#depthref <- get_true_depth(data_sensors)

