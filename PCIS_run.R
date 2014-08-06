source("PMI_IVS.R")

### Comment the below out if not using shell script to run code and just enter
### values for 'fn' and 'dat.dir' manually..
### ------
# Read in of a string argument to indicate segment
args <- commandArgs(TRUE)
fn <- args[1]   #   filename - e.g. 'bank_8_24fh_1.csv'
dat.dir <- args[2]   #  dataset - e.g. 'Bank'
### ------

# Read in data
# ----
out.dir <- paste("~/", dat.dir, "/", sep = "")

inp.file <- fn
fn <- unlist(strsplit(fn, split = "/"))
fn <- fn[length(fn)]

names.1 <- scan(file = inp.file, what = character(), sep = ",",
                nlines = 1)
dat <- read.table(inp.file, skip = 1, sep = ",", fill = TRUE)

id <- dat[,1]
dat <- dat[,-1]
names.1 <- names.1[-1]
# store output variable in 'y'
y <- dat[,ncol(dat)]
# store input variables in 'x'
x <- dat[,1:(ncol(dat) - 1)]

names(y) <- names.1[ncol(dat)]
names(x) <- names.1[1:(ncol(dat) - 1)]

ndata <- length(y)

# Call function 'PCIS_framework' from file "PMI.R" to compute best inputs
# based on PCIS. Output is written to file within this function
PCIS_framework(x, y, iter = min(ncol(x), 30), 
               outfile = paste(out.dir, fn, "_PCIS_out.csv", sep = ""), 
               mifile = paste(out.dir, fn, "_PCIS_mi.txt", sep = ""),
               silent = FALSE)



