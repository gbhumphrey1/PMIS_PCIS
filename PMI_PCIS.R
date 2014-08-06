library(KernSmooth)
library(nnet)
#-------------------------------------------------------------------------------
standardize <- function(Z) {
  Z.mean <- mean(Z)
  Z.sd <- sd(Z)
  Z.stand <- (Z - Z.mean) / Z.sd
  return(Z.stand)
}
#-------------------------------------------------------------------------------
dot.prod <- function(s, d) {
  s1 <- s[1:d]
  s2 <- s[(d + 1):(2*d)]
  return(sum(s1 * s2))
}
#-------------------------------------------------------------------------------
density.estimation.d1 <- function(Z) {

  dens <- 1 / (1 + exp(-Z))
  return(dens)
}
#-------------------------------------------------------------------------------
get.norm.sigma <- function(Z) {

  Z <- as.matrix(Z)
  d <- ncol(Z)
  sd <- sqrt(apply(Z, 2, var))
  sigma <- sd * (4/(d + 2))^(1/(d + 4)) * N^(-1/(d + 4))
}
#-------------------------------------------------------------------------------
kernel.est.uvn <- function(Z) {

  N <- length(Z)
  d <- 1
# compute sigma & constant
  sigma <- bw.nrd(Z)
  constant <- sqrt(2*pi) * sigma * N

# Commence main loop
  dens <- vector()
  for(h in 1:N) {
    dis.Z <- (Z - Z[h])^2
    exp.dis <- exp(-dis.Z / (2*sigma^2))
    dens[h] <- sum(exp.dis) / constant
  }
  return(dens)
}
#-------------------------------------------------------------------------------
density.estimation.d2 <- function(Z.1, Z.2) {

  d.1 <- density.estimation.d1(Z.1)
  d.2 <- density.estimation.c1(Z.2)

  dens <- d.1 * d.2

  return(dens)
}
#-------------------------------------------------------------------------------
kernel.est.mvn <- function(Z) {

# Compute covariance and determinant of cov
  N <- nrow(Z)
  d <- ncol(Z)
  Cov <- cov(Z)
  det.Cov <- det(Cov)

# Compute sigma & constant
  sigma <- (4/(d + 2))^(1/(d + 4)) * N^(-1/(d + 4))
  constant <- (sqrt(2*pi)*sigma)^d * sqrt(det.Cov) * N

# Commence main loop
  dens <- vector()
  for(h in 1:N) {
    dist.val <- mahalanobis(Z, center = Z[h,], cov = Cov)
    exp.dis <- exp(-dist.val / (2*sigma^2))
    dens[h] <- sum(exp.dis) / constant
  }

  return(dens)
}
#-------------------------------------------------------------------------------
pmi.calc <- function(Y, X, Y.bin = FALSE) {
# Use Gaussian kernel to compute marginal densities F(Y) & F(X) and
# joint density F(X,Y)
  N <- length(Y)
  pdf.X <- kernel.est.uvn(X)
  if(Y.bin) {
    pdf.Y <- density.estimation.d1(Y)
    pdf.XY <- density.estimation.d2(Y, X)
  } else {
    pdf.Y <- kernel.est.uvn(Y)
    pdf.XY <- kernel.est.mvn(cbind(Y, X))
  }
  calc <- log(pdf.XY / (pdf.Y * pdf.X))
  return(sum(calc) / N)
}
#-------------------------------------------------------------------------------
enp <- function(Z, linear = FALSE) {
# Compute effective number of parameters

  Z <- as.matrix(Z)
  if(linear) {
    A.mat <- Z %*% t(Z)
    A.inv <- chol2inv(A.mat)
    H.mat <- t(Z) %*% A.inv %*% Z
    return(sum(diag(H.mat)))
  } else {
    N <- nrow(Z)
    d <- ncol(Z)
    sigma <- (4/(d + 2))^(1/(d + 4)) * N^(-1/(d + 4))
    Cov <- cov(Z)
    det.Cov <- det(Cov)
    inv.Cov <- chol2inv(Cov)
    constant <- (sqrt(2*pi)*sigma)^d * sqrt(det.Cov)

    H.mat <- array(dim = c(N, N))
    for(h in 1:N) {
      dis.Z <- vector()
      for(j in 1:d) {
        dis.Z <- cbind(dis.Z, Z[h,j] - Z[,j])
      }
      if(d == 1) {
        exp.dis <- exp(-(dis.Z)^2 / (2*sigma^2)) / constant
      } else {
        T.dis.Z <- apply(dis.Z, 1, "%*%", inv.Cov)
        M.dis.Z <- apply(cbind(t(T.dis.Z), dis.Z), 1, dot.prod, d = d)
        exp.dis <- exp(-M.dis.Z / (2*sigma^2)) / constant
      }
      H.mat[h,] <- exp.dis / sum(exp.dis)
    }
  }
  return(sum(diag(H.mat)))
}
#-------------------------------------------------------------------------------
GRNN.mah <- function(Z.out, Z.in) {

  Z.in <- as.matrix(Z.in)
  N <- nrow(Z.in)
  d <- ncol(Z.in)
  sigma <- (4 / (d + 2))^(1 / (d + 4)) * N^(-1 / (d + 4))
  Cov <- cov(Z.in)
  det.Cov <- det(Cov)
  constant <- (sqrt(2*pi)*sigma)^d * sqrt(det.Cov)

# Commence main loop
  Z.hat <- vector()
  for(i in 1:N) {
    dis <- vector()
    for(j in 1:d) {
      dis <- cbind(dis, Z.in[-i,j] - Z.in[i,j])
    }
    if(d == 1) {
      exp.dis <- exp(-(dis)^2 / (2*sigma^2)) / constant
    } else {
      dist.val <- mahalanobis(Z.in[-i,], center = Z.in[i,], cov = Cov)
      exp.dis <- exp(-dist.val / (2*sigma^2)) / constant
    }


    dot.A <- sum(Z.out[-i] * exp.dis)      # Sum numerator
    dot.B <- sum(exp.dis)              # Sum denominator
    dot.B <- max(dot.B, 1e-6)
    Z.hat[i] <- dot.A / dot.B
  }
  return(Z.hat)
}
#-------------------------------------------------------------------------------
GRNN.gauss <- function(Z.out, Z.in) {

  Z.in <- as.matrix(Z.in)
  N <- nrow(Z.in)
  d <- ncol(Z.in)
  sigma <- (4 / (d + 2))^(1 / (d + 4)) * N^(-1 / (d + 4))
  Cov <- cov(Z.in)
  det.Cov <- det(Cov)
  constant <- (sqrt(2*pi)*sigma)^d * sqrt(det.Cov)

# Commence main loop
  Z.hat <- vector()
  for(i in 1:N) {
    dist.val <- vector()
    for(j in 1:d) {
      dist.val <- cbind(dist.val, (Z.in[-i,j] - Z.in[i,j])^2)
    }
    if(d > 1) {
      dist.val <- rowSums(dist.val)   # Sum Euclidean distances
    }
    exp.dis <- exp(-dist.val / (2*sigma^2)) / constant


    dot.A <- sum(Z.out[-i] * exp.dis)      # Sum numerator
    dot.B <- sum(exp.dis)              # Sum denominator
    dot.B <- max(dot.B, 1e-6)
    Z.hat[i] <- dot.A / dot.B
  }
  return(Z.hat)
}
#-------------------------------------------------------------------------------
PNN <- function(Z.out, Z.in) {

  Z.in <- as.matrix(Z.in)
  N <- nrow(Z.in)
  d <- ncol(Z.in)
  sigma <- (4/(d + 2))^(1 /(d + 4)) * N^(-1/(d + 4))

# Commence main loop
  Z.hat <- vector()
  class.vals <- unique(Z.out)
  class.A <- which(Z.out == class.vals[1])
  class.B <- which(Z.out == class.vals[2])
  for(i in 1:N) {
    if(d > 1) {
      dist.A <- rowSums((Z.in[class.A[class.A != i],] - matrix(Z.in[i,], nrow = 1))^2)
      dist.B <- rowSums((Z.in[class.B[class.B != i],] - matrix(Z.in[i,], nrow = 1))^2)
    } else {
      dist.A <- (Z.in[i,] - Z.in[class.A,])^2
      dist.B <- (Z.in[i,] - Z.in[class.B,])^2
    }
    dist.A <- exp(-dist.A / (2*sigma^2))
    dist.B <- exp(-dist.B / (2*sigma^2))
    dot.A <- sum(dist.A)/length(class.A[class.A != i])
    dot.B <- sum(dist.B)/length(class.B[class.B != i])
    Z.hat[i] <- ifelse(dot.A > dot.B, class.vals[1], class.vals[2])
  }
  return(Z.hat)
}
#-------------------------------------------------------------------------------
MAD <- function(Z) {
# To compute median absolute deviation
#--------------------------------------

# Compute median PMI values
  median.Z <- median(Z)

# Compute median absolute deviation of PMI values
  MAD.Z <- abs(Z - median.Z)

# Compute S
  S <- 1.4826 * median(MAD.Z)

# Compute MAD Z.score
  Z.score <- MAD.Z / S

  return(Z.score)

}
#-------------------------------------------------------------------------------
PCIS_framework <- function(x, y, iter = ncol(x), outfile = "out.txt", 
                           mifile = "mi.txt", silent = FALSE) {

#** Stepwise PCIS**
# Author : Greer Humphrey
# School of Civil and Environmental Engineering
# University of Adelaide, SA, 5005, Australia

  cpu.time <- proc.time()
  n.inputs <- ncol(x)
  n.data <- nrow(x)
  inp.names <- names(x)

# record number of operations performed on data
  nfevals.1 <- 0 # number of linear model calls
  nfevals.2 <- 0 # number of correlation calls

# Standardize the output data with zero mean and unit variance
  y.stand <- standardize(y)
  y <- y.stand

# Standardize the input data with zero mean and unit variance
  z.stand <- apply(x, 2, standardize)
  x <- z.stand

# Mock up an array to keep track of unselected inputs
  input.tracker <- 1:n.inputs
  counter <- n.inputs
  n.selected <- 0
  z.in <- vector()
  scores <- NULL
  max.iter <- (min(iter, ncol(x)) + 1)


  for(iter.1 in 1:max.iter) {
    if(n.selected > 0) {
    # Reset the array with the output in it back to the original values
      y <- y.stand
    # Compute output with already selected inputs and output residuals
      y.hat <- lm(y ~ z.in)$fitted.values
      nfevals.1 <- nfevals.1 + 1

    # Calculate residuals u  Note: the y array is reset back to
    # the original values at start of next iteration (see above)
      u <- y - y.hat

    # Standardize output residuals
      u.stand <- standardize(u)

      mean.y <- mean(y)
      AIC.k <- n.data * log(sum(u^2) / n.data) + 2*(n.selected + 1)
      BIC.score <- n.data * log(sum(u^2) / n.data) + log(n.data)*(n.selected + 1)
      RMSE <- sqrt(sum(u^2) / n.data)
      R2 <- 1 - sum(u^2)/sum((y - mean.y)^2)

      scores.0 <- data.frame(Input = format(inp.names[input.tracker[tag]]), 
                             PC = signif(PC[tag], 4),
                             AIC_k = signif(AIC.k, 4),
                             BIC = signif(BIC.score, 4),
                             RMSE = signif(RMSE, 4),
                             R2 = signif(R2, 4),
                             Hampel = signif(Z.score[tag], 4),
                             ModEvals = nfevals.1,
                             PCEvals = nfevals.2,
                             CPUtime = signif(as.numeric(proc.time()[2] - cpu.time[2]), 4),
                             ElapsedTime = signif(as.numeric(proc.time()[3] - cpu.time[3]), 4))

      scores <- rbind(scores, scores.0)
      if(iter.1 == 2) {
        write.table(scores.0, file = outfile, row.names = FALSE, col.names = TRUE, sep = ",")
      } else {
        write.table(scores.0, file = outfile, row.names = FALSE, col.names = FALSE, sep = ",", append = TRUE)
      }


    # Output results to screen and file
      if(!silent) print(scores, quote = FALSE)
      if(iter.1 == max.iter) break # break if reached maximum iterations

    # Remove selected input from input data set and add to selected data set
      n.inputs <- n.inputs - 1
      input.tracker <- input.tracker[-tag]
      x <- as.matrix(x[,-tag])

    # Go through each input in turn to compute PMI
    # due to the already selected inputs
      z.res <- vector()
      Input <- vector()
      PC <- vector()
      for(current.inp in 1:n.inputs) {
        z <- x[,current.inp]

      # Now compute input residuals
        z.hat <- lm(z ~ z.in)$fitted.values
        nfevals.1 <- nfevals.1 + 1
        v <- z - z.hat

      # Standardize input residuals
        v.stand <- standardize(v)

      # Store input residuals to use for bootstrapping
        z.res <- cbind(z.res, v.stand)
      # Compute correlation between input and output
        PC[current.inp] <- cor(u.stand, v.stand)^2
        nfevals.2 <- nfevals.2 + 1
        Input[current.inp] <- input.tracker[current.inp]
      }
    } else {
    # Go through each input in turn to compute MI
      Input <- vector()
      PC <- vector()
      for(current.inp in 1:n.inputs) {
        z <- x[,current.inp]
      # Compute correlation between input and output
        PC[current.inp] <- cor(y, z)^2
        nfevals.2 <- nfevals.2 + 1
        Input[current.inp] <- input.tracker[current.inp]
      }
    }

    if(n.selected == 0) {
      write.table(data.frame(Input = format(inp.names[Input]), 
                             PC = format(PC, digits = 4)), file = mifile, 
                  row.names = FALSE, quote = FALSE)
    }

  # Find input with highest PC score
    best.PC <- max(PC)
    tag <- which.max(PC)

  # Output the input that has been selected
    if(!silent) cat("\nSelect input", inp.names[input.tracker[tag]], "....\n")

  # Calculate various scores
    Z.score <- MAD(PC)

    n.selected <- n.selected + 1
    z.in <- cbind(z.in, as.vector(x[,tag]))
}

# The program will only reach this point upon termination
  cat("\nPCIS ROUTINE COMPLETED\n")

  cpu.time <- proc.time() - cpu.time
  write("", file = outfile, append = TRUE)
  out.dat <- as.data.frame(t(as.numeric(cpu.time)[1:3]))
  names(out.dat) <- names(cpu.time)[1:3]
  write.table(out.dat, file = outfile, append = TRUE, sep = ",",
              na = "", row.names = FALSE, col.names = TRUE)
}
#-------------------------------------------------------------------------------
PMI_grnn_framework <- function(x, y, iter = ncol(x), outfile = "out.txt", 
                               mifile = "mi.txt", ybin = FALSE,
                               silent = FALSE) {
# Stepwise PMI
# Use tabulated critical values to identify significant inputs

# Used Gaussian kernel for both PDF estimation & GRNN predictions
# Average MI/PMI computed without considering effect of negative values of MI
  cpu.time <- proc.time()
  n.inputs <- ncol(x)
  n.data <- nrow(x)
  inp.names <- names(x)

  I.crit <- c(50, 60, 70, 80, 90, 100, 120, 140, 160, 180, 200, 220, 240, 260, 
              280, 300, 400, 500, 600, 700, 800, 900, 1000, 2000, 3000, 4000, 5000, 
              0.2224, 0.2031, 0.1879, 0.1756, 0.1657, 0.1572, 0.1434, 0.1321, 
              0.1237, 0.1166, 0.1103, 0.1055, 0.1005, 0.0965, 0.0928, 0.0896, 
              0.0775, 0.0689, 0.0627, 0.0578, 0.0539, 0.0507, 0.0481, 0.0333, 
              0.0268, 0.0230, 0.0204)
  dim(I.crit) <- c(27, 2)
  I.crit <- as.data.frame(I.crit)
  names(I.crit) <- c("n", "val")

# Fit model to tabulated values
  I.crit.mod <- nls(val ~ a*n^b, data = I.crit, 
                    start = list(a = 1, b = 1))
  new.data <- data.frame(n = n.data, val = NA)
  I.crit.val <- predict(I.crit.mod, newdata = new.data)

# record number of operations performed on data
  nfevals.1 <- 0 # number of grnn calls
  nfevals.2 <- 0 # number of pmi calls

# Standardize the output data with zero mean and unit variance
  if(!ybin) {
    y.stand <- standardize(y)
    y <- y.stand
  } else {
    y.stand <- y
  }

# Standardize the input data with zero mean and unit variance
  z.stand <- apply(x, 2, standardize)
  x <- z.stand

# Mock up an array to keep track of unselected inputs
  input.tracker <- 1:n.inputs
  counter <- n.inputs
  n.selected <- 0
  z.in <- vector()
  scores <- NULL
  max.iter <- (min(iter, ncol(x)) + 1)

  for(iter.1 in 1:max.iter) {
    if(n.selected > 0) {

    # Reset the array with the output in it back to the original values
      y <- y.stand
    # Compute output with already selected inputs and output residuals
      if(!ybin) {
        y.hat <- GRNN.gauss(y, z.in)
        nfevals.1 <- nfevals.1 + 1
      } else {
        y.hat <- PNN(y, z.in)
        nfevals.1 <- nfevals.1 + 1
      }
    # Calculate residuals u  Note: the y array is reset back to
    # the original values at start of next iteration (see above)
      u <- y - y.hat

    # Standardize output residuals
      if(!ybin) {
        u.stand <- standardize(u)
      } else {
        u.stand <- u
      }

      mean.y <- mean(y)
      p <- enp(z.in)
      AIC.score <- n.data * log(sum(u^2) / n.data) + 2*p
      AIC.k <- n.data * log(sum(u^2) / n.data) + 2*(n.selected + 1)
      BIC.score <- n.data * log(sum(u^2) / n.data) + log(n.data)*p
      RMSE <- sqrt(sum(u^2) / n.data)
      R2 <- 1 - sum(u^2)/sum((y - mean.y)^2)

      scores.0 <- data.frame(Input = format(inp.names[input.tracker[tag]]), 
                             PMI = signif(PMI[tag], 4),
                             MC_I_95 = signif(I.crit.val, 4),
                             AIC_p = signif(AIC.score, 4),
                             AIC_k = signif(AIC.k, 4),
                             BIC = signif(BIC.score, 4),
                             RMSE = signif(RMSE, 4),
                             R2 = signif(R2, 4),
                             Hampel = signif(Z.score[tag], 4),
                             ModEvals = nfevals.1,
                             PMIEvals = nfevals.2,
                             CPUtime = signif(as.numeric(proc.time()[2] - cpu.time[2]), 4),
                             ElapsedTime = signif(as.numeric(proc.time()[3] - cpu.time[3]), 4))

      scores <- rbind(scores, scores.0)
      if(iter.1 == 2) {
        write.table(scores.0, file = outfile, row.names = FALSE, col.names = TRUE, sep = ",")
      } else {
        write.table(scores.0, file = outfile, row.names = FALSE, col.names = FALSE, sep = ",", append = TRUE)
      }
#      write("", file = outfile, append = TRUE)

    # Output results to screen and file
      if(!silent) print(scores, quote = FALSE)
      if(iter.1 == max.iter) break # break if reached maximum iterations

    # Remove selected input from input data set and add to selected data set
      n.inputs <- n.inputs - 1
      input.tracker <- input.tracker[-tag]
      x <- as.matrix(x[,-tag])

    # Go through each input in turn to compute PMI
    # due to the already selected inputs
      Input <- vector()
      PMI <- vector()
      for(current.inp in 1:n.inputs) {
      # Compute and standardise input residuals
        z <- x[,current.inp]
        z.hat <- GRNN.gauss(z, z.in)
        nfevals.1 <- nfevals.1 + 1
        v <- z - z.hat
        v.stand <- standardize(v)
      # Compute PMI
        PMI[current.inp] <- pmi.calc(u.stand, v.stand, ybin)
        nfevals.2 <- nfevals.2 + 1
        Input[current.inp] <- input.tracker[current.inp]
      }
    } else {
    # Go through each input in turn to compute MI
      Input <- vector()
      PMI <- vector()
      for(current.inp in 1:n.inputs) {
        z <- x[,current.inp]
        PMI[current.inp] <- pmi.calc(y, z, ybin)
        nfevals.2 <- nfevals.2 + 1
        Input[current.inp] <- input.tracker[current.inp]
      }
    }

    if(!silent) print(cbind(inp.names[Input], format(PMI, digits = 4)), quote = FALSE)
    if(n.selected == 0) {
      write.table(data.frame(Input = format(inp.names[Input]), 
                             PMI = format(PMI, digits = 4)), file = mifile, 
                  row.names = FALSE, quote = FALSE)
    }

  # Find input with highest PMI score
    best.PMI <- max(PMI)
    tag <- which.max(PMI)
    if(!silent) cat("\nSelect input", inp.names[input.tracker[tag]], "....\n")

# Calculate various scores
    Z.score <- MAD(PMI)

    z.in <- cbind(z.in, as.vector(x[,tag]))
    n.selected <- n.selected + 1
  }

# The program will only reach this point upon termination
  cat("\nPMI ROUTINE COMPLETED\n")

  cpu.time <- proc.time() - cpu.time
#  write.table(scores, file = outfile, row.names = FALSE, col.names = TRUE, sep = ",")
  write("", file = outfile, append = TRUE)
  out.dat <- as.data.frame(t(as.numeric(cpu.time)[1:3]))
  names(out.dat) <- names(cpu.time)[1:3]
  write.table(out.dat, file = outfile, append = TRUE, sep = ",",
              na = "", row.names = FALSE, col.names = TRUE)
  return()
}
#-------------------------------------------------------------------------------