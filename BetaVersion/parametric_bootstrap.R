##############################
##############################
# parametric bootstrap
##############################
##############################

#input file needed: tip_info_test.txt, test_nshare.csv
#install.packages('ape') #if ape is not installed yet
#Replace “path_to_AAF” with the directory where the uncompressed AAF package is
#setwd('path_to_AAF/AAF/')

# set parameters
nboot <- 10 # number of replicates
k <- 21
filter <- 1
info.file <- 'data/tip_info_test.txt'
n.table.file <- 'data/test_nshare.csv'

# load the two input files. info.file is with coverage, read length and
# sequencing error, same file inputed for tip-correction. n.table.file is a
# matrix with numbers of shared kmers between any pair of species.
info <- read.table(info.file, header=T)
spname <- info$spname
nspp <- length(info$spname)
coverage <- info$coverage
E <- info$errorrate
r <- info$readlength
CORRECTIONFACTOR <- 2 #correction for estimation.
kmernames <- as.character(spname)
truen.table <- read.csv(n.table.file, header=T)

# create true dist file
truedist <- array(0,c(nspp,nspp))
for(i in 1:nspp){
  for(j in 1:nspp) {
    n.total <- min(truen.table[i,i],truen.table[j,j])
    truedist[i,j] <- -(1/k)*log(truen.table[i,j]/n.total)
  }
}

library(ape)

system("rm outfile outtree")
system("rm consensus_trees_read_parametric")
system("rm consensus_trees_total_parametric")

system("touch consensus_trees_read_parametric")
system("touch consensus_trees_total_parametric")

for(iboot in 1:nboot){

  # add approximate variation in truedist
  dist <- truedist
  for(ii in 1:(nspp-1))
    for(jj in (ii+1):nspp){
    
      # make sp i have the smallest n.total
      index <- c(ii,jj)
      index <- index[order(c(truen.table[ii,ii], truen.table[jj,jj]))]
      i <- index[1]
      j <- index[2]

      nn <- min(truen.table[i,i], truen.table[j,j])
      D <- truedist[i,j]

      # not filtered
      if(filter==1){
        Li <- coverage[i] * (r[i] - k + 1)/r[i]
        Lj <- coverage[j] * (r[j] - k + 1)/r[j]

        pri <- 1 - exp(-Li)
        prj <- 1 - exp(-Lj)

        pei <- 1 - (exp(-Li*(1-E[i])^k)-exp(-Li))/(1-exp(-Li))
        pej <- 1 - (exp(-Lj*(1-E[j])^k)-exp(-Lj))/(1-exp(-Lj))

        #NOTE: these are calculated for the sp with the smallest n.total
        pta <- Li*(1-(1-E[i])^k)
        psa <- Li*(1-(1-E[i]/3)^(k*D))

        ms <- nn*exp(-k*D)*(pei*pej*pri*prj + psa)
        vs <- nn*exp(-k*D)*(pei*pej*pri*prj*(1 - pei*pej*pri*prj) + psa)

        mt <- nn*(pei*pri + pta)
        # with covariances
        w <- 0
        for(ik in 1:(k-1))	w <- w + (ik/k)*max(0,2*(r[i]-2*k+ik+1)/(r[i]-k+1))
        vt <- nn*(pei*pri*(1 - pei*pri)*(1+(1/Li)*(r[i]-k)) + pta*(1+w))

        Vs <- log(1+vs/ms^2)
        Vt <- log(1+vt/mt^2)

        approxSE <- CORRECTIONFACTOR*((Vs + Vt)^.5)/k
        dist[i,j] <- truedist[i,j] + approxSE*rnorm(1)
        dist[j,i] <- dist[i,j]
      }
      # filtered
      if(filter==2){
        Li <- coverage[i] * (r[i] - k + 1)/r[i]
        Lj <- coverage[j] * (r[j] - k + 1)/r[j]
      
        pr2i <- 1 - exp(-Li) - Li*exp(-Li)
        pr2j <- 1 - exp(-Lj) - Lj*exp(-Lj)

        EE <- (1 - E[i])^k
        q0i <- exp(-Li)
        q1i <- Li*exp(-Li)
        G1i <- (exp(-Li*EE) - q0i - q1i)/(1 - q0i - q1i)
        G2i <- EE*Li*exp(-Li*EE)/(1 - q0i - q1i)
        pe2i <- 1 - G1i - G2i

        EE <- (1 - E[j])^k
        q0j <- exp(-Lj)
        q1j <- Lj*exp(-Lj)
        G1j <- (exp(-Lj*EE) - q0j - q1j)/(1 - q0j - q1j)
        G2j <- EE*Lj*exp(-Lj*EE)/(1 - q0j - q1j)
        pe2j <- 1 - G1j - G2j

        ms <- nn*exp(-k*D)*(pe2i*pe2j*pr2i*pr2j)
        vs <- nn*exp(-k*D)*(pe2i*pe2j*pr2i*pr2j*(1 - pe2i*pe2j*pr2i*pr2j))

        mt <- nn*(pe2i*pr2i)
        vt <- nn*(pe2i*pr2i*(1 - pe2i*pr2i))

        Vs <- log(1+vs/ms^2)
        Vt <- log(1+vt/mt^2)

        approxSE <- CORRECTIONFACTOR*((Vs + Vt)^.5)/k
        dist[i,j] <- truedist[i,j] + approxSE*rnorm(1)
        dist[j,i] <- dist[i,j]
      }
    }

  # write phylip file for fitch_kmerX
  write(nspp, file="infile")
  for(i in 1:nspp){
    spnames10 <- paste(c(kmernames[i],character(length=10-nchar(kmernames[i]))),collapse=" ")
    writeline <- paste(c(spnames10, dist[i,]), collapse=" ")
    write(writeline, file="infile", append=T)
  }

  # fitch_kmerX
  system2('./fitch_kmerX', input = c('K',k,'Y'), stdout = NULL)

  # import tree and convert to a covariance matrix
  Ptree <- read.tree(file="outtree")

  W <- cophenetic.phylo(Ptree)
  W <- W[sort(rownames(W)),sort(rownames(W))]

  # write phylip file for PHYLIP consense
  system("cat outtree >> consensus_trees_read_parametric")
  system("rm outfile outtree")

  ##############################
  # STAGE 2: bootstrap for evolutionary variation
  ##############################

  distread <- dist
  dist <- W

  # create covariance matrix for distances
  distv <- array(0,c((nspp*(nspp-1)/2)^2,1))
  distvcounter <- 0
  for(i1 in 1:(nspp-1))
    for(j1 in (i1+1):nspp){
      for(i2 in 1:(nspp-1))
        for(j2 in (i2+1):nspp){
          distvcounter <- 1 + distvcounter
        
          d1 <- dist[i1,j1]
          d2 <- dist[i2,j2]
        
          da <- dist[i1,j2]
          db <- dist[i1,i2]
          dc <- dist[i2,j1]
          dd <- dist[j1,j2]
        
          jointdist <- max(0,.5*(2*(d1 + d2) - (da + db + dc + dd)))
          nt1 <- min(truen.table[i1,i1], truen.table[j1,j1])
          nt2 <- min(truen.table[i2,i2], truen.table[j2,j2])
        
          eD <- exp(jointdist)
          w <- (eD^k - 1)
          for(i in 1:(k-1)) w <- w + 2*(eD^i - 1)
        
          approxCov <- (w/(k^2*nt1^.5*nt2^.5))
          distv[distvcounter] <- approxCov
        }
    }

  distV <- array(0,c(nspp*(nspp-1)/2,nspp*(nspp-1)/2))
  vcounter <- 0
  for(ii in 1:(nspp*(nspp-1)/2))
    for(jj in 1:(nspp*(nspp-1)/2)){
        vcounter <- 1 + vcounter
        distV[ii,jj] <- distv[vcounter]
        distV[jj,ii] <- distv[vcounter]
    }
  iD <- chol(distV)

  approxdistbootlist <- NULL
  
  dist <- distread
  
  # add approximate variation in dist
  e <- iD %*% rnorm(n=nspp*(nspp-1)/2)
  ecounter <- 0
  for(ii in 1:(nspp-1))
    for(jj in (ii+1):nspp){
        ecounter <- ecounter + 1
  
        dist[ii,jj] <- dist[ii,jj] + e[ecounter]
        dist[jj,ii] <- dist[ii,jj]
    }


  # write phylip file for fitch_kmer
  write(nspp, file="infile")
  for(i in 1:nspp){
      spnames10 <- paste(c(kmernames[i],character(length=10-nchar(kmernames[i]))),collapse=" ")
      writeline <- paste(c(spnames10, dist[i,]), collapse=" ")
      write(writeline, file="infile", append=T)
  }
  # fitch_kmerX
  system2('./fitch_kmerX', input = c('K',k,'Y'), stdout = NULL)
  # write phylip file for PHYLIP consense
  system("cat outtree >> consensus_trees_total_parametric")
  system("rm outfile")
  system("rm outtree")
}
system("./consense", input = c('consensus_trees_read_parametric','Y'))
system("mv outfile consense_outfile_read_parametric")
system("mv outtree consensus_read_parametric.tre")

system("./consense", input = c('consensus_trees_total_parametric','Y'))
system("mv outfile consense_outfile_total_parametric")
system("mv outtree consensus_total_parametric.tre")
system("rm infile")
