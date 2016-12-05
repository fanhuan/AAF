# Figure3R
#path to AAF

klist <- c(15,17,19,21,23,25,27,29,31,33,35,37,39)

splist = c("CE","CI","CS","FA","FF","FL","FM","FR","FT","FV","GB","IB","IP1","IP2","LB","LC","LF","LG","LH","LR","LX","TD")

gslist <- c(1117239483,1106896350,509025482,469640962,249777258,353529542,350740303,1008785924,687870506,825221352,513827382,1670344889,1700338357,2068951903,1086176530,1277481710,1060165924,1170269382,1162554067,1407387220,1127149933,911834516)
gclist <- c(0.397,0.419,0.387,0.385,0.384,0.382,0.358,0.326,0.345,0.3596,0.3602,0.3234,0.3223,0.3977,0.3991,0.3512,0.382,0.356,0.3447,0.3553,0.3665,0.3662)
coverage <- c(1.84,1.8,4.03,4.75,3.58,5.81,4.24,8.99,7.02,11.68,4.41,5.49,5.99,4.79,1.62,1.64,9.22,1.71,2.79,6.39,2.39,3.92)
# for D=.1 (based on primate, since they have similar divergence time)
D <- 0.1
kk <- 15:39
plot(exp(-kk*D) ~ kk, type='l', lty=5, ylim=c(0,1), xlab='k', ylab=expression(p[h])) #theoretical line without homoplasy
#text(x=30, y=.95, labels='D', cex=1.5)
#text(x=25, y=.8, labels='d = 0.02', cex=1.1)
text(x=25, y=.25, labels='d = 0.1', cex=1.1)

for(i.sp in 1:length(splist)){ #loop for species
  sp <- splist[i.sp]
  g <- gslist[i.sp]
  h <- gclist[i.sp]
  phlist <- NULL
  for(ik in 1:length(klist)){ #loop for k 
    k <- klist[ik]
    filename <- paste("path_to_AAF/AAF/TropicalTrees/Q_k/k", k,"/",splist[i.sp],".hist", sep="")
    Q <- read.table(filename, header=FALSE) #the kmer abundance file for i.sp with ik 
    q <- as.numeric(Q[,2]) #the second column of the kmer abundance file, which is the frequency of each abundance
    nn <- sum(q)
    q <- q/nn # Percentage of each abundance
    
    q.breaks <- round(coverage[i.sp]*(1:floor(length(q)/coverage[i.sp])))
    
    # if want to exclude singletons as caused by sequencing error
    q.breaks <- c(2,q.breaks)
    q.corrected <- array(0,c(1,length(q.breaks)-1))
    for(i in 1:(length(q.breaks)-1))
      q.corrected[i] <- sum(q[q.breaks[i]:(q.breaks[i+1]-1)])
    q <- q.corrected/sum(q.corrected)
    
    ekD <- exp(-k*D)
    gg <- (g-k+1)    
    pp <- 2*(.5 - h + h^2)^(k)
    temp <- 0
    for(i in 1:length(q)) temp <- temp + (1-ekD)^i*q[i]
    ph1 <- 1-temp
    
    ph2 <- (1-ekD)*(1-(1-pp)^gg)
    
    temp <- 0
    for(i in 1:length(q)) temp <- temp + (1-(1-pp)^gg)^i*q[i]
    ph3 <- .5*ekD*temp
    
    ph <- ph1 + (1-ph1)*ph2 + (1-ph1)*(1-ph2)*ph3
    phlist <- c(phlist,ph)
  }
  lines(phlist ~ klist, col="red")
  print(splist[i.sp])
}
