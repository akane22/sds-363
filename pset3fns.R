library(ks)

default.gridsize <- function(d)
{
  if (d==1)      gridsize <- 401
  else if (d==2) gridsize <- rep(151,d)
  else if (d==3) gridsize <- rep(51, d)
  else if (d>=4) gridsize <- rep(21, d)
  ##else gridsize <- NA
  
  return(gridsize)
}

default.bgridsize <- function(d)
{
  if (d==1)      gridsize <- 401
  else if (d==2) gridsize <- rep(151,d)
  else if (d==3) gridsize <- rep(31, d)
  else if (d==4) gridsize <- rep(15, d)
  else gridsize <- NA
  
  return(gridsize)
}

default.bflag <- function(d, n)
{
  if (d==1) thr <- 1
  else if (d==2) thr <- 500
  else if (d>2) thr <- 1000
  bf <- n>thr
  
  return(bf)
}
ks.defaults <- function(x) {
  d <- ncol(x)
  n <- nrow(x)
  w <- rep(1,n)
  binned <- default.bflag(d=d, n=n)
  bgridsize <- default.bgridsize(d)
  gridsize <- default.gridsize(d)
  if (length(gridsize)==1) {gridsize <- rep(gridsize, d)}
  if (length(bgridsize)==1) {bgridsize <- rep(bgridsize, d)}
  
  return(list(d=d, n=n, w=w, binned=binned, bgridsize=bgridsize, gridsize=gridsize))
}

parse.name <- function(x)
{
  if (is.vector(x))
  {
    d <- 1
    x.names <- deparse(substitute(x))
  }
  else
  {  
    d <- ncol(x)
    x.names <- colnames(x)
    if (is.null(x.names))
    {
      x.names <- strsplit(deparse(substitute(x)), "\\[")[[1]][1]
      x.names <- paste(x.names, "[, ", 1:d,"]",sep="") 
    }
  }
  return(x.names)
}

kda.nd <- function(x, x.group, Hs, prior.prob, gridsize, supp, eval.points, binned, bgridsize, w, compute.cont, approx.cont)
{
  if (is.data.frame(x)) {x <- as.matrix(x)}
  gr <- levels(x.group)
  ##gr <- sort(unique(x.group))
  m <- length(gr)
  d <- ncol(x)
  
  ## find largest bandwidth matrix to initialise grid
  detH <- vector() 
  for (j in 1:m){
    detH[j] <- det(Hs[((j-1)*d+1) : (j*d),])  
  }
  Hmax.ind <- which.max(detH)
  Hmax <- Hs[((Hmax.ind-1)*d+1) : (Hmax.ind*d),]
  ##Hmax <- Hns(x) 
  xmin <- apply(x, 2, min) - supp*max(sqrt(diag(Hmax)))
  xmax <- apply(x, 2, max) + supp*max(sqrt(diag(Hmax)))
  
  
  if (binned & d > 4){ stop("Binning only available for 1- to 4-d data")}
  
  fhat.list <- list()
  for (j in 1:m)
  {
    xx <- x[x.group==gr[j],]
    ww <- w[x.group==gr[j]]   
    H <- Hs[((j-1)*d+1) : (j*d),]     
    if (binned){
      fhat.temp <- kdde(x=xx, bgridsize=bgridsize, H=H, xmin=xmin, xmax=xmax, w=ww, deriv.order=0, binned=TRUE)
    } else {
      fhat.temp <- kde(x=xx, H=H, eval.points=eval.points, w=ww)
    }
    ## compute individual density estimate
    fhat.list$estimate <- c(fhat.list$estimate, list(fhat.temp$estimate))
    fhat.list$eval.points <- fhat.temp$eval.points
    fhat.list$x <- c(fhat.list$x, list(xx))
    fhat.list$H <- c(fhat.list$H, list(H))
    fhat.list$w <- c(fhat.list$w, list(ww))
    
    ## compute prob contour levels
  }
  
  fhat.list$names <- parse.name(x)  ## add variable names
  fhat.list$binned <- binned
  fhat.list$gridded <- fhat.temp$gridded
  
  if (is.null(prior.prob)){
    pr <- rep(0, length(gr))
    for (j in 1:length(gr)){
      pr[j] <- length(which(x.group==gr[j]))
    }
    pr <- pr/nrow(x)
    fhat.list$prior.prob <- pr
  } else {
    fhat.list$prior.prob <- prior.prob
  }
  fhat.list$x.group <- x.group
  
  class(fhat.list) <- "kda"
  return (fhat.list)
}

kda_rogue <- function(x, x.group, prior.prob=NULL,kde.flag,Hs)
{
  supp=3.7
  compute.cont <- TRUE
  approx.cont <- TRUE
  
  eval.points <- x
  if (!is.factor(x.group)) {x.group <- factor(x.group)}
  ##gr <- sort(unique(x.group))
  gr <- levels(x.group)
  m <- length(gr)
  
  ## default values 
  ksd <- ks.defaults(x=x)
  d <- ksd$d; n <- ksd$n; w <- ksd$w
  binned <- ksd$binned
  bgridsize <- ksd$bgridsize
  gridsize <- ksd$gridsize
  
  if (!is.null(prior.prob)){
    if (!(identical(all.equal(sum(prior.prob), 1), TRUE))){
      stop("Sum of weights not equal to 1")
    }
  }
  if(missing(Hs)){
    Hs <- Hkda(x=x, x.group=x.group, bw="plugin", binned=default.bflag(d=d, n=n))
  }
  ## Compute KDA on grid
  if (kde.flag){
    fhat.list <- kda.nd(x=x, x.group=x.group, Hs=Hs, prior.prob=prior.prob, gridsize=gridsize, supp=supp, binned=binned, bgridsize=bgridsize,compute.cont=compute.cont, approx.cont=approx.cont,eval.points=eval.points,w=w)
  }## Compute KDA at eval.points
  fhat <- kda.nd(x=x, x.group=x.group, Hs=Hs, prior.prob=prior.prob, gridsize=gridsize, supp=supp, binned=FALSE, bgridsize=bgridsize,eval.points=eval.points, compute.cont=compute.cont, approx.cont=approx.cont,w=w)
  fhat.wt <- matrix(0, ncol=m, nrow=nrow(eval.points))  
  
  
  for (j in 1:m){
    fhat.wt[,j] <- fhat$estimate[[j]]* fhat$prior.prob[j]
  }
  ## Assign y according largest weighted density value 
  disc.gr.temp <- apply(fhat.wt, 1, which.max)
  disc.gr <- as.factor(gr[disc.gr.temp])
  if (is.numeric(gr)) {disc.gr <- as.numeric(levels(disc.gr))[disc.gr]}
  
  if (kde.flag) {
    fhat.list$x.group.estimate <- disc.gr
  } else {
    fhat.list <- disc.gr
  }
  fhat.list$type <- "kda"
  
  return(fhat.list)
}
elem <- function(i, d)
{
  elem.vec <- rep(0, d)
  elem.vec[i] <- 1
  
  return(elem.vec)
}      
grid.interp <- function(x, gridx, f)
{
  d <- ncol(x)
  n <- nrow(x)
  
  gridsize <- sapply(gridx,length)
  gind <- matrix(0, nrow=n, ncol=d)
  
  for (i in 1:n){
    for (j in 1:d)
    {
      tsum <- sum(x[i,j] >= gridx[[j]])
      if (tsum==0) {
        gind[i,j] <- 1
      } else {
        gind[i,j] <- tsum
      }
    }
  }
  for (j in 1:d) {gind[gind[,j]>=gridsize[j],j] <- gridsize[j]-1}
  
  bperm <- list()
  for (j in 1:d) {bperm[[j]] <- elem(1,2)}
  binary.perm <- as.matrix(expand.grid(bperm))
  colnames(binary.perm) <- NULL
  
  gind.list <- list()
  fx <- rep(0, length=n)
  for (i in 1:n)
  {
    gind.list[[i]] <- matrix(gind[i,], nrow=2^d, ncol=d, byrow=TRUE) + binary.perm
    w <- matrix(0, nrow=2^d, ncol=d)
    gridw <- matrix(0, nrow=2^d, ncol=d)
    for (j in 1:d)
    {
      gind.list[[i]][,j][gind.list[[i]][,j]>=gridsize[j]] <- gridsize[j]
      gridw[,j] <- gridx[[j]][gind.list[[i]][,j]]
    }
    w <- abs(matrix(as.numeric(x[i,]), nrow=2^d, ncol=d, byrow=TRUE) - gridw)
    w <- apply(w, 1, prod)
    ##w <- apply(abs(matrix(as.numeric(x[i,]), nrow=2^d, ncol=d, byrow=TRUE) - gridw), 1, prod)
    ##w <- 1/apply(abs(sweep(gridw, 2, x[i,])), 1, prod)
    ##w[w>1e5] <- 1e5
    w <- w/sum(w)
    ##fx[i] <- sum(w*f[gind.list[[i]]])
    fx[i] <- sum(w*f[gind.list[[i]][2^d:1,]])
  }
  
  
  return(fx)
}

predict.kda <- function(object, x)
{
  fhat <- object
  m <- length(fhat$prior.prob)
  if (is.vector(fhat$x[[1]])) {
    n <- length(x) 
  } else if (is.vector(x)){
    n <- 1
  } else {
    n <- nrow(x)
  }
  fhat.temp <- matrix(0, ncol=m, nrow=n)
  for (j in 1:m)
  {    
    fhat.temp[,j] <- fhat$prior.prob[j]*grid.interp(x=x, gridx=fhat$eval.points, f=fhat$estimate[[j]])
  }
  est.group <- apply(fhat.temp, 1, which.max)
  est.group <- unlist(est.group)
  ##est.group <- unique(fhat$x.group)[est.group]
  est.group <- as.factor(sort(unique(fhat$x.group))[est.group])
  return(est.group)
}

compare.kda.cv <- function(x, x.group)
{ 
  # multi-dimensional   
  n <- nrow(x)
  d <- ncol(x)
  H <- Hkda(x, x.group)
  
  ## classify data x using KDA rules based on x itself
  ## kda.group <- kda(x, x.group, Hs=H, y=x, prior.prob=prior.prob)
  ## comp <- compare(x.group, kda.group)
  
  gr <- sort(unique(x.group)) 
  kda.cv.gr <- x.group
  
  for (i in 1:n)
  {
    H.mod <- H
    ### find group that x[i] belongs to 
    ind <- which(x.group[i]==gr)
    indx <- x.group==gr[ind]
    indx[i] <- FALSE
    
    
    kda.cv.gr[i] <- kda_rogue(x[-i,], x.group[-i], kde.flag=FALSE, Hs = H)[i]
  }
  kda.cv.gr <- unlist(kda.cv.gr)
  return(compare(x.group, as.numeric(kda.cv.gr))) 
}