library(audio)
library(fields)
# FFT analysis of single time series
freq_trim <- function(x, n)
{
y<-fft(x)
y[n:(length(y)-(n-1))] <- 0
y <- Re(fft(y, inverse=T)/length(y))
op <- par(mfrow=c(3,1), mar=c(3,4,2,2)+.1)
plot(x, type='l', main="FFT removed higher frequencies")
lines(y, col='red', lwd=3)
plot(Mod(fft(x)[1: ceiling((length(x)+1)/2) ]), type='l', ylab="FFT")
plot(Mod(fft(y)[1: ceiling((length(y)+1)/2) ]), type='l', ylab="Truncated FFT")
par(op)
}

map_cor <- function(x)
{
  image(t(cor(x))[ncol(cor(x)):1,])
}
# A function to look at a time series
eda.ts <- function (x, bands=FALSE) {
  op <- par(no.readonly = TRUE)
  par(mar=c(0,0,0,0), oma=c(1,4,2,1))
  # Compute the Ljung-Box p-values
  # (we only display them if needed, i.e., 
  # if we have any reason of
  # thinking it is white noise).
  p.min <- .05
  k <- 15
  p <- rep(NA, k)
  for (i in 1:k) {
    p[i] <- Box.test(
      x, i, type = "Ljung-Box"
    )$p.value
  }
  if( max(p)>p.min ) {
    par(mfrow=c(5,1))
  } else {
    par(mfrow=c(4,1))
  }
  if(!is.ts(x))
    x <- ts(x)
  plot(x, axes=FALSE);
  axis(2); axis(3); box(lwd=2)
  if(bands) {
    a <- time(x)
    i1 <- floor(min(a))
    i2 <- ceiling(max(a))
    y1 <- par('usr')[3]
    y2 <- par('usr')[4]
    if( par("ylog") ){
      y1 <- 10^y1
      y2 <- 10^y2
    }
    for (i in seq(from=i1, to=i2-1, by=2)) {
      polygon( c(i,i+1,i+1,i), c(y1,y1,y2,y2), 
               col='grey', border=NA )
    }
    lines(x)
  }
  acf(x, axes=FALSE)
  axis(2, las=2)
  box(lwd=2)
  mtext("ACF", side=2, line=2.5)
  pacf(x, axes=FALSE)
  axis(2, las=2)
  box(lwd=2)
  mtext("ACF", side=2, line=2.5)
  spectrum(x, col=par('fg'), log="dB", 
           main="", axes=FALSE )
  axis(2, las=2)
  box(lwd=2)
  mtext("Spectrum", side=2, line=2.5)
  abline(v=1, lty=2, lwd=2)
  abline(v=2:10, lty=3)
  abline(v=1/2:5, lty=3)
  if( max(p)>p.min ) {
    main <- 
    plot(p, type='h', ylim=c(0,1), 
         lwd=3, main="", axes=F)
    axis(2, las=2)
    box(lwd=2)
    mtext("Ljung-Box p-value", side=2, line=2.5)
    abline(h=c(0,.05),lty=3)
  }
  par(op)
}
plot.ma <- function (x, k=20, ...) {
  plot(abs(x), ...)
  a <- time(x)
  b <- predict(loess(abs(x) ~ a))
  lines(b ~ as.vector(a), col='red', lwd=3)
  k <- 20
  lines(filter(abs(x), rep(1/k,k)), col='blue', lwd=3)
}

recurrence_plot <- function (x, ...) {
  image(outer(x, x, function (a, b) abs(a-b)), ...)
  box()
}

phase_plane_plot <- function (
  x, 
  col=rainbow(length(x)-1), 
  xlab = "x", ylab = "dx/dt", 
  ...) {
  plot( x[-1], diff(x), col = col, 
        xlab = xlab, ylab = ylab, ... )
}

data_driven_time_warp <- function (y) {
  cbind(
    x = cumsum(c(0, abs(diff(y)))),
    y = y
  )
}

see.ts <- function (name, ma=NULL, ar=NULL, d=0, n=2000) {
  order=c(length(ar), d, length(ma))
  x <- arima.sim(list(ma=ma, ar=ar, order=order), n)
  op <- par(mfrow=c(4,1))
  plot(x, main=name)
  acf(x)
  pacf(x)
  spectrum(x, spans=10, col=par('fg'))
  par(op)
}

# create a movie of the EEG
map_cor <- function(x)
{
  image(t(cor(x))[ncol(cor(x)):1,])
}

movie <- function(data, name)
{
  dm <- data.matrix(data)
  for ( i in 1:nrow(data))
  {
    filename <- paste(name,i,".jpg",sep="")
    jpeg(file=filename)
    m <- matrix(dm[i,], nrow=8, ncol=8)
    image(m,zlim=c(-200,200));
    dev.off()
  }
}

plot.movie <- function(data)
{
  dm <- data.matrix(data)
  for ( i in 1:nrow(data))
  {
    m <- matrix(dm[i,], nrow=8, ncol=8)
    image(m, zlim=c(-300, 300), col=heat.colors(32));
    wait(0.1);
  }
}

# calculate dtw distance matrix 
dtw.mat <- function(data)
{
  res <- matrix(0, nrow=(ncol(data)), ncol=(ncol(data)))
  rowCount = min(4000, nrow(data))
  for( i in 1:ncol(data))
  {
    for( j in i:ncol(data))
    {
      d <- dtw(data[1:rowCount,i], data[1:rowCount,j])
      res[i,j] = d$normalizedDistance
    }
  }
  res
}

draw.maps <- function(name, sleepmap, awakemap, seizuremap)
{
  op <- par(mfrow=c(3,1))
  image.plot(x=1:64, y=1:64,z=sleepmap,xlab="", ylab="",main=paste(name, " sleep"))
  image.plot(x=1:64, y=1:64,z=awakemap,xlab="", ylab="",main=paste(name, " awake"))
  image.plot(x=1:64, y=1:64,z=seizuremap,xlab="", ylab="",main=paste(name, " seizure"))
  par(op)
}
draw.spectra <- function(sleepv, awakev, seizurev)
{
  op<-par(mfrow=c(3,1))
  spectrum(sleepv, main="Sleep")
  spectrum(awakev, main="Awake")
  spectrum(seizurev, main="Seizure")
  par(op)
}

zero.crossing <- function(v)
{
  l <- length(v);
  v[1:(l-1)] <=0 & v[2:l] >0
}

freq.analysis <- function(v)
{
  d <- diff(which(zero.crossing(v)))
  hist(420/d)
  res <- c(mean(420/d), sd(420/d))
  names(res) <- c("mean", "sd")
  res
}
# split vector into list of vectors, each of which is one wave, cluster waves according to similarity. Classify new waves according to nearest neighbours. 
# check whether classification of new waves matches their source (sleep, wake, seizure)
 
# matrix -> window_size -> list of most unusual signals in each window
# matrix -> window_size -> clusters by similarity in each window






