dat <- data.frame(x0=c(1:100), x1 = col_a + mean(pred), x2 = rep(mean(pred), 100) )
data <- dat


intersects <- function(x1, x2) {
  seg1 <- which(!!diff(x1 > x2))     # location of first point in crossing segments
  above <- x2[seg1] > x1[seg1]       # which curve is above prior to crossing
  slope1 <- x1[seg1+1] - x1[seg1]
  slope2 <- x2[seg1+1] - x2[seg1]
  x <- seg1 + ((x2[seg1] - x1[seg1]) / (slope1 - slope2))
  y <- x1[seg1] + slope1*(x - seg1)
  data.frame(x=x, y=y, pindex=seg1, pabove=(1:2)[above+1L]) 
  # pabove is greater curve prior to crossing
}

fillColor <- function(data, addLines=TRUE) {
  ## Find points of intersections
  ints <- intersects(data[,2], data[,3]) # because the first column is for Dates
  intervals <- findInterval(1:nrow(data), c(0, ints$x))
  
  ## 1) Make plot for horsepower
  par(mar = c(4.2, 2.2, 0.7, 0.4)) # for overleaf 620x498 
  matplot(data, type="n", col=2:3, lty=1, lwd=4, xaxt='n', ylim=c(0,80),  xlab="Horsepower [hp]",  cex.axis=1.5, font.main = 1.5, cex.lab=1.5)

  #for horsepower plot:
  faktor = (grid_a[100] - grid_a[1])/100
   axis(1, at=c(100-250/faktor , 100-200/faktor , 100-150/faktor , 100-100/faktor, 100-50/faktor, 100), labels = c(150, 200, 250, 300, 350, 400),
        cex.axis=1.5, font.main = 1.5, cex.lab=1.5)
 
   ## 2) Make plot for weight
   par(mar = c(4.2, 2.2, 0.7, 0.4)) # für overleaf 620x498 
   matplot(data, type="n", col=2:3, lty=1, lwd=4, xaxt='n', ylim=c(0,80),  xlab="Weight [lbs]",  cex.axis=1.5, font.main = 1.5, cex.lab=1.5)
   
   faktor = (grid_b[100] - grid_b[1])/100
   axis(1, at=c(1 , 100-3000/faktor , 100-2000/faktor, 100-1000/faktor, 100), labels = c(2000, 3000, 4000, 5000, 6000),
        cex.axis=1.5, font.main = 1.5, cex.lab=1.5)
   
   ## 3) Make plot for length
   par(mar = c(4.2, 2.2, 0.7, 0.4)) # für overleaf 620x498 
   matplot(data, type="n", col=2:3, lty=1, lwd=4, xaxt='n', ylim=c(0,80),  xlab="Length [in]",  cex.axis=1.5, font.main = 1.5, cex.lab=1.5)
   
   faktor = (grid_c[100] - grid_c[1])/100
   axis(1, at=c(100-80/faktor, 100-60/faktor , 100-40/faktor, 100-20/faktor, 100), labels = c(180, 200, 220, 240, 260),
        cex.axis=1.5, font.main = 1.5, cex.lab=1.5)
   
   


  ## Draw the polygons
  for (i in (seq_along(table(intervals))-1) ) {
    xstart <- ifelse(i == 1, 1, ints$x[i-1])
    ystart <- ifelse(i == 1, data[1,2], ints$y[i-1])
    xend <- ints$x[i]
    yend <- ints$y[i]
    x <- seq(nrow(data))[intervals == i]
    polygon(c(xstart, x, xend, rev(x)), c(ystart, data[x,2], yend, rev(data[x,3])),
            col=ints$pabove[i]%%2+2)
  }
  
  # add end of plot
  
  xstart <- ints[dim(ints)[1],1]
  ystart <- ints[dim(ints)[1],2]
  xend <- nrow(data)
  yend <- data[dim(data)[1],2]
  x <- seq(nrow(data))[intervals == max(intervals)]
  polygon(c(xstart, x, xend, rev(x)), c(ystart, data[x,2], yend, rev(data[x,3])),
          col=ints[dim(ints)[1]-1,4]%%2+2)
  
  ## Add lines for curves
  #if (addLines)
  #  invisible(lapply(1:2, function(x) lines(seq(nrow(data)), data[,x], col=x%%2+2, lwd=2)))
  lines(y = data[,2], x = c(1:100), col="blue", lwd=3)
  lines(y = data[,3], x = c(1:100), lwd=3)
  
}

## Plot the data
fillColor(dat,FALSE)
