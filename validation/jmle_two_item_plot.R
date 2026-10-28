# Plot already-retained enumeration. This script never fits or generates data.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1)
out <- normalizePath(args[1],mustWork=TRUE)
rows <- jsonlite::read_json(file.path(out,'enumeration.json'))$rows
positive <- Filter(function(r) r$true_contrast==.5,rows)
null <- Filter(function(r) r$true_contrast==0,rows)
x <- vapply(positive,function(r) r$informative_persons,numeric(1))
values <- function(data,method,key) vapply(data,function(r) r$methods[[method]][[key]],numeric(1))
colours <- c(CML='#117C73',JMLE='#B85F1D')
draw <- function() {
  par(mfrow=c(1,3),mar=c(5.1,4.4,2.7,.8),oma=c(1,0,2.5,0),cex=.88)
  for (panel in 1:3) {
    ylim <- list(c(.4,1.27),c(.62,1.25),c(0,1.05))[[panel]]
    plot(range(x),ylim,type='n',log='x',xaxt='n',xlab='Informative persons M\n(two items per person)',
         ylab=c('Mean contrast given finite MLE','RMS profile SE / sampling SD','Coverage given interval')[panel],
         main=c('Item contrast: truth = 0.5','SE calibration at null','Nominal 95% Wald coverage')[panel])
    axis(1,at=c(8,32,128,512,2048,16384),labels=c('8','32','128','512','2048','16384'))
    abline(h=pretty(ylim),col='#EEEEEE',lwd=.7)
    for (method in names(colours)) {
      data <- if (panel==2) null else positive
      key <- c('mean_given_finite','rms_se_over_sampling_sd','coverage_given_interval')[panel]
      lines(x,values(data,method,key),col=colours[method],lwd=2)
    }
    if (panel==1) {
      abline(h=.5,lty=2);legend('topright',c('CML','JMLE','True contrast'),
          col=c(colours,'black'),lty=c(1,1,2),lwd=c(2,2,1),bty='n',cex=.9)
    } else if (panel==2) {
      abline(h=1,lty=2);abline(h=1/sqrt(2),lty=3,col=colours['JMLE'])
      legend('topright',c('CML','JMLE'),col=colours,lty=1,lwd=2,bty='n',cex=.9)
    } else {
      abline(h=.95,lty=2)
      lines(x,values(null,'JMLE','coverage_given_interval'),lty=3,col=colours['JMLE'],lwd=2)
      legend('bottomleft',c('CML, delta=0.5','JMLE, delta=0.5','JMLE, delta=0'),
          col=c(colours,colours['JMLE']),lty=c(1,1,3),lwd=2,bty='n',cex=.85)
    }
  }
  mtext('Two-item Rasch oracle: exact binomial enumeration',outer=TRUE,side=3,line=.6,cex=1.1)
  mtext('Analytic reference, not application performance; finite-outcome conditioning is explicit.',outer=TRUE,side=1,line=0,cex=.8)
}
for (extension in c('pdf','svg','png')) {
  filename <- file.path(out,paste0('incidental_parameter_diagnostic.',extension))
  if (file.exists(filename)) stop('Retained plot already exists')
  if (extension=='pdf') pdf(filename,width=12,height=4.5)
  if (extension=='svg') svg(filename,width=12,height=4.5)
  if (extension=='png') png(filename,width=2040,height=765,res=170,type='cairo')
  draw();dev.off()
}
