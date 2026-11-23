# Replot the corrected sensory records and compare dose-response assumptions.
# Run from this directory with: Rscript replot_sensory.R
# Requires only base R, stats, graphics and grDevices; use a CJK font for Chinese.
# The normalized fixed-plateau model is a sensitivity analysis. The theses do
# not specify their normalization rule or fitting equation, so it is not claimed
# to reproduce their undisclosed fitting method.

options(stringsAsFactors=FALSE, scipen=5)
set.seed(20261123)
dir.create('figures',showWarnings=FALSE)
dir.create('tables',showWarnings=FALSE)
read_input <- function(f) read.csv(file.path('inputs',f),fileEncoding='UTF-8',check.names=FALSE,na.strings=c('','NA'))
write_table <- function(d,f) write.csv(d,file.path('tables',f),row.names=FALSE,na='',fileEncoding='UTF-8')
comp <- read_input('compounds.csv')
attrs <- read_input('attributes.csv')
s <- read_input('single_intensity_records.csv')
sa <- read_input('single_attribute_records.csv')
mix <- read_input('mixture_attribute_records.csv')
codes <- read_input('mixture_codes.csv')
old <- read_input('previous_4PL_parameters.csv')
summary_ref <- read_input('single_workbook_summary_reference.csv')
stopifnot(nrow(s)==1674,length(unique(s$sample))==93,length(unique(s$compound))==15,
          all(s$intensity>=0 & s$intensity<=12), all(table(s$sample)==18),
          nrow(mix)==1132,length(unique(mix$code))==65,
          s$intensity[s$sample=='正丁醇(C6)' & s$subject=='S12']==12)

# Resample complete assessor trajectories, preserving paired concentrations.
B <- 2000L
means <- list()
for(cmp in comp$compound) {
    d <- s[s$compound==cmp,]
    samples <- unique(d$sample[order(d$level)])
    subjects <- sort(unique(d$subject))
    mat <- sapply(samples,function(k) {z<-d[d$sample==k,];z$intensity[match(subjects,z$subject)]})
    stopifnot(nrow(mat)==18,all(is.finite(mat)))
    bmeans <- replicate(B,colMeans(mat[sample.int(18,18,replace=TRUE),,drop=FALSE]))
    ci <- t(apply(bmeans,1,quantile,c(.025,.975)))
    m <- d[match(samples,d$sample),c('sample','compound','level','ppb','mapping_ok')]
    m$n <- 18L
    m$mean <- colMeans(mat)
    m$sd <- apply(mat,2,sd)
    m$lower <- ci[,1];m$upper <- ci[,2]
    means[[cmp]] <- m
}
means <- do.call(rbind,means);rownames(means)<-NULL
write_table(means,'single_means_and_95CI.csv')
mean_check <- merge(means[,c('sample','mean','n')],summary_ref,by='sample',all.y=TRUE)
mean_check$difference <- mean_check$mean-mean_check$workbook_mean
mean_check$unique_reference <- !duplicated(mean_check$sample) & !duplicated(mean_check$sample,fromLast=TRUE)
write_table(mean_check,'individual_means_vs_workbook_summary.csv')

freq_rows <- function(d,ids,attr_ids) {
    rows <- list();i<-0L
    for(k in ids) for(a in attr_ids) {
        z<-d[d$sample==k,];v<-z[[a]];i<-i+1L
        rows[[i]]<-data.frame(sample=k,attribute=a,n_valid=sum(!is.na(v)),count=sum(v,na.rm=TRUE),frequency=mean(v,na.rm=TRUE))
    }
    do.call(rbind,rows)
}
single_freq <- freq_rows(sa,unique(means$sample),attrs$attribute)
single_freq <- merge(single_freq,unique(sa[,c('sample','compound','level','ppb','mapping_ok')]),by='sample',sort=FALSE)
single_freq <- merge(single_freq,attrs,by='attribute',sort=FALSE)
mix_freq <- freq_rows(mix,codes$sample,attrs$attribute)
mix_freq <- merge(mix_freq,codes,by='sample',sort=FALSE)
mix_freq <- merge(mix_freq,attrs,by='attribute',sort=FALSE)
stopifnot(nrow(single_freq)==93*16,nrow(mix_freq)==65*16,
          all(single_freq$frequency>=0 & single_freq$frequency<=1),
          all(mix_freq$frequency>=0 & mix_freq$frequency<=1))
write_table(single_freq,'single_attribute_frequencies.csv')
write_table(mix_freq,'mixture_attribute_frequencies.csv')
ns <- aggregate(subject~code,mix,length);names(ns)[2]<-'n_recorded'
ns$n_planned <- 18L;ns$n_missing <- 18L-ns$n_recorded
write_table(ns,'mixture_record_coverage.csv')

# Exact previous model: bottom + (top-bottom)*logistic(s*(log10(C)-log10(EC50))).
# The conventional Hill exponent is s/log(10), not s.
p_old <- function(cmp) {
    a<-old[old$compound==cmp,]
    if(nrow(a)!=1 || !is.finite(a$EC50_ppb)) return(NULL)
    c(a$bottom,a$top,a$slope,log10(a$EC50_ppb))
}
pred4 <- function(p,x) p[1]+(p[2]-p[1])*plogis(p[3]*(x-p[4]))
pred2 <- function(p,x) 12*plogis(p[1]*(x-p[2]))
fit2 <- function(x,y) {
    starts<-expand.grid(s=c(.3,1,3),e=as.numeric(quantile(x,c(.25,.5,.75))))
    fs<-lapply(seq_len(nrow(starts)),function(i) tryCatch(nlminb(as.numeric(starts[i,]),function(p) sum((pred2(p,x)-y)^2),
          lower=c(.01,min(x)-6),upper=c(30,max(x)+6),control=list(iter.max=600,eval.max=2000,rel.tol=1e-11)),error=function(e) NULL))
    fs<-Filter(function(f)!is.null(f)&&is.finite(f$objective),fs)
    if(!length(fs))return(NULL)
    f<-fs[[which.min(sapply(fs,function(z)z$objective))]]
    list(p=f$par,sse=f$objective,convergence=f$convergence)
}
fit4 <- function(x,y,slope_max=10,start=NULL,full=TRUE) {
    lo<-c(0,0,.05,min(x)-3);hi<-c(12,12,slope_max,max(x)+3)
    starts<-list(c(min(y)*.5,11.8,1,median(x)))
    if(!is.null(start)) starts<-c(list(pmax(lo,pmin(hi,start))),starts)
    if(full) for(h in c(.4,2,7)) for(e in as.numeric(quantile(x,c(.25,.5,.75)))) starts[[length(starts)+1]]<-c(min(y)*.5,11.9,h,e)
    fs<-lapply(starts,function(p)tryCatch(nlminb(p,function(q) {
        if(q[2]<=q[1]) return(1e6+sum(q*q))
        sum((pred4(q,x)-y)^2)
    },lower=lo,upper=hi,control=list(iter.max=800,eval.max=2500,rel.tol=1e-11)),error=function(e)NULL))
    fs<-Filter(function(f)!is.null(f)&&is.finite(f$objective)&&f$par[2]>f$par[1],fs)
    if(!length(fs))return(NULL)
    f<-fs[[which.min(sapply(fs,function(z)z$objective))]]
    list(p=f$par,sse=f$objective,convergence=f$convergence)
}
fits2<-list();fits4new<-list();audit<-list()
for(cmp in comp$compound) {
    m<-means[means$compound==cmp,];m<-m[order(m$level),]
    p<-p_old(cmp)
    if(is.null(p)) {
        audit[[cmp]]<-data.frame(compound=cmp,levels=nrow(m),old_EC50=NA,old_Hill_n=NA,old_bottom=NA,old_top=NA,
          old_boot_boundary_fraction=NA,old_R2=NA,old_mean_RMSE=NA,fixed_0_12_EC50=NA,fixed_0_12_Hill_n=NA,
          fixed_0_12_R2=NA,fixed_0_12_mean_RMSE=NA,expanded_slope_4PL_EC50=NA,expanded_slope_4PL_slope=NA,
          expanded_slope_4PL_mean_RMSE=NA,R_refit_4PL_EC50=NA,R_refit_4PL_mean_RMSE=NA,old_vs_refit_max_prediction_change=NA,
          LOCO_RMSE_2PL=NA,LOCO_RMSE_4PL=NA,normalization_invariance_error=NA,
          status='Concentration mapping unresolved: display all C1-C8 observed values; withhold concentration fit')
        next
    }
    x<-log10(m$ppb);y<-m$mean;f2<-fit2(x,y);fits2[[cmp]]<-f2
    f10<-fit4(x,y,10,p);fits4new[[cmp]]<-f10
    f30<-fit4(x,y,30,f10$p)
    old_a<-old[old$compound==cmp,];sse4<-sum((pred4(p,x)-y)^2)
    stopifnot(abs((1-sse4/sum((y-mean(y))^2))-old_a$R2_mean)<1e-7)
    # Merely dividing the response and both plateaus by 12 cannot move EC50.
    pn<-p;pn[1:2]<-pn[1:2]/12
    invariance<-max(abs(pred4(p,x)/12-pred4(pn,x)))
    stopifnot(invariance<1e-12)
    cv2<-rep(NA_real_,length(y));cv4<-cv2
    for(j in seq_along(y)) {
        fj2<-fit2(x[-j],y[-j]);cv2[j]<-pred2(fj2$p,x[j])-y[j]
        if(length(y)-1>=5) {
            fj4<-fit4(x[-j],y[-j],10,p,full=TRUE)
            if(!is.null(fj4))cv4[j]<-pred4(fj4$p,x[j])-y[j]
        }
    }
    audit[[cmp]]<-data.frame(compound=cmp,levels=length(y),old_EC50=10^p[4],old_Hill_n=p[3]/log(10),
        old_bottom=p[1],old_top=p[2],old_boot_boundary_fraction=old_a$boundary_fraction,old_R2=old_a$R2_mean,
        old_mean_RMSE=sqrt(sse4/length(y)),fixed_0_12_EC50=10^f2$p[2],fixed_0_12_Hill_n=f2$p[1]/log(10),
        fixed_0_12_R2=1-f2$sse/sum((y-mean(y))^2),fixed_0_12_mean_RMSE=sqrt(f2$sse/length(y)),
        expanded_slope_4PL_EC50=10^f30$p[4],expanded_slope_4PL_slope=f30$p[3],expanded_slope_4PL_mean_RMSE=sqrt(f30$sse/length(y)),
        R_refit_4PL_EC50=10^f10$p[4],R_refit_4PL_mean_RMSE=sqrt(f10$sse/length(y)),
        old_vs_refit_max_prediction_change=max(abs(pred4(p,x)-pred4(f10$p,x))),
        LOCO_RMSE_2PL=sqrt(mean(cv2^2)),LOCO_RMSE_4PL=if(all(is.na(cv4)))NA_real_ else sqrt(mean(cv4^2,na.rm=TRUE)),
        normalization_invariance_error=invariance,status=if(old_a$boundary_fraction>.2)'Frequent bootstrap boundary contact: EC50 conditional on bounds' else 'Fewer bootstrap boundary contacts: still model dependent')
}
audit<-do.call(rbind,audit);rownames(audit)<-NULL
write_table(audit,'model_assumption_audit.csv')
saveRDS(list(means=means,audit=audit,fits2=fits2,fits4=fits4new,single_freq=single_freq,mix_freq=mix_freq),file.path('tables','R_analysis.rds'))
capture.output(sessionInfo(),file='tables/R_sessionInfo.txt')

blue<-'#246F91';orange<-'#C77837';ink<-'#283842';gray<-'#7E8B93'
palette<-c('#BE4948','#267FA6','#2D9272','#8267AB','#C79A2B','#52A5B1','#9F7156','#687784')
font_family<-Sys.getenv('SENSORY_FONT',unset='Noto Sans CJK SC')
open_pdf<-function(name,w,h) cairo_pdf(file.path('figures',paste0(name,'.pdf')),width=w,height=h,family=font_family,pointsize=10,bg='white',onefile=TRUE)
axis_conc<-function(m) {
    if(!all(m$mapping_ok==1)) return(list(x=m$level,lims=range(m$level)+c(-.3,.3),at=m$level,labels=paste0('C',m$level),mapped=FALSE))
    x<-log10(m$ppb); rr<-range(x);at<-pretty(rr,n=4);at<-at[at>=rr[1]-.05 & at<=rr[2]+.05]
    list(x=x,lims=rr+c(-.22,.22),at=at,labels=parse(text=paste0('10^',at)),mapped=TRUE)
}
curve_panel<-function(cmp,i,mode='fitted') {
    m<-means[means$compound==cmp,];m<-m[order(m$level),];ax<-axis_conc(m)
    normalize<-mode=='compare';fac<-if(normalize)12 else 1
    plot(NA,xlim=ax$lims,ylim=if(normalize)c(0,1.04) else c(0,12.4),xaxs='i',yaxs='i',axes=FALSE,xlab='',ylab='')
    abline(h=if(normalize)c(.25,.5,.75,1) else c(3,6,9,12),col='#EBEEF0',lwd=.6)
    axis(1,at=ax$at,labels=ax$labels,col=ink,col.axis=ink,lwd=.65,tck=-.025,cex.axis=.82)
    axis(2,at=if(normalize)c(0,.25,.5,.75,1) else c(0,3,6,9,12),las=1,col=ink,col.axis=ink,lwd=.65,tck=-.025,cex.axis=.85)
    box(bty='l',lwd=.7,col=ink)
    p<-if(cmp%in%names(fits4new))fits4new[[cmp]]$p else NULL
    if(ax$mapped && !is.null(p) && mode!='observed') {
        xx<-seq(min(ax$x),max(ax$x),length.out=400)
        a<-audit[audit$compound==cmp,]
        lt<-if(mode=='compare')1 else if(a$old_boot_boundary_fraction>.2)2 else 1
        lines(xx,pred4(p,xx)/fac,col=blue,lwd=1.6,lty=lt)
        if(mode=='compare') lines(xx,pred2(fits2[[cmp]]$p,xx)/12,col=orange,lwd=1.45,lty=2)
    } else lines(ax$x,m$mean/fac,col=gray,lwd=1.1)
    arrows(ax$x,m$lower/fac,ax$x,m$upper/fac,angle=90,code=3,length=.035,col=gray,lwd=.85)
    points(ax$x,m$mean/fac,pch=21,bg='white',col=ink,lwd=1.0,cex=.85)
    title(main=paste0('(',letters[i],') ',comp$cn[match(cmp,comp$compound)]),adj=0,line=.5,cex.main=1.02,font.main=1)
    if(!ax$mapped) mtext('浓度对应待核实',side=1,line=1.45,cex=.68,col=gray)
}
for(mode in c('fitted','observed','compare')) {
    name<-switch(mode,fitted='01_All15_Intensity_4PL_R',observed='02_All15_Observed_Response_R',compare='03_All15_Model_Comparison_R')
    open_pdf(name,9.8,13.4)
    par(mfrow=c(5,3),mar=c(2.25,2.6,1.85,.65),oma=c(4.2,2.3,.4,.2),mgp=c(1.3,.35,0),tcl=-.2,cex=.98)
    for(i in seq_len(nrow(comp)))curve_panel(comp$compound[i],i,mode)
    mtext('浓度（ppb，对数刻度；异丁酸乙酯按C级显示）',side=1,outer=TRUE,line=.35,cex=.95)
    mtext(if(mode=='compare')'归一化强度（平均评分 / 12）' else '气味强度（0–12）',side=2,outer=TRUE,line=.7,cex=1.0)
    if(mode=='compare')mtext('蓝实线：四参数模型；橙虚线：固定上下平台的两参数模型',side=1,outer=TRUE,line=2.1,cex=.92,col=gray)
    if(mode=='fitted')mtext('点与误差线：均值及95%区间；蓝线：4PL；虚线：边界接触较多的拟合',side=1,outer=TRUE,line=2.1,cex=.88,col=gray)
    if(mode=='observed')mtext('点与误差线：均值及95%区间；连线连接实测均值',side=1,outer=TRUE,line=2.1,cex=.95,col=gray)
    dev.off()
}

# Radar geometries and labels are drawn entirely in R.
# Values are selection frequencies, not the 0-12 intensity scores. The main
# single-compound album follows the thesis's panel-specific radial limits;
# a second album uses a common 0-1 radial limit. Neither rescales the data.
radar_frame<-function(labels,title,subtitle=NULL,cex_label=.85,rmax=1) {
    n<-length(labels);theta<-pi/2-2*pi*(0:(n-1))/n
    ux<-cos(theta);uy<-sin(theta)
    plot.new();plot.window(xlim=c(-1.48,1.48),ylim=c(-1.61,1.43),asp=1)
    ticks<-if(rmax==1)c(.25,.5,.75,1) else seq(rmax/4,rmax,length.out=4)
    for(r in ticks)polygon(r/rmax*ux,r/rmax*uy,border='#C9D2D8',col=NA,lty=3,lwd=.75)
    segments(0,0,ux,uy,col='#DEE4E8',lwd=.7)
    for(i in seq_len(n)) {
        adjx<-if(abs(ux[i])<.15).5 else if(ux[i]>0)0 else 1
        text(ux[i]*1.10,uy[i]*1.10,labels[i],cex=cex_label,col=ink,adj=c(adjx,.5),xpd=NA)
    }
    # Radial tick labels are placed along one spoke, with a white background.
    for(r in ticks) text(.035,r/rmax,sprintf('%.2g',r),pos=4,cex=.68,col=gray)
    text(0,1.39,title,cex=1.13,col=ink,font=1,xpd=NA)
    if(!is.null(subtitle))text(0,-1.53,subtitle,cex=.72,col=gray,xpd=NA)
    invisible(list(x=ux,y=uy))
}
single_radar<-function(cmp,i,common=FALSE) {
    f<-single_freq[single_freq$compound==cmp,]
    lev<-sort(unique(f$level))
    rmax<-if(common)1 else max(.4,min(1,ceiling(max(f$frequency)*5)/5))
    xy<-radar_frame(attrs$cn,paste0('(',letters[i],') ',comp$cn[match(cmp,comp$compound)]),cex_label=.87,rmax=rmax)
    for(k in lev) {
        z<-f[f$level==k,];v<-z$frequency[match(attrs$attribute,z$attribute)]
        xx<-v/rmax*xy$x;yy<-v/rmax*xy$y
        lines(c(xx,xx[1]),c(yy,yy[1]),col=palette[k],lwd=1.12)
        points(xx,yy,col=palette[k],pch=c(15,16,17,25,18,4,8,1)[k],cex=.54,lwd=.75)
    }
    legend(x=0,y=-1.26,legend=paste0('C',lev),col=palette[lev],lty=1,pch=c(15,16,17,25,18,4,8,1)[lev],
           ncol=4,bty='n',cex=.8,xjust=.5,yjust=.5,x.intersp=.55,y.intersp=1.05,seg.len=1.0)
}
single_radar_pages<-split(seq_len(nrow(comp)),ceiling(seq_len(nrow(comp))/6))
open_pdf('04_Single_Compound_Radar_All15_R',11.7,16.54)
for(g in single_radar_pages) {
    par(mfrow=c(3,2),mar=c(.2,.3,.2,.3),oma=c(1.6,.5,.5,.5),xpd=NA,cex=1.12)
    for(i in g)single_radar(comp$compound[i],i)
    if(length(g)<6) for(k in seq_len(6-length(g)))plot.new()
    mtext('径向刻度：属性选择比例；各面板上限见刻度；颜色区分浓度级',side=1,outer=TRUE,line=.35,cex=.95,col=gray)
}
dev.off()

mix_attr_order<-c(1,15:2)
mix_attrs<-attrs[mix_attr_order,]
base<-mix_freq[mix_freq$code=='M-0.5',]
base_v<-base$frequency[match(mix_attrs$attribute,base$attribute)]
mixture_radar<-function(code,large=FALSE) {
    f<-mix_freq[mix_freq$code==code,]
    n<-ns$n_recorded[match(code,ns$code)]
    xy<-radar_frame(mix_attrs$cn,code,subtitle=paste0('n = ',n),cex_label=if(large).91 else .86)
    bx<-base_v*xy$x;by<-base_v*xy$y
    polygon(bx,by,border=NA,col=adjustcolor('#8E9BA4',alpha.f=.13))
    lines(c(bx,bx[1]),c(by,by[1]),col='#8C9AA4',lwd=1.15,lty=2)
    v<-f$frequency[match(mix_attrs$attribute,f$attribute)]
    xx<-v*xy$x;yy<-v*xy$y
    polygon(xx,yy,border=NA,col=adjustcolor('#CC627A',alpha.f=.16))
    lines(c(xx,xx[1]),c(yy,yy[1]),col='#C65371',lwd=1.6)
    points(xx,yy,pch=16,col='#C65371',cex=.36)
}
mix_pages<-split(seq_len(nrow(codes)),ceiling(seq_len(nrow(codes))/12))
open_pdf('05_Mixture_Radar_All65_R',11.7,16.54)
for(g in mix_pages) {
    par(mfrow=c(4,3),mar=c(.0,.15,.05,.15),oma=c(2.0,.5,.5,.5),xpd=NA,cex=1.02)
    for(i in g)mixture_radar(codes$code[i])
    if(length(g)<12) for(k in seq_len(12-length(g)))plot.new()
    mtext('径向刻度：属性选择比例（0–1）   实线：样品   灰色虚线：R=0.5基底',side=1,outer=TRUE,line=.4,cex=.9,col=gray)
}
dev.off()
single_high<-paste0('M1-',c('L','Z1','Z2','Q','S','C'))
stopifnot(all(single_high%in%codes$code))
open_pdf('06_Mixture_Radar_Six_SingleHigh_R',11.7,16.54)
par(mfrow=c(3,2),mar=c(.2,.3,.2,.3),oma=c(1.8,.5,.5,.5),xpd=NA,cex=1.12)
for(k in single_high)mixture_radar(k,large=TRUE)
mtext('径向刻度：属性选择比例（0–1）   实线：样品   灰色虚线：R=0.5基底',side=1,outer=TRUE,line=.4,cex=.95,col=gray)
dev.off()

open_pdf('07_Single_Compound_Radar_CommonScale_R',11.7,16.54)
for(g in single_radar_pages) {
    par(mfrow=c(3,2),mar=c(.2,.3,.2,.3),oma=c(1.6,.5,.5,.5),xpd=NA,cex=1.12)
    for(i in g)single_radar(comp$compound[i],i,common=TRUE)
    if(length(g)<6)for(k in seq_len(6-length(g)))plot.new()
    mtext('径向刻度：属性选择比例；统一上限为1；颜色区分浓度级',side=1,outer=TRUE,line=.35,cex=.95,col=gray)
}
dev.off()

# Each single radar is also available as a separate editable vector figure.
dir.create('figures/single_radar_panels',showWarnings=FALSE)
for(i in seq_len(nrow(comp))) {
    cairo_pdf(file.path('figures/single_radar_panels',sprintf('single_radar_%02d.pdf',i)),width=5.6,height=5.8,family=font_family,pointsize=10,bg='white')
    par(mar=c(.1,.1,.1,.1),xpd=NA,cex=1.12)
    single_radar(comp$compound[i],i)
    dev.off()
}
write_table(data.frame(panel=letters[seq_len(nrow(comp))],compound=comp$compound,english=comp$en),'panel_index.csv')
cat('Completed in R:',R.version.string,'\n')
cat('Single compounds:',nrow(comp),'dose labels:',nrow(means),'ratings:',nrow(s),'\n')
cat('Mixture samples:',nrow(codes),'available records:',nrow(mix),'missing records:',65*18-nrow(mix),'\n')
cat('Models with frequent bootstrap boundary contact:',sum(audit$old_boot_boundary_fraction>.2,na.rm=TRUE),'\n')
cat('Maximum normalization invariance error:',max(audit$normalization_invariance_error,na.rm=TRUE),'\n')
cat('All plots drawn by R graphics/grDevices.\n')
