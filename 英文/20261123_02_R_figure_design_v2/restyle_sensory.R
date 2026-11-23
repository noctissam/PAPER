# R figure design v2. Run in this directory: Rscript restyle_sensory.R
# Rendering changes only: the cached R analysis and all observations are reused.
# Recalculate the source analysis with analysis/replot_sensory.R if inputs change.
options(stringsAsFactors=FALSE, scipen=5)
dir.create('figures',showWarnings=FALSE)
for(k in c('single_radar_pages','mixture_radar_pages','single_radar_panels','single_response_panels'))
    dir.create(file.path('figures',k),showWarnings=FALSE)
inp<-function(f)read.csv(file.path('inputs',f),fileEncoding='UTF-8',check.names=FALSE,na.strings=c('','NA'))
comp<-inp('compounds.csv');attrs<-inp('attributes.csv');codes<-inp('mixture_codes.csv')
cache<-readRDS('tables/R_analysis.rds')
means<-cache$means;audit<-cache$audit;fits4<-cache$fits4;fits2<-cache$fits2
sf<-cache$single_freq;mf<-cache$mix_freq
coverage<-read.csv('tables/mixture_record_coverage.csv',fileEncoding='UTF-8')
stopifnot(nrow(comp)==15,nrow(means)==93,nrow(sf)==93*16,nrow(mf)==65*16,
          all(means$n==18),nrow(codes)==65,length(fits4)==14,length(fits2)==14)
reference_sf<-read.csv('tables/single_attribute_frequencies.csv',fileEncoding='UTF-8')
reference_mf<-read.csv('tables/mixture_attribute_frequencies.csv',fileEncoding='UTF-8')
stopifnot(max(abs(reference_sf$frequency-sf$frequency))<1e-13,
          max(abs(reference_mf$frequency-mf$frequency))<1e-13)

# Warm, ordered concentration colours replace the original categorical rainbow.
CLEVEL<-c('#DDA978','#CF8C65','#BD705A','#A95657','#91465A','#763C5A','#5A3354','#3D2946')
INK<-'#263A3D';MUTED<-'#788889';TEAL<-'#286E65';COPPER<-'#B87D4E'
GRID<-'#E3E9E6';RAY<-'#F0F3F1';OUTER<-'#CAD5CF'
FAMILY<-Sys.getenv('SENSORY_FONT',unset='Noto Sans CJK SC')
openpdf<-function(name,w,h,pointsize=10) cairo_pdf(file.path('figures',paste0(name,'.pdf')),
    width=w,height=h,family=FAMILY,pointsize=pointsize,bg='white',onefile=TRUE)
pred4<-function(p,x)p[1]+(p[2]-p[1])*plogis(p[3]*(x-p[4]))
pred2<-function(p,x)12*plogis(p[1]*(x-p[2]))
cmp_title<-function(cmp)if(cmp=='反,顺-2,6-壬二烯醛')'反,顺-2,6-\n壬二烯醛' else cmp

# Complete response grid: four columns, three full rows and a fourth row with
# three data panels plus a legend. All fifteen compounds remain visible.
axis_info<-function(m) {
    if(!all(m$mapping_ok==1))return(list(x=m$level,lims=range(m$level)+c(-.25,.25),
        ticks=m$level,labels=paste0('C',m$level),mapped=FALSE))
    x<-log10(m$ppb);rr<-range(x);step<-if(diff(rr)>5)2 else 1
    at<-seq(ceiling(rr[1]),floor(rr[2]),by=step)
    list(x=x,lims=rr+c(-.18,.18),ticks=at,labels=parse(text=paste0('10^',at)),mapped=TRUE)
}
response_panel<-function(cmp,i,mode='fitted') {
    m<-means[means$compound==cmp,];m<-m[order(m$level),]
    ax<-axis_info(m);fac<-if(mode=='compare')12 else 1
    plot(NA,xlim=ax$lims,ylim=c(0,12.4)/fac,axes=FALSE,xlab='',ylab='',xaxs='i',yaxs='i')
    abline(h=6/fac,col=GRID,lwd=.55)
    axis(1,at=ax$ticks,labels=ax$labels,lwd=.6,lwd.ticks=.55,col='#869690',col.axis=INK,tck=-.024,cex.axis=.79)
    axis(2,at=c(0,6,12)/fac,labels=c(0,6,12)/fac,lwd=.6,lwd.ticks=.55,col='#869690',col.axis=INK,
         las=1,tck=-.025,cex.axis=.83)
    p<-if(cmp%in%names(fits4))fits4[[cmp]]$p else NULL
    if(!is.null(p) && ax$mapped && mode!='observed') {
        xx<-seq(min(ax$x),max(ax$x),length.out=400)
        frequent<-audit$old_boot_boundary_fraction[match(cmp,audit$compound)]>.2
        lines(xx,pred4(p,xx)/fac,col=TEAL,lwd=1.55,lty=if(mode=='fitted' && frequent)2 else 1)
        if(mode=='compare')lines(xx,pred2(fits2[[cmp]]$p,xx)/12,col=COPPER,lwd=1.3,lty=2)
    }else lines(ax$x,m$mean/fac,col=TEAL,lwd=1.0)
    arrows(ax$x,m$lower/fac,ax$x,m$upper/fac,angle=90,code=3,length=.027,col='#9AAEAA',lwd=.75)
    points(ax$x,m$mean/fac,pch=21,bg='white',col=TEAL,cex=.78,lwd=.95)
    label<-paste0(letters[i],'  ',cmp_title(cmp))
    mtext(label,side=3,line=if(grepl('\n',label)) .05 else .45,adj=0,cex=.88,col=INK)
    if(!ax$mapped)mtext('浓度级',side=1,line=1.4,cex=.68,col=MUTED)
}
response_legend<-function(mode) {
    par(mar=c(1,1,1,0))
    plot.new();plot.window(xlim=c(0,1),ylim=c(0,1))
    if(mode=='fitted') {
        legend('center',legend=c('均值及95%区间','四参数拟合','边界接触较多'),
               col=c(TEAL,TEAL,TEAL),pch=c(21,NA,NA),pt.bg='white',lty=c(NA,1,2),
               lwd=c(.8,1.55,1.55),bty='n',cex=.8,seg.len=1.65,y.intersp=1.8)
    }else if(mode=='compare') {
        legend('center',legend=c('四参数模型','固定平台两参数模型','均值及95%区间'),
               col=c(TEAL,COPPER,TEAL),pch=c(NA,NA,21),pt.bg='white',lty=c(1,2,NA),
               lwd=c(1.5,1.3,.8),bty='n',cex=.78,seg.len=1.65,y.intersp=1.8)
    }else legend('center',legend=c('实测均值','95%均值区间'),col=c(TEAL,'#9AAEAA'),pch=c(21,NA),
                pt.bg='white',lty=c(1,1),lwd=c(1,.8),bty='n',cex=.82,y.intersp=1.8)
}
for(mode in c('fitted','observed','compare')) {
    name<-switch(mode,fitted='V2_01_All15_Intensity',observed='V2_02_All15_Observed',compare='V2_03_Model_Comparison')
    openpdf(name,7.6,8.0,pointsize=9.6)
    par(mfrow=c(4,4),mar=c(2.0,2.4,2.3,.45),oma=c(2,2.3,.3,.2),mgp=c(1.3,.35,0),cex=1)
    for(i in seq_len(nrow(comp)))response_panel(comp$compound[i],i,mode)
    response_legend(mode)
    mtext('浓度 / ppb（对数刻度；异丁酸乙酯按C级显示）',side=1,outer=TRUE,line=.55,cex=.84,col=INK)
    mtext(if(mode=='compare')'归一化强度（平均评分 / 12）' else '气味强度 / 0–12分',side=2,outer=TRUE,line=.7,cex=.9,col=INK)
    dev.off()
}

# Radar: clean circular rings, restrained spokes, percentage tick labels and
# a shared legend outside the panels. Values and label identities are unchanged.
radar_frame<-function(labels,label,rmax=1,n=NULL,small=FALSE) {
    nattr<-length(labels);theta<-pi/2-2*pi*(0:(nattr-1))/nattr
    ux<-cos(theta);uy<-sin(theta)
    plot.new();plot.window(xlim=c(-1.43,1.43),ylim=c(-1.28,1.52),asp=1,xaxs='i',yaxs='i')
    circle<-seq(0,2*pi,length.out=241)
    ticks<-if(rmax==1)c(.25,.5,.75,1) else pretty(c(0,rmax),n=4)
    ticks<-ticks[ticks>0 & ticks<=rmax+1e-12]
    for(r in ticks)lines(r/rmax*cos(circle),r/rmax*sin(circle),col=if(abs(r-rmax)<1e-8)OUTER else GRID,lwd=.6)
    segments(0,0,ux,uy,col=RAY,lwd=.45)
    for(j in seq_len(nattr)) {
        a<-if(abs(ux[j])<.13).5 else if(ux[j]>0)0 else 1
        text(ux[j]*1.075,uy[j]*1.075,labels[j],adj=c(a,.5),cex=if(small).80 else .86,col=INK,xpd=NA)
    }
    # Place ring values between the two upper spokes rather than on the fruit axis.
    ta<-pi/2-pi/(nattr*1.5)
    for(r in ticks)text(r/rmax*cos(ta),r/rmax*sin(ta),paste0(round(r*100),'%'),
        pos=4,offset=.10,cex=.65,col=MUTED,xpd=NA)
    text(-1.34,1.43,label,adj=c(0,.5),cex=if(small)1.0 else 1.05,col=INK,xpd=NA)
    if(!is.null(n))text(1.34,1.43,paste0('n=',n),adj=c(1,.5),cex=.70,col=MUTED,xpd=NA)
    invisible(list(x=ux,y=uy))
}
single_radar<-function(cmp,i,common=FALSE,small=FALSE) {
    d<-sf[sf$compound==cmp,];lev<-sort(unique(d$level))
    rmax<-if(common)1 else max(.4,min(1,ceiling(max(d$frequency)*5)/5))
    xy<-radar_frame(attrs$cn,paste0(letters[i],'  ',cmp),rmax=rmax,small=small)
    for(k in lev) {
        z<-d[d$level==k,];v<-z$frequency[match(attrs$attribute,z$attribute)]
        x<-v/rmax*xy$x;y<-v/rmax*xy$y
        lines(c(x,x[1]),c(y,y[1]),col=CLEVEL[k],lwd=1.1)
        ok<-v>0
        points(x[ok],y[ok],pch=21,col=CLEVEL[k],bg='white',cex=.36,lwd=.7)
    }
}
level_legend<-function(levels=1:8) {
    par(mar=c(0,0,0,0),cex=1)
    plot.new();plot.window(xlim=c(0,1),ylim=c(0,1))
    legend('center',legend=paste0('C',levels),col=CLEVEL[levels],lty=1,pch=21,pt.bg='white',
           horiz=TRUE,bty='n',cex=.9,lwd=1.2,pt.cex=.65,seg.len=1.6,x.intersp=.5)
}
mix_legend<-function() {
    par(mar=c(0,0,0,0),cex=1)
    plot.new();plot.window(xlim=c(0,1),ylim=c(0,1))
    legend('center',legend=c('混合样品','R=0.5基底'),col=c(TEAL,COPPER),lty=c(1,2),
           lwd=c(1.7,1.2),horiz=TRUE,bty='n',cex=.9,seg.len=2.6,x.intersp=.7)
}
footer<-function(label) {
    par(mar=c(0,0,0,0),cex=1);plot.new()
    text(.985,.5,label,adj=1,cex=.72,col=MUTED)
}
panel_layout<-function(ncol,nrow,panel_height) {
    mat<-rbind(rep(1,ncol),matrix(2:(nrow*ncol+1),nrow=nrow,byrow=TRUE),rep(nrow*ncol+2,ncol))
    layout(mat,heights=c(.42,rep(panel_height,nrow),.24))
}
single_groups<-split(seq_len(nrow(comp)),ceiling(seq_len(nrow(comp))/6))
for(common in c(FALSE,TRUE)) {
    key<-if(common)'V2_07_Single_CommonScale' else 'V2_04_Single_Radar'
    for(page in seq_along(single_groups)) {
        g<-single_groups[[page]];nr<-ceiling(length(g)/3)
        openpdf(sprintf('single_radar_pages/%s_p%02d',key,page),15.6,.66+nr*4.7)
        panel_layout(3,nr,4.7);level_legend()
        par(mar=c(.05,.15,.05,.15),cex=1)
        for(i in g)single_radar(comp$compound[i],i,common=common)
        if(length(g)<3*nr)for(k in seq_len(3*nr-length(g)))plot.new()
        footer(if(common)'属性选择率（统一100%刻度）' else '属性选择率（各面板上限见刻度）')
        dev.off()
    }
}

# Use the original workbook's clockwise attribute order for both single and
# mixed samples. This replaces the reversed mixture order of the thesis.
ma<-attrs[1:15,]
b<-mf[mf$code=='M-0.5',]
base_v<-b$frequency[match(ma$attribute,b$attribute)]
mixture_radar<-function(code,small=FALSE) {
    d<-mf[mf$code==code,]
    xy<-radar_frame(ma$cn,code,n=coverage$n_recorded[match(code,coverage$code)],small=small)
    bx<-base_v*xy$x;by<-base_v*xy$y
    v<-d$frequency[match(ma$attribute,d$attribute)]
    x<-v*xy$x;y<-v*xy$y
    polygon(x,y,border=NA,col=adjustcolor(TEAL,alpha.f=.11))
    lines(c(bx,bx[1]),c(by,by[1]),col=COPPER,lwd=1.15,lty=2)
    lines(c(x,x[1]),c(y,y[1]),col=TEAL,lwd=1.65)
    points(x[v>0],y[v>0],pch=16,col=TEAL,cex=.30)
}
mix_groups<-split(seq_len(nrow(codes)),ceiling(seq_len(nrow(codes))/12))
for(page in seq_along(mix_groups)) {
    g<-mix_groups[[page]];nr<-ceiling(length(g)/4)
    openpdf(sprintf('mixture_radar_pages/V2_05_Mixture_Radar_p%02d',page),16.2,.66+nr*3.65)
    panel_layout(4,nr,3.65);mix_legend()
    par(mar=c(.0,.13,.0,.13),cex=1)
    for(i in g)mixture_radar(codes$code[i],small=TRUE)
    if(length(g)<4*nr)for(k in seq_len(4*nr-length(g)))plot.new()
    footer('属性选择率（统一100%刻度）')
    dev.off()
}
single_high<-paste0('M1-',c('L','Z1','Z2','Q','S','C'))
openpdf('V2_06_Mixture_Six_SingleHigh',15.6,.66+2*4.7)
panel_layout(3,2,4.7);mix_legend()
par(mar=c(.05,.15,.05,.15),cex=1)
for(code in single_high)mixture_radar(code)
footer('属性选择率（统一100%刻度）');dev.off()

# Separate vectors preserve practical editability for future manuscript layouts.
for(i in seq_len(nrow(comp))) {
    cmp<-comp$compound[i];lev<-sort(unique(sf$level[sf$compound==cmp]))
    openpdf(sprintf('single_radar_panels/V2_single_radar_%02d',i),5.2,5.30)
    layout(matrix(1:2,ncol=1),heights=c(4.75,.55))
    par(mar=c(.05,.13,.05,.13),cex=1)
    single_radar(cmp,i)
    level_legend(lev);dev.off()
    openpdf(sprintf('single_response_panels/V2_single_response_%02d',i),3.3,3.0)
    par(mar=c(3.0,3.15,2.6,.6),mgp=c(1.65,.45,0),cex=1)
    response_panel(cmp,i)
    mtext(if(cmp=='异丁酸乙酯')'浓度级' else '浓度 / ppb',side=1,line=1.9,cex=.9,col=INK)
    mtext('气味强度',side=2,line=2,cex=.9,col=INK)
    dev.off()
}

write.csv(sf[,c('sample','attribute','frequency','n_valid','count')],
          'tables/V2_plotted_single_frequencies.csv',row.names=FALSE,fileEncoding='UTF-8')
write.csv(mf[,c('code','attribute','frequency','n_valid','count')],
          'tables/V2_plotted_mixture_frequencies.csv',row.names=FALSE,fileEncoding='UTF-8')
write.csv(data.frame(level=paste0('C',1:8),colour=CLEVEL),
          'tables/V2_concentration_palette.csv',row.names=FALSE)
capture.output(sessionInfo(),file='tables/V2_R_sessionInfo.txt')
cat('R figure design v2 complete.\n')
cat('Cached statistics reused unchanged: 15 compounds, 93 dose labels, 65 mixture samples.\n')
cat('All contours, circular grids, labels and legends drawn in R.\n')
cat('Concentration palette:',paste(CLEVEL,collapse=', '),'\n')
