

########################################################################
######################### basic functions
####################################################################
####################################################################

inversevector<-function(x)
{
y<-x
n<-length(x)
y[n:1]<-x
y
}


#################################

which.is.max<-function(x)
{
	y <- seq(length(x))[x == max(x)]
	if(length(y) > 1)
		y <- sample(y, 1)
	y
}

which.is.max.matrix<-function(x)
{
       temp<-max(x)
       lines<-apply(x, 1,is.element,el=temp)
       y.row<-seq(length(lines))[lines]
       if(length(y.row)>1)
             y.row<-sample(y.row,1)
       y.col<-which.is.max(x[y.row,])
       c(y.row,y.col)
} 

####################

which.is.min<-function(x)
{
	y <- seq(length(x))[x == min(x)]
	if(length(y) > 1)
		y <- sample(y, 1)
	y
}

###################################
rmna<-function(x)
{
x[!is.na(x)]
}

###################################

min.vec<-function(x,y)
{
n<-length(x)
z<-1:n
for(i in 1:n)
 z[i]<-min(x[i], y[i])
z
}

CenterVector<-function(x)
{
  result<-x
  y<-mean(x[x!=999], na.rm=T)
  result[x!=999]<-result[x!=999]-y
  result
 }



###################################################################################################
########## part 1: cluster
############################## return the relative difference (vector of length n-1) for vector x
##################################### need to track all the merge, height, size information,
################ also the meanvalue for each tree

clust.one.array.c<-function(array, nuc, chr, cen)
{
  chset<-sort(unique(chr))
  k<-length(chset)
  N<-length(array)
  clust.result<-1:(N*7)
  M<-0  

storage.mode(N)="integer"
storage.mode(array)="double"
storage.mode(clust.result)="double"
storage.mode(chr)="integer"
storage.mode(k)="integer"
storage.mode(chset)="integer"
storage.mode(nuc)="integer"
storage.mode(cen)="integer"
storage.mode(M)="integer"  

junk=.C("clacarray", N, array, chr, k, chset, nuc, cen, result=clust.result, M=M, PACKAGE="clac")

result<-junk$result
M<-junk$M
result<-matrix(result, nrow=N, byrow=F)
clust.final<-result[1:M, ]
} ## column 1: size
  ## column 2: height
  ## column 3: meanvalue
  ## column 4,5: merge
  ## column 6: chrom index
  ## column 7: arm index, 1---left, 2---right

###########################################################################################
############# Part II: pick nodes
###########################################################################################

#######################################
chrom.arm.length<-function(ch, nucpos, centro)
{
result<-matrix(0, nrow=max(ch), ncol=2)
chset<-unique(ch)
for(i in sort(chset))
 {
   result[i,1]<-sum(ch==i & nucpos<=centro[i])
   result[i,2]<-sum(ch==i & nucpos>centro[i])
 }
result
}


########################################### estimate the background noise sd
essd<-function(clustfinal, array, chromarml,sizecut=10)
{
clustindex<-NULL
chset<-unique(clustfinal[,6])
for(i in sort(chset))
for(j in 1:2)
{
  temp<-clustfinal[,6]==i & clustfinal[,7]==j
  if(sum(temp) >1)
    {
       temptree<-list(height=clustfinal[temp,2], merge=clustfinal[temp,4:5])
       tempresult<-cutree(temptree, h=0.5)
       t<-max(tempresult)
       for(k in 1:t)
         {
           if(sum(tempresult==k)>sizecut)
                 tempresult[tempresult==k]<-0
          }
       clustindex<-c(clustindex,tempresult)
    }
}
c(sd(array[clustindex!=0]),mean(array[clustindex!=0]))
}

################################################# the cut line for lsize--vs--meanvalue
normalline<-function(clust.final)
{
size<-clust.final[,1]
lsize<-log(size)
lsize<-lsize/max(lsize)
meanv<-clust.final[,3]

xtemp<-seq(0.3, 1, 0.1)
updownquan<-NULL
xindex<-NULL
for(i in 1: (length(xtemp)-1))
 {
   flag<- lsize>=xtemp[i]-0.05 & lsize<=xtemp[i]+0.05
   if(sum(flag)>0)
    {
   pickmeanv<-meanv[flag]
   temp<-quantile(abs(pickmeanv), 0.99)
   updownquan<-cbind(updownquan, c(-temp, temp))
   xindex<-c(xindex, xtemp[i])
     }
 }

xnew<-xindex
test1<-lm(updownquan[1,]~xnew)
test2<-lm(updownquan[2,]~xnew)

return(rbind(test1$coefficients, test2$coefficients))
}

################################################# pick points based on lsize--vs--meanvalue
meansizepick<-function(clust.final, esd, basesd, linepara, thresh=1)
{
size<-clust.final[,1]
lsize<-log(size)
lsize<-lsize/max(lsize)
meanv<-clust.final[,3]

newline<-linepara*esd/basesd

type1<- meanv> (newline[2,1]+newline[2,2]*lsize) | meanv<(newline[1,1]+newline[1,2]*lsize)
cbind(abs(meanv)>thresh, type1)
}


################################################## pick points based on height--size
heisiz.score<-function(clust.final, lsi, k=2)
{
size<-clust.final[,1]
lsize<-log(size)
lsize<-lsize/max(lsize)
height<-clust.final[,2]
height+0.55/(lsi-1.001)^k*(lsize-1.001)^k
}


################################################## pick points for lsilo< alpha < lsihi
finalpick<-function(clust.final, lsilo, lsihi, esd, basesd, linepara, thresh=1, cutheight=0.55, kk=2)
{
temp<-meansizepick(clust.final, esd, basesd, linepara, thresh)
lo<-heisiz.score(clust.final, lsilo,k=kk) < cutheight 
hi<-heisiz.score(clust.final, lsihi,k=kk) < cutheight 

((lo & !hi )& temp[,2]) | temp[,1]
}

########################################################################
######################### part 3: get cluster from the picked notes
########################################################################

select.clust.c<-function(clustfinal, selectnode, ch.s, nuc.s, centro)
{
clactree.s<-clustfinal[, 4:7]
M.s<-nrow(clactree.s)
N.s<-length(ch.s)
chset<-sort(unique(ch.s))
k<-length(chset)
seleaf<-rep(0, N.s)
senode<-selectnode
cc<-0

storage.mode(M.s)="integer"
storage.mode(clactree.s)="integer"
storage.mode(N.s)="integer"
storage.mode(ch.s)="integer"
storage.mode(k)="integer"
storage.mode(chset)="integer"
storage.mode(nuc.s)="integer"
storage.mode(centro)="integer"
storage.mode(senode)="integer"
storage.mode(seleaf)="integer"
storage.mode(cc)="integer"

junk<-.C("selectleaf", M.s, clactree.s, N.s, ch.s, k, chset, nuc.s, centro, senode, seleaf=seleaf, cc=cc, PACKAGE="clac")

junk$seleaf
}

######################################################### calculate the region mean
makemean.f<-function(sample, clustindex, ch.s, nuc.s ,centro)
{
N<-length(sample)
result<-rep(0,N)
seleafs<-clustindex

chset<-sort(unique(ch.s))
k<-length(chset)

storage.mode(sample)="double"
storage.mode(result)="double"

storage.mode(seleafs)="integer"
storage.mode(N)="integer"
storage.mode(ch.s)="integer"
storage.mode(k)="integer"
storage.mode(chset)="integer"
storage.mode(nuc.s)="integer"
storage.mode(centro)="integer"

junk.mean<-.C("makemean", N, sample, seleafs, ch.s, chset, k, nuc.s, centro, result=result, PACKAGE="clac")
junk.mean$result
}

#########################################################################
################# prepare fdr
interval.indicator<-function(clust.total)
{
k<-ncol(clust.total)
result<-clust.total
for(i in (k-1):1)
  {
   temp<-apply(as.matrix(result[,(i+1):k]), 1, sum)
   result[temp!=0, i]<- 0
   }
result
}
 

########################################################################
######################### part 4:  deal the data
########################################################################
##### the .Fortran function take "999" as missing value
clac.smootharray.f<-function(NORMAL, CANCER, ch, size=5)
{
NORMAL<-as.matrix(NORMAL)
WHOLEdata<-cbind(NORMAL, CANCER)
m<-ncol(NORMAL)
chset<-unique(ch)
WHOLEdata[is.na(WHOLEdata)]<-999
WHOLEdata.sm<-NULL
p<-ncol(WHOLEdata)
storage.mode(p) <- "integer"
storage.mode(size)<-"integer"
#dyn.load("/home/wp57/boss/winter2004/CLACcode/CLAC.so")
for(j in sort(chset)){
   ntemp<-sum(ch==j)
   if(ntemp>1)
     {
      temp<-WHOLEdata[ch==j, ]
      result<-temp
      storage.mode(ntemp)<-"integer"
      storage.mode(temp) <- "double" 
      storage.mode(result) <- "double"
      junk <- .Fortran("avesmooth",
			ntemp,
			p,
                        size, temp, result=result, PACKAGE="clac")
      y<-junk$result
      }
    if(ntemp==1)
      y<-WHOLEdata[ch==j,]

    WHOLEdata.sm<-rbind(WHOLEdata.sm, y)
  }
CANCER.sm<-WHOLEdata.sm[,-(1:m)]
NORMAL.sm<-WHOLEdata.sm[, 1:m]
return(list(CANCER.sm=as.matrix(CANCER.sm), NORMAL.sm=as.matrix(NORMAL.sm)))
}

####################### work on NORMAL to prepare the fdr
####### index[i]=0:  female/female; index[i]=1: male/female
clac.normal.fdrprepare<-function(NORMAL.sm, index, kcut=c(seq(0.27, 0.42, 0.03),0.47,0.50,1), chr, nucposi, centro)
{
array<-NORMAL.sm
array[array==999]<-0
array[is.na(array)]<-0
m<-ncol(array)
chromarml<-chrom.arm.length(chr, nucposi, centro)

for(i in 1:m)
  if(index[i]==1)
  {
   normme<-mean(array[chr<23,i])
   array[,i]<-array[,i]-normme
   }
 
#### do clust to each array
normal.total1<-clust.one.array.c(array[,1], nucposi, chr, centro)

n<-nrow(normal.total1)
L<-n*m*7
clust.total<-1:L
dim(clust.total)<-c(n,7,m) 
clust.total[,,1]<-normal.total1

if(m>1)
 for(i in 2:m)
   clust.total[,,i]<-clust.one.array.c(array[,i], nucposi, chr, centro)

#### combine all normal array
dist<-NULL
for(i in 1:m)
{
temp<-clust.total[,,i]
if(index[i]==1)
  { 
   temp<-clust.total[clust.total[,6,i]<23, ,i]
  }
dist<-rbind(dist,temp)
}

normsd<-sd(as.vector(array[chr<23,]), na.rm=T)
norm.cutline<-normalline(dist)

cutheight=0.55

#### pick on each array
alpha<-NULL
N<-length(chr)
for(j in 1:m)
{
normal.clust.result1<-matrix(0, nrow=length(chr), ncol=length(kcut)-1)
normal.sd<-sd(array[chr<23, j])
for(i in 1:(length(kcut)-1))
  {
    normal.pick<- finalpick(clust.total[,,j], kcut[i], kcut[i+1], normal.sd, normsd, norm.cutline, thresh=min(1, 3*normal.sd), kk=2)
    if(index[j]==1)
       normal.pick<-normal.pick & (clust.total[,6,j]<23)    
    normal.clust.result1[,i]<-select.clust.c(clust.total[,,j], normal.pick, chr, nucposi, centro)
  }
normal.interval.clust<-interval.indicator(normal.clust.result1)
alpha.interval<-apply(normal.interval.clust!=0, 2, sum)
alpha.final<-cumsum(inversevector(alpha.interval))
alpha<-cbind(alpha,alpha.final)
}
alpha<-as.matrix(alpha)
alpharesult<-apply(alpha, 1, sum)/(N*m-sum(chr>=23)*sum(index==1))
return(list(normsd=normsd, norm.cutline=norm.cutline, alpha=inversevector(alpharesult), kcut=kcut, chromarml=chromarml))
}

###############################################################################
########################### for one cancer array
clac.array<-function(array, fdr0, normfdr, chr, nucposi, centro)
{
normsd<-normfdr$normsd
alpha<-normfdr$alpha
kcut<-normfdr$kcut
norm.cutline<-normfdr$norm.cutline
#chromarml<-normfdr$chromarml

kk<-length(kcut)

########## throw away the missing gene
index<-is.na(array) | array==999
currentarray<-array[!index]
rawsd<-sd(currentarray)
ch<-chr[!index]
nucpos<-nucposi[!index]
N<-length(currentarray)

chromarml<-chrom.arm.length(ch, nucpos, centro)

##### do clustering on each chrom arm
currentclust<-clust.one.array.c(currentarray, nucpos, ch, centro)

##### fix bias
esd.array<-essd(currentclust, currentarray, chromarml,sizecut=31)
currentclust[,3]<-currentclust[,3]+esd.array[2] ######## not straitforward to use "+" here...think!
               
CANCER.clust.interval<-1:(N*(kk-1))
dim(CANCER.clust.interval)<-c(N, kk-1)  ### clust indicator

############ pick out nodes
for(t in 1:(kk-1))
  {
   array.pick<-finalpick(currentclust, kcut[t], 1, esd.array[1], normsd, norm.cutline, thresh=min(1, 3*rawsd))
   CANCER.clust.interval[, t]<-select.clust.c(currentclust, array.pick, ch, nucpos, centro)
  }

############### calculate FDR
      beta<-apply(CANCER.clust.interval!=0, 2, sum)
      beta.final<-beta/N
      fdrm<-alpha/beta.final
      fdrm[is.na(fdrm)]<-0
 ################# select the cut with FDR closest to fdr0
  t<-which.is.min(abs(fdrm-fdr0))
  fdr.final<-fdrm[t]
  alpha.final<-alpha[t]
  temp<-makemean.f(currentarray, CANCER.clust.interval[, t], ch, nucpos,centro)

  result<-array
  result[]<-0
  result[!index]<-temp
  return(c(fdr.final, alpha.final, result))
}

#########################################################
#################   part 5: plot
###########################################################

####################################################################### fix fdr, only color the called genes
##################################    while gray the left       
clac.plarray<-function(sample, clustindex, chr,nucposi,centr, graylevel=0.5)
  {
    n=length(chr)
    Ch=max(chr)
    ran=range(nucposi)
    plot(nucposi,rep(1,length(nucposi)),type="n",axes=F,ylim=c(0,max(chr)+1),
    xlim=c(ran[1],ran[2]),ylab="",xlab="")

upcol<-gray(graylevel)
downcol<-gray(graylevel)

chset<-unique(chr)
for(j in sort(chset)){
   jp=Ch-j+1
   nuc=nucposi[chr==j]
   y=sample[chr==j]
   indexy<-clustindex[chr==j]

   y[is.na(y)]<-0
   y[y==999]<-0
   yposi<-y
   yposi[y<0]<-0
   ynega<-y
   ynega[y>0]<-0

   pick<-(1:length(indexy))[indexy!=0]
   npick<-(1:length(indexy))[indexy==0]

  ### plot unpicked one
  if(length(npick)>0)
  {
   segments(nuc[npick],jp,nuc[npick],jp+yposi[npick],col=upcol)
   segments(nuc[npick],jp,nuc[npick],jp+ynega[npick],col=downcol)
  }
  ### plot picked one
  if(length(pick)>0)
  {
   segments(nuc[pick],jp,nuc[pick],jp+yposi[pick],col=2)
   segments(nuc[pick],jp,nuc[pick],jp+ynega[pick],col=3)
   }
   xlim=range(0,max(nuc))
   segments(xlim[1],jp,xlim[2],jp)
  

   text(-5000000,jp,labels=j,cex=.7)
   segments(centr[j],jp+.2,centr[j],jp-.2,col=1)
   }
}


#######################################
clac.plarray.lines<-function(sample, clustindex, chr,nucposi,centr, graylevel=0.9)
  {
    n=length(chr)
    Ch=max(chr)
    ran=range(nucposi)
    plot(nucposi,rep(1,length(nucposi)),type="n",axes=F,ylim=c(0,max(chr)+1),
    xlim=c(ran[1],ran[2]),ylab="",xlab="")

upcol<-gray(0.5)
downcol<-gray(0.5)

con<-log(10)
sample<-sample/con *log(2)

chset<-unique(chr)

###########
## First, draw background scales
for(j in c(sort(chset), Ch+1)){
   jp=Ch-j+1
   nuc<-nucposi[chr==j | chr==j-1]
   xlim=range(0,max(nuc))
   for(copy in 2:10)
   segments(xlim[1], jp+log(copy)/con, xlim[2], jp+log(copy)/con, col=gray(graylevel))
}
##########
## Second, draw data
for(j in sort(chset)){
   jp=Ch-j+1
   nuc=nucposi[chr==j]
   y=sample[chr==j]
   indexy<-clustindex[chr==j]

   y[is.na(y)]<-0
   y[y==999]<-0
   yposi<-y
   yposi[y<0]<-0
   ynega<-y
   ynega[y>0]<-0

   pick<-(1:length(indexy))[indexy!=0]
   npick<-(1:length(indexy))[indexy==0]

  ### plot unpicked one
   if(length(npick)>0)
    {
    segments(nuc[npick],jp,nuc[npick],jp+yposi[npick],col=upcol)
    segments(nuc[npick],jp,nuc[npick],jp+ynega[npick],col=downcol)
    }
  ### plot picked one
   if(length(pick)>0)
     {
    segments(nuc[pick],jp,nuc[pick],jp+yposi[pick],col=2)
    segments(nuc[pick],jp,nuc[pick],jp+ynega[pick],col=3)
     }
   xlim=range(0,max(nuc))
   segments(xlim[1],jp,xlim[2],jp)
   if(j<23)
      text(-5000000,jp,labels=j,cex=.7)
   if(j==23)
      text(-5000000, jp, labels="X", cex=0.7)
   if(j==24)
      text(-5000000, jp, labels="Y", cex=0.7)
   segments(centr[j],jp+0.25,centr[j],jp-0.25,col=6)
   }
}


######################################
clac.consensus<-function(sampleM, chr, nucposi,centr,graylevel=0.9, cut=0)
  {
    n<-length(chr)
    ran<-range(nucposi)
 
   M<-ncol(sampleM)    
   color<-rainbow(6*M)
   height<-M   
   chset<-unique(chr)    

####### for legend location, if there are more than 20 chromosome(most probably from human genome), put legend in the middle
####### otherwise, put the legend at the bottom of the picture    
if(length(chset)>20)
{
plot(nucposi,rep(1,length(nucposi)),type="n",axes=F,ylim=c(0,max(chr)*2+1),
    xlim=c(ran[1],ran[2]),ylab="",xlab="")
}
if(length(chset)<20)
{
plot(nucposi,rep(1,length(nucposi)),type="n",axes=F,ylim=c(-7,max(chr)*2+1),
    xlim=c(ran[1],ran[2]),ylab="",xlab="")
}

Ch=max(chr)  
############
## first plot background scale lines
for(j in sort(chset)){
   jp=Ch*2-j*2+1
   nuc<-nucposi[chr==j]
   xlim=range(0,max(nuc))
   for(copy in seq(0.2, 0.8, 0.2))
   {
   segments(xlim[1], jp+copy, xlim[2], jp+copy, col=gray(graylevel))
   segments(xlim[1], jp-copy, xlim[2], jp-copy, col=gray(graylevel))
   }
   segments(xlim[1], jp+1, xlim[2], jp+1, col=gray(graylevel-0.1))
   segments(xlim[1], jp-1, xlim[2], jp-1, col=gray(graylevel-0.1))
}
############
## second plot data
for(j in sort(chset))
   {
   jp=Ch*2-j*2+1
   nuc=nucposi[chr==j]
   xlim=range(0,max(nuc))
   y=sampleM[chr==j, ]
   y[is.na(y)]=0  
   posicount<-apply(y>0, 1, sum)+1
   negacount<-apply(y<0, 1, sum)+1
   for(i in 1:length(nuc)){
      if(posicount[i]-1+negacount[i]-1>=cut)
        {
         segments(nuc[i],jp,nuc[i],jp+(posicount[i]-1)/height,col=color[M-posicount[i]])
      ######if(negacount[i]-1 >=cut)
         segments(nuc[i],jp,nuc[i],jp-(negacount[i]-1)/height,col=color[M+negacount[i]])
         }
    }
   if(j<23)
     text(-5000000,jp,labels=j,cex=.7)
   if(j==23)
     text(-5000000, jp, labels="X", cex=0.7)
   if(j==24)
     text(-5000000, jp, labels="Y", cex=0.7)
   segments(centr[j],jp+.5,centr[j],jp-.5,col=6)
   segments(xlim[1],jp,xlim[2],jp)                  
   }

###############
## Third. Plot legend

ends<-range(nucposi)
start<-ends[1]+(ends[2]-ends[1])*1/2
finish<-ends[1]+(ends[2]-ends[1])*6/7

if(length(chset)>20)
  level<-5*2
if(length(chset)<20)
  level<-(-4)

segments(start, level-2, start,level+2)
segments(start, level-2, finish, level-2)
segments(finish,level-2, finish, level+2)
segments(finish,level+2, start, level+2)

interval<-(finish-start)/40
M<-100
color<-rainbow(6*M)
height<-M   
xpoints<-seq(start+interval, finish-interval, length=2*M)

for(h in seq(0.2, 1, 0.2))
  {
   segments(xpoints[1], level+h, xpoints[2*M], level+h, col=gray(0.9))
   segments(xpoints[1], level-h, xpoints[2*M], level-h, col=gray(0.9))
  }
for(i in 1:M)
 {
  segments(xpoints[i+M], level, xpoints[i+M],level+i/height, col=color[M-i])
  segments(xpoints[M+1-i], level, xpoints[M+1-i],level-i/height, col=color[M+i])
 } 
segments(xpoints[1], level, xpoints[2*M], level, col=1)
gap<-M%/%2
yh<-level-3
text(xpoints[1], yh, 100)
text(xpoints[gap], yh, 50)
text((xpoints[M]+xpoints[M+1])/2, yh, 0)
text(xpoints[M+gap], yh, 50)
text(xpoints[2*M], yh, 100)
text(start-100, level+3, "Percent of Samples")
}

#####################################################
plot.allarray<-function(arrayraw, arrayregion, fdr, arrayname, chr=ch, nucposi=nucpos, centr=centro)
{
    n<-length(chr)
    m<-ncol(arrayraw)
    plot(1:n,rep(1,n),type="n",axes=F,ylim=c(0,2*m+1), xlim=c(1,n),ylab="",xlab="")

upcol<-"#FFD600"
downcol<-"#EBFF00"
chcol<-"#00FFFF"

    for(j in 1:m)
     {
    jp<-2*m+1-2*j
   segments(0, jp+0.5,n, jp+0.5, col=upcol)
   segments(0, jp+1,n, jp+1, col=upcol)
   segments(0, jp-0.5,n, jp-0.5, col=downcol)
    y<-arrayraw[,j]
    y[y==999]<-0
    y[is.na(y)]<-0
    segments((1:n)[y>0],jp,(1:n)[y>0],jp+y[y>0],col=upcol)
    segments((1:n)[y<0],jp,(1:n)[y<0],jp+y[y<0],col=downcol)
     
    interval<-arrayregion[,j]
    pick<-(1:n)[interval>0]
    if(sum(pick)>0)
       segments((1:n)[pick], jp, (1:n)[pick], jp+y[pick], col=2)    
    pick<-(1:n)[interval<0]
    if(sum(pick)>0)
       segments((1:n)[pick], jp, (1:n)[pick], jp+y[pick], col=3)    

   text(-n/50,jp,labels=j,cex=.7)
   text(n/50, jp+1, labels=arrayname[j], cex=0.5)
   text(n/50, jp+0.5,labels=paste("FDR=", round(fdr[j],3)), cex=0.5)
     }

   chset<-unique(chr)
   count<-chset
   for(j in sort(chset))
      count[j]<-sum(chr==j)
   count<-cumsum(count)
   segments(count, -0, count, 2*m+1, col=chcol)
   count<-c(0, count)
   for(j in sort(chset))
      text((count[j]+count[j+1])/2, 2*m+1, j, cex=0.5, col=chcol)
  }
##############################################################################
##############################################################################
######################## part 6: final: out put functions for Excel
##############################################################################
##############################################################################

clac.from.Excel<-function(column.name, SampleID, Data, windowsize, FDR, log2trans, centerColumn, chromosomeOption)
{
ch<-Data[,column.name=="CH"]
nucpos<-Data[,column.name=="NUCPOS"]
CANCER<-as.matrix(Data[,column.name=="Disease"])
NORMAL<-as.matrix(Data[,column.name=="NN" | column.name=="NMF"])
index<-column.name[column.name=="NN" | column.name=="NMF"]
index<-(index=="NMF")+0

Cancer.SampleID<-SampleID[column.name=="Disease"]

picknm<-!((ch==999) | (nucpos==999))
ch.total<-ch[picknm]
nucpos.total<-nucpos[picknm]
CANCER<-CANCER[picknm,]
NORMAL<-NORMAL[picknm,]

CANCER<-as.matrix(CANCER)
NORMAL<-as.matrix(NORMAL)

if(log2trans)
  {
    CANCER<-log(CANCER)/log(2)
    CANCER[CANCER==(log(999)/log(2))]<-999
    NORMAL<-log(NORMAL)/log(2)
    NORMAL[NORMAL==(log(999)/log(2))]<-999
   }
if(centerColumn)
  {
     CANCER<-apply(CANCER, 2, CenterVector)
     NORMAL<-apply(NORMAL, 2, CenterVector)
   }
return(list(ch=ch.total, nucpos=nucpos.total, CANCER=CANCER, NORMAL=NORMAL, index=index, fdr=FDR, picknm=picknm, windowsize=windowsize, Cancer.SampleID=Cancer.SampleID, ChromosomeOption=chromosomeOption))
}

########## ch.total, nucpos.total don't have missing value; picknm is the index of 
#################################################
clac.preparenormal.Excel<-function(inputdata)
{
NORMAL<-inputdata$NORMAL
CANCER<-inputdata$CANCER
ch.total<-inputdata$ch
nucpos.total<-inputdata$nucpos
windowsize<-inputdata$windowsize
index<-inputdata$index
targetFDR<-inputdata$fdr
picknm<-inputdata$picknm
chromosomeOption<-inputdata$ChromosomeOption

MissData<-apply(CANCER==999, 2, sum)/nrow(CANCER)
N.MissData<-apply(NORMAL==999, 2, sum)/nrow(NORMAL)

if(chromosomeOption)
 centro<-rep(1e+9, 24)
if(!chromosomeOption)
 centro<-c(123000000, 93500000, 91500000, 51000000, 47700000, 60500000, 58900000, 45000000, 49000000, 40000000,53000000, 35400000, 15000000, 15600000, 17000000, 39000000, 24000000, 16000000, 28500000, 27800000, 12200000, 11900000, 58500000, 10000000)

Smoothdata<-clac.smootharray.f(NORMAL, CANCER, ch.total, size=windowsize)
CANCER.sm<-Smoothdata$CANCER.sm
NORMAL.sm<-Smoothdata$NORMAL.sm

NORMAL.sm[NORMAL.sm==999]<-0
normal.result<-clac.normal.fdrprepare(NORMAL.sm, index=index, kcut=c(seq(0.27, 0.43, 0.02),0.47,0.50,0.55,0.6,0.65,0.7, 1), chr=ch.total, nucposi=nucpos.total, centro)

m<-ncol(CANCER.sm)
Nnum<-ncol(NORMAL.sm)
return(list(normal.result=normal.result, ch.total=ch.total, nucpos.total=nucpos.total, centro=centro, m=m, CANCER=CANCER, CANCER.sm=CANCER.sm, targetFDR=targetFDR, NORMAL.sm=NORMAL.sm, picknm=picknm, Nnum=Nnum,MissData=MissData, N.MissData=N.MissData))
}

clac.onetumorarray.Excel<-function(NormalResult, i)
{
normal.result<-NormalResult$normal.result
ch.total<-NormalResult$ch.total
nucpos.total<-NormalResult$nucpos.total
centro<-NormalResult$centro
CANCER.sm<-NormalResult$CANCER.sm
targetFDR<-NormalResult$targetFDR

clac.array(CANCER.sm[,i], fdr0=targetFDR, normfdr=normal.result, chr=ch.total, nucposi=nucpos.total, centro=centro)
}

### CLAC.result<-NormalResult$CANCER.sm
### m<-NormalResult$m
### CLAC.result[,i]<-clac.onetumorarray.Excel(NormalResult,i)

clac.total.sd<-function(i, CANCER.sm, CANCER.result)
{
x<-CANCER.sm[,i]
y<-CANCER.result[,i]
sd(x[y==0], na.rm=T)
}

##### calculate consensus FDR
consensusFDR<-function(fdr, ConsensusCount)
{
N<-length(ConsensusCount)
m<-length(fdr)
p<-median(fdr)
ccset<-sort(unique(ConsensusCount))

tempresult<-1:(max(ccset)+1)
for(i in ccset)
{
curcc<-i
nu<-1-pbinom(curcc-1, size=m, prob=p)
den<-sum(ConsensusCount>=curcc)/N
tempresult[curcc+1]<-nu/den
}
tempresult[ConsensusCount+1]
}
 
clac.finalsummary.Excel<-function(NormalResult, CLAC.result)
{
NORMAL.sm<-NormalResult$NORMAL.sm
CANCER.sm<-NormalResult$CANCER.sm
ch.total<-NormalResult$ch.total
picknm<-NormalResult$picknm
Nnum<-NormalResult$Nnum
MissData<-NormalResult$MissData
N.MissData<-NormalResult$N.MissData

fdr.final<-CLAC.result[1,]
alpha.final<-CLAC.result[2,]
CANCER.result<-CLAC.result[-c(1,2),]
CANCER.result<-as.matrix(CANCER.result)
m<-ncol(CANCER.result)

CANCER.sm[CANCER.sm==999]<-0
Nvar<-apply(as.matrix(1:m), 1, clac.total.sd, CANCER.sm=CANCER.sm, CANCER.result=CANCER.result)
arraySD<-apply(CANCER.sm, 2, sd)
N.arraySD<-apply(as.matrix(NORMAL.sm[ch.total<23, ]), 2, sd)

rownum<-length(picknm)
CANCER.total<-matrix(0, rownum, m)
CANCER.total[picknm,]<-CANCER.result

CANCER.smooth<-matrix(0, rownum, m)
CANCER.smooth[picknm, ]<-CANCER.sm

Consensus<-apply(CANCER.total!=0, 1, sum)
arrayCount<-apply(CANCER.total!=0, 2, sum)

conFDR<-consensusFDR(alpha.final, Consensus)

return(list(Region=CANCER.total, Smooth=CANCER.smooth, Consensus=cbind(Consensus,conFDR), fdr=fdr.final, Nnum=Nnum, MissData=MissData,N.MissData=N.MissData,  Nvar=Nvar, arrayCount=arrayCount, arraySD=arraySD, N.arraySD=N.arraySD))
}


clac.PlotOneArray<-function(inputdata, clac.result, Sample) 
{
ch.total<-inputdata$ch
nucpos.total<-inputdata$nucpos
cSampleID<-inputdata$Cancer.SampleID
picknm<-inputdata$picknm

Cancer.sm<-(clac.result$Smooth)[picknm, ]
Cancer.sm<-as.matrix(Cancer.sm)
Cancer.Region<-(clac.result$Region)[picknm,]
Cancer.Region<-as.matrix(Cancer.Region)
fdr<-clac.result$fdr

pick<-(1:length(cSampleID))[cSampleID==Sample]

picksample<-Cancer.sm[, pick]
# picksample[picksample==999]<-0
#if(log10scale)
#  picksample<-picksample*log(2)/log(10)

clustindex<-Cancer.Region[, pick]
curFDR<-fdr[pick]

centro<-c(123000000, 93500000, 91500000, 51000000, 47700000, 60500000, 58900000, 45000000, 49000000, 40000000,53000000, 35400000, 15000000, 15600000, 17000000, 39000000, 24000000, 16000000, 28500000, 27800000, 12200000, 11900000, 58500000, 10000000)
clac.plarray.lines(picksample, clustindex, ch.total,nucpos.total,centro, graylevel=0.5)
title(main=paste("CLAC Plot for Sample: ", Sample, "; (FDR=", round(curFDR, 3), ")", sep=""))
}

clac.PlotConsensus<-function(inputdata, clac.result, SampleIDPick)
{
ch.total<-inputdata$ch
nucpos.total<-inputdata$nucpos
cSampleID<-inputdata$Cancer.SampleID

#Cancer.sm<-clac.result$Smooth
Cancer.Region<-clac.result$Region
#fdr<-clac.result$fdr

n<-length(cSampleID)
pickIndi<-rep(0,n)
for(i in 1:n)
  if(is.element(cSampleID[i], SampleIDPick))
      pickIndi[i]<-1

if(sum(pickIndi)<1)
         return(FALSE)
  
sampleM<-as.matrix(Cancer.Region[, pickIndi==1])
centro<-c(123000000, 93500000, 91500000, 51000000, 47700000, 60500000, 58900000, 45000000, 49000000, 40000000,53000000, 35400000, 15000000, 15600000, 17000000, 39000000, 24000000, 16000000, 28500000, 27800000, 12200000, 11900000, 58500000, 10000000)
clac.consensus(sampleM, ch.total, nucpos.total,centro)
title(main=paste("CLAC Consensus Plot for ", sum(pickIndi), " Samples"))
return(TRUE)
}

#######################################################################################################################
#######################################################################################################################
################################### output function in R
#######################################################################################################################
#######################################################################################################################

clac.preparenormal.R<-function(CANCER, NORMAL, Normal.Type, chromosome.number, nucleotide.position, windowsize=5, targetFDR=0.01, chromosomeOption=FALSE, centromere=NULL)
{
ch.total<-chromosome.number
nucpos.total<-nucleotide.position

NORMAL<-as.matrix(NORMAL)
CANCER<-as.matrix(CANCER)
CANCER[is.na(CANCER)]=999
NORMAL[is.na(NORMAL)]=999

MissData<-apply(CANCER==999, 2, sum)/nrow(CANCER)
N.MissData<-apply(NORMAL==999, 2, sum)/nrow(NORMAL)

if(chromosomeOption)
 centro<-rep(1e+9, 24)
if(!chromosomeOption)
 {
 if(is.null(centromere))
   centro<-c(123000000, 93500000, 91500000, 51000000, 47700000, 60500000, 58900000, 45000000, 49000000, 40000000,53000000, 35400000, 15000000, 15600000, 17000000, 39000000, 24000000, 16000000, 28500000, 27800000, 12200000, 11900000, 58500000, 10000000)
 if(!is.null(centromere))
   centro<-centromere
  }
Smoothdata<-clac.smootharray.f(NORMAL, CANCER, ch.total, size=windowsize)
CANCER.sm<-Smoothdata$CANCER.sm
NORMAL.sm<-Smoothdata$NORMAL.sm

NORMAL.sm[NORMAL.sm==999]<-0
normal.result<-clac.normal.fdrprepare(NORMAL.sm, index=Normal.Type, kcut=c(seq(0.27, 0.43, 0.02),0.47,0.50,0.55,0.6,0.65,0.7, 1), chr=ch.total, nucposi=nucpos.total, centro)

m<-ncol(CANCER.sm)
Nnum<-ncol(NORMAL.sm)

n<-nrow(CANCER.sm)
picknm<-rep(TRUE, n)

return(list(normal.result=normal.result, ch.total=ch.total, nucpos.total=nucpos.total, centro=centro, m=m, CANCER=CANCER, CANCER.sm=CANCER.sm, targetFDR=targetFDR, NORMAL.sm=NORMAL.sm, picknm=picknm, Nnum=Nnum,MissData=MissData, N.MissData=N.MissData))
}

##############################
clac.onetumorarray.R<-function(NormalResult,i)
{
normal.result<-NormalResult$normal.result
ch.total<-NormalResult$ch.total
nucpos.total<-NormalResult$nucpos.total
centro<-NormalResult$centro
CANCER.sm<-NormalResult$CANCER.sm
targetFDR<-NormalResult$targetFDR

clac.array(CANCER.sm[,i], fdr0=targetFDR, normfdr=normal.result, chr=ch.total, nucposi=nucpos.total, centro=centro)
}

###############################
clac.tumorarray.R<-function(NormalResult, tumorarrayIndex)
{
temp<-dim(NormalResult$CANCER.sm)
CLAC.result<-matrix(0, nrow=temp[1]+2, ncol=temp[2])
for(i in tumorarrayIndex)
  CLAC.result[,i]<-clac.onetumorarray.R(NormalResult,i)
return(fdr=CLAC.result[1,], RegionMean=CLAC.result[-c(1,2), ], alpha=CLAC.result[2,], tumorarrayIndex=tumorarrayIndex)
}

################################

clac.PlotSingleArray.R<-function(i, NormalResult, clac.result, centromere=NULL, graylevel=0.9)
{
ch.total<-NormalResult$ch.total
nucpos.total<-NormalResult$nucpos.total

Cancer.sm<-NormalResult$CANCER.sm
Cancer.sm<-as.matrix(Cancer.sm)

Cancer.Region<-clac.result$RegionMean
Cancer.Region<-as.matrix(Cancer.Region)

fdr<-clac.result$fdr
tumorarrayIndex<-clac.result$tumorarrayIndex

if(is.null(centromere))
  {
  centro<-c(123000000, 93500000, 91500000, 51000000, 47700000, 60500000, 58900000, 45000000, 49000000, 40000000,53000000, 35400000, 15000000, 15600000, 17000000, 39000000, 24000000, 16000000, 28500000, 27800000, 12200000, 11900000, 58500000, 10000000)
   }
if(!is.null(centromere))
   {
    centro<-centromere
    }
clac.plarray.lines(Cancer.sm[,i], Cancer.Region[,i], ch.total,nucpos.total,centro, graylevel)
}
#####################################
clac.PlotConsensus.R<-function(clac.result, chromosome.number, nucleotide.position, sample.index, centromere=NULL,graylevel=0.9)
{
sampleM<-as.matrix((clac.result$RegionMean)[,sample.index]) 
if(is.null(centromere))
  {
  centro<-c(123000000, 93500000, 91500000, 51000000, 47700000, 60500000, 58900000, 45000000, 49000000, 40000000,53000000, 35400000, 15000000, 15600000, 17000000, 39000000, 24000000, 16000000, 28500000, 27800000, 12200000, 11900000, 58500000, 10000000)
   }
if(!is.null(centromere))
   {
    centro<-centromere
    }
clac.consensus(sampleM, chromosome.number, nucleotide.position,centro, graylevel)
}
#######################################
clac.PlotAllArray.R<-function(NormalResult, clac.result, ArrayName=NULL,centromere=NULL)
{
ch.total<-NormalResult$ch.total
nucpos.total<-NormalResult$nucpos.total

Cancer.sm<-NormalResult$CANCER.sm
Cancer.sm<-as.matrix(Cancer.sm)

Cancer.Region<-clac.result$RegionMean
Cancer.Region<-as.matrix(Cancer.Region)

fdr<-clac.result$fdr
tumorarrayIndex<-clac.result$tumorarrayIndex

if(is.null(centromere))
  {
  centro<-c(123000000, 93500000, 91500000, 51000000, 47700000, 60500000, 58900000, 45000000, 49000000, 40000000,53000000, 35400000, 15000000, 15600000, 17000000, 39000000, 24000000, 16000000, 28500000, 27800000, 12200000, 11900000, 58500000, 10000000)
   }
if(!is.null(centromere))
   {
    centro<-centromere
    }

if(is.null(ArrayName))
   {
    arrayname<-paste("Array ", 1:ncol(Cancer.sm))
    }
if(!is.null(centromere))
   {
    arrayname<-ArrayName
    }

plot.allarray(Cancer.sm, Cancer.Region, fdr, arrayname, chr=ch.total, nucposi=nucpos.total, centr=centro)
}











