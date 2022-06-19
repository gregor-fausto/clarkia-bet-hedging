####
####
## Code for conceptual figure describing seed bank data
####
####

# set font sizes
pt12 = 1
pt10 = 10/12
pt9 = 9/12
pt8 = 8/12
pt7 = 7/12
pt6 = 6/12

# scale figure dimensions
scale = .9

# - Panel for seed bag trials ----

adj.y = 5
adj.y2 = 10

# two one panel plots
tiff("products/figures/seed-bag-trials.tif",
     units='px',height = 800*2*scale, width = 800*3*scale,res=800,compression='lzw')

par(mfrow=c(1,1),mar=c(0,.25,0,0),
    oma=c(1.5,1.5,0,0), mgp=c(3,.1,0))

# times at which bags were collected
t.sample = c(0,4,4,12,16,16,24,28,28,36)/36

# plotting frame
plot(NA,
     xlim=c(0,40),ylim=c(6,40-adj.y*2),
     type='n',frame=FALSE,
     axes=FALSE)

# x axis (time)
axis(1, c(0,10,20,30,40), 
     cex.axis = pt6, pos = 5,
     las = 1,  tck=-0.01, padj = -1,
     col = 'black', col.ticks = 1, 
     cex.axis = pt6)

# y axis (experimental round)
axis(2, c(35-adj.y*2,22.5-adj.y,10), tck=-0.01, 
     labels = c("Age 0","Age 1","Age 2"), 
     col = "black", col.ticks = c(1,1,1,1,1,1,0), las =2 ,
     cex.axis=pt9)

# full time of each set of seed bags in the round
arrows(0, 35-adj.y*2, 12, 35-adj.y*2, length=0.05, angle=90, code=1,lwd=1.5)
arrows(0, 22.5-adj.y, 24, 22.5-adj.y, length=0.05, angle=90, code=1,lwd=1.5)
arrows(0, 10, 36, 10, length=0.05, angle=90, code=1,lwd=1.5)

# add observations of seedlings
segments(x0=c(4,16,28),
         y0=c(35-adj.y*2,22.5-adj.y,10),
         y1=c(39-adj.y*2,26.5-adj.y,14),
         lwd=1.5,lty='dotted')

# add points for observations at each time
points(x=c(4,4,12,16,16,24,28,28,36),
       y=c(c(35,39,35)-adj.y*2,c(22.5,26.5,22.5)-adj.y,10,14,10),
       pch=c(19,21,19,19,21,19,19,21,19),cex=1,
       bg=c("#0077bb","white","#0077bb","#0077bb","white","#0077bb","#0077bb","white","#0077bb"),
       col=c("#0077bb","#e7298a","#0077bb","#0077bb","#e7298a","#0077bb","#0077bb","#e7298a","#0077bb"))

# observations of age 0 seeds
text(x=4,34-adj.y*2-.5,expression(paste(y)[11]),cex=pt7)
text(x=12,34-adj.y*2-.5,expression(paste(y)[12]),cex=pt7)
text(x=6.5,39-adj.y*2,expression(paste(y["g,1"])),cex=pt7)

# observations of age 1 seeds
text(x=(16),21.5-adj.y-.5,expression(paste(y)[13]),cex=pt7)
text(x=(24),21.5-adj.y-.5,expression(paste(y)[14]),cex=pt7)
text(x=18.5,26.5-adj.y,expression(paste(y["g,2"])),cex=pt7)

# observations of age 2 seeds
text(x=(28),9-.5,expression(paste(y)[15]),cex=pt7)
text(x=(36),9-.5,expression(paste(y)[16]),cex=pt7)
text(x=30.5,14,expression(paste(y["g,3"])),cex=pt7)

# legend
legend("topright",
       pch = c(19,21),
       col = c('#0077bb','#e7298a'),#,'black','black'),
       lty = c(FALSE,FALSE,1,3),
       legend = c("Intact seed count", "Seedling count"),#,"Survival",'Germination'),
       cex=pt8,
       box.lty=0,
       bg=NA)

# time (months)
rect(c(0,4,12,16,24,28,36),5.5,c(4,12,16,24,28,36,40),6,
     col=c('gray50','gray90'),border=0,lwd=0)
# add time of year
text(c(0,4,12,16,24,28,36,40),6.75,
     c("O","J"),cex=pt7)
mtext("Time (months)",side=1,line=.25,cex=pt9)

dev.off()



# - Panel for lab germination and viability trials ----

tiff("products/figures/lab-trials.tif",
     units='px',height = 800*1.5*scale, width = 800*2*scale,res=800,compression='lzw')

par(mfrow=c(1,1),mar=c(0,0,0,0),
    oma=c(0,0,0,0)+.1, mgp=c(3,.1,0))

# plotting frame
plot(NA,
     xlim=c(9,21),ylim=c(.4,.64),
     type='n',frame=FALSE,
     axes=FALSE)

# germination trails
segments(x0=c(10),
         x1=c(20),
         y0=c(.6))
points(x=c(10,20),y=c(.6,.6),pch=19)

# viability trials
segments(x0=c(10),
         x1=c(20),
         y0=c(.4))
points(x=c(10,20),y=c(.4,.4),pch=19)

# connect germination and viability trials
segments(x0=20,x1=10.25,y0=.6,y1=.41,lty='dotted')
shape::Arrows(x0=20,x1=10.25,y0=.6,y1=.41, arr.type="triangle", 
              arr.width=.2,arr.length=.2,lty=0)

# label observations
text(x=9.5,.625,expression(paste(n)["g,1"]^"viab"),cex=pt7)
text(x=20.5,.625,expression(paste(y)["g,1"]^"viab"),cex=pt7)

text(x=9.5,.425,expression(paste(n)["v,1"]^"viab"),cex=pt7)
text(x=20.5,.425,expression(paste(y)["v,1"]^"viab"),cex=pt7)

dev.off()


# - Functions ----
## Exponential survival function

f.exp = function(t,lam=lambda){
  return(exp(-(t/beta)^alpha))
}

## Discrete survival function

f.discrete = function(t,t_min,t_max){
  t_step = t_max-.01
  if(t<t_step){
    val = exp(-(t_min))
  } else {
    val = exp(-(t_max))
  }
  return(val)
}

# - Simulated data ----
# this is sort of messy
# it started out that I simulated data from a seed bank with constant survivorship
# so used an exponential survival function
# then I moved to plot the survival in discrete terms, so wrote a function to translate the
# continuous to exponential
# I'm leaving it like this because the continuous can be converted to the discrete
# and this is a conceptual illustration but if I were starting over I would 
# start with a discrete simulation

# number of seeds starting the trials
n = 100
# scale parameter
lambda = 1
# translate scale parameter to terms for Weibull
beta=1/lambda
# shape parameter
alpha=1

# conditional probability of germination (all *intact* seeds have x% germination probability)
g = c(.25,.25,.25)

# viability at the end of the year for age 0, 1, and 2 seeds
v = c(.9,.8,.7)

# - Calculate age-specific probabilities and survivorship schedule ----

## sample times
t.sample = t = c(0,4,4,12,16,16,24,28,28,36)/36

# Calculate age-specific survival probability
Oct_0 = 1

Janpre_1 = f.exp(t=t.sample[2])
Jangerm_1 = g[1]
Janpost_1 = (1-g[1])
Oct_1 = f.exp(t=t.sample[4])/f.exp(t=t.sample[3])

Janpre_2 = f.exp(t=t.sample[5])/f.exp(t=t.sample[4])
Jangerm_2 = g[2]
Janpost_2 = (1-g[2])
Oct_2 = f.exp(t=t.sample[7])/f.exp(t=t.sample[6])

Janpre_3 = f.exp(t=t.sample[8])/f.exp(t=t.sample[7])
Jangerm_3 = g[3]
Janpost_3 = (1-g[3])
Oct_3 = f.exp(t=t.sample[10])/f.exp(t=t.sample[9])

# Calculate survivorship schedule
l.x=c()
l.x[1]=1
l.x[2]=Janpre_1
l.x[3]=Janpost_1*Janpre_1
l.x[4]=Oct_1*Janpost_1*Janpre_1
l.x[5]=Janpre_2*Oct_1*Janpost_1*Janpre_1
l.x[6]=Janpost_2*Janpre_2*Oct_1*Janpost_1*Janpre_1
l.x[7]=Oct_2*Janpost_2*Janpre_2*Oct_1*Janpost_1*Janpre_1
l.x[8]=Janpre_3*Oct_2*Janpost_2*Janpre_2*Oct_1*Janpost_1*Janpre_1
l.x[9]=Janpost_3*Janpre_3*Oct_2*Janpost_2*Janpre_2*Oct_1*Janpost_1*Janpre_1
l.x[10]=Oct_3*Janpost_3*Janpre_3*Oct_2*Janpost_2*Janpre_2*Oct_1*Janpost_1*Janpre_1
l.x[11]=NA

# conditional probability of germination

g.v=c()
g.v[1] = g[1]/(1-(1-v[1]^(1/3))*(1-g[1]))
g.v[2] = g[2]/(1-(1-(v[1]*(v[2]/v[1])^(1/3)))*(1-g[2]))
g.v[3] = g[3]/(1-(1-(v[2]*(v[3]/v[2])^(1/3)))*(1-g[3]))

# Calculate age-specific survival probability
Oct_0.v = 1

Janpre_1.v = f.exp(t=t.sample[2])*(g[1]+(1-g[1])*v[1]^(1/3))
Jangerm_1.v = g.v[1]
Janpost_1.v = (1-g.v[1])
Oct_1.v = (f.exp(t=t.sample[4])*v[1])/(f.exp(t=t.sample[3])*(v[1])^(1/3))

Janpre_2.v =( f.exp(t=t.sample[5])*(g[2]+(1-g[2])*v[1]*(v[2]/v[1])^(1/3)))/(f.exp(t=t.sample[4])*v[1])
Jangerm_2.v = g.v[2]
Janpost_2.v = (1-g.v[2])
Oct_2.v = (f.exp(t=t.sample[7])*v[2])/(f.exp(t=t.sample[6])*v[1]*(v[2]/v[1])^(1/3))

Janpre_3.v = (f.exp(t=t.sample[8])*(g[3]+(1-g[3])*v[2]*(v[3]/v[2])^(1/3)))/(f.exp(t=t.sample[7])*v[2])
Jangerm_3.v = g.v[3]
Janpost_3.v = (1-g.v[3])
Oct_3.v = (f.exp(t=t.sample[10])*v[3])/(f.exp(t=t.sample[9])*v[2]*(v[3]/v[2])^(1/3))

# Calculate survivorship schedule
l.x.viability=c()
l.x.viability[1]=1
l.x.viability[2]=Janpre_1.v
l.x.viability[3]=Janpost_1.v*Janpre_1.v
l.x.viability[4]=Oct_1.v*Janpost_1.v*Janpre_1.v
l.x.viability[5]=Janpre_2.v*Oct_1.v*Janpost_1.v*Janpre_1.v
l.x.viability[6]=Janpost_2.v*Janpre_2.v*Oct_1.v*Janpost_1.v*Janpre_1.v
l.x.viability[7]=Oct_2.v*Janpost_2.v*Janpre_2.v*Oct_1.v*Janpost_1.v*Janpre_1.v
l.x.viability[8]=Janpre_3.v*Oct_2.v*Janpost_2.v*Janpre_2.v*Oct_1.v*Janpost_1.v*Janpre_1.v
l.x.viability[9]=Janpost_3.v*Janpre_3.v*Oct_2.v*Janpost_2.v*Janpre_2.v*Oct_1.v*Janpost_1.v*Janpre_1.v
l.x.viability[10]=Oct_3.v*Janpost_3.v*Janpre_3.v*Oct_2.v*Janpost_2.v*Janpre_2.v*Oct_1.v*Janpost_1.v*Janpre_1.v
l.x.viability[11]=NA


# - Panel for lab trial data ----
# figure shows viability estimated in the lab in October
# and viability inferred to the intervening times, with January marked in open circle

tiff("products/figures/viability-data.tif",
     units='px',height = 800*2*scale, width = 800*2*scale,res=800,compression='lzw')

par(mfrow=c(1,1),mar=c(0,.25,0,0),
    oma=c(1.5,1.5,0,0), mgp=c(3,.1,0))
t.sample = c(0,4,4,12,16,16,24,28,28,36)/36

## Panel B
## Viability and inferred viability
plot(NA,
     xlim=c(0,40),ylim=c(.6,1),
     type='n',frame=FALSE,
     axes=FALSE)

axis(1, c(0,10,20,30,40), 
     cex.axis = pt6,
     las = 1,  tck=-0.01, padj = -1,
     col = 'black', col.ticks = 1, 
     cex.axis = pt6)

axis(2, seq(.5,1,by=.1),
     labels = seq(.5,1,by=.1),
     tck=-0.01, las=1, hadj=1.1,
     col = 'black', col.ticks = 1, 
     cex.axis = pt6)
mtext('Pr(seed is viable)',side=2,line=1,adj=.5,col='black',cex=pt9)
mtext("Time (months)",side=1,line=.5,cex=pt9)

x0=seq(0,1,by=1/3)
lines(x=x0*12,v[1]^(x0),
      lty='dotdash',col='#1b9e77')

x1=seq(1,2,by=.01)
lines(x=x1*12,v[1]*(v[2]/v[1])^(x1-1),
      lty='dotdash',col='#1b9e77')

x2=seq(2,3,by=.01)
lines(x=x2*12,v[2]*(v[3]/v[2])^(x2-2),
      lty='dotdash',col='#1b9e77')

points(c(0,12,24,36),
       c(1,v),
       pch=19,
       col='#1b9e77')
points(c(4,16,28), 
       c(v[1]^(1/3), v[1]*(v[2]/v[1])^(1/3), v[2]*(v[3]/v[2])^(1/3)),
       pch=21,bg='white',
       col='#1b9e77')

legend("topright",
       pch = c(19,21),
       col = c('#1b9e77','#1b9e77'),
       bg = c(NA,'white'),
       legend = c("Estimated viability", "Inferred viability"),
       cex=pt8,
       box.lty=0)

rect(c(0,4,12,16,24,28,36),.6,c(4,12,16,24,28,36,40),.61,
     col=c('gray50','gray90'),lwd=0,border=0)
text(c(0,4,12,16,24,28,36,40),.625,
     c("O","J"),cex=pt7)

dev.off()

# - Panel for combined data ----

#1b9e77 purple
#e7298a green
tiff("products/figures/survival.tif",
     units='px',height = 800*2*scale, width = 800*3*scale,res=800,compression='lzw')

par(mfrow=c(1,1),mar=c(0,.25,0,0),
    oma=c(1.5,1.5,0,0), mgp=c(3,.1,0))

t.sample = c(0,4,4,12,16,16,24,28,28,36)/36

# plot survival function
plot(NA,
     xlim=c(0,40),ylim=c(0,1),
     type='n',frame=FALSE,
     axes=FALSE)

# x axis
axis(1, c(0,10,20,30,40), 
     cex.axis = pt6,
     las = 1,  tck=-0.01, padj = -1,
     col = 'black', col.ticks = 1, 
     cex.axis = pt6)

# y axis
axis(2, c(0,.2,.4,.6,.8,1),
     labels = c(0,.2,.4,.6,.8,1),
     tck=-0.01, las=1, hadj=1.1,
     col = 'black', col.ticks = 1, 
     cex.axis = pt6)

mtext('Pr(remaining in seed bank)',
      side=2,line=1,adj=1.5,col='black',cex=pt9)

l.x.s=l.x[c(3,6,9)]+l.x[c(2,5,8)]*c(Jangerm_1,Jangerm_2,Jangerm_3)

# plot survivorship points on to curve
segments(x0=c(4,16,28),
         y0=c(l.x.s),y1=c(l.x[c(3,6,9)]),
         lty=c('dotted'), col = '#e7298a',lwd=1.5)

## PERSISTENCE CURVE
t1 = seq(0/36,4/36,by=.001)
vals=unlist(lapply(t1,f.discrete,t_min=min(t1),t_max=max(t1)))
lines(t1*36, vals, col = '#0077bb',lwd=1.5)

t2 = seq(4/36,12/36,by=.001)
vals=unlist(lapply(t2,f.discrete,t_min=min(t2),t_max=max(t2)))
lines(t2*36, (1-g[1])*vals, col = '#0077bb',lwd=1.5)

t3 = seq(12/36,16/36,by=.001)
vals=unlist(lapply(t3,f.discrete,t_min=min(t3),t_max=max(t3)))
lines(t3*36, (1-g[1])*vals, col = '#0077bb',lwd=1.5)

t4 = seq(16/36,24/36,by=.001)
vals=unlist(lapply(t4,f.discrete,t_min=min(t4),t_max=max(t4)))
lines(t4*36, (1-g[1])*(1-g[2])*vals, col = '#0077bb',lwd=1.5)

t5 = seq(24/36,28/36,by=.001)
vals=unlist(lapply(t5,f.discrete,t_min=min(t5),t_max=max(t5)))
lines(t5*36, (1-g[1])*(1-g[2])*vals, col = '#0077bb',lwd=1.5)

t6 = seq(28/36,36/36,by=.001)
vals=unlist(lapply(t6,f.discrete,t_min=min(t6),t_max=max(t6)))
lines(t6*36, (1-g[1])*(1-g[2])*(1-g[3])*vals, col = '#0077bb',lwd=1.5)

## VIABILITY CURVE
l.x.s=l.x.viability[c(2,5,8)]

points(t.sample[1:10]*36,l.x.viability[1:10],col='#1b9e77',pch=19,cex=.5)


legend("topright",
       col = c('#0077bb',"#e7298a",'#1b9e77'),
       lty = c(1,3,NA),
       pch = c(NA,NA,19),
       legend = c("Seed mortality",
                  "Seed germination",
                  "Probability after\nadjusting for viability"),
       cex=pt8,
       box.lty=0,inset=c(0,-.1))

mtext("Time (months)",side=1,line=.5,cex=pt9)


dev.off()



