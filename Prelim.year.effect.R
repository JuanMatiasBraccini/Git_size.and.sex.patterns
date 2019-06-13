#Need to run after line 521

This.sp=names(Table.species)
Select.fn=function(dat)
{
  Use="NO"
  
  #Minimum number of observations overall
  N=nrow(dat)
  
  #Minimum number of observations per block
  Table.blk=with(dat,table(BLOCK))
  Table.blk=Table.blk[Table.blk>=min.block]
  
  #Minimum number of blocks
  if(N>=overall.n & length(Table.blk)>=block.n)Use="YES"
  
  return(Use)
  
}

Selected=vector('list',length(This.sp))
names(Selected)=This.sp
for(i in 1:length(This.sp))Selected[[i]]=Select.fn(subset(DATA,SPECIES==This.sp[[i]] & FL>0 & !is.na(FL) & FL<800 & !(is.na(BLOCK) | BLOCK==0 )))

ID=which(Selected=="YES")

This.sp=names(ID)


This.sp=This.sp[match(c("BW","TK","GM","WH","HZ","LG","MI","SO","TG","BT"),This.sp)]
Species.labels=c("Dusky shark","Sandbar shark","Gummy shark","Whiskery shark","Smooth hammerhead",
                 "Spinner shark","Milk shark","Spot-tail shark","Tiger shark","Blacktip sharks")

Comm.Species.labels=c("Dusky shark","Sandbar shark","Gummy shark","Whiskery shark")



#2.2. Run GLM
fn.size=function(dat)
{
  
  #remove lat-long combos with less than 10 observations
  dat$PASTED=paste(dat$LAT,dat$LONG)
  These.Lats.Longs=table(dat$PASTED)
  These.Lats.Longs=These.Lats.Longs[These.Lats.Longs>10]
  dat=subset(dat,PASTED%in%names(These.Lats.Longs))
  
  if(This.sp[i]=="PE")dat=subset(dat,!LAT==(-15))
  if(This.sp[i]=="BW")dat=subset(dat,!LONG==(121))
  if(This.sp[i]%in%c("WH","HZ","LG"))dat=subset(dat,!LONG%in%c(116,117))
  
  #put as factors
  dat$BOTDEPTH.bin=as.factor(round(dat$BOTDEPTH/10)*10)
  dat$LONG=as.factor(dat$LONG)
  dat$LAT=as.factor(dat$LAT)
  dat$year=as.factor(dat$year)
  
  #model  
  if(!i%in%simpler.model)model<- glm(log(FL)~LAT+LONG+BOTDEPTH.bin+year, data=dat, family=gaussian, maxit=500)
  if(i%in%simpler.model)model<- glm(log(FL)~LAT+BOTDEPTH.bin+year, data=dat, family=gaussian, maxit=500)
  #model<- glm(log(FL)~Mid.Lat*Mid.Long+BOTDEPTH+Method, data=dat, family=gaussian, maxit=500)
  
  #Anova
  Signifcance1=anova(model,test="Chisq")
  
  #Deviance explained
  Dev.exp=Dsquared(model,adjust=F)$d3
  
  return(list(Anova=Signifcance1,Dev.exp=Dev.exp,data=dat,model=model)) 
}


Output.GLM=vector('list',length(This.sp))
names(Output.GLM)=This.sp
simpler.model=match(c("BT","MI","PE","SO","TG","PN"),This.sp)
for(i in 1:length(Output.GLM))Output.GLM[[i]]=fn.size(subset(DATA,SPECIES==This.sp[i] & FL>0 & !is.na(FL) & FL<800 & !(is.na(BLOCK) | BLOCK==0 )))


#note: change order of terms to see if d.f. change, if so, then confounding!!


#Check fit diagnostics
#fit plots
for(i in 1:length(Output.GLM)) fn.plot.diag(Output.GLM[[i]]$model,This.sp[i])

#predictions
pred.fn=function(dat,Model,Thisvar,NAME)
{
  newdat=dat[,match(Thisvar,names(dat))]
  
  biasCorr <- Model$deviance/Model$df.residual/2
  dat$pred=exp(predict(Model,newdat)+biasCorr)
  
  plot(dat$FL,dat$pred,main=NAME,pch=19,col=2,ylab="Predicted",xlab="Observed")
  lines(dat$FL,dat$FL,lwd=2)
  return(dat)
}

ThisVar=vector('list',length(Output.GLM))
ThisVar[[1]]=ThisVar[[2]]=ThisVar[[3]]=ThisVar[[4]]=ThisVar[[5]]=ThisVar[[6]]=c("LAT","LONG","BOTDEPTH.bin","year")
ThisVar[[7]]=ThisVar[[8]]=ThisVar[[9]]=ThisVar[[10]]=ThisVar[[11]]=ThisVar[[12]]=c("LAT","BOTDEPTH.bin","year")

Store.preds=vector('list',length(Output.GLM))
par(mfcol=c(4,3),mai=c(0.4,0.4,0.12,0.1),oma=c(1,1.2,0.15,0.01))
for(i in 1:length(Output.GLM))Store.preds[[i]]=pred.fn(Output.GLM[[i]]$data,Output.GLM[[i]]$model,ThisVar[[i]],This.sp[i])


#2.3. Get term significance 
for(i in 1:length(Output.GLM))
{
  print("---------------------------------------------")
  print("Species")
  print(names(Output.GLM)[i])
  print("---------------------------------------------")
  print("Anova table")
  print(Output.GLM[[i]]$Anova)
  print("---------------------------------------------")
  print("Deviance explained by terms")
  print(round(fn.dev.exp.term(Output.GLM[[i]]),1))
  print("---------------------------------------------")
  print("Total deviance explained")
  print(round(Output.GLM[[i]]$Dev.exp),1)
}




#Coefficients
par(mfcol=c(4,3),mai=c(0.4,0.4,0.12,0.1),oma=c(1,1.2,0.15,0.01),las=1,mgp=c(1,0.7,0))
for(i in 1:length(Output.GLM))
{
  fn.plot.coef(Output.GLM[[i]]$model,unique(Output.GLM[[i]]$data$year),"year")
  mtext(Species.labels[i],3)
  
}
mtext("Year",1,outer=T)
mtext("Coefficient",2,outer=T,las=3)


