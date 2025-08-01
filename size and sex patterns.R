#SCRIPT FOR ANALYSING SIZE FREQUENCY AND SEX RATIO PATTERNS


require(lubridate)  #for dates manipulation
library(PBSmapping)
data(worldLLhigh)
library(plotrix)
library(lme4) #mixed models
library(vegan)
#library(ReporteRs)
library(dplyr)
library(Hmisc)
library(tidyverse)
library(patchwork)
library(ggrepel)

# DATA SECTION-------------------------------------------------------------------------
fn.user=function(x1,x2)paste(x1,Sys.getenv("USERNAME"),x2,sep='/')
if(!exists('handl_OneDrive')) source(fn.user(x1='C:/Users',
                                             x2='OneDrive - Department of Primary Industries and Regional Development/Matias/Analyses/SOURCE_SCRIPTS/Git_other/handl_OneDrive.R'))

#Sharks data base  
User="Matias"
source(handl_OneDrive('Analyses/SOURCE_SCRIPTS/Git_other/Source_Shark_bio.R'))

#controls
Export.dat="YES"   #exporting size data for population dynamics
do.paper=FALSE

if(do.paper)
{
  #Southern Oscillation Index
  SOI=read.csv(handl_OneDrive("Data/SOI.1975_2013.csv"))
  
  #Mean Freo sea level
  Freo=read.csv(handl_OneDrive("Data/Freo_mean_sea_level.csv"))
  names(Freo)[3]="Freo"
  
  #Sea Surface Temperature
  Temperature=read.csv(handl_OneDrive("Data/SST.nice.format.csv"))
  
  
  #Shark zones
  library(rgdal)
  JA_Northern_Shark=readOGR(handl_OneDrive("Data/Mapping/Shark_shape_files/JA_Northern_Shark.shp"), layer="JA_Northern_Shark") 
  WA_Northern_Shark=readOGR(handl_OneDrive("Data/Mapping/Shark_shape_files/NorthCoastShark_s43.shp"), layer="NorthCoastShark_s43") 
  WA_Northern_Shark_2=readOGR(handl_OneDrive("Data/Mapping/Shark_shape_files/NorthWestCoastShark_s43.shp"), layer="NorthWestCoastShark_s43") 
  SDGDLL_zone1=readOGR(handl_OneDrive("Data/Mapping/Shark_shape_files/SDGDLL_zone1.shp"), layer="SDGDLL_zone1") 
  SDGDLL_zone2=readOGR(handl_OneDrive("Data/Mapping/Shark_shape_files/SDGDLL_zone2.shp"), layer="SDGDLL_zone2") 
  WCDGDLL=readOGR(handl_OneDrive("Data/Mapping/Shark_shape_files/WCDGDLL.shp"), layer="WCDGDLL") 
  
  Current.year=2016
}

source(handl_OneDrive('Analyses/SOURCE_SCRIPTS/Git_other/ggplot.themes.R'))  #my themes

Bathymetry_120=read.table(handl_OneDrive("Data/Mapping/get_data112_120.cgi"))
Bathymetry_138=read.table(handl_OneDrive("Data/Mapping/get_data120.05_138.cgi"))
Bathymetry=rbind(Bathymetry_120,Bathymetry_138)

# PROCEDURE SECTION-------------------------------------------------------------------------
fn.subs=function(YEAR) substr(YEAR,start=3,stop=4)


#Export data for bycatch study
write.csv(DATA,handl_OneDrive("Analyses/Ecosystem indices and multivariate/Shark-bycatch/WA.csv"),row.names=F)



#Add Freo, SOI and Temperature
if(do.paper)
{
  DATA=merge(DATA,SOI, by.x=c("Month","year"),by.y=c("Month.soi","Year.soi"),all.x=T)
  DATA=merge(DATA,Freo, by.x=c("Month","year"),by.y=c("Month","Year"),all.x=T)
  DATA=merge(DATA,Temperature, by.x=c("Month","year","Lat.round","Long.round"),by.y=c("Month","Year","Lat","Long"),all.x=T)
}

#Select Species
if(do.paper)
{
  Sks=c("AA","BSG","BT","BW","CA","CP","ES","GB","GG","GM","GN","HH","HS","HZ","LE","LG",
        "MI","MS","OF","PE","PJ","PN","SD","SI","SO","SW","TG","TK","TN","TS","WB",
        "WC","WD","WH","WP","WS","WW","BN","BU","BX","GF","GR","HG","HW","LP","NS","PC",
        "RB","RW","SC","SF","SN","WE","WG","ZE")
  DATA=subset(DATA,SPECIES%in%Sks)
}
id.scalies=grep(".T",DATA$SPECIES)
scalies=DATA[id.scalies,]
DATA=DATA[-id.scalies,]


if(do.paper)
{
  Commercial.Sks=c("BW","GM","TK","WH")
  Sks=unique(DATA$SPECIES)
  write.csv(DATA,handl_OneDrive("Analyses/Size and sex patterns/DATA.csv"),row.names=F)
}



# Predict NA FL if TL available-------------------------------------------------------------------------
LH=read.csv(handl_OneDrive('Data/Life history parameters/Life_History.csv'))
All.species.names=read.csv(handl_OneDrive("Data/Species_names_shark.only.csv"))
All.species.names=All.species.names%>%
  mutate(Name=tolower(Name))%>%
  rename(SNAME=Name)

Res.vess=c('FLIN','NAT',"HAM","HOU","RV BREAKSEA","RV Gannet","RV GANNET","RV SNIPE 2")

N.species.with.length=DATA%>%
                        filter(!BOAT%in%Res.vess)%>%
                        mutate(Length=ifelse(!is.na(FL)|!is.na(TL)|!is.na(PL),'Yes','No'))%>%
                        filter(Length=='Yes')%>%
                        group_by(COMMON_NAME,SCIENTIFIC_NAME)%>%
                        tally()%>%
                        arrange(COMMON_NAME)%>%
                        data.frame%>%
                        filter(n>50)%>%
  mutate(COMMON_NAME=ifelse(COMMON_NAME=="Bronze whaler","Copper shark",COMMON_NAME))
# Keep.species=c("angel sharks","copper shark","dusky shark","great hammerhead","grey nurse shark","gummy shark",
#                "lemon shark","milk shark","pigeye shark","sandbar shark","sawsharks",
#                "scalloped hammerhead","shortfin mako","smooth hammerhead","spinner shark",
#                "spurdogs","tiger shark","whiskery shark","wobbegongs")
Keep.species=sort(tolower(N.species.with.length$COMMON_NAME))
Keep.species=ifelse(Keep.species=='southern eagle ray',"eagle ray",
             ifelse(Keep.species=='port jackson',"port jackson shark",
             ifelse(Keep.species=='fiddler ray',"southern fiddler ray",
             ifelse(Keep.species=='eastern school shark',"school shark",
             ifelse(Keep.species=='sawsharks',"common sawshark",
             Keep.species)))))

Species.with.no.fork.but.TL=c('ZE','FR','SH')  #species where FL cannot be measured
Species.with.no.fork.but.DW=c('ER','SR')  #species where FL cannot be measured

DATA=DATA%>%
      rename(SP=SPECIES,
             LAT=Mid.Lat,
             LONG=Mid.Long)%>%
      dplyr::select(SHEET_NO,SP,FL,TL,PL,SEX,Month,year,BOAT,MESH_SIZE,Method,LAT,LONG,zone,BOTDEPTH)%>%
      left_join(All.species.names%>%
                  dplyr::select(SPECIES,SNAME,SP),by='SP')%>%
      left_join(LH%>%dplyr::select(SPECIES,a_FL.to.TL,b_FL.to.TL,Max.TL,LF_o),
                by="SPECIES")%>%
      mutate(FINYEAR=ifelse(Month>6,paste(year,"-",fn.subs(year+1),sep=""),
                            paste(year-1,"-",fn.subs(year),sep="")),
             FL=ifelse(is.na(FL),(TL-b_FL.to.TL)/a_FL.to.TL,FL),
             FL=ifelse(SP%in%Species.with.no.fork.but.TL,TL,FL),
             FL=ifelse(SP%in%Species.with.no.fork.but.DW,PL,FL))%>%
      filter(!is.na(FL))%>%
      mutate(dummy=ifelse(!is.na(Max.TL),Max.TL,1e4))%>%
      filter(FL<dummy)%>%
      mutate(SNAME=tolower(SNAME),
             SNAME=ifelse(SNAME=="common sawshark","sawsharks",SNAME),
             MESH_SIZE=ifelse(MESH_SIZE=="10\"","10",
                       ifelse(MESH_SIZE=="6\"","6",
                       ifelse(MESH_SIZE=="5\r\n5","5",
                       ifelse(MESH_SIZE=="7\"","7",
                       ifelse(MESH_SIZE=="5\"","5",
                       ifelse(MESH_SIZE=="4\"","4",
                       ifelse(MESH_SIZE=="8\"","8",
                       MESH_SIZE))))))),
             MESH_SIZE=as.numeric(MESH_SIZE))%>%
      filter(SNAME%in%Keep.species)%>%
      mutate(Size.type=ifelse(SP%in%Species.with.no.fork.but.TL,'TL',
                       ifelse(SP%in%Species.with.no.fork.but.DW,'PL',
                       "FL")))

# Extract data for pop din model -------------------------------------------------------------------------
DATA.pop.din=DATA%>%
  filter(!BOAT%in%Res.vess)

dis.species=sort(unique(DATA.pop.din$SNAME))

# Explore spatial structure -------------------------------------------------------------------------
# Spatial distribution of pop dyn data 
do.dis=FALSE
if(do.dis)
{
  Check.this.sp=with(DATA.pop.din%>%filter(!is.na(LONG))%>%
                       filter(!is.na(BOAT))%>%
                       filter(Method=='GN')%>%        #select gillnet only
                       filter(MESH_SIZE%in%c(6.5,7)),table(SPECIES))
  Check.this.sp=Check.this.sp[Check.this.sp>200]
  
  #by zone and mesh
  fn.spatio.temp.for.SS=function(SPe)
  {
    d=DATA.pop.din%>%
      filter(SPECIES==SPe)%>%
      filter(!is.na(LONG))%>%
      filter(!is.na(BOAT))%>%
      filter(Method=='GN')%>%        #select gillnet only
      filter(MESH_SIZE%in%c(6.5,7))  #select only 6.5 and 7 inch mesh data (same as fleet)
    Zn=table(d$zone)
    Zn=Zn[Zn>50]
    if(length(Zn)>0)
    {
      Zn=names(Zn)
      NM=unique(d$SNAME)
      
      for(z in 1:length(Zn))
      {
        d1=d%>%
          filter(zone==Zn[z])%>%
          mutate(BOAT_mesh=paste(BOAT,MESH_SIZE,sep='-'))%>%
          group_by(LAT,LONG,year,BOAT_mesh)%>%
          tally()
        Limx=range(d1$LONG,na.rm=T)
        Limy= range(d1$LAT,na.rm=T) 
        
        
        p=d1%>%
          ggplot(aes(LONG,LAT,size=n,col=BOAT_mesh))+
          facet_wrap(~year)+
          xlim(Limx)+
          ggtitle(Zn[z])+
          theme_PA()+theme(legend.position = 'top')+
          geom_point(alpha=0.7)+
          geom_contour(data = Bathymetry%>%filter(V1>=Limx[1] & V1<=Limx[2] & V2>=Limy[1] & V2<=Limy[2])%>%
                         rename(LONG=V1,LAT=V2)%>%mutate(n=1,BOAT=unique(d1$BOAT_mesh[1])),
                       aes(LONG, LAT, z=V3),breaks=c(-50,-100,-200,-500),linetype="solid",colour="grey70",alpha=0.6)
        print(p)
        ggsave(handl_OneDrive(paste0("Analyses/Size and sex patterns/spatial_length_comp/Pop din_spatial_",
                                     paste(NM,Zn[z],sep='_'),".tiff")),
               width = 10,height = 8,compression = "lzw")
      }
    }
  }
  for(s in 1:length(Check.this.sp)) fn.spatio.temp.for.SS(SPe=names(Check.this.sp)[s])
  
  #by zone and mesh -histogram
  fn.spatio.temp_hist.for.SS=function(SPe)
  {
    d=DATA.pop.din%>%
      filter(SPECIES==SPe)%>%
      filter(!is.na(LONG))%>%
      filter(!is.na(BOAT))%>%
      filter(!is.na(SEX))%>%
      filter(zone%in%c('West','Zone1','Zone2'))%>%
      filter(Method=='GN')%>%        #select gillnet only
      filter(MESH_SIZE%in%c(6.5,7))  #select only 6.5 and 7 inch mesh data (same as fleet)
    
    d.NSF=DATA.pop.din%>%
      filter(SPECIES==SPe)%>%
      filter(!is.na(LONG))%>%
      filter(!is.na(BOAT))%>%
      filter(!is.na(SEX))%>%
      filter(zone%in%c('North'))%>%
      filter(Method=='LL')%>%
      mutate(MESH_SIZE='LL')
    
    d=rbind(d,d.NSF)%>%
      mutate(TL=ifelse(is.na(TL),a_FL.to.TL*FL+b_FL.to.TL,TL))
    Zn=table(d$zone)
    Zn=Zn[Zn>50]
    if(length(Zn)>0)
    {
      NM=unique(d$SNAME)
      p=d%>%
        mutate(year=as.numeric(substr(FINYEAR,1,4)),
               Mesh.sex=paste(MESH_SIZE,SEX),
               TL.bin=10*floor(TL/10))%>%
        group_by(Mesh.sex,TL.bin,zone,year)%>%
        tally()%>%ungroup()
      DRoP=p%>%group_by(zone,year)%>%summarise(n=sum(n))%>%filter(n<10)%>%mutate(drop=paste(zone,year))
      p=p%>%
        mutate(drop=paste(zone,year))%>%
        filter(!drop%in%DRoP$drop)%>%
        ggplot(aes(TL.bin,n,color=Mesh.sex))+
        geom_line()+
        facet_grid(year~zone)+
        theme_PA(str.siz=9,axs.t.siz=7)+theme(legend.title=element_blank(),legend.position = 'top')+
        xlab('TL (cm)')+ guides(colour = guide_legend(nrow = 1))
      
      print(p)
      ggsave(handl_OneDrive(paste0("Analyses/Size and sex patterns/spatial_length_comp/Pop din_spatial_",NM,"_histogram.tiff")),
             width = 6,height = 7,compression = "lzw")
      
      
    }
  }
  for(s in 1:length(Check.this.sp)) fn.spatio.temp_hist.for.SS(SPe=names(Check.this.sp)[s])
  
  #by Sex
  colfunc <- colorRampPalette(c("pink",'grey', "blue"))
  fn.spatio.temp_sex.for.SS=function(SPe)
  {
    d=DATA.pop.din%>%
      filter(SPECIES==SPe)%>%
      filter(!is.na(LONG))%>%
     # filter(!is.na(BOAT))%>%
     # filter(Method=='GN')%>%        #select gillnet only
    #  filter(MESH_SIZE%in%c(6.5,7))%>%  #select only 6.5 and 7 inch mesh data (same as fleet)
      filter(!is.na(SEX))
    
    NM=unique(d$SNAME)
    d1=d%>%
      mutate(SEX=ifelse(SEX=='F','Female','Male'),
             LAT=round(LAT),
             LONG=round(LONG))%>%
      group_by(LAT,LONG,SEX)%>%
      tally()%>%
      ungroup()%>%
      spread(SEX,n,fill=0)%>%
      mutate(N=Female+Male,
             Prop.male=Male/N)%>%
      filter(N>=5)%>%
      mutate(Prop.male1=as.character(round(Prop.male,1)),
             Prop.male1=ifelse(Prop.male1==0,0.01,Prop.male1))
    
    if(nrow(d1)>10)
    {
      Limx=range(d1$LONG,na.rm=T)
      Limy= range(d1$LAT,na.rm=T) 
      Kl.rnge=sort(unique(d1$Prop.male1))
      Kls=colfunc(length(Kl.rnge))
      names(Kls)=Kl.rnge
      p=d1%>%
        mutate(Prop.male1=factor(Prop.male1,levels=Kl.rnge))%>%
        ggplot(aes(LONG,LAT,color=Prop.male1))+
        xlim(Limx)+ylim(Limy)+
        ggtitle('Proportion males')+
        theme_PA()+theme(legend.position = 'none')+
        geom_contour(data = Bathymetry%>%filter(V1>=Limx[1] & V1<=Limx[2] & V2>=Limy[1] & V2<=Limy[2])%>%
                       rename(LONG=V1,LAT=V2)%>%mutate(N=0.01,Prop.male=0.001),
                     aes(LONG, LAT, z=V3),breaks=c(-50,-100,-200,-500),linetype="solid",colour="grey70",alpha=0.6)+
        geom_point(size=5)+
        scale_color_manual(values = Kls)+
        geom_text_repel(aes(label = Prop.male1),box.padding = 0.5)
      print(p)
      ggsave(handl_OneDrive(paste0("Analyses/Size and sex patterns/spatial_sex_comp/Pop din_spatial_sex_",
                                   NM,".tiff")),
             width = 8,height = 8,compression = "lzw")
      do.zone=FALSE
      if(do.zone)
      {
        Zn=table(d$zone)
        Zn=Zn[Zn>50]
        if(length(Zn)>0)
        {
          Zn=names(Zn)
          NM=unique(d$SNAME)
          
          for(z in 1:length(Zn))
          {
            d1=d%>%
              filter(zone==Zn[z])%>%
              mutate(SEX=ifelse(SEX=='F','Female','Male'))%>%
              group_by(LAT,LONG,SEX)%>%
              tally()%>%
              ungroup()%>%
              spread(SEX,n)%>%
              mutate(N=Female+Male,
                     Prop.male=Male/N)%>%
              filter(N>=10)%>%
              mutate(Prop.male1=as.character(round(Prop.male,1)))
            
            Limx=range(d1$LONG,na.rm=T)
            Limy= range(d1$LAT,na.rm=T) 
            
            Kl.rnge=sort(unique(d1$Prop.male1))
            Kls=colfunc(length(Kl.rnge))
            names(Kls)=Kl.rnge
            p=d1%>%
              mutate(Prop.male1=factor(Prop.male1,levels=Kl.rnge))%>%
              ggplot(aes(LONG,LAT,size=Prop.male,color=Prop.male1))+
              xlim(Limx)+
              ggtitle(paste(Zn[z],'=Proportion males'))+
              theme_PA()+theme(legend.position = 'top')+
              geom_point(alpha=0.9)+
              geom_contour(data = Bathymetry%>%filter(V1>=Limx[1] & V1<=Limx[2] & V2>=Limy[1] & V2<=Limy[2])%>%
                             rename(LONG=V1,LAT=V2)%>%mutate(N=0.1,Prop.male=0.01),
                           aes(LONG, LAT, z=V3),breaks=c(-50,-100,-200,-500),linetype="solid",colour="grey70",alpha=0.6)+
              scale_color_manual(values = Kls)
            print(p)
            ggsave(handl_OneDrive(paste0("Analyses/Size and sex patterns/spatial_sex_comp/Pop din_spatial_sex_",
                                         paste(NM,Zn[z],sep='_'),".tiff")),
                   width = 10,height = 8,compression = "lzw")
          }
        }
      }
      
    }
    
  }
  for(s in 1:length(Check.this.sp)) fn.spatio.temp_sex.for.SS(SPe=names(Check.this.sp)[s])
  
  
  #Spatial patterns by month
  colfunc <- colorRampPalette(c("yellow", "blue"))
  fun.spatial.month=function(sp)
  {
    d=DATA%>%
      filter(SP==sp)%>%
      filter(FL>LF_o)%>%
      filter(!is.na(LAT))%>%
      mutate(Size=25*round(FL/25),
             LAT=round(LAT,1),
             LONG=round(LONG,1))%>%
      group_by(Month,SEX,Size,LAT,LONG)%>%
      tally()
    unik.size=sort(unique(d$Size))
    CLs=colfunc(length(unik.size))
    names(CLs)=unik.size
    p=d%>%
      mutate(Size=factor(Size,levels=unik.size))%>%
      ggplot(aes(LONG,LAT,color=Size,size=n))+
      geom_point(alpha=0.8)+
      facet_wrap(~Month)+
      scale_color_manual(values = CLs)
    return(p)
  }
  check.this=c('TK','BW')
  for(i in 1:length(check.this))
  {
    NM=check.this[i]
    fun.spatial.month(sp=NM)
    ggsave(handl_OneDrive(paste0("Analyses/Size and sex patterns/spatial_length_comp/spatial_length_by month_",NM,".tiff")),
           width = 8,height = 8,compression = "lzw")
  }
  
}


do.dis=FALSE
if(do.dis)
{
  Min.length.display=90
  max.nrow=6
  chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) 
  for(s in 1:length(dis.species))
  {
    NN=capitalize(dis.species[s])
    dd1=DATA.pop.din%>%
      filter(SNAME==dis.species[s] & !is.na(LAT))%>%
      mutate(BLK.lat=as.numeric(substr(LAT,1,3)),
             BLK.long=as.numeric(substr(LONG,1,3)),
             BLK=paste0(abs(BLK.lat),(BLK.long-100)),
             FL.group=10*(floor(FL/10)))
    Levls=expand.grid(abs(max(dd1$BLK.lat)):abs(min(dd1$BLK.lat)),((min(dd1$BLK.long)-100):(max(dd1$BLK.long)-100)))%>%
      mutate(BLK=paste0(Var1,Var2))%>%
      arrange(Var1,Var2)
    
    
    LAts=sort(unique(Levls$Var1))
    Nrow=length(LAts)
    if(length(LAts)<=max.nrow)
    {
      lat.list=list(LAts)
    }else
    {
      lat.list=chunk2(LAts,ceiling(Nrow/max.nrow))
    }
    
    for(xx in 1:length(lat.list))
    {
      print(paste('spatial length comp----',NN,"---latitudes ----",paste(lat.list[[xx]],collapse='_')))
      dd2=dd1%>%
        filter(BLK.lat%in%-lat.list[[xx]])
      if(nrow(dd2)>Min.length.display)
      {
          Levls=expand.grid(abs(max(dd2$BLK.lat)):abs(min(dd2$BLK.lat)),((min(dd2$BLK.long)-100):(max(dd2$BLK.long)-100)))%>%
          mutate(BLK=paste0(Var1,Var2))%>%
          arrange(Var1,Var2)
        Ncol=length(unique(Levls$Var2))
        dd2=dd2%>%
          filter(BLK%in%Levls$BLK)%>%
          mutate(BLK=factor(BLK,levels=Levls$BLK))
        
        p=dd2%>%
          group_by(BLK,FL.group)%>%
          tally()%>%
          ggplot(aes(FL.group,n))+
          geom_bar(stat='identity')+
          facet_wrap(~BLK,drop=FALSE,ncol=Ncol)+
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size=8,hjust=1))+
          ggtitle(NN)
        print(p)
        ggsave(handl_OneDrive(paste0("Analyses/Size and sex patterns/spatial_length_comp/",
                                     NN,"_",paste(lat.list[[xx]],collapse = '_'),".tiff")),
               width = 8,height = 8,compression = "lzw")
      }
     }
    
    
  }
}
do.dis2=FALSE
if(do.dis2)
{
  scale_values <- function(x,MIN,MAX){(x-MIN)/(MAX-MIN)}
  for(s in 1:length(dis.species))
  {
    NN=capitalize(dis.species[s])
    dd1=DATA.pop.din%>%
      filter(SNAME==dis.species[s] & !is.na(LAT))%>%
      mutate(BLK.lat=as.numeric(substr(LAT,1,3)),
             BLK.long=as.numeric(substr(LONG,1,3)),
             BLK.lat.bottom=BLK.lat-0.5,
             BLK.lat.top=BLK.lat+0.5,
             BLK.long.left=BLK.long-0.5,
             BLK.long.right=BLK.long+0.5,
             BLK=paste0(abs(BLK.lat),(BLK.long-100)),
             FL.group=10*(floor(FL/10)),
             scaled.lat.bottom=scale_values(BLK.lat.bottom,MIN=min(BLK.lat),MAX=max(BLK.lat)),
             scaled.long.left=scale_values(BLK.long.left,MIN=min(BLK.long),MAX=max(BLK.long)),
             scaled.lat.top=scale_values(BLK.lat.top,MIN=min(BLK.lat),MAX=max(BLK.lat)),
             scaled.long.right=scale_values(BLK.long.right,MIN=min(BLK.long),MAX=max(BLK.long)))
    
    p=dd1%>%ggplot(aes(LONG,LAT))
    hists=dd1%>%group_by(BLK.lat,BLK.long,BLK,scaled.long.left,scaled.lat.bottom,
                         scaled.long.right,scaled.lat.top,FL.group)%>%tally()
    
    unik.blok=sort(unique(hists$BLK))
    for(x in 1:length(unik.blok))
    {
      b=hists%>%filter(BLK==unik.blok[x])
      p=p+
        inset_element(p=b%>%ggplot(aes(FL.group,n))+geom_bar(stat='identity'),
                      left=unique(b$scaled.long.left),
                      bottom=unique(b$scaled.lat.bottom),
                      right=unique(b$scaled.long.right),
                      top=unique(b$scaled.lat.top))
    }  
    
  }
}


#Export size frequency other gears data for population dynamics modelling
if (Export.dat=="YES")
{
  hndl=handl_OneDrive('Analyses/Data_outs/')
  for(s in 1:length(dis.species))
  {
    dd1=DATA.pop.din%>%filter(SNAME==dis.species[s])
    NN=capitalize(dis.species[s])
    if(NN=="Angel sharks") NN="Australian angelshark"
    if(nrow(dd1)>0)
    {
      #Observations (TDGDLF & NSF)
      Observations=dd1%>%
        filter(!is.na(FL))%>%
        mutate(Keep=ifelse((Method=='LL' & LAT>=(-25)) | (Method=='GN'& MESH_SIZE%in%c(6.5,7) & LAT<(-25)),'Yes','No'))%>%
        filter(Keep=='Yes')
      N.obs=Observations%>%
        group_by(FINYEAR,Method,zone,SPECIES)%>%
        tally()%>%
        rename(N.observations=n)
      N.shots=Observations%>%
        distinct(SHEET_NO,.keep_all = T)%>%
        group_by(FINYEAR,Method,zone,SPECIES)%>%
        tally()%>%
        rename(N.shots=n)
      write.csv(full_join(N.shots,N.obs,by=c('FINYEAR','Method','zone','SPECIES')),
                paste(hndl,NN,'/',NN,"_Size_composition_Observations.csv",sep=''),row.names=F)
      
      
       #Data
      gn=dd1%>%filter(Method=="GN" & LAT<(-25))
      if(nrow(gn)>0)
      {
        #raw by zone
        zn=unique(gn$zone)
        for(x in 1:length(zn))
        {
          nm=zn[x]
          a=gn%>%filter(zone==zn[x])
          a_6.5=subset(a,!is.na(FL) & MESH_SIZE=="6.5",select=c(Month,FINYEAR,year,FL,SEX,Size.type))
          a_7=subset(a,!is.na(FL) & MESH_SIZE=="7",select=c(Month,FINYEAR,year,FL,SEX,Size.type))
          if(nrow(a_6.5)>0)write.csv(a_6.5,paste(hndl,'/',NN,'/',NN,"_Size_composition_",nm,".6.5.inch.raw.csv",sep=""),row.names=F)
          if(nrow(a_7)>0)write.csv(a_7,paste(hndl,'/',NN,'/',NN,"_Size_composition_",nm,".7.inch.raw.csv",sep=""),row.names=F)
        }
        
        #table of observations
        # fn.table.shots=function(dat,SP)
        # {
        #   zn=unique(dat$zone)
        #   b=vector('list',length(zn))
        #   for(x in 1:length(zn))
        #   {
        #     a=gn%>%filter(zone==zn[x] & !is.na(FL))
        #     a$Number=1
        #     Obs=aggregate(Number~FINYEAR,a,sum)
        #     a$Dup=paste(a$year,a$Month,a$SHEET_NO)
        #     bb=a[!duplicated(a$Dup),]
        #     bb$Number=1
        #     Shots=aggregate(Number~FINYEAR,bb,sum)
        #     this=merge(Obs,Shots,by="FINYEAR")
        #     names(this)[2:3]=c("N.observations","N.shots")
        #     this$Species=unique(a$SPECIES)
        #     this$Fishery="TDGDLF"
        #     this$zone=zn[x]
        #     b[[x]]=this
        #   }
        #   write.csv(do.call(rbind,b),paste(hndl,'/',SP,'/',SP,"_Size_composition_Numb_obs_size.freq.TDGDLF.csv",sep=""),row.names=F)
        # }
        # fn.table.shots(dat=gn,SP=NN)
      }
      ll=dd1%>%filter(Method=="LL" & LAT>=(-25)) 
      if(nrow(ll)>0)
      {
        ll=ll%>%dplyr::select(Month,FINYEAR,year,FL,SEX,Size.type)
        write.csv(ll,paste(hndl,NN,'/',NN,
                           "_Size_composition_NSF.LONGLINE.csv",sep=''),row.names=F)
        
      }
      dl=dd1%>%filter(Method=="DL") 
      if(nrow(dl)>0)
      {
        dl=dl%>%dplyr::select(Month,FINYEAR,year,FL,SEX,Size.type)
        write.csv(dl,paste(hndl,NN,'/',NN,
                           "_Size_composition_dropline.csv",sep=''),row.names=F)
      }
      tr=dd1%>%filter(Method=="TW") 
      if(nrow(tr)>0)
      {
        tr=tr%>%dplyr::select(Month,FINYEAR,year,FL,SEX,Size.type)
        write.csv(tr,paste(hndl,NN,'/',NN,
                           "_Size_composition_Pilbara_Trawl.csv",sep=''),row.names=F)
      }
      rm(gn,ll,tr,dl)
      
      print(paste('----Exporting length composition for ----------',NN))
    }
  }

  #Export gillnet data for mean weight analysis
  Dat.weight.an=subset(DATA,Method=="GN" & MESH_SIZE%in%c("6","6.5","7"))
  write.csv(Dat.weight.an,handl_OneDrive("Analyses/Catch and effort/Survey.weight.csv"),row.names=F)
  
}

# Do analyses for spatial size paper-------------------------------------------------------------------------
if(do.paper)
{
  Table.soak=aggregate(SOAK.TIME~Method,DATA,mean,na.rm=T)
  Table.methods=table(DATA$Method,useNA='ifany')
  
  #Add effort accordingly (hook-hours for LL, km-gn-hour for GN)  INCOMPLETE! fix some netlengths, etc
  #DATA$Effort=with(DATA,ifelse(Method=="LL",SOAK.TIME*,ifelse(Method=="GN",SOAK.TIME*NET_LENGTH,NA)))
  
  #Get mid point of block
  DATA$LAT=-as.numeric(with(DATA,ifelse(!is.na(Mid.Lat),substr(Mid.Lat,2,3),substr(BLOCK,1,2))))
  DATA$LONG=as.numeric(with(DATA,ifelse(!is.na(Mid.Long),substr(Mid.Long,1,3),substr(BLOCK,3,4))))
  
  
  #Average block depth
  Aver.depth=aggregate(BOTDEPTH~BLOCK,DATA,mean)
  names(Aver.depth)[2]="BlockMeanDepth"
  DATA=merge(DATA,Aver.depth,by="BLOCK",all.x=T)
  
  
  #Number of males and females by zone for commercial sp
  a=subset(DATA,SPECIES%in%Commercial.Sks)
  table(a$zone,a$SEX,as.character(a$SPECIES))
  rm(a)
  DATA=subset(DATA,year<=Current.year)
  
}

# Proportion of neonate duskies in commercial gillnet catch-------------------------------------------------------------------------
FL.dusky.neonate=82.5  #mid-point between max 0+ FL (93cm) and min 1+ FL (72cm) (Simpfendorfer et al 2002)
Data.dusky=subset(DATA,Method=="GN" & MESH_SIZE%in%c("6.5","7") & SPECIES=="BW")
Data.dusky$Finyear=factor(with(Data.dusky,ifelse(Month%in%1:6,
            paste(year-1,"-",year,sep=""),ifelse(Month%in%7:12,paste(year,"-",year+1,sep=""),NA))))
TTa=with(subset(Data.dusky,!is.na(FL)),table(Finyear))
names(TTa[which(TTa>20)])
Data.dusky=subset(Data.dusky,Finyear%in%names(TTa[which(TTa>20)]))
FINYEARs=levels(Data.dusky$Finyear)
fn.prop.dus=function(FL.neo, min.obs, min.yrs)
{
  a=subset(Data.dusky,!BLOCK==0)
  Data.dusky.neonates=subset(a,FL<=FL.neo)
  
  BL.yr=table(Data.dusky.neonates$BLOCK,Data.dusky.neonates$Finyear)
  BL.yr[BL.yr<min.obs]=0
  BL.yr[BL.yr>=min.obs]=1   #select blocks with at least 5 records per year with observations
  Use.blk=rowSums(BL.yr)
  Use.blk=names(Use.blk[Use.blk>=min.yrs])
  Table.prop.dusky=aggregate(Number~Finyear+BLOCK,subset(a,BLOCK%in%Use.blk),sum)
  names(Table.prop.dusky)[match("Number",names(Table.prop.dusky))]="Total.n"
  Table.prop.dusky.neo=aggregate(Number~Finyear+BLOCK,subset(Data.dusky.neonates,BLOCK%in%Use.blk),sum)
  names(Table.prop.dusky.neo)[match("Number",names(Table.prop.dusky.neo))]="Neonate.n"
  Prop.neo=merge(Table.prop.dusky,Table.prop.dusky.neo,by=c("Finyear","BLOCK"))
  Prop.neo$prop=Prop.neo$Neonate.n/Prop.neo$Total.n
  ADD=FINYEARs[which(!FINYEARs%in%Prop.neo$Finyear)]
  ADD=data.frame(Finyear=ADD,BLOCK=NA,Total.n=NA,Neonate.n=NA,prop=NA)
  Prop.neo=rbind(Prop.neo,ADD)
  Prop.neo=Prop.neo[order(Prop.neo$Finyear,Prop.neo$BLOCK),]
  Uni=sort(as.numeric(Use.blk))
  Prop.neo$Index= as.numeric(substr(Prop.neo$Finyear,start=1,stop=4))
  YRS=unique( Prop.neo$Finyear)
  Y=as.numeric(substr(YRS,start=3,stop=4))
  Y=ifelse(nchar(Y)==1,paste(0,Y,sep=''),Y)
  Y1=as.numeric(substr(YRS,start=8,stop=9))
  Y1=ifelse(nchar(Y1)==1,paste(0,Y1,sep=''),Y1)
  NMS=paste(Y,'-',Y1,sep="")
  
  for(q in 1:length(Uni))
  {
    dd=subset(Prop.neo,BLOCK==Uni[q] & Total.n>=min.obs)
    LE=unique(dd$BLOCK)
    ADD=FINYEARs[which(!FINYEARs%in%dd$Finyear)]
    ADD=data.frame(Finyear=ADD,BLOCK=NA,Total.n=NA,Neonate.n=NA,prop=NA)
    ADD$Index= as.numeric(substr(ADD$Finyear,start=1,stop=4))
    dd=rbind(dd,ADD)    
    dd=dd[order(dd$Index,dd$BLOCK),]    
    mp=barplot(dd$prop,col="cadetblue",xaxt='n',ylim=c(0,1),main=paste('block',LE),cex.axis=1.25)
    text(mp,dd$prop+.05,dd$Total.n,cex=.8)
    axis(1,at=mp,F,tck=-0.03)
    box()
    axis(1,at=mp[seq(1,nrow(mp),2)],F,tck=-0.06,cex.axis=1)
    if(q%in%c(4,7))axis(1,at=mp[seq(1,nrow(mp),2)],NMS[seq(1,nrow(mp),2)],tck=-0.06,cex.axis=1.25)
  }
}
tiff(handl_OneDrive("Analyses/Population dynamics/Prop.neonate.dusky.tiff"),width = 2400, height = 2000,units = "px", res = 300,compression = "lzw")
par(mfcol=c(4,2),mar=c(1,2,1,1),oma=c(3,3,.1,.1),las=1,mgp=c(2.5,.75,0),cex.axis=.8,cex.lab=1.1)
fn.prop.dus(FL.dusky.neonate,min.obs=5,min.yrs=5)
mtext("Proportion of neonates in catch",2,line=1,outer=T,las=3,cex=1.5)
mtext("Financial year",1,line=1.75,outer=T,cex=1.5) 
dev.off()



# Export outputs from spatial size distribution paper-------------------------------------------------------------------------
if(do.paper)
{
  setwd(handl_OneDrive("Analyses/Size and sex patterns"))
  
  
  #Table 1. Summary of numbers, size and sex ratios for all shark and ray species
  Tab1.fun=function(Spec)
  {
    dat=subset(DATA,SPECIES==Spec)
    
    #Numbers                
    N.rel=length(dat$SPECIES)  
    
    #Size
    FL.rel.min=FL.rel.max=NA
    FL.rel.min=min(dat$FL,na.rm=T) 
    FL.rel.max=max(dat$FL,na.rm=T) 
    
    #sex ratio 
    Sex.rel=table(dat$SEX,useNA='ifany')
    id.M=match("M",names(Sex.rel))
    id.F=match("F",names(Sex.rel))
    id.U=match("U",names(Sex.rel))
    Rel.M=Rel.F=NA
    Rel.M=Sex.rel[id.M]
    Rel.F=Sex.rel[id.F]
    Rel.U=Sex.rel[id.U]
    
    #Proportion by gear
    Prop.gear=table(dat$Method)
    Prop.gear.GN.LL=Prop.gear[match(c("GN","LL"),names(Prop.gear))]
    Prop.gear.GN.LL=100*Prop.gear.GN.LL/sum(Prop.gear)
    Prop.gear.GN=Prop.gear.GN.LL[1]
    Prop.gear.LL=Prop.gear.GN.LL[2]
    
    return(list(N.rel=N.rel,FL.rel.min=FL.rel.min,FL.rel.max=FL.rel.max,
                Rel.M=Rel.M,Rel.F=Rel.F,Rel.U=Rel.U,Prop.gear.LL=Prop.gear.LL,
                Prop.gear.GN=Prop.gear.GN))
    
  }
  Tabl1.list=vector('list',length(Sks))
  names(Tabl1.list)=Sks
  for (i in 1:length(Sks)) Tabl1.list[[i]]=Tab1.fun(Sks[i])
  
  Tabl1.matrix=matrix(nrow=length(Sks),ncol=length(Tabl1.list[[1]]))
  
  for (i in 1:length(Sks))Tabl1.matrix[i,]=do.call(cbind,Tabl1.list[[i]])
  colnames(Tabl1.matrix)=names(Tabl1.list[[1]])
  rownames(Tabl1.matrix)=Sks
  
  Tabl1.matrix=as.data.frame(Tabl1.matrix)
  Tabl1.matrix$Species=rownames(Tabl1.matrix)
  Tabl1.matrix=merge(Tabl1.matrix,SPECIES.names,by="Species")
  Tabl1.matrix=Tabl1.matrix[order(-Tabl1.matrix$N.rel),]
  names(Tabl1.matrix)=c("Code","N","FL.min","FL.max","N.male","N.female","N.unknown","LL","GN","Species","Scient.Name")
  Tabl1.matrix=Tabl1.matrix[,match(c("Species","Scient.Name","Code",
                                     "FL.min","FL.max","N.male","N.female","N.unknown","N","LL","GN"),names(Tabl1.matrix))]
  Tabl1.matrix$N=Tabl1.matrix$N.male+Tabl1.matrix$N.female+Tabl1.matrix$N.unknown
  
  write.csv(Tabl1.matrix,"Table1.csv",row.names=F)
  
  
  #Overall Sex ratio ratio
  All.sp=Tabl1.matrix$Code 
  
  Chi.squar=function(datos)
  {
    X2=NA
    TABLA=table(datos$SEX)
    if(length(TABLA)>1) X2=chisq.test(TABLA)
    return(X2)
  }
  Store.X2=data.frame(Code=All.sp,X2.p=NA)
  for(i in 1:length(All.sp))Store.X2[i,2]=Chi.squar(subset(DATA,SPECIES%in%All.sp[i]& SEX%in%c("M","F")))
  Store.X2$X2.p=ifelse(!is.na(Store.X2$X2.p),round(Store.X2$X2.p,3),Store.X2$X2.p)
  
  
  This.sp=names(Table.species)
  
  
  #Range analysis
  DATA=subset(DATA,!(SPECIES=="GM" & Mid.Lat>(-26)))
  DATA=subset(DATA,!BLOCK%in%c(2415,2315,2519))  #remove land blocks
  
  Lat.range=c(-36,-10)
  Long.range=c(112,129)
  fn.plot.sp=function(dat)plot(dat$LONG,dat$LAT,xlim=Long.range,ylim=Lat.range,main=This.sp[i])
  for(i in 1:length(This.sp))fn.plot.sp(subset(DATA,SPECIES==This.sp[i]))
  
  Table.Blok.Sp=table(DATA$SPECIES,DATA$BLOCK)
  Table.Blok.Sp=ifelse(Table.Blok.Sp>0,1,0)
  Table.Blok.Sp=data.frame(SPECIES=rownames(Table.Blok.Sp),N.blok=rowSums(Table.Blok.Sp))
  Table.Blok.Sp=subset(Table.Blok.Sp,SPECIES%in%names(Table.species[Table.species>=100]))
  Table.Blok.Sp=merge(Table.Blok.Sp,SPECIES.names,by.x="SPECIES",by.y="Species",all.x=T)
  
  
  
  
  #Keep only gillnets of 6 to 7 inch mesh and LL
  Meshes=c("6","6.5","7")
  DATA=subset(DATA,Method%in%c("GN","LL"))
  DATA$Keep=with(DATA,ifelse(Method=="LL","YES",ifelse(MESH_SIZE %in%Meshes,"YES","NO")))
  DATA=subset(DATA,Keep=="YES")
  
  Table.mesh=table(DATA$MESH_SIZE)
  Table.mesh=round(Table.mesh/sum(Table.mesh),2)
  
  YRS=sort(unique(DATA$year))
  n.yrs=length(YRS)
  
  
  
  
  #2. Output commercial species FL data for 6.5 and 7 inch gillnet for population dynamics
  #2.1 Separate data by zone
  fn.byzone=function(a,coef.1,coef.2)
  {
    a$ZONE=as.character(a$zone)
    a$ZONE=with(a,ifelse(ZONE=="2","Zone2",ifelse(ZONE=="1","Zone1",ZONE)))
    a$ZONE=with(a,ifelse(is.na(ZONE) & LONG>=116.5 & LAT<=(-26),"Zone2",
                         ifelse(is.na(ZONE) & LONG<116.5 & LAT<=(-33),"Zone1",
                                ifelse(is.na(ZONE) & LAT>(-33) & LAT<=(-26) & LONG<116.5,"WC",
                                       ifelse(is.na(ZONE) & LAT>(-26) & LONG<114,"Closed",
                                              ifelse(is.na(ZONE) & LAT>(-26) & LONG>=114 & LONG<123.75,"North",
                                                     ifelse(is.na(ZONE) & LAT>(-26) & LONG>=123.75,"Joint",ZONE)))))))
    a=subset(a,!ZONE=="North")
    
    a$FL=with(a,ifelse(is.na(FL) & !is.na(TL),(TL-coef.1)/coef.2,FL))
    zones=unique(a$ZONE)
    Zones.size=vector('list',length(zones))
    names(Zones.size)=zones
    for(z in 1:length(Zones.size)) Zones.size[[z]]=subset(a,ZONE==zones[z])
    return(Zones.size)  
  }
  #note: coef.1 and .2 derived in De Wysiecki and Braccini 2017
  WH.size=fn.byzone(subset(DATA,SPECIES=="WH" & Method=="GN" & MESH_SIZE%in%c("6.5","7")),
                    coef.1=6.267,coef.2=1.063)
  GM.size=fn.byzone(subset(DATA,SPECIES=="GM" & Method=="GN" & MESH_SIZE%in%c("6.5","7")),
                    coef.1=2.007,coef.2=1.019)
  BW.size=fn.byzone(subset(DATA,SPECIES=="BW" & Method=="GN" & MESH_SIZE%in%c("6.5","7")),
                    coef.1=1.486,coef.2=1.202)
  TK.size=fn.byzone(subset(DATA,SPECIES=="TK" & Method=="GN" & MESH_SIZE%in%c("6.5","7")),
                    coef.1=5.87,coef.2=1.129)
  
  
  #Export sex ratio of commercial species for population dynamcis
  Sex.list=list(WH.size,GM.size,BW.size,TK.size)
  names(Sex.list)=c("WH","GM","BW","TK")
  Zne.sx.Ratio=list(WH=NA,GM=NA,BW=NA,TK=NA)
  fn.Zne.sx.Ratio=function(dat)
  {
    names(dat)[match(c("Zone2","West","Zone1"),names(dat))] =c("Zn2","WC","Zn1")
    n=length(dat)
    Stor=vector('list',n)
    names(Stor)=names(dat)
    for(i in 1:n)
    {
      aa=subset(dat[[i]],SEX%in%c("M","F"))
      Stor[[i]]=table(aa$SEX)
    }
    return(Stor)
  }
  
  #Export proportion of males in catch data for population dynamics modelling
  hndl=handl_OneDrive("Data/Population dynamics/Prop.males.in.catch/")
  if (Export.dat=="YES")
  {
    for(q in 1:length(Zne.sx.Ratio))
    {
      Zne.sx.Ratio[[q]]= fn.Zne.sx.Ratio(Sex.list[[q]])
      tab=do.call(rbind,Zne.sx.Ratio[[q]])
      
      tab.all=sum(tab[,2])/(sum(tab[,1])+sum(tab[,2]))
      tab1=tab[,2]/(tab[,1]+tab[,2])
      id=match(c("Zn1","Zn2","WC"),names(tab1))
      tab1=matrix(tab1[id],ncol=3)
      colnames(tab1)=c("Zn1","Zn2","WC")
      write.csv(tab1,paste(hndl,"prop.males.",names(Zne.sx.Ratio)[q],".csv",sep=""),row.names=F)
      write.csv(tab.all,paste(hndl,"prop.males.All.",names(Zne.sx.Ratio)[q],".csv",sep=""),row.names=F)
    }
  }
  
  fn.Zne.sx.Ratio.yr=function(dat,NM)
  {
    n=length(dat)
    for(i in 1:n)
    {
      aa=subset(dat[[i]],SEX%in%c("M","F"))
      dd=aggregate(Number~SEX+year,aa,sum)
      wide <- reshape(dd, v.names = "Number", idvar = "SEX",
                      timevar = "year", direction = "wide")
      p.male=unlist(wide[2,2:ncol(wide)]/colSums(wide[,2:ncol(wide)]))
      Yr=as.numeric(substr(names(p.male),8,18))
      Mod=lm(p.male~Yr)
      plot(Yr,p.male,pch=19,col=2,main=paste(NM,names(dat)[i]))
      legend("topright",paste("slope signif=",round(anova(Mod)$"Pr(>F)"[1],4)),bty='n')
      
    }
    
  }
  for(q in 1:length(Zne.sx.Ratio))fn.Zne.sx.Ratio.yr(Sex.list[[q]],names(Sex.list)[q])
  
  rm(Sex.list) 
  
  
  
  
  #2.2 create composition bins
  
  fn.size.comp=function(dat,Min,Max,interval)
  {
    dat=dat[order(dat$year),]
    
    dat$FINYEAR=with(dat,ifelse(Month>6,paste(year,"-",fn.subs(year+1),sep=""),
                                paste(year-1,"-",fn.subs(year),sep="")))
    
    Rango=c(Min,Max)
    SEQ=seq(Rango[1],Rango[2],interval)
    dat$FL.bin=floor(dat$FL/interval)*interval
    dat$FINYEAR=factor(dat$FINYEAR)
    
    #dat1=subset(dat,Mid.Lat<(-26) & Method=="GN" & MESH_SIZE=="7")
    #dat2=subset(dat,Mid.Lat<(-26) & Method=="GN" & MESH_SIZE=="6.5")  
    
    #Table.7.inch=with(dat1,table(FINYEAR,FL.bin))
    #Table.6.5.inch=with(dat2,table(FINYEAR,FL.bin))
    Table.6.5.and.7.inch=with(dat,table(FINYEAR,FL.bin))
    
    return(Table.6.5.and.7.inch)
    #   return(list(Table.7.inch=Table.7.inch,Table.6.5.inch=Table.6.5.inch,
    #               Table.6.5.and.7.inch=Table.6.5.and.7.inch))
  }
  
  Out.size.WH=WH.size
  Out.size.GM=GM.size
  Out.size.BW=BW.size
  Out.size.TK=TK.size
  
  Comm.sp=c("BW","TK","GM","WH")
  Min.FL=c(56,42,23,15)  #FL at birth
  Max.FL= c(310,200,170,150)  #max FL   
  
  #By zone
  for(i in 1:length(Out.size.WH)) Out.size.WH[[i]]=fn.size.comp(subset(Out.size.WH[[i]],
                                                                       FL<=Max.FL[4]& FL>Min.FL[4]),Min.FL[4],Max.FL[4],2)
  
  for(i in 1:length(Out.size.GM)) Out.size.GM[[i]]=fn.size.comp(subset(Out.size.GM[[i]],
                                                                       FL<=Max.FL[3] & FL>Min.FL[3]),Min.FL[3],Max.FL[3],2)
  
  for(i in 1:length(Out.size.BW)) Out.size.BW[[i]]=fn.size.comp(subset(Out.size.BW[[i]],
                                                                       FL<=Max.FL[1] & FL>Min.FL[1]),Min.FL[1],Max.FL[1],2)
  
  for(i in 1:length(Out.size.TK)) Out.size.TK[[i]]=fn.size.comp(subset(Out.size.TK[[i]],
                                                                       FL<=Max.FL[2] & FL>Min.FL[2]),Min.FL[2],Max.FL[2],2)
  
  # #By zone and sex
  # Out.size.WH.m=Out.size.WH.f=WH.size
  # Out.size.GM.m=Out.size.GM.f=GM.size
  # Out.size.BW.m=Out.size.BW.f=BW.size
  # Out.size.TK.m=Out.size.TK.f=TK.size
  # 
  # for(i in 1:length(Out.size.WH)) 
  # {
  #   Out.size.WH.m[[i]]=fn.size.comp(subset(Out.size.WH[[i]],SEX=='M' & FL<=Max.FL[4] & FL>Min.FL[4]),Min.FL[4],Max.FL[4],2)
  #   Out.size.WH.f[[i]]=fn.size.comp(subset(Out.size.WH[[i]],SEX=='F' & FL<=Max.FL[4] & FL>Min.FL[4]),Min.FL[4],Max.FL[4],2)
  # }
  # 
  # 
  # for(i in 1:length(Out.size.GM))
  # {
  #   Out.size.GM.m[[i]]=fn.size.comp(subset(Out.size.GM[[i]],SEX=='M' & FL<=Max.FL[3]& FL>Min.FL[3]),Min.FL[3],Max.FL[3],2)
  #   Out.size.GM.f[[i]]=fn.size.comp(subset(Out.size.GM[[i]],SEX=='F' & FL<=Max.FL[3]& FL>Min.FL[3]),Min.FL[3],Max.FL[3],2)
  # }
  # 
  # for(i in 1:length(Out.size.BW))
  # {
  #   Out.size.BW.m[[i]]=fn.size.comp(subset(Out.size.BW[[i]],SEX=='M' & FL<=Max.FL[1]& FL>Min.FL[1]),Min.FL[1],Max.FL[1],2)  
  #   Out.size.BW.f[[i]]=fn.size.comp(subset(Out.size.BW[[i]],SEX=='F' & FL<=Max.FL[1]& FL>Min.FL[1]),Min.FL[1],Max.FL[1],2)  
  # }
  # 
  # 
  # for(i in 1:length(Out.size.TK))
  # {
  #   Out.size.TK.m[[i]]=fn.size.comp(subset(Out.size.TK[[i]],SEX=='M' & FL<=Max.FL[2] & FL>Min.FL[2]),Min.FL[2],Max.FL[2],2)
  #   Out.size.TK.f[[i]]=fn.size.comp(subset(Out.size.TK[[i]],SEX=='F' & FL<=Max.FL[2] & FL>Min.FL[2]),Min.FL[2],Max.FL[2],2)
  # }
  
  #Fix nonsense species
  DATA$SPECIES=with(DATA,ifelse(SPECIES=="CP" & LAT>(-25),"BW",SPECIES))
  
  
  
  
  setwd(handl_OneDrive("Analyses/Size and sex patterns"))
  
  #Species diversity
  fn.map=function()
  {
    plotMap(worldLLhigh, xlim=XAXIS,ylim=YAXIS,
            plt = NULL,col="grey80",tck = 0.025, tckMinor = 0.0125,
            xlab="",ylab="",axes=F)
    lines(rbind(129,129),rbind(-15,-31.5),lty=2,col="grey50")
    axis(2,seq(YAXIS[1],YAXIS[2],1),F,tck=-0.03)
    axis(2,seq(YAXIS[1],YAXIS[2],4),F,tck=-0.08)
    axis(1,seq(XAXIS[1],XAXIS[2],1),F,tck=-0.03)
    axis(1,seq(XAXIS[1],XAXIS[2],4),F,tck=-0.08)
    box()
  }
  
  AxisY=function() axis(2,seq(YAXIS[1],YAXIS[2],4),-seq(YAXIS[1],YAXIS[2],4),cex.axis=1.25,las=2,tck=-0.08)
  AxisX=function() axis(1,seq(XAXIS[1],XAXIS[2],4),seq(XAXIS[1],XAXIS[2],4),cex.axis=1.25,tck=-0.08)
  
  fn.diversity=function(dat,LEG)
  {
    dat=subset(dat,!BLOCK==0)
    a=dat[,match(c("BLOCK","LAT","LONG"),names(dat))]
    a=a[!duplicated(a$BLOCK),]
    
    #Species richness (Number of species)
    Table.Blok.Sp=table(dat$SPECIES,dat$BLOCK)
    Table.Blok.Sp=ifelse(Table.Blok.Sp>0,1,0)
    Table.Blok.Sp=data.frame(BLOCK=colnames(Table.Blok.Sp),N.species=colSums(Table.Blok.Sp))
    Table.Blok.Sp=merge(Table.Blok.Sp,a,by="BLOCK",all.x=T)
    MAX=max(Table.Blok.Sp$N.species)
    fn.map()
    points(Table.Blok.Sp$LONG+0.5,Table.Blok.Sp$LAT-0.5,pch=19,col="black",
           cex=(Table.Blok.Sp$N.species/MAX)*2)
    mtext(LEG,3)
    if(LEG=="Gillnet")
    {
      legend("topleft","Richness",bty='n',cex=1.5)
      legend("right",c("1","10","20"),pch=19,col=1,pt.cex=c(1/MAX,10/MAX,20/MAX)*2,bty='n',cex=1.3)
      AxisY()
    }
    
    
    #Species evenness
    Table.Blok.Sp=table(dat$BLOCK,dat$SPECIES)
    H=diversity(Table.Blok.Sp)  #Shannon-Weaver
    J=H/log(specnumber(Table.Blok.Sp))  #Pielou's evenness
    Table.eve=data.frame(BLOCK=names(J),Evenness=J)
    Table.eve=merge(Table.eve,a,by="BLOCK",all.x=T)
    fn.map()
    points(Table.eve$LONG+0.5,Table.eve$LAT-0.5,pch=19,col="black",
           cex=Table.eve$Evenness*2)
    if(LEG=="Gillnet")
    {
      legend("topleft","Evenness",bty='n',cex=1.5)
      legend("right",c("0.25","0.5","1"),pch=19,col=1,pt.cex=c(0.25,0.5,1)*2,bty='n',cex=1.3)
      AxisY()
    }
    
    if(LEG=="Longline")
    {
      legend(116,-22,"Western",bty='n',cex=1.3)
      legend(116,-24,"Australia",bty='n',cex=1.3)
    }
    
    #Species diversity
    Table.div=data.frame(BLOCK=names(H),Diversity=H)
    Table.div=merge(Table.div,a,by="BLOCK",all.x=T)
    fn.map()
    points(Table.div$LONG+0.5,Table.div$LAT-0.5,pch=19,col="black",
           cex=Table.div$Diversity)
    if(LEG=="Gillnet")
    {
      legend("topleft","Diversity",bty='n',cex=1.5)
      legend("right",c("0.5","1","2"),pch=19,col=1,pt.cex=c(0.5,1,2),bty='n',cex=1.3)
      AxisY()
    }
    AxisX()
  }
  
  
  XAXIS=c(110,130)
  YAXIS=c(-37,-13)
  
  tiff("Figure Richness.tiff",width = 1300, height = 2200,units = "px", res = 300,compression = "lzw")
  par(mfcol=c(3,2),mai=c(0.1,.1,0.15,.01),oma=c(3,3,0.15,0.01))
  fn.diversity(subset(DATA,Method=="GN"),"Gillnet")
  
  fn.diversity(subset(DATA,Method=="LL"),"Longline")
  mtext("Longitude (?E)",side=1,line=1.6,font=1,las=0,cex=1.4,outer=T)
  mtext("Latitude (?S)",side=2,line=1,font=1,las=0,cex=1.4,outer=T)
  dev.off()
  
  
  #Percent deviance explained
  #Overall
  Dsquared <- function(model, adjust = F)
  {
    if(!is.null(model$deviance))
    {
      d2 <- (model$null.deviance - model$deviance) / model$null.deviance
      if (adjust)
      {
        n <- length(model$fitted.values)
        p <- length(model$coefficients)
        d2 <- 1 - ((n - 1) / (n - p)) * (1 - d2)
      }
      d3=d2*100
      d1=model$deviance
      d2=model$null.deviance
    }
    
    if(is.null(model$deviance))
    {
      d1=NA
      d2=NA
      d3=NA
    }
    return(list(d1=d1,d2=d2,d3=d3))
  }
  
  #Each term
  fn.dev.exp.term=function(dat)
  {
    100*(dat$Anova$Deviance[2:length(dat$Anova$Deviance)]/dat$model$null.deviance)
  }
  
  
  #Model fit diagnostics
  fn.plot.diag=function(MODEL,SPECIES)
  {
    RES=MODEL$residuals   #residuals
    Std.RES=RES/sd(RES)   #standardised residuals (res/SD(res))
    PREDS=predict(MODEL)
    
    par(mfcol=c(2,2),las=1,mar=c(3,3,2,1),oma=c(2.5,.1,.1,.1),las=1,mgp=c(2,.5,0),cex.axis=.8,cex.lab=1.1)
    qqnorm(RES,main="",ylim=c(-5,5),xlim=c(-5,5),ylab="Residuals",xlab="Quantiles of standard normal distribution")
    qqline(RES, col = 'grey40',lwd=1.5,lty=2)
    
    hist(Std.RES,xlim=c(-5,5),ylab="Frequency",xlab="Stan. residuals",main="",col="grey",breaks=50)
    box()
    
    plot(PREDS,Std.RES,ylim=c(-4,12),ylab="Stan. residuals",xlab="Expected values")
    abline(0,0,lwd=1.5,lty=2,col='grey40')
    
    plot(PREDS,sqrt(abs(Std.RES)),ylim=c(0,2.6),ylab="Square root of stan. residuals",xlab="Expected values")
    mtext(paste(SPECIES),3,outer=T,lin=-1)
  }
  
  
  
  #Get blocks
  Blk.pos=subset(DATA,select=c("BLOCK","LAT","LONG"))
  Blk.pos=Blk.pos[!duplicated(Blk.pos$BLOCK),]
  Blk.pos=subset(Blk.pos,!(is.na(BLOCK) | BLOCK==0))
  
  
  
  #Exploratory analyses
  cfac=function(x,breaks=NULL)
  {
    x=round(x,2)  
    if(is.null(breaks)) breaks=unique(quantile(x,na.rm=T))
    x=cut(x,breaks,include.lowest=T,right=F)
    levels(x)=paste(breaks[-length(breaks)],ifelse(diff(breaks)>1,
                                                   c(paste('-',breaks[-c(1,length(breaks))]-1,sep=''),'+'),''),sep='')
    return(x)
  }
  fn.explore=function(dat)
  {
    dat=subset(dat,!is.na(Mid.Long)  & Mid.Long>0)
    dat=subset(dat,!is.na(Mid.Lat) & Mid.Lat<0)
    dat$Sex.bin=with(dat,ifelse(SEX=="M",1,0))
    dat$Soi=cfac(dat$SOI)
    dat$FREO=cfac(dat$Freo)
    dat$Lat=cfac(dat$Mid.Lat)
    dat$Lon=cfac(dat$Mid.Long)
    
    # Main Effect of factors
    par(mfcol=c(1,1),mai=c(.3,1,1,.1),oma=c(1,.1,.1,.1),las=1)
    plot.design(Sex.bin~as.factor(BLOCK)+as.factor(year)+as.factor(Month)+Soi+FREO+Lat+Lon,data=dat,
                cex.lab=1.5,cex.axis=.75,main=paste("Sex ratio",unique(dat$COMMON_NAME)))
    
    par(mfcol=c(4,2),mai=c(.7,.7,.1,.1),oma=c(2,2,.1,.1),las=1)
    plot(Sex.bin~as.factor(Month),data=dat,ylab="",xlab="Month") 
    plot(Sex.bin~as.factor(year),data=dat,ylab="",xlab="year")
    plot(Sex.bin~as.factor(BLOCK),data=dat,ylab="",xlab="block")
    plot(Sex.bin~Soi,data=dat,ylab="",xlab="SOI")
    plot(Sex.bin~FREO,data=dat,ylab="",xlab="Freo")
    plot(Sex.bin~Lat,data=dat,ylab="",xlab="Lat")
    plot(Sex.bin~Lon,data=dat,ylab="",xlab="Long")
    mtext(paste("Sex ratio",unique(dat$COMMON_NAME)),2,0,outer=T,las=3,cex=1.5)
    
    #Interactions
    par(mfcol=c(1,1))
    interaction.plot(as.factor(dat$year),as.factor(dat$BLOCK),dat$Sex.bin,ylab="response",xlab="block")
    
    
    
    #Block
    agg.Male=aggregate(Number~BLOCK,subset(dat,SEX=="M"),sum)
    names(agg.Male)[2]="N.male"
    agg.Fem=aggregate(Number~BLOCK,subset(dat,SEX=="F"),sum)
    names(agg.Fem)[2]="N.female"
    agg.sex=merge(agg.Male,agg.Fem,by=c("BLOCK"),all.x=T,all.y=T)
    agg.sex$N.male=with(agg.sex,ifelse(is.na(N.male),0,N.male))
    agg.sex$N.female=with(agg.sex,ifelse(is.na(N.female),0,N.female))
    agg.sex$ratio=round(agg.sex$N.male/agg.sex$N.female,2)
    agg.sex.block=agg.sex
    
    #Year
    agg.Male=aggregate(Number~year,subset(dat,SEX=="M"),sum)
    names(agg.Male)[2]="N.male"
    agg.Fem=aggregate(Number~year,subset(dat,SEX=="F"),sum)
    names(agg.Fem)[2]="N.female"
    agg.sex=merge(agg.Male,agg.Fem,by=c("year"),all.x=T,all.y=T)
    agg.sex$N.male=with(agg.sex,ifelse(is.na(N.male),0,N.male))
    agg.sex$N.female=with(agg.sex,ifelse(is.na(N.female),0,N.female))
    agg.sex$ratio=round(agg.sex$N.male/agg.sex$N.female,2)
    agg.sex.year=agg.sex
    
    #Month
    agg.Male=aggregate(Number~Month,subset(dat,SEX=="M"),sum)
    names(agg.Male)[2]="N.male"
    agg.Fem=aggregate(Number~Month,subset(dat,SEX=="F"),sum)
    names(agg.Fem)[2]="N.female"
    agg.sex=merge(agg.Male,agg.Fem,by=c("Month"),all.x=T,all.y=T)
    agg.sex$N.male=with(agg.sex,ifelse(is.na(N.male),0,N.male))
    agg.sex$N.female=with(agg.sex,ifelse(is.na(N.female),0,N.female))
    agg.sex$ratio=round(agg.sex$N.male/agg.sex$N.female,2)
    agg.sex.Month=agg.sex
    
    #Block and Year 
    agg.Male=aggregate(Number~BLOCK+year,subset(dat,SEX=="M"),sum)
    names(agg.Male)[3]="N.male"
    agg.Fem=aggregate(Number~BLOCK+year,subset(dat,SEX=="F"),sum)
    names(agg.Fem)[3]="N.female"
    agg.sex=merge(agg.Male,agg.Fem,by=c("BLOCK","year"),all.x=T,all.y=T)
    agg.sex$N.male=with(agg.sex,ifelse(is.na(N.male),0,N.male))
    agg.sex$N.female=with(agg.sex,ifelse(is.na(N.female),0,N.female))
    agg.sex$ratio=round(agg.sex$N.male/agg.sex$N.female,2)
    agg.sex.Blk.yr=agg.sex
    
    #Block, Year and Month
    agg.Male=aggregate(Number~BLOCK+year+Month,subset(dat,SEX=="M"),sum)
    names(agg.Male)[4]="N.male"
    
    agg.Fem=aggregate(Number~BLOCK+year+Month,subset(dat,SEX=="F"),sum)
    names(agg.Fem)[4]="N.female"
    agg.sex=merge(agg.Male,agg.Fem,by=c("BLOCK","year","Month"),all.x=T,all.y=T)
    agg.sex$N.male=with(agg.sex,ifelse(is.na(N.male),0,N.male))
    agg.sex$N.female=with(agg.sex,ifelse(is.na(N.female),0,N.female))
    agg.sex$ratio=round(agg.sex$N.male/agg.sex$N.female,2)
    agg.sex.Blk.yr.mn=agg.sex
    
    return(list(BLK=agg.sex.block,YR=agg.sex.year,MN=agg.sex.Month,
                BLK.Yr=agg.sex.Blk.yr,Blk.Yr.Mn=agg.sex.Blk.yr.mn))
  }
  
  depth=function(dat)
  {
    a=aggregate(Number~BOTDEPTH,subset(dat,SEX=="M"),sum)
    b=aggregate(Number~BOTDEPTH,subset(dat,SEX=="F"),sum)
    agg.sex=merge(a,b,by=c("BOTDEPTH"),all.x=T,all.y=T)
    plot(agg.sex$BOTDEPTH,agg.sex$Number.x/agg.sex$Number.y,xlim=c(0,150),ylab="male:female",xlab="depth (m)",
         main=unique(dat$SPECIES))
  }
  
  for(i in 1:length(Commercial.Sks))depth(subset(DATA,SPECIES%in%Commercial.Sks[i]& SEX%in%c("M","F")))
  
  
  #Colinearity
  fn.colinearity=function(dat)
  {
    
    dat$BOTDEPTH.bin=as.factor(round(dat$BOTDEPTH/10)*10)
    dat$LONG=as.factor(dat$LONG)
    dat$LAT=as.factor(dat$LAT)
    
    Covars=dat[,match(c("BOTDEPTH.bin","LAT","LONG"),names(dat))]
    #Covars=dat[,match(c("year","Month","BOTDEPTH","Moon","SOI","Freo","Mid.Lat","Mid.Long","Temperature"),names(dat))]
    Covars=na.omit(Covars)
    
    panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
    {
      usr <- par("usr"); on.exit(par(usr))
      par(usr = c(0, 1, 0, 1))
      r <- abs(cor(x, y))
      txt <- format(c(r, 0.123456789), digits=digits)[1]
      txt <- paste(prefix, txt, sep="")
      if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
      text(0.5, 0.5, txt, cex = cex.cor * r)
    }
    pairs(Covars, lower.panel=panel.smooth, upper.panel=panel.cor,main=This.sp[i])
    
  }
  
  
  
  #Explore raw data for patterns in sandbar FL and Whiskery sex ratio
  fn.plot.all.FL=function(a,YLIM,XLIM,N.min)
  {
    a$Number=1
    This=aggregate(FL~Mid.Lat+Mid.Long,a,mean,na.rm=T)
    Number=aggregate(Number~Mid.Lat+Mid.Long,a,sum,na.rm=T)
    This=merge(This,Number,by=c("Mid.Lat","Mid.Long"),all.x=T)
    This=subset(This,Number>N.min)
    This$Prop=This$FL/max(This$FL)
    
    par(mfcol=c(2,1),mai=c(0.7,1,.25,.1),mgp=c(1,.5,0))
    plot(This$Mid.Long,This$Mid.Lat,pch=19,cex=(This$Prop),ylim=YLIM,xlim=XLIM,col=4,
         main=paste("Mean FL (as % of max mean FL) for",unique(a$SPECIES)),xlab="Longitude",ylab="Latitude")
    legend("topright",c("10%","50%","100%"),bty="n",col=4,pch=c(19,19,19),pt.cex=c(0.1,0.5,1),cex=1.1)
    
    b=aggregate(FL~BOTDEPTH,a,mean,na.rm=T)
    plot(b$BOTDEPTH,b$FL,xlab="Depth",ylab="FL",main=paste("All records,",unique(a$SPECIES)),xlim=c(0,250))
  }
  
  fn.plot.all.FL(subset(DATA,SPECIES=="TK"),c(-36,-14),c(113,120),N.min=6)
  fn.plot.all.FL(subset(DATA,SPECIES=="BW"),c(-36,-14),c(113,128),N.min=6)
  fn.plot.all.FL(subset(DATA,SPECIES=="GM"),c(-36,-14),c(113,128),N.min=6)
  
  
  #Sex ratio by shot for shots with more than n observations
  fn.plot.all.sex=function(a,YLIM,XLIM,N.min)
  {
    a$Number=1
    mal=aggregate(Number~Mid.Lat+Mid.Long,subset(a,SEX=="M"),sum)
    names(mal)[3]="N.male"
    fem=aggregate(Number~Mid.Lat+Mid.Long,subset(a,SEX=="F"),sum)
    names(fem)[3]="N.female"
    both=merge(mal,fem,by=c("Mid.Lat","Mid.Long"),all.x=T,all.y=T)
    both$N.male=with(both,ifelse(is.na(N.male),0,N.male))
    both$N.female=with(both,ifelse(is.na(N.female),0,N.female))
    both$Sum=both$N.male+both$N.female
    both$Prop.male=both$N.male/both$Sum
    both$Prop.male=ifelse(both$Prop.male==0,0.01,both$Prop.male)
    both$Prop.female=both$N.female/both$Sum
    both$Prop.female=ifelse(both$Prop.female==0,0.01,both$Prop.female)
    both=subset(both,Sum>N.min)
    
    par(mfcol=c(2,1),mai=c(.7,1,.25,.1),mgp=c(1,.5,0))
    plot(both$Mid.Long,both$Mid.Lat,pch=19,cex=(both$Prop.male),ylim=YLIM,xlim=XLIM,col=4,
         main=paste("Percentage of males per shot for",unique(a$SPECIES)),xlab="",ylab="Latitude")
    legend("topright",c("10%","50%","100%"),bty="n",col=4,pch=c(19,19,19),pt.cex=c(0.1,0.5,1),cex=1.1)
    
    plot(both$Mid.Long,both$Mid.Lat,pch=19,cex=(both$Prop.female),ylim=YLIM,xlim=c(113,128),col="pink",
         main=paste("Percentage of females per shot for",unique(a$SPECIES)),xlab="Longitude",ylab="Latitude")
    
  }
  
  
  fn.plot.all.sex(subset(DATA,SPECIES=="WH"),c(-36,-26),c(113,128),N.min=9)
  fn.plot.all.sex(subset(DATA,SPECIES=="GM"),c(-36,-26),c(113,128),N.min=9)
  fn.plot.all.sex(subset(DATA,SPECIES=="BW"),c(-36,-25),c(113,128),N.min=9)
  fn.plot.all.sex(subset(DATA,SPECIES=="TK"),c(-36,-17),c(113,120),N.min=9)
  
  
  
  #Schooling behaviour by shot
  fn.school.sex=function(a,N.min)
  {
    A="NO SCHOOLING"  
    test=subset(a,SEX%in%c("M","F"))
    if(nrow(test)>0)
    {
      a$Number=1
      mal=aggregate(Number~Mid.Lat+Mid.Long,subset(a,SEX=="M"),sum)
      names(mal)[3]="N.male"
      fem=aggregate(Number~Mid.Lat+Mid.Long,subset(a,SEX=="F"),sum)
      names(fem)[3]="N.female"
      both=merge(mal,fem,by=c("Mid.Lat","Mid.Long"),all.x=T,all.y=T)
      both$N.male=with(both,ifelse(is.na(N.male),0,N.male))
      both$N.female=with(both,ifelse(is.na(N.female),0,N.female))
      both$Sum=both$N.male+both$N.female
      both$Prop.male=both$N.male/both$Sum
      both$Prop.female=both$N.female/both$Sum
      both=subset(both,Sum>N.min)
      
      if(nrow(both)>0)
      {
        both$School=with(both,ifelse(Prop.male>0.9| Prop.female>0.9,"Schooling","No.schooling"))
        A=table(both$School)
      }
    }
    return(A)
  }
  Schooling=vector('list',length(This.sp))
  names(Schooling)=This.sp
  for(i in 1:length(This.sp))Schooling[[i]]=fn.school.sex(subset(DATA,SPECIES==This.sp[[i]]),N.min=19)
  
  
  #Export data for royal society open science repository
  RSOS.export=subset(DATA,!is.na(LAT),select=c(SPECIES,SCIENTIFIC_NAME,SEX,FL,zone,Method,LAT,BOTDEPTH,Month))
  RSOS.export=subset(RSOS.export,LAT<0)
  colnames(RSOS.export)[match("BOTDEPTH",names(RSOS.export))]="DEPTH"
  RSOS.export$SCIENTIFIC_NAME=as.character(RSOS.export$SCIENTIFIC_NAME)
  RSOS.export$SCIENTIFIC_NAME=with(RSOS.export,ifelse(SPECIES=="CA","Scyliorhinidae",
                                                      ifelse(SPECIES=="MI","Rhizoprionodon acutus",
                                                             ifelse(SPECIES=="BN","Carcharhinus altimanus",
                                                                    ifelse(SPECIES=="BU","Carcharhinus leucas",
                                                                           ifelse(SPECIES=="BX","Hexanchus spp",
                                                                                  ifelse(SPECIES=="GR","Carcharhinus amblyrhynchos",
                                                                                         ifelse(SPECIES=="HG","Sphyrna mokarran",
                                                                                                ifelse(SPECIES=="HW","Eusphyra blochii",
                                                                                                       ifelse(SPECIES=="LP","",
                                                                                                              ifelse(SPECIES=="NS","Carcharhinus cautus",
                                                                                                                     ifelse(SPECIES=="PC","Pristis clavata",
                                                                                                                            ifelse(SPECIES=="RB","Carcharhinus melanopterus",
                                                                                                                                   ifelse(SPECIES=="RW","Triaenodon obesus",
                                                                                                                                          ifelse(SPECIES=="SC","Pristiophorus cirratus",
                                                                                                                                                 ifelse(SPECIES=="SF","Hemitriakis falcata",
                                                                                                                                                        ifelse(SPECIES=="SN","Hemipristis elongata",
                                                                                                                                                               ifelse(SPECIES=="WE","Hemigaleidae",
                                                                                                                                                                      ifelse(SPECIES=="WG","Mustelus Sp.B",
                                                                                                                                                                             SCIENTIFIC_NAME)))))))))))))))))))
  write.csv(RSOS.export,"Data.for.RSOS.csv",row.names=F)
  
  
  
  #Select species for analysis (at least 100 observations overall, & > 10 observation in at least 2 zones)
  overall.n=100
  block.n=10
  min.block=10
  min.zn=50
  zone.n=2
  min.yr.zone=10
  
  Select.fn=function(dat)
  {
    Use="NO"
    
    #Minimum number of observations overall
    N=nrow(dat)
    
    #Minimum number of observations per block
    Table.blk=with(dat,table(BLOCK))
    Table.blk=Table.blk[Table.blk>=min.block]
    
    #Minimum number of observations per zone
    Table.zn=with(dat,table(zone))
    Table.zn=Table.zn[Table.zn>=min.zn]
    
    #Table YRs
    Table.Yrs.Zn=with(dat,table(zone,year))
    Table.Yrs.Zn[Table.Yrs.Zn>0]=1
    Zone_yr=rowSums(Table.Yrs.Zn)
    Zone_yr=Zone_yr[Zone_yr>=min.yr.zone]
    
    #Minimum number of blocks
    #  if(N>=overall.n & length(Table.blk)>=block.n)Use="YES"
    
    #Selection criteria
    if(N>=overall.n & length(Table.zn)>=zone.n & length(Zone_yr)>0)Use="YES"
    
    return(Use)
    
  }
  
  Selected=vector('list',length(This.sp))
  names(Selected)=This.sp
  for(i in 1:length(This.sp))Selected[[i]]=Select.fn(subset(DATA,SPECIES==This.sp[[i]] & year<=Current.year & FL>0 & !is.na(FL) & FL<800 & !(is.na(BLOCK) | BLOCK==0 )))
  ID=which(Selected=="YES")
  
  This.sp=names(ID)
  
  #sort species by importance
  a=subset(DATA,SPECIES%in%This.sp)
  SP.sort=rev(sort(table(a$SPECIES)))
  rm(a)
  This.sp=This.sp[match(names(SP.sort),This.sp)]
  
  Species.labels=c("Dusky shark","Sandbar shark","Gummy shark","Whiskery shark","Smooth hammerhead",
                   "Spinner shark","Spot-tail shark","Milk shark","Blacktip sharks","Tiger shark",
                   "Pencil Shark","Scalloped Hammerhead")
  
  Comm.Species.labels=c("Dusky shark","Sandbar shark","Gummy shark","Whiskery shark")
  
  names(This.sp)=Species.labels
  Species.Scient=c("C. obscurus","C. plumbeus","M. antarcticus","F. macki",
                   "S. zygaena","C. brevipinna","C. sorrah","R. acutus","C. limbatus & C. tilstoni",
                   "G. cuvier","H. hyugaensis","S. lewini")
  names(Species.Scient)=This.sp
  
  Species.labels.Scientific=c(expression(italic("C. obscurus")),expression(italic("C. plumbeus")),
                              expression(italic("M. antarcticus")),expression(italic("F. macki")),
                              expression(italic("S. zygaena")),expression(italic("C. brevipinna")),
                              expression(italic("C. sorrah")),expression(italic("R. acutus")),
                              expression(italic("C. limbatus & C. tilstoni")),expression(italic("G. cuvier")),        
                              expression(italic("H. hyugaensis")),expression(italic("S. lewini")))
  
  
  Store.explore=vector('list',length(This.sp))
  names(Store.explore)=This.sp
  #for(i in 1:length(This.sp))Store.explore[[i]]=fn.explore(subset(DATA,SPECIES%in%This.sp[i]& SEX%in%c("M","F")))
  
  
  
  #Small-scale patterns
  a=subset(DATA,SPECIES%in%This.sp & FL>0 & !is.na(FL) & FL<800 & !(is.na(BLOCK) | BLOCK==0 ))
  #show zone with max N obsv
  fun.max.N=function(a)
  {
    d=subset(DATA,SPECIES==a & FL>0 & !is.na(FL) & FL<800 & !(is.na(BLOCK) | BLOCK==0 ))
    Table.blk=table(d$zone,d$Month)
    d=rowSums(Table.blk)
    return(names(d[which(d==max(d))]))
  }
  Blk.sp=rep(NA,length(This.sp))
  names(Blk.sp)=This.sp
  for(i in 1:length(This.sp)) Blk.sp[i]= fun.max.N(This.sp[[i]])
  
  
  #Table.blk=table(a$BLOCK,a$Month,as.character(a$SPECIES))
  #Blk.sp=c(2114,3415,3415,3415,3115,2114,2114,2114,3215,3415)
  #names(Blk.sp)=sort(unique(as.character(a$SPECIES)))
  Blk.sp=Blk.sp[match(This.sp,names(Blk.sp))]
  
  Small.scale=function(dat,sP)
  {
    boxplot(FL~Month,dat,main="",col='grey80',cex.axis=0.8)
    # legend("topleft",sP,bty='n',cex=1.35,adj=c(0.1,0))
    mtext(sP,3,line=0,cex=1.1)
  }
  
  
  
  #Variablity by month
  tiff("Figure S5.tiff",width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
  par(mfcol=c(4,3),mai=c(.2,0.25,.175,.1),oma=c(2,2,.2,.1),mgp=c(1,.5,0),las=1)
  for(i in 1:length(This.sp))Small.scale(subset(DATA,SPECIES==This.sp[[i]] & !is.na(FL) &
                                                  !(is.na(BLOCK) | BLOCK==0 )),Species.labels.Scientific[i])
  mtext("Fork length (cm)",2,line=0.25,las=3,cex=1.25,outer=T)
  mtext("Month",1,line=0.5,cex=1.25,outer=T)
  dev.off() 
  
  #Example variablity at latitude
  Small.scale.LAT=function(dat,sP)
  {
    boxplot(FL~LAT,dat,main="",col='grey80',cex.axis=1)
    # legend("topleft",sP,bty='n',cex=1.35,adj=c(0.1,0))
    mtext(sP,3,line=0,cex=1)
  }
  
  tiff("Figure S4.tiff",width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
  par(mfcol=c(4,3),mai=c(.2,0.25,.175,.1),oma=c(2,2,.2,.1),mgp=c(1,.5,0),las=1)
  for(i in 1:length(This.sp))Small.scale.LAT(subset(DATA,SPECIES==This.sp[[i]] & !is.na(FL) &
                                                      !(is.na(BLOCK) | BLOCK==0 )),Species.labels.Scientific[i])
  mtext("Fork length (cm)",2,line=0.25,las=3,cex=1.25,outer=T)
  mtext("Latitude (?S)",1,line=0.5,cex=1.25,outer=T)
  dev.off() 
  
  
  
  
  Previous.approach="NO"  #(i.e. use lat and long as predictors)
  
  if(Previous.approach=="YES")
  {
    Small.scale=function(dat,COL,PCH)
    {
      dat$Mon=factor(dat$Month,levels=1:12)
      Tab=table(dat$Mon)
      Tab=Tab[Tab>9]
      dat=subset(dat,Mon%in% names(Tab))
      ag=aggregate(FL~Mon,dat,mean)
      ID=which(!levels(dat$Mon)%in%ag$Mon)
      if(length(ID)>0)
      {
        ADD=data.frame(Mon=ID,x1=NA)
        colnames(ADD)=colnames(ag)
        ag=rbind(ag,ADD)    
      }
      ag=ag[order(ag$Mon),]
      lines(ag$Mon,ag[,2],col=COL,lwd=2)
      points(ag$Mon,ag[,2],col=COL,cex=1.25,pch=PCH,bg="white")
      
    }
    Small.scale.Sex=function(dat,COL,PCH)
    {
      dat$Mon=factor(dat$Month,levels=1:12)
      Tab=table(dat$Mon)
      Tab=Tab[Tab>9]
      dat=subset(dat,Mon%in% names(Tab))
      dat$Number=1
      
      agg.Male=aggregate(Number~Mon,subset(dat,SEX=="M"),sum)
      names(agg.Male)[2]="N.male"
      agg.Fem=aggregate(Number~Mon,subset(dat,SEX=="F"),sum)
      names(agg.Fem)[2]="N.female"
      agg.sex=merge(agg.Male,agg.Fem,by=c("Mon"),all.x=T,all.y=T)
      agg.sex$N.male=with(agg.sex,ifelse(is.na(N.male),0,N.male))
      agg.sex$N.female=with(agg.sex,ifelse(is.na(N.female),0,N.female))
      agg.sex$Prop=agg.sex$N.male/(agg.sex$N.male+agg.sex$N.female)
      agg.sex=agg.sex[,c(1,4)]
      
      ID=which(!levels(dat$Mon)%in% agg.sex$Mon)
      if(length(ID)>0)
      {
        ADD=data.frame(Mon=ID,x1=0)
        colnames(ADD)=colnames(agg.sex)
        agg.sex=rbind(agg.sex,ADD)    
      }
      agg.sex=agg.sex[order(agg.sex$Mon),]
      lines(agg.sex$Mon,agg.sex[,2],col=COL,lwd=2)
      points(agg.sex$Mon,agg.sex[,2],col=COL,cex=1.25,pch=PCH,bg="white")
    }
    
    tiff("Figure Month.tiff",width = 1300, height = 2200,units = "px", res = 300,compression = "lzw")
    par(mfcol=c(2,1),mai=c(.2,0.25,.1,.1),oma=c(2,2,.1,.1),mgp=c(1,.5,0),las=1)
    EXP=0.95
    #FL
    plot(1:12,ylim=c(70,170),ylab="",xlab="",col="transparent",xaxt='n')
    COLS=c("black","grey80","grey45","grey65","grey30","black","grey80","grey45","grey55","grey80")
    PCHs=c(rep(17,5),21,rep(19,4))
    IDS=c(1,5,6,9,10)
    for(i in IDS)Small.scale(subset(DATA,SPECIES==This.sp[[i]] & zone==Blk.sp[i] & FL>0 & !is.na(FL) & FL<800 &
                                      !(is.na(BLOCK) | BLOCK==0 )),COLS[i],PCHs[i])
    legend("top",Species.labels.Scientific[IDS],pch=PCHs[IDS],col=COLS[IDS],bty='n',cex=EXP,pt.cex=EXP)  
    mtext("Mean fork length (cm)",2,line=2,las=3,cex=1.5)
    axis(1,1:12,F,tck=-0.02)
    axis(1,seq(2,12,2),F,tck=-0.05)
    
    #sex ratio
    plot(1:12,ylim=c(0,1.1),ylab="",xlab="",col="transparent")
    COLS=c("black","grey80","grey45","black","grey30","grey45","grey80","grey45","black","black")
    PCHs=c(17,17,17,17,17,19,19,17,21,19)
    IDS=c(3:4,7:10)
    for(i in IDS)Small.scale.Sex(subset(DATA,SPECIES==This.sp[i]& SEX%in%c("M","F")),COLS[i],PCHs[i])
    legend("top",Species.labels.Scientific[IDS],pch=PCHs[IDS],col=COLS[IDS],bty='n',cex=EXP,pt.cex=EXP)  
    mtext("Proportion of males",2,line=2,las=3,cex=1.5)
    mtext("Month",1,line=1.65,cex=1.5)
    axis(1,1:12,F,tck=-0.02)
    axis(1,seq(2,12,2),F,tck=-0.05)
    dev.off()
    
  }
  
  
  
  #2. Spatial patterns in size
  Use.Lat.Long="NO"
  
  if(Use.Lat.Long=="YES")
  {
    
    #2.1. Prelim analyses
    #Colinearity
    for(i in 1:length(This.sp))fn.colinearity(subset(DATA,SPECIES==This.sp[i] & FL>0 & !is.na(FL) & FL<800 & !(is.na(BLOCK) | BLOCK==0 )))
    
    #Confounding
    tabl.terms=function(dat)
    {
      dat$LONG=as.factor(dat$LONG)
      dat$LAT=as.factor(dat$LAT)
      table(dat$LAT,dat$LONG)
    }
    Table.Conf=vector('list',length(This.sp))
    names(Table.Conf)=This.sp
    for(i in 1:length(This.sp))Table.Conf[[i]]=tabl.terms(subset(DATA,SPECIES==This.sp[i] & FL>0 & !is.na(FL) & FL<800 & !(is.na(BLOCK) | BLOCK==0 )))
    
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
      
      #model  
      if(!i%in%simpler.model)model<- glm(log(FL)~LAT+LONG+BOTDEPTH.bin, data=dat, family=gaussian, maxit=500)
      if(i%in%simpler.model)model<- glm(log(FL)~LAT+BOTDEPTH.bin, data=dat, family=gaussian, maxit=500)
      #model<- glm(log(FL)~Mid.Lat*Mid.Long+BOTDEPTH+Method, data=dat, family=gaussian, maxit=500)
      
      #Anova
      Signifcance1=anova(model,test="Chisq")
      
      #Deviance explained
      Dev.exp=Dsquared(model,adjust=F)$d3
      
      return(list(Anova=Signifcance1,Dev.exp=Dev.exp,data=dat,model=model)) 
    }
    simpler.model=match(c("BT","MI","PE","SO","TG","PN"),This.sp)
  }
  
  
  #compare size distribution for GN and LL used in the same latitudinal bands
  comp.meth.fn=function(dat,SNAME)
  {
    YMAX=max(dat$FL,na.rm=T)
    Tab=table(dat$LAT,dat$Method)
    Tab=ifelse(Tab>5,1,0)
    Tab=subset(Tab,rowSums(Tab)==2) #get lats where both gears used
    Common.lat=rownames(Tab)
    fn.dens=function(a)density(a,adjust=2,na.rm=T,from=10,to=YMAX)
    if(length(Tab)>0)
    {
      GN.Same.lat=fn.dens(subset(dat, Method=="GN" & as.numeric(as.character(LAT))%in%Common.lat)$FL)
      LL.Same.lat=fn.dens(subset(dat, Method=="LL" & as.numeric(as.character(LAT))%in%Common.lat)$FL)
      plot(GN.Same.lat,main="",ylim=c(0,max(c(GN.Same.lat$y,LL.Same.lat$y))),
           xlim=c(0,YMAX),lwd=3,xlab="",ylab="",cex.axis=1.25)
      lines(LL.Same.lat,col="grey60",lwd=3)
      mtext(SNAME,3,cex=1.15)    
    }
    #   if(length(Tab)==0)
    #   {
    #     plot(1:10,1:10,col="transparent",yaxt='n',xaxt='n',ann=F,main=This.sp[i])
    #     text(5,5,"no data for both gears")
    #   }
    
  }
  tiff("Size distr. BN LL same Lats.tiff",width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
  par(mfcol=c(3,3),mai=c(0.4,0.6,0.2,0.1),oma=c(1,1.2,0.15,0.01),las=1,mgp=c(1,.6,0))
  for(i in 1:length(This.sp))comp.meth.fn(subset(DATA,SPECIES==This.sp[i] & LAT <0),Species.labels.Scientific[i])
  plot(1,xaxt='n',yaxt='n',ylab='',xlab='',col="transparent")
  box(col="white")
  legend("center",c("Gillnet","Longline"),lty=1,lwd=3,col=c(1,"grey60"),bty="n",cex=2.5)
  mtext("Density",2,-.75,outer=T,cex=1.6,las=3)
  mtext("Fork length (cm)",1,-.5,outer=T,cex=1.6)
  dev.off()
  
  #explore data used in GLM
  fn.explr=function(dat)
  {
    #remove lat-long combos with less than 10 observations
    if(This.sp[i]=="PE")dat=subset(dat,!LAT==(-15))
    if(This.sp[i]=="BW")dat=subset(dat,!LONG==(121))
    if(This.sp[i]%in%c("WH","HZ","LG"))dat=subset(dat,!LONG%in%c(116,117))
    
    #put as factors
    par(mfcol=c(3,2))
    Cl=1:length(levels(as.factor(dat$zone)))
    with(dat,interaction.plot(as.factor(round(BOTDEPTH/10)*10),as.factor(zone),  log(FL),fixed = TRUE,col=Cl,main=This.sp[i]))
    plot(log(FL)~as.factor(round(BOTDEPTH/10)*10),dat,main="10")
    plot(log(FL)~as.factor(round(BOTDEPTH/25)*25),dat,main="25")
    with(dat,interaction.plot(as.factor(year),as.factor(zone),  log(FL),fixed = TRUE,col=Cl,main="annual trends by zone"))
    plot(log(FL)~as.factor(round(BOTDEPTH/25)*25),subset(dat,SEX=="F"),main="25 female",col="pink")
    plot(log(FL)~as.factor(round(BOTDEPTH/25)*25),subset(dat,SEX=="M"),main="25 male",col="blue")
    
  }
  for(i in 1:length(This.sp))fn.explr(subset(DATA,SPECIES==This.sp[i] & FL>0 & !is.na(FL) & FL<800 & !(is.na(BLOCK) | BLOCK==0 )))
  
  
  #2.2. Run GLM
  #note: change order of terms to see if d.f. change, if so, then confounding!!
  
  #set zone levels (then relevel within GLM function)
  Zone.order=c("Joint","North","Closed","West","Zone1","Zone2")
  DATA$zone=factor(DATA$zone,levels=Zone.order)
  
  #check number of observations for each combo of factors
  fn.chk.combo=function(d,depth.bin)
  {
    D=subset(DATA,SPECIES==d,select=c(BOTDEPTH,FL,zone,Method,SEX))
    D$Depth.bin=round(D$BOTDEPTH/depth.bin)*depth.bin
    datFL=subset(D,BOTDEPTH<=250 & FL>0 & !is.na(FL) & FL<800 & !is.na(zone))
    datSx=subset(D,SEX%in%c("F","M") & !is.na(zone))
    
    Tab.FL=with(datFL,table(Method,Depth.bin,zone))
    a=as.data.frame.table(Tab.FL, responseName = "value")
    wide=reshape(a,v.names ="value" , idvar = c("Method","zone"),
                 timevar = "Depth.bin", direction = "wide")
    colnames(wide)[3:ncol(wide)]=substr(colnames(wide)[3:ncol(wide)],7,15)
    wide.FL=cbind(data.frame(Species=d),wide)
    
    Tab.Sx=with(datSx,table(Method,Depth.bin,zone))
    a=as.data.frame.table(Tab.Sx, responseName = "value")
    wide=reshape(a,v.names ="value" , idvar = c("Method","zone"),
                 timevar = "Depth.bin", direction = "wide")
    colnames(wide)[3:ncol(wide)]=substr(colnames(wide)[3:ncol(wide)],7,15)
    wide.Sx=cbind(data.frame(Species=d),wide)
    
    return(list(FL=wide.FL,Sx=wide.Sx))
  }
  Store.chek.FL=Store.chek.Sx=vector('list',length(This.sp))
  for(i in 1:length(This.sp))
  {
    dummy=fn.chk.combo(This.sp[[i]],25)
    Store.chek.FL[[i]]=dummy$FL
    Store.chek.Sx[[i]]=dummy$Sx
  }
  
  
  
  if(Use.Lat.Long=="NO")
  {
    fn.size=function(dat,depth.bin)
    {
      dat=subset(dat,BOTDEPTH<=250 & FL>0 & !is.na(FL) & FL<800 & !is.na(zone))
      
      #put as factors
      dat$zone=factor(dat$zone)
      dat$BOTDEPTH.bin=as.factor(round(dat$BOTDEPTH/depth.bin)*depth.bin)
      dat$year=as.factor(dat$year)
      
      
      #model  
      model<- glm(log(FL)~zone+BOTDEPTH.bin, data=dat, family=gaussian, maxit=500)
      
      #Anova
      Signifcance1=anova(model,test="Chisq")
      
      #Deviance explained
      Dev.exp=Dsquared(model,adjust=F)$d3
      
      return(list(Anova=Signifcance1,Dev.exp=Dev.exp,data=dat,model=model)) 
    }
  }
  
  Output.GLM=vector('list',length(This.sp))
  names(Output.GLM)=This.sp
  Output.GLM.25=Output.GLM
  for(i in 1:length(Output.GLM))Output.GLM[[i]]=fn.size(subset(DATA,SPECIES==This.sp[i]),10)
  for(i in 1:length(Output.GLM))Output.GLM.25[[i]]=fn.size(subset(DATA,SPECIES==This.sp[i]),25)
  
  
  
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
  if(Use.Lat.Long=="YES")
  {
    ThisVar[[1]]=ThisVar[[2]]=ThisVar[[3]]=ThisVar[[4]]=ThisVar[[5]]=ThisVar[[6]]=c("LAT","LONG","BOTDEPTH.bin")
    ThisVar[[7]]=ThisVar[[8]]=ThisVar[[9]]=ThisVar[[10]]=ThisVar[[11]]=ThisVar[[12]]=c("LAT","BOTDEPTH.bin")  
  }
  if(Use.Lat.Long=="NO") for (i in 1:length(ThisVar)) ThisVar[[i]]=c("zone","BOTDEPTH.bin")
  
  
  Store.preds=vector('list',length(Output.GLM))
  par(mfcol=c(5,3),mai=c(0.4,0.4,0.12,0.1),oma=c(1,1.2,0.15,0.01))
  for(i in 1:length(Output.GLM))Store.preds[[i]]=pred.fn(Output.GLM[[i]]$data,Output.GLM[[i]]$model,ThisVar[[i]],This.sp[i])
  
  
  #2.3. Get term significance 
  ColNam=c("Species","P_Zone","P_Depth","Dev.exp_Zone","Dev.exp_Depth","Dev.exp_Tot")
  Term.Mat=matrix(nrow=length(Output.GLM),ncol=length(ColNam))
  colnames(Term.Mat)=ColNam
  Term.Mat=as.data.frame(Term.Mat)
  for(i in 1:length(Output.GLM))
  {
    Term.Mat$Species[i]=names(Output.GLM)[i]
    Anov=round(Output.GLM[[i]]$Anova$"Pr(>Chi)",3)
    Term.Mat$P_Zone[i]=Anov[2]
    Term.Mat$P_Depth[i]=Anov[3]  
    Dev=round(fn.dev.exp.term(Output.GLM[[i]]))
    Term.Mat$Dev.exp_Zone[i]=Dev[1]
    Term.Mat$Dev.exp_Depth[i]=Dev[2]
    Term.Mat$Dev.exp_Tot[i]=Dev[1]+Dev[2]
    #Term.Mat$Dev.exp_Tot[i]=round(Output.GLM[[i]]$Dev.exp)
  }
  Term.Mat$P_Zone=ifelse(Term.Mat$P_Zone==0,'<0.001',Term.Mat$P_Zone)
  Term.Mat$P_Depth=ifelse(Term.Mat$P_Depth==0,'<0.001',Term.Mat$P_Depth)
  
  Term.Mat=merge(Term.Mat,data.frame(Scientific=Species.Scient,Species=names(Species.Scient)),by="Species")
  Term.Mat=Term.Mat[ match(This.sp,Term.Mat$Species),]
  Term.Mat=Term.Mat[,match(c("Scientific","P_Zone","P_Depth","Dev.exp_Zone","Dev.exp_Depth","Dev.exp_Tot"),names(Term.Mat))]
  names(Term.Mat)[1]="Species"
  
  #Export as table
  Model.tbl=function(WD,Tbl,Doc.nm,caption,paragph,HdR.col,HdR.bg,Hdr.fnt.sze,Hdr.bld,
                     body.fnt.sze,Grid.col,Fnt.hdr,Fnt.body,
                     HDR.names,HDR.span,HDR.2nd)
  {
    mydoc = docx(Doc.nm)  #create r object
    mydoc = addSection( mydoc, landscape = T )   #landscape table
    # add title
    if(!is.na(caption))mydoc = addParagraph(mydoc, caption, stylename = "TitleDoc" )
    
    # add a paragraph
    if(!is.na(paragph))mydoc = addParagraph(mydoc , paragph, stylename="Citationintense")
    
    #add table
    MyFTable=FlexTable(Tbl,header.column=F,add.rownames =F,
                       header.cell.props = cellProperties(background.color=HdR.bg), 
                       header.text.props = textProperties(color=HdR.col,font.size=Hdr.fnt.sze,
                                                          font.weight="bold",font.family =Fnt.hdr), 
                       body.text.props = textProperties(font.size=body.fnt.sze,font.family =Fnt.body))
    
    #Add header
    MyFTable = addHeaderRow(MyFTable,text.properties=textBold(),value=HDR.names,colspan=HDR.span)
    
    #Add second header
    MyFTable = addHeaderRow(MyFTable, text.properties = textBold(),value =HDR.2nd)
    
    
    # table borders
    MyFTable = setFlexTableBorders(MyFTable,
                                   inner.vertical = borderNone(),inner.horizontal = borderNone(),
                                   outer.vertical = borderNone(),
                                   outer.horizontal = borderProperties(color=Grid.col, style="solid", width=4))
    
    # set columns widths (in inches)
    #MyFTable = setFlexTableWidths( MyFTable, widths = Col.width)
    
    mydoc = addFlexTable( mydoc, MyFTable)   
    mydoc = addSection( mydoc, landscape = F ) 
    
    # write the doc 
    writeDoc( mydoc, file = paste(Doc.nm,".docx",sep=''))
  }
  
  Model.tbl(WD=getwd(),Tbl=Term.Mat,Doc.nm="Model FL",caption=NA,paragph=NA,
            HdR.col='black',HdR.bg='white',Hdr.fnt.sze=12,Hdr.bld='normal',body.fnt.sze=12,
            Grid.col='black',Fnt.hdr= "Arial",Fnt.body= "Arial",
            HDR.names=c('Species','P','Deviance explained (%)'),
            HDR.span=c(1,2,3),HDR.2nd=c("","Zone","Depth","Zone","Depth","total"))
  
  
  
  #2.4. Show Spatial patterns
  
  Zone.mid.point="NO"  #show by block
  Zone.mid.Lat=c(-13.71, -16.97, -21.80, -29, -34, -34)
  Zone.mid.Long= c(125.63, 118.76, 111.83, 112.34, 113.43, 122.49) 
  Zone.mids=data.frame(zone=c("Joint","North","Closed","West","Zone1","Zone2"),
                       Zone.mid.Lat=Zone.mid.Lat,
                       Zone.mid.Long=Zone.mid.Long)
  
  fn.show.spatial.FL=function(dat,pred,WHAT)
  {
    if(Zone.mid.point=="YES")
    {
      dat=merge(dat,Zone.mids,by="zone",all.x=T)
      dat$LAT=dat$Zone.mid.Lat
      dat$LONG=dat$Zone.mid.Long
    }
    
    dat=subset(dat,!(Mid.Long>114.5342 & Mid.Long<115.3386 & Mid.Lat< (-27.13656) & Mid.Lat> (-28.42635)))
    dat=subset(dat,!(Mid.Long>114.7817 & Mid.Long<119.4223 & Mid.Lat< (-23.65974) & Mid.Lat> (-26.35147)))
    dat=subset(dat,!(Mid.Long>114.6579 & Mid.Long<115.5242 & Mid.Lat< (-22.70642) & Mid.Lat> (-23.37935)))
    
    dat$LONG=as.numeric(as.character((dat$LONG)))
    dat$LAT=as.numeric(as.character((dat$LAT)))
    
    pred$LONG=as.numeric(as.character((pred$LONG)))
    pred$LAT=as.numeric(as.character((pred$LAT)))
    
    
    #Observations
    if(WHAT=="obs") meanFL=aggregate(FL~LONG+LAT,dat,mean) #plot observations
    if(WHAT=="preds") meanFL=aggregate(FL~LONG+LAT,pred,mean)  #plot predictions
    
    plotMap(worldLLhigh, xlim=XAXIS,ylim=YAXIS,
            plt = NULL,col="grey95",tck = 0.025, tckMinor = 0.0125,
            xlab="",ylab="",axes=F)
    
    #add zones
    plot(WA_Northern_Shark,add=T,col="grey85")
    plot(JA_Northern_Shark,ylim=c(-39,-11),xlim=c(108,130),add=T,col="grey40")
    plot(WA_Northern_Shark_2,add=T,col=1,angle=0,density=seq(5,35,2))
    plot(WA_Northern_Shark_2,add=T,col=1,angle=90,density=seq(5,35,2))
    plot(WCDGDLL,add=T,col="grey70")
    plot(SDGDLL_zone1,add=T,col="white")
    plot(SDGDLL_zone2,add=T,col="grey55")
    
    lines(rbind(129,129),rbind(-15,-31.5),lty=2,col="black")
    points(meanFL$LONG,meanFL$LAT-0.5,cex=meanFL$FL/150,pch=19,col=1)
    
    text(109,-14,paste("N=",nrow(dat),sep=""),cex=.8,pos=4)
    axis(2,seq(YAXIS[1],YAXIS[2],1),F,tck=-0.03)
    axis(2,seq(YAXIS[1],YAXIS[2],4),F,tck=-0.065)
    axis(1,seq(XAXIS[1],XAXIS[2],1),F,tck=-0.03)
    axis(1,seq(XAXIS[1],XAXIS[2],4),F,tck=-0.065)  
    
    box() 
    
  }
  
  tiff("Figure 1.tiff",width = 1100, height = 1900,units = "px", res = 300,compression = "lzw")
  par(mfcol=c(4,3),mai=c(0.1,.15,0.125,.05),oma=c(2.25,2,0.1,0.05),las=1,mgp=c(1,.55,0))
  for(i in 1:length(Output.GLM))
  {
    fn.show.spatial.FL(Output.GLM[[i]]$data,Store.preds[[i]],'obs')
    
    if(names(Output.GLM)[i]=="BT")mtext(Species.labels.Scientific[i],3,cex=0.6)
    if(!names(Output.GLM)[i]=="BT")mtext(Species.labels.Scientific[i],3,cex=0.7)
    if(i%in%c(1:4))  axis(2,seq(YAXIS[1],YAXIS[2],4),-seq(YAXIS[1],YAXIS[2],4),cex.axis=1.15,las=2,tck=-0.08)
    if(i%in%c(4,8,12)) axis(1,seq(XAXIS[1],XAXIS[2],4),seq(XAXIS[1],XAXIS[2],4),cex.axis=1.15,tck=-0.08)
    
    if(i==1)
    {
      legend(111,-21,"Western",bty='n',cex=1.25)
      legend(111,-24,"Australia",bty='n',cex=1.25)
    }
    
    # add size reference
    if(i==9)
    {
      points(rep(117,3),c(-23,-26,-29),cex=c(100,200,300)/150,pch=19,col=1)
      text(rep(117,3),c(-23,-26,-29),c("100 cm","200 cm","300 cm"),pos=4)
    }
    
    # add zone reference
    if(i==11)
    {
      PoSi=seq(-32.25,-21.5,length.out=6)
      Zn.leg=rev(c("JANSF","WANCSF","Closure","WCDGDLF","Zone 1","Zone 2"))
      Zn.col=rev(c("grey40","grey85",NA,"grey70","white","grey55"))
      pCH=rev(c(22,22,12,22,22,22))
      for(t in 1:length(Zn.leg))
      {
        points(117,PoSi[t],pch=pCH[t],bg=Zn.col[t],cex=1.35)
        text(117,PoSi[t],Zn.leg[t],pos=4,cex=0.9)
      }
    }
  }
  
  mtext("Longitude (?E) ",side=1,line=0.95,font=1,las=0,cex=1.15,outer=T)
  mtext("Latitude (?S)",side=2,line=0.725,font=1,las=0,cex=1.15,outer=T)
  dev.off()
  
  
  
  #coefs
  fn.plot.coef=function(Mod,dat,what)
  {
    all.coef=summary(Mod)$coefficients
    id=grep(what,row.names(all.coef))
    
    if(length(id)>0)
    {
      x=1:length(c(0,all.coef[id,1]))
      y=c(0,all.coef[id,1])
      CI=2*c(0,all.coef[id,2])
      YLIM=c(min(y-CI)*1.25,max(y+CI)*1.25)
      plot(x,y,ylab="",xlab="",xaxt='n',ylim=YLIM,pch=19)
      segments(x,y-CI,x,y+CI)
      axis(1,x,dat)
    }
    
  }
  
  tiff("Coef_Zone_FL.tiff",width = 1500, height = 2000,units = "px", res = 300,compression = "lzw")
  par(mfcol=c(4,3),mai=c(0.2,.3,0.18,.1),oma=c(2.25,2,0.15,0.1),las=1,mgp=c(1,.5,0))
  
  for(i in 1:length(Output.GLM))
  {
    fn.plot.coef(Output.GLM[[i]]$model,levels(Output.GLM[[i]]$data$zone),"zone")
    mtext(Species.labels[i],3,cex=0.9)
    
  }
  mtext("Zone",1,line=1,outer=T,cex=1.5)
  mtext("Coefficient",2,line=0,outer=T,las=3,cex=1.5)
  dev.off()
  
  
  #show variability in data
  Zone.names=Zone.order
  names(Zone.names)=c("JANSF","WANCSF","Closure","WCDGDLF","Zone 1","Zone 2")
  
  fn.show.spatial.FL_obs=function(dat,sP)
  {
    #NMs=names(Zone.names[match(levels(dat$zone),Zone.names)])
    NMs=levels(dat$zone)
    boxplot(FL~zone,dat,main="",col='grey80',names=NMs,cex.axis=.9)
    mtext(sP,3,line=0,cex=1.1)
  }
  
  tiff("Figure S1.tiff",width = 2400, height = 2000,units = "px", res = 300,compression = "lzw")
  par(mfcol=c(4,3),mai=c(.2,0.2,.175,.04),oma=c(2,1.5,.2,.01),mgp=c(1,.4,0),las=1)
  for(i in 1:length(This.sp))fn.show.spatial.FL_obs(Output.GLM[[i]]$data,Species.labels.Scientific[i])
  mtext("Fork length (cm)",2,line=0.15,las=3,cex=1.25,outer=T)
  mtext("Zone",1,line=0.5,cex=1.25,outer=T)
  dev.off() 
  
  
  #2.5 Show Depth effects
  
  #show variability in data
  fn.show.spatial.FL.depth_obs=function(dat,sP)
  {
    dat=subset(dat,!is.na(BOTDEPTH) & BOTDEPTH<300 & BOTDEPTH>0)
    dat$BOTDEPTH.bin=as.factor(round(dat$BOTDEPTH/10)*10)
    
    boxplot(FL~BOTDEPTH.bin,dat,main="",col='grey80')
    mtext(sP,3,line=0,cex=1.1)
  }
  
  tiff("Figure S2.tiff",width = 2000, height = 2400,units = "px", res = 300,compression = "lzw")
  par(mfcol=c(4,3),mai=c(.2,0.2,.175,.04),oma=c(2,1.5,.2,.01),mgp=c(1,.4,0),las=1)
  for(i in 1:length(This.sp))fn.show.spatial.FL.depth_obs(Output.GLM[[i]]$data,Species.labels.Scientific[i])
  mtext("Fork length (cm)",2,line=0.15,las=3,cex=1.25,outer=T)
  mtext("Depth (m)",1,line=0.5,cex=1.25,outer=T)
  dev.off() 
  
  CI.approach="NO"
  if(CI.approach=="YES")
  {
    fn.show.spatial.FL.depth=function(dat,pred,WHAT)
    {
      dat=subset(dat,!is.na(BOTDEPTH) & BOTDEPTH<300 & BOTDEPTH>0)
      dat$BOTDEPTH.bin=as.factor(round(dat$BOTDEPTH/10)*10)
      
      #Observations
      if(WHAT=="obs") meanFL=aggregate(FL~BOTDEPTH.bin,dat,mean) #plot observations
      if(WHAT=="preds") meanFL=aggregate(FL~BOTDEPTH.bin,pred,mean)  #plot predictions
      
      x=as.numeric(as.character(meanFL$BOTDEPTH.bin))
      y=meanFL$FL
      
      Unic=as.character(sort(unique(meanFL$BOTDEPTH.bin)))
      CI=data.frame(Low=rep(NA,length(Unic)),UP=NA)
      sub=function(d)quantile(d$FL,c(0.025,0.975),na.rm=T)
      for (t in 1:length(Unic))CI[t,]=sub(subset(dat,BOTDEPTH.bin==Unic[t]))
      Ymax=max(CI)*1.2
      plot(x,y,ylim=c(0,Ymax),xlim=c(0,max(dat$BOTDEPTH)*1.1),ylab="",xlab="",xaxt='n',yaxt='n',pch=19,col=1)
      segments(x,CI$Low,x,CI$UP)
      axis(2,seq(0,Ymax,50),seq(0,Ymax,50),cex.axis=1.25)
      axis(1,seq(0,300,50),seq(0,300,50),cex.axis=1.25)
      axis(2,seq(0,Ymax,10),F,tck=-0.02)
      axis(1,seq(0,300,10),F,tck=-0.02)
      
    }
    tiff("Figure S2.tiff",width = 2000, height = 2400,units = "px", res = 300,compression = "lzw")
    par(mfcol=c(5,3),mai=c(.2,0.25,.175,.1),oma=c(2,2,.2,.1),mgp=c(1,.5,0),las=1)
    for(i in 1:length(Output.GLM))
    {
      fn.show.spatial.FL.depth(Output.GLM[[i]]$data,Store.preds[[i]],'obs')
      mtext(Species.labels.Scientific[i],3,line=0,cex=1.1)
    }
    mtext("Fork length (cm)",2,line=0.45,las=3,cex=1.25,outer=T)
    mtext("Depth (m)",1,line=0.5,cex=1.25,outer=T)
    dev.off()
    
  }
  
  
  #coefs
  tiff("Coef_Depth_FL.tiff",width = 1500, height = 2000,units = "px", res = 300,compression = "lzw")
  par(mfcol=c(4,3),mai=c(0.2,.3,0.18,.1),oma=c(2.25,2,0.15,0.1),las=1,mgp=c(1,.5,0))
  
  for(i in 1:length(Output.GLM))
  {
    fn.plot.coef(Output.GLM[[i]]$model,levels(Output.GLM[[i]]$data$BOTDEPTH.bin),"BOTDEPTH.bin")
    mtext(Species.labels[i],3,cex=0.9)
    
  }
  mtext("Depth (m)",1,line=1,outer=T,cex=1.5)
  mtext("Coefficient",2,line=0,outer=T,las=3,cex=1.5)
  dev.off()
  
  
  
  #3. Spatial patterns in sex ratio
  
  #3.1 Prelim analyses
  
  #Select species for analysis (at least 100 observations overall, & > 10 observation in at least 2 zones)
  This.sp=names(Table.species)
  Selected=vector('list',length(This.sp))
  names(Selected)=This.sp
  for(i in 1:length(This.sp))Selected[[i]]=Select.fn(subset(DATA,SPECIES==This.sp[[i]] & year<=Current.year & SEX%in%c("M","F") & !(is.na(BLOCK) | BLOCK==0 )))
  ID=which(Selected=="YES")
  
  This.sp=names(ID)
  
  #sort species by importance
  a=subset(DATA,SPECIES%in%This.sp)
  SP.sort=rev(sort(table(a$SPECIES)))
  rm(a)
  This.sp=This.sp[match(names(SP.sort),This.sp)]
  
  Species.labels=c("Dusky shark","Sandbar shark","Gummy shark","Whiskery shark","Smooth hammerhead",
                   "Spinner shark","Spot-tail shark","Milk shark","Blacktip sharks","Tiger shark",
                   "Western Wobbegong","Pencil Shark","Banded Wobbegong","Scalloped Hammerhead")
  
  Comm.Species.labels=c("Dusky shark","Sandbar shark","Gummy shark","Whiskery shark")
  
  names(This.sp)=Species.labels
  Species.Scient=c("C. obscurus","C. plumbeus","M. antarcticus","F. macki",
                   "S. zygaena","C. brevipinna","C. sorrah","R. acutus","C. limbatus & C. tilstoni",
                   "G. cuvier","Orectolobus sp.","H. hyugaensis","Orectolobus ornatus","S. lewini")
  names(Species.Scient)=This.sp
  
  Species.labels.Scientific=c(expression(italic("C. obscurus")),expression(italic("C. plumbeus")),
                              expression(italic("M. antarcticus")),expression(italic("F. macki")),
                              expression(italic("S. zygaena")),expression(italic("C. brevipinna")),
                              expression(italic("C. sorrah")),expression(italic("R. acutus")),
                              expression(italic("C. limbatus & C. tilstoni")),expression(italic("G. cuvier")),
                              expression(italic("Orectolobus sp.")),expression(italic("H. hyugaensis")),
                              expression(italic("O. ornatus")),expression(italic("S. lewini")))
  
  
  Small.scale.sex=function(dat,sP)
  {
    dat=subset(dat,!SEX=="U")
    dat$Month=as.factor(dat$Month)
    dat$SEX=as.factor(dat$SEX)
    Ag=aggregate(Number~Month+SEX,dat,sum)
    wide <- reshape(Ag, v.names = "Number", idvar = "Month",
                    timevar = "SEX", direction = "wide")
    wide$Number.M=with(wide,ifelse(is.na(Number.M),0,Number.M))
    wide$prop=(wide$Number.M)/(wide$Number.F+wide$Number.M)
    # plot(wide$prop,ylim=c(0,1),xlim=c(1,12),type='h',main=sP,cex.main=1.25,ylab="",xlab="",lwd=4)
    plot(wide$prop,ylim=c(0,1),xlim=c(1,12),type='h',ylab="",xlab="",lwd=2.5)
    abline(h=0.5,lty=2,col="grey50")
    #legend("topleft",sP,bty='n',cex=1.35,adj=c(0.1,0))
    mtext(sP,3,line=0,cex=1.1)
    box()
  }
  
  tiff("Figure S6.tiff",width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
  par(mfcol=c(5,3),mai=c(.2,0.25,.175,.1),oma=c(2,2,.2,.1),mgp=c(1,.5,0),las=1)
  
  for(i in 1:length(This.sp))Small.scale.sex(subset(DATA,SPECIES==This.sp[[i]]  &
                                                      !(is.na(BLOCK) | BLOCK==0 )),Species.labels.Scientific[i])
  mtext("Proportion of males",2,line=0.25,las=3,cex=1.25,outer=T)
  mtext("Month",1,line=0.5,cex=1.25,outer=T)
  dev.off() 
  
  
  #Colinearity
  if(Use.Lat.Long=="YES") for(i in 1:length(This.sp))fn.colinearity(subset(DATA,SPECIES==This.sp[i] & SEX%in%c("M","F")))
  
  #3.2. Run GLM
  if(Use.Lat.Long=="YES")
  {
    simpler.model=match(c("LG","BT","MI","PE","SO","TG"),This.sp)
    fn.sex=function(dat)
    {
      dat$Number=1
      dat$BOTDEPTH.bin=round(dat$BOTDEPTH/10)*10
      
      #remove lat-long combos with less than 10 observations
      dat$PASTED=paste(dat$LAT,dat$LONG)
      These.Lats.Longs=table(dat$PASTED)
      These.Lats.Longs=These.Lats.Longs[These.Lats.Longs>10]
      dat=subset(dat,PASTED%in%names(These.Lats.Longs))
      
      #get number of males and females
      if(!i%in%simpler.model)
      {
        agg.Male=aggregate(Number~LONG+LAT+BOTDEPTH.bin,subset(dat,SEX=="M"),sum)
        names(agg.Male)[4]="N.male"
        agg.Fem=aggregate(Number~LONG+LAT+BOTDEPTH.bin,subset(dat,SEX=="F"),sum)
        names(agg.Fem)[4]="N.female"
        agg.sex=merge(agg.Male,agg.Fem,by=c("LONG","LAT","BOTDEPTH.bin"),all.x=T,all.y=T)
        agg.sex$N.male=with(agg.sex,ifelse(is.na(N.male),0,N.male))
        agg.sex$N.female=with(agg.sex,ifelse(is.na(N.female),0,N.female))
        agg.sex$Sum=agg.sex$N.male+agg.sex$N.female
        agg.sex=subset(agg.sex,Sum>4)
        
        #put as factor
        agg.sex$BOTDEPTH.bin=as.factor(agg.sex$BOTDEPTH)
        agg.sex$LONG=as.factor(agg.sex$LONG)
        agg.sex$LAT=as.factor(agg.sex$LAT)
        
        #model
        model<- glm(cbind(N.male,N.female)~LAT+LONG+BOTDEPTH.bin, data=agg.sex, family="binomial", maxit=500)
      }
      
      if(i%in%simpler.model)
      {
        agg.Male=aggregate(Number~LAT+BOTDEPTH.bin,subset(dat,SEX=="M"),sum)
        names(agg.Male)[3]="N.male"
        agg.Fem=aggregate(Number~LAT+BOTDEPTH.bin,subset(dat,SEX=="F"),sum)
        names(agg.Fem)[3]="N.female"
        agg.sex=merge(agg.Male,agg.Fem,by=c("LAT","BOTDEPTH.bin"),all.x=T,all.y=T)
        agg.sex$N.male=with(agg.sex,ifelse(is.na(N.male),0,N.male))
        agg.sex$N.female=with(agg.sex,ifelse(is.na(N.female),0,N.female))
        agg.sex$Sum=agg.sex$N.male+agg.sex$N.female
        agg.sex=subset(agg.sex,Sum>4)
        
        #put as factor
        agg.sex$BOTDEPTH.bin=as.factor(agg.sex$BOTDEPTH)
        agg.sex$LAT=as.factor(agg.sex$LAT)
        
        #model
        model<- glm(cbind(N.male,N.female)~LAT+BOTDEPTH.bin, data=agg.sex, family="binomial", maxit=500)
      }
      
      #Anova
      Signifcance1=anova(model,test="Chisq")
      
      #Deviance explained
      Dev.exp=Dsquared(model,adjust=F)$d3
      
      return(list(Anova=Signifcance1,Dev.exp=Dev.exp,data=dat,model=model)) 
    }
  }
  
  if(Use.Lat.Long=="NO")
  {
    fn.sex=function(dat,depth.bin)
    {
      dat=subset(dat,BOTDEPTH<=250 & !is.na(zone))
      dat$Number=1
      dat$BOTDEPTH.bin=round(dat$BOTDEPTH/depth.bin)*depth.bin
      
      agg.Male=aggregate(Number~zone+BOTDEPTH.bin,subset(dat,SEX=="M"),sum)
      names(agg.Male)[match("Number",names(agg.Male))]="N.male"
      agg.Fem=aggregate(Number~zone+BOTDEPTH.bin,subset(dat,SEX=="F"),sum)
      names(agg.Fem)[match("Number",names(agg.Fem))]="N.female"
      agg.sex=merge(agg.Male,agg.Fem,by=c("zone","BOTDEPTH.bin"),all.x=T,all.y=T)
      agg.sex$N.male=with(agg.sex,ifelse(is.na(N.male),0,N.male))
      agg.sex$N.female=with(agg.sex,ifelse(is.na(N.female),0,N.female))
      agg.sex$Sum=agg.sex$N.male+agg.sex$N.female
      agg.sex=subset(agg.sex,Sum>4)
      
      #put as factor
      agg.sex$BOTDEPTH.bin=as.factor(agg.sex$BOTDEPTH)
      agg.sex$zone=factor(agg.sex$zone)
      
      
      #model
      model<- glm(cbind(N.male,N.female)~zone+BOTDEPTH.bin, data=agg.sex, family="binomial", maxit=500)
      
      #Anova
      Signifcance1=anova(model,test="Chisq")
      
      #Deviance explained
      Dev.exp=Dsquared(model,adjust=F)$d3
      
      return(list(Anova=Signifcance1,Dev.exp=Dev.exp,data=dat,model=model))     
    }
  }
  
  Output.GLM.sex=vector('list',length(This.sp))
  names(Output.GLM.sex)=This.sp
  
  for(i in 1:length(Output.GLM.sex))Output.GLM.sex[[i]]=fn.sex(subset(DATA,SPECIES==This.sp[i] & BOTDEPTH>0 & SEX%in%c("M","F")),10)
  
  
  #Get term significance 
  ColNam=c("Species","P_Zone","P_Depth","Dev.exp_Zone","Dev.exp_Depth","Dev.exp_Tot")
  Term.Mat=matrix(nrow=length(Output.GLM.sex),ncol=length(ColNam))
  colnames(Term.Mat)=ColNam
  Term.Mat=as.data.frame(Term.Mat)
  
  for(i in 1:length(Output.GLM.sex))
  {
    Term.Mat$Species[i]=names(Output.GLM.sex)[i]
    Anov=round(Output.GLM.sex[[i]]$Anova$"Pr(>Chi)",3)
    Term.Mat$P_Zone[i]=Anov[2]
    Term.Mat$P_Depth[i]=Anov[3]  
    Dev=round(fn.dev.exp.term(Output.GLM.sex[[i]]))
    Term.Mat$Dev.exp_Zone[i]=Dev[1]
    Term.Mat$Dev.exp_Depth[i]=Dev[2]
    Term.Mat$Dev.exp_Tot[i]=Dev[1]+Dev[2]
  }
  Term.Mat$P_Zone=ifelse(Term.Mat$P_Zone==0,'<0.001',Term.Mat$P_Zone)
  Term.Mat$P_Depth=ifelse(Term.Mat$P_Depth==0,'<0.001',Term.Mat$P_Depth)
  
  Term.Mat=merge(Term.Mat,data.frame(Scientific=Species.Scient,Species=names(Species.Scient)),by="Species")
  Term.Mat=Term.Mat[ match(This.sp,Term.Mat$Species),]
  Term.Mat=Term.Mat[,match(c("Scientific","P_Zone","P_Depth","Dev.exp_Zone","Dev.exp_Depth","Dev.exp_Tot"),names(Term.Mat))]
  names(Term.Mat)[1]="Species"
  Model.tbl(WD=getwd(),Tbl=Term.Mat,Doc.nm="Model Sex",caption=NA,paragph=NA,
            HdR.col='black',HdR.bg='white',Hdr.fnt.sze=12,Hdr.bld='normal',body.fnt.sze=12,
            Grid.col='black',Fnt.hdr= "Arial",Fnt.body= "Arial",
            HDR.names=c('Species','P','Deviance explained (%)'),
            HDR.span=c(1,2,3),HDR.2nd=c("","Zone","Depth","Zone","Depth","total"))
  
  
  #3.3. Show Latitude and longitude effect
  Zone.mid.point="NO"   #show only mid point of zone instead of by block
  #Zone.mid.point="YES"  
  fn.show.spatial.sex=function(dat)
  {
    Rad=0.7
    if(Zone.mid.point=="YES")
    {
      dat=merge(dat,Zone.mids,by="zone",all.x=T)
      dat$LAT=dat$Zone.mid.Lat
      dat$LONG=dat$Zone.mid.Long
      Rad=1.75
    }
    agg.Male=aggregate(Number~LONG+LAT+BOTDEPTH.bin,subset(dat,SEX=="M"),sum)
    names(agg.Male)[4]="N.male"
    agg.Fem=aggregate(Number~LONG+LAT+BOTDEPTH.bin,subset(dat,SEX=="F"),sum)
    names(agg.Fem)[4]="N.female"
    agg.sex=merge(agg.Male,agg.Fem,by=c("LONG","LAT","BOTDEPTH.bin"),all.x=T,all.y=T)
    agg.sex$N.male=with(agg.sex,ifelse(is.na(N.male),0,N.male))
    agg.sex$N.female=with(agg.sex,ifelse(is.na(N.female),0,N.female))
    agg.sex$Sum=agg.sex$N.male+agg.sex$N.female
    agg.sex=subset(agg.sex,Sum>4)
    
    agg.sex$LONG=as.numeric(as.character(agg.sex$LONG))
    agg.sex$LAT=as.numeric(as.character(agg.sex$LAT))
    
    #Observations  
    plotMap(worldLLhigh, xlim=XAXIS,ylim=YAXIS,
            plt = NULL,col="grey80",tck = 0.025, tckMinor = 0.0125,
            xlab="",ylab="",axes=F)
    lines(rbind(129,129),rbind(-15,-31.5),lty=2,col="black")
    
    #add zones
    plot(WA_Northern_Shark,add=T,col="grey85")
    plot(JA_Northern_Shark,ylim=c(-39,-11),xlim=c(108,130),add=T,col="grey40")
    plot(WA_Northern_Shark_2,add=T,col=1,angle=0,density=seq(5,35,2))
    plot(WA_Northern_Shark_2,add=T,col=1,angle=90,density=seq(5,35,2))
    plot(WCDGDLL,add=T,col="grey70")
    plot(SDGDLL_zone1,add=T,col="white")
    plot(SDGDLL_zone2,add=T,col="grey55")
    
    
    for(j in 1:nrow(agg.sex))
    {
      floating.pie(agg.sex$LONG[j],agg.sex$LAT[j]-1,c(agg.sex$N.male[j]+.1,agg.sex$N.female[j]+.1),
                   radius=Rad,col=c("white","black"))
    }
    axis(2,seq(YAXIS[1],YAXIS[2],1),F,tck=-0.03)
    axis(2,seq(YAXIS[1],YAXIS[2],4),F,tck=-0.08)
    axis(1,seq(XAXIS[1],XAXIS[2],1),F,tck=-0.03)
    axis(1,seq(XAXIS[1],XAXIS[2],4),F,tck=-0.08)
    
    text(109,-14,paste("N=",nrow(dat),sep=""),cex=.75,pos=4)
    box() 
  } 
  
  tiff("Figure 2.tiff",width = 950, height = 2000,units = "px", res = 300,compression = "lzw")
  par(mfcol=c(5,3),mai=c(0.1,.1,0.18,.01),oma=c(2.25,2,0.15,0.01),las=1,mgp=c(1,.5,0))
  for(i in 1:length(Output.GLM.sex))
  {
    fn.show.spatial.sex(Output.GLM.sex[[i]]$data)
    
    mtext(Species.labels.Scientific[i],3,cex=0.75)
    if(i%in%c(1:5))  axis(2,seq(YAXIS[1],YAXIS[2],4),-seq(YAXIS[1],YAXIS[2],4),cex.axis=1.15,las=2,tck=-0.08)
    if(i%in%c(5,10,14)) axis(1,seq(XAXIS[1],XAXIS[2],4),seq(XAXIS[1],XAXIS[2],4),cex.axis=1.15,tck=-0.08)
    
    if(i==1)
    {
      legend(111,-20,"Western",bty='n',cex=1)
      legend(111,-23,"Australia",bty='n',cex=1)
    }
    
  }
  
  #add scale
  plot(1:10,col="transparent",xaxt='n',yaxt='n',ann=F)
  Zn.leg=rev(c("JANSF","WANCSF","Closure","WCDGDLF","Zone 1","Zone 2"))
  Zn.col=rev(c("grey40","grey85",NA,"grey70","white","grey55"))
  pCH=rev(c(22,22,12,22,22,22))
  for(t in 1:length(Zn.leg))
  {
    points(2,3.5+t,pch=pCH[t],bg=Zn.col[t],cex=1.85)
    text(2.25,3.5+t,Zn.leg[t],pos=4,cex=0.9)
  }
  box(col="white")
  
  mtext("Longitude (?E) ",side=1,line=0.95,font=1,las=0,cex=1.15,outer=T)
  mtext("Latitude (?S)",side=2,line=0.55,font=1,las=0,cex=1.15,outer=T)
  dev.off()
  
  
  
  
  
  #3.4. Depth
  fn.show.spatial.sex.depth=function(dat,sP)
  {
    agg.Male=aggregate(Number~BOTDEPTH.bin,subset(dat,SEX=="M"),sum)
    names(agg.Male)[2]="N.male"
    agg.Fem=aggregate(Number~BOTDEPTH.bin,subset(dat,SEX=="F"),sum)
    names(agg.Fem)[2]="N.female"
    agg.sex=merge(agg.Male,agg.Fem,by=c("BOTDEPTH.bin"),all.x=T,all.y=T)
    agg.sex$N.male=with(agg.sex,ifelse(is.na(N.male),0,N.male))
    agg.sex$N.female=with(agg.sex,ifelse(is.na(N.female),0,N.female))
    agg.sex$Sum=agg.sex$N.male+agg.sex$N.female
    agg.sex=subset(agg.sex,Sum>4)
    agg.sex=agg.sex[order(agg.sex$BOTDEPTH.bin),]
    agg.sex$BOTDEPTH.bin=as.numeric(as.character(agg.sex$BOTDEPTH.bin))
    agg.sex$prop=agg.sex$N.male/agg.sex$Sum
    x=agg.sex$BOTDEPTH.bin
    y=agg.sex$prop
    plot(x,y,ylim=c(0,1),xlim=c(0,max(x)*1.1),ylab="",xlab="",type='h',xaxt='n',yaxt='n',lwd=2.5,pch=19,col=1,cex=1.5)
    abline(h=0.5,lty=2,col="grey50")
    axis(2,seq(0,1,.2),F,cex.axis=1.25)
    axis(1,seq(0,300,50),seq(0,300,50),cex.axis=1.25)
    axis(1,seq(0,300,10),F,tck=-0.02)
    mtext(sP,3,line=0,cex=0.95)
  }
  
  tiff("Figure S3.tiff",width = 1300, height = 2200,units = "px", res = 300,compression = "lzw")
  par(mfcol=c(5,3),mai=c(0.2,.1,0.175,.01),oma=c(2,3,0.15,0.01),las=1,mgp=c(2.5,.6,0))
  for(i in 1:length(Output.GLM.sex))
  {
    fn.show.spatial.sex.depth(Output.GLM.sex[[i]]$data,Species.labels.Scientific[i])
    if(i %in% 1:5)   axis(2,seq(0,1,.2),seq(0,1,.2),cex.axis=1.25)
    
  }
  mtext("Depth (m)",side=1,line=0.7,font=1,las=0,cex=1.25,outer=T)
  mtext("Proportion of males",side=2,line=1.45,font=1,las=0,cex=1.25,outer=T)
  dev.off()
  
  
  
  
  
  
  
  #NOTE USED #####
  
  #1.2 Test Block effect
  
  #Select appropriate data for comparisons
  fn.subset=function(dat)
  {
    Bk.yr=table(dat$BLOCK,dat$year)
    Bk.yr=ifelse(Bk.yr>1,1,0)
    Sum.Yr=colSums(Bk.yr)
    ID=which(Sum.Yr>=5)  #drop years with less than 5 occurrences
    dat=subset(dat,year%in%names(ID))
    if(nrow(dat)>1)
    {
      Bk.Mn=table(dat$BLOCK,dat$Month)
      Bk.Mn=ifelse(Bk.Mn>1,1,0)
      Sum.Mn=colSums(Bk.Mn)
      ID=which(Sum.Mn>=5)
      dat=subset(dat,Month%in%names(ID))
    }
    
    return(dat)
  }
  
  
  This.sp=This.sp[-match(c("LE","PJ","SD","WE","ES"),This.sp)] #Drop species with no data after criteria
  
  Data.GLM=vector('list',length(This.sp))
  names(Data.GLM)=This.sp
  for(i in 1:length(Data.GLM))Data.GLM[[i]]=fn.subset(subset(DATA,SPECIES%in%This.sp[i]& SEX%in%c("M","F")))
  
  MaxDepth=250
  Sex.fn=function(dat)
  {
    Mn=sort(unique(dat$Month))
    DaT=vector('list',length(Mn))
    names(DaT)=Mn
    ANOV=data.frame(Month=Mn,p=NA,N.blk=NA)
    for (x in 1:length(Mn))
    {
      dat1=subset(dat,!(is.na(BLOCK) | BLOCK==0 ) & Month==Mn[x])
      Table=table(dat1$BLOCK)
      ID=which(Table>=5)
      dat1=subset(dat1,BLOCK%in%names(Table[ID]))
      dat1$BLOCK=as.factor(dat1$BLOCK)
      
      agg.Male=aggregate(Number~BLOCK,subset(dat1,SEX=="M"),sum)
      names(agg.Male)[2]="N.male"
      agg.Fem=aggregate(Number~BLOCK,subset(dat1,SEX=="F"),sum)
      names(agg.Fem)[2]="N.female"
      agg.sex=merge(agg.Male,agg.Fem,by=c("BLOCK"),all.x=T,all.y=T)
      agg.sex$N.male=with(agg.sex,ifelse(is.na(N.male),0,N.male))
      agg.sex$N.female=with(agg.sex,ifelse(is.na(N.female),0,N.female))
      agg.sex=merge(agg.sex,Blk.pos,by="BLOCK",all.x=T)
      DaT[[x]]=agg.sex
      ANOV$N.blk[x]=nrow(agg.sex)
      #glm
      if(nrow(agg.sex)>1)
      {
        model<- glm(cbind(N.male,N.female)~BLOCK, data=agg.sex, family="binomial", maxit=500)
        
        #Anova table
        ANOV$p[x]=round(anova(model,test="Chisq")$"Pr(>Chi)"[2],3)
        
      }
      
    }
    return(list(ANOV=ANOV,DaT=DaT))
  }
  
  SEX.ratio.GN=vector('list',length(This.sp))
  names(SEX.ratio.GN)=This.sp
  for(i in 1:length(This.sp))SEX.ratio.GN[[i]]=Sex.fn(Data.GLM[[i]])
  
  #get p values for each species
  Months=1:12
  see.p=function(dat)
  {
    dat$Species=names(SEX.ratio.GN)[i]
    
    ID=which(!Months%in%dat$Month)
    if(length(ID)>0)
    {
      Add=data.frame(Month=Months[ID],p=NA, N.blk=NA, Species=names(SEX.ratio.GN)[i])
      dat=rbind(dat,Add)
      dat=dat[order(dat$Month),]
    }
    return(dat)
  }
  Store.p=SEX.ratio.GN
  for(i in 1:length(SEX.ratio.GN))Store.p[[i]]=see.p(SEX.ratio.GN[[i]]$ANOV)
  Store.p=do.call(rbind,Store.p)
  Store.p.wide.p=reshape(Store.p[,c(1,2,4)],idvar ="Species",v.names=c("p"),timevar="Month",direction="wide")
  Store.p.wide.N=reshape(Store.p[,c(1,3,4)],idvar ="Species",v.names=c("N.blk"),timevar="Month",direction="wide")
  
  Store.p.wide.p=merge(Store.p.wide.p,Tabl1.matrix[,c(3,9)],by.x="Species",by.y="Code",all.x=T)
  Store.p.wide.p=Store.p.wide.p[order(-Store.p.wide.p$N),]
  
  write.csv(Store.p.wide.p,"Sex.ratio.p.values.csv",row.names=F)
  
  #ID for which species model is important and significant
  #note: keep only those explaining more than 20% of the deviance
  #Important=rep(NA,length(SEX.ratio.GN))
  #for(i in 1:length(SEX.ratio.GN))Important[i]=(SEX.ratio.GN[[i]]$Dev.exp)
  
  
  
  
  
  #Show term effects 
  
  #1. Latitude and longitude effect
  Species.labels=c("Blacktip shark", "Dusky shark","Gummy shark","Smooth hammerhead","Spinner shark",
                   "Milk shark","Pencil shark","Tiger shark","Sandbar shark","Whiskery shark")
  Comm.Species.labels=c("Dusky shark","Sandbar shark","Gummy shark","Whiskery shark")
  ID=match(c("BW","TK","GM","WH"),names(SEX.ratio.GN))
  Mnth=c(11,5,3,3)
  
  MNTH=c("November","May","March","March")
  
  fn.show.block=function(dat)
  {
    
    #Observations
    dat1=subset(dat,!(is.na(BLOCK) | BLOCK==0 ))
    Table=table(dat1$BLOCK)
    ID=which(Table>=5)
    dat1=subset(dat1,BLOCK%in%names(Table[ID]))
    dat1$BLOCK=as.factor(dat1$BLOCK)
    
    agg.Male=aggregate(Number~BLOCK,subset(dat1,SEX=="M"),sum)
    names(agg.Male)[2]="N.male"
    agg.Fem=aggregate(Number~BLOCK,subset(dat1,SEX=="F"),sum)
    names(agg.Fem)[2]="N.female"
    agg.sex=merge(agg.Male,agg.Fem,by=c("BLOCK"),all.x=T,all.y=T)
    agg.sex$N.male=with(agg.sex,ifelse(is.na(N.male),0,N.male))
    agg.sex$N.female=with(agg.sex,ifelse(is.na(N.female),0,N.female))
    agg.sex=merge(agg.sex,Blk.pos,by="BLOCK",all.x=T)
    agg.sex$Prop.male=(agg.sex$N.male+0.1)/(agg.sex$N.male+agg.sex$N.female)
    
    plotMap(worldLLhigh, xlim=XAXIS,ylim=YAXIS,
            plt = NULL,col="grey70",tck = 0.025, tckMinor = 0.0125,
            xlab="",ylab="",axes=F)
    lines(rbind(129,129),rbind(-15,-31.5),lty=2,col="grey50")
    #points((agg.sex$LONG),(agg.sex$LAT-1),col=1,pch=19,cex=agg.sex$Prop.male*EXP)
    #pie chart
    for(j in 1:nrow(agg.sex))
    {
      floating.pie(agg.sex$LONG[j],agg.sex$LAT[j]-1,c(agg.sex$N.male[j]+.1,agg.sex$N.female[j]+.1),
                   radius=0.45,col=c("grey95","black"))
    }
    axis(2,seq(YAXIS[1],YAXIS[2],2),F)
    axis(1,seq(XAXIS[1],XAXIS[2],2),F)
    box() 
  }
  
  
  XAXIS=c(110,134.5)
  YAXIS=c(-37,-13)
  
  tiff("Figure 1.tiff",width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
  par(mfcol=c(2,2),mai=c(0.6,0.5,0.1,0),oma=c(0.1,0.1,0.1,0.1))
  for (i in 1:length(ID))
  {
    fn.show.block(subset(Data.GLM[[ID[i]]],SPECIES==names(SEX.ratio.GN)[ID[i]] & SEX%in%c("M","F") & Month==Mnth[i]))
    #fn.show.block(SEX.ratio.GN[[ID[i]]]$DaT[[Mnth[i]]])
    
    legend(108,-12,Comm.Species.labels[i],bty='n',cex=1.65)
    legend(108.5,-14,paste("(",MNTH[i],")",sep=""),bty='n',cex=1.3)
    if(i%in%c(1:2))  axis(2,seq(YAXIS[1],YAXIS[2],4),-seq(YAXIS[1],YAXIS[2],4),cex.axis=1.25,las=2)
    if(i%in%c(2,4)) axis(1,seq(XAXIS[1],XAXIS[2],4),seq(XAXIS[1],XAXIS[2],4),cex.axis=1.25)
    
    if(i==1)
    {
      legend(114,-22,"Western",bty='n',cex=1.75)
      legend(114,-25,"Australia",bty='n',cex=1.75)
    }
    
  }
  
  mtext("Longitude (?E)",side=1,line=-1.125,font=1,las=0,cex=1.7,outer=T)
  mtext("Latitude (?S)",side=2,line=-1.5,font=1,las=0,cex=1.7,outer=T)
  dev.off()
  
  
  
  
  #1.3 Test Month effect
  
  #Select appropriate data for comparisons
  fn.subset.Mn=function(dat)
  {
    Bk.yr=table(dat$BLOCK,dat$year)
    Bk.yr=ifelse(Bk.yr>1,1,0)
    Sum.Yr=colSums(Bk.yr)
    ID=which(Sum.Yr>=5)  #drop years with less than 5 occurrences
    dat=subset(dat,year%in%names(ID))
    if(nrow(dat)>1)
    {
      Bk.Mn=table(dat$BLOCK,dat$Month)
      Bk.Mn=ifelse(Bk.Mn>1,1,0)
      Sum.Blk=rowSums(Bk.Mn)
      ID=names(which(Sum.Blk==max(Sum.Blk)))
      dat=subset(dat,BLOCK==ID[1])
    }
    
    return(dat)
  }
  
  Data.GLM.Mn=vector('list',length(This.sp))
  names(Data.GLM.Mn)=This.sp
  for(i in 1:length(Data.GLM.Mn))Data.GLM.Mn[[i]]=fn.subset.Mn(subset(DATA,SPECIES%in%This.sp[i]& SEX%in%c("M","F")))
  
  
  Sex.Mn.fn=function(dat)
  {
    Mn=sort(unique(dat$BLOCK))
    ANOV=data.frame(BLOCK=Mn,p=NA,N.months=NA)
    for (x in 1:length(Mn))
    {
      dat1=subset(dat,!(is.na(BLOCK) | BLOCK==0 ) & BLOCK==Mn[x])
      Table=table(dat1$Month)
      ID=which(Table>=5)
      dat1=subset(dat1,Month%in%names(Table[ID]))
      dat1$Month=as.factor(dat1$Month)
      
      agg.Male=aggregate(Number~Month,subset(dat1,SEX=="M"),sum)
      names(agg.Male)[2]="N.male"
      agg.Fem=aggregate(Number~Month,subset(dat1,SEX=="F"),sum)
      names(agg.Fem)[2]="N.female"
      agg.sex=merge(agg.Male,agg.Fem,by=c("Month"),all.x=T,all.y=T)
      agg.sex$N.male=with(agg.sex,ifelse(is.na(N.male),0,N.male))
      agg.sex$N.female=with(agg.sex,ifelse(is.na(N.female),0,N.female))
      
      ANOV$N.months[x]=nrow(agg.sex)
      #glm
      if(ncol(agg.sex)>1)
      {
        model<- glm(cbind(N.male,N.female)~Month, data=agg.sex, family="binomial", maxit=500)
        
        #Anova table
        ANOV$p[x]=round(anova(model,test="Chisq")$"Pr(>Chi)"[2],3)
        
      }
      
    }
    return(ANOV)
  }
  
  SEX.ratio.GN=vector('list',length(This.sp))
  names(SEX.ratio.GN)=This.sp
  for(i in 1:length(This.sp))SEX.ratio.GN[[i]]=Sex.Mn.fn(Data.GLM.Mn[[i]])
  
  #get p values for each species
  Store.p=do.call(rbind,SEX.ratio.GN)
  Store.p$Species=rownames(Store.p)
  write.csv(Store.p,"Sex.ratio.Month.p.values.csv",row.names=F)
  
  
  #Show term effects 
  fn.show.block=function(dat)
  {
    
    #Observations
    agg.Male=aggregate(Number~Month,subset(dat,SEX=="M"),sum)
    names(agg.Male)[2]="N.male"
    agg.Fem=aggregate(Number~Month,subset(dat,SEX=="F"),sum)
    names(agg.Fem)[2]="N.female"
    agg.sex=merge(agg.Male,agg.Fem,by=c("Month"),all.x=T,all.y=T)
    agg.sex$N.male=with(agg.sex,ifelse(is.na(N.male),0,N.male))
    agg.sex$N.female=with(agg.sex,ifelse(is.na(N.female),0,N.female))
    agg.sex$Prop.male=(agg.sex$N.male)/(agg.sex$N.male+agg.sex$N.female)
    
    
    ID=which(!Months%in%agg.sex$Month)
    if(length(ID)>0)
    {
      Add=data.frame(Month=Months[ID],N.male=NA, N.female=NA, Prop.male=NA)
      agg.sex=rbind(agg.sex,Add)
      agg.sex=agg.sex[order(agg.sex$Month),]
    }
    return(agg.sex)
  }
  
  Store.Mn.out=Data.GLM.Mn
  for(i in 1:length(Data.GLM.Mn)) Store.Mn.out[[i]]=fn.show.block(Data.GLM.Mn[[i]])
  
  #COL=grey.colors(length(Data.GLM.Mn), start = 0.3, end = 0.9, gamma = 2.2, alpha = NULL)
  COL=c("black","grey60","black","grey60")
  LTY=c(1,2,2,1)
  
  tiff("Figure 2.tiff",width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
  par(mfcol=c(1,1),mai=c(0.8,0.8,0.1,0),oma=c(0.1,0.1,0.1,0.1),mgp=c(2.5,.6,0))
  plot(Months,ylab="Proportion male",xlab="Month",cex.lab=1.75,cex.axis=1.5,ylim=c(0,1),col="transparent",las=1)
  for(i in 10:length(Data.GLM.Mn)) lines(Store.Mn.out[[i]]$Month,Store.Mn.out[[i]]$Prop.male,col=COL[i-9],lty=LTY[i-9],lwd=2)
  legend("topright",c("Whiskery shark","Gummy shark","Sandbar shark","Dusky shark"),bty='n',lty=LTY,col=COL,lwd=2,cex=1.5)
  
  dev.off()
  
  
  
  #1.4 Test Depth effect
  MaxDepth=250
  Sex.fn.Depth=function(dat)
  {
    dat=subset(dat,!(is.na(BLOCK) | BLOCK==0 ))
    Table.Mn.Blk=table(dat$BLOCK,dat$Month)
    Mn=sort(unique(dat$Month))
    ANOV=data.frame(Month=Mn,p=NA,coef=NA)
    for (x in 1:length(Mn))
    {
      dat1=subset(dat,Month==Mn[x])
      Table=table(dat1$BOTDEPTH,dat1$BLOCK)
      Table=ifelse(Table>=1,1,0)
      ID=colSums(Table)
      ID=names(which(ID==max(ID)))
      dat1=subset(dat1,BLOCK==ID[1])
      
      if(nrow(subset(dat1,SEX=="M"))>0)
      {
        agg.Male=aggregate(Number~BOTDEPTH,subset(dat1,SEX=="M"),sum)
        names(agg.Male)[2]="N.male"
        agg.Fem=aggregate(Number~BOTDEPTH,subset(dat1,SEX=="F"),sum)
        names(agg.Fem)[2]="N.female"
        agg.sex=merge(agg.Male,agg.Fem,by=c("BOTDEPTH"),all.x=T,all.y=T)
        agg.sex$N.male=with(agg.sex,ifelse(is.na(N.male),0,N.male))
        agg.sex$N.female=with(agg.sex,ifelse(is.na(N.female),0,N.female))
        agg.sex$Tot=agg.sex$N.male+agg.sex$N.female
        agg.sex=subset(agg.sex,Tot>=4)
        
        #glm
        if(nrow(agg.sex)>1)
        {
          model<- glm(cbind(N.male,N.female)~BOTDEPTH, data=agg.sex, family="binomial", maxit=500)
          
          #Anova table
          ANOV$p[x]=round(anova(model,test="Chisq")$"Pr(>Chi)"[2],3)
          ANOV$coef[x]=coef(model)[2]
        }
        
      }
      
    }
    return(ANOV)
  }
  
  SEX.ratio.GN=vector('list',length(This.sp))
  names(SEX.ratio.GN)=This.sp
  for(i in 1:length(This.sp))SEX.ratio.GN[[i]]=Sex.fn.Depth(Data.GLM[[i]])
  
  #get p values for each species
  see.p=function(dat)
  {
    dat$Species=names(SEX.ratio.GN)[i]
    
    ID=which(!Months%in%dat$Month)
    if(length(ID)>0)
    {
      Add=data.frame(Month=Months[ID],p=NA, coef=NA, Species=names(SEX.ratio.GN)[i])
      dat=rbind(dat,Add)
      dat=dat[order(dat$Month),]
    }
    return(dat)
  }
  Store.p=SEX.ratio.GN
  for(i in 1:length(SEX.ratio.GN))Store.p[[i]]=see.p(SEX.ratio.GN[[i]])
  Store.p=do.call(rbind,Store.p)
  Store.p.wide.p=reshape(Store.p[,c(1,2,4)],idvar ="Species",v.names=c("p"),timevar="Month",direction="wide")
  Store.p.wide.coef=reshape(Store.p[,c(1,3,4)],idvar ="Species",v.names=c("coef"),timevar="Month",direction="wide")
  
  Store.p.wide.p=merge(Store.p.wide.p,Tabl1.matrix[,c(3,9)],by.x="Species",by.y="Code",all.x=T)
  Store.p.wide.p=Store.p.wide.p[order(-Store.p.wide.p$N),]
  
  write.csv(Store.p.wide.p,"Sex.ratio.p.values.Depth.csv",row.names=F)
  
  
  #Show term effects 
  PCH=21
  COL1=c("black","grey80","grey45","white")
  fn.show.block=function(dat)
  {
    #Observations
    dat=subset(dat,!(is.na(BLOCK) | BLOCK==0 ))
    
    SP=names(table(dat$SPECIES))
    Mn=sort(unique(dat$Month))
    
    par(mfcol=c(3,4))
    for (x in 1:length(Mn))
    {
      dat.M=subset(dat,Month==Mn[x])
      
      SP1=names(table(dat.M$SPECIES)) 
      plot(1:max(dat.M$BOTDEPTH,na.rm=T),ylim=c(0,1),col='transparent',las=1,ylab='',xlab="",main=Mn[x])
      for(n in 1:length(SP1))
      {
        dat1=subset(dat.M,SPECIES==SP1[n])
        Table=table(dat1$BOTDEPTH,dat1$BLOCK)
        Table=ifelse(Table>=1,1,0)
        ID=colSums(Table)
        ID=names(which(ID==max(ID)))
        dat1=subset(dat1,BLOCK==ID[1])
        
        if(nrow(subset(dat1,SEX=="M"))>0)
        {
          agg.Male=aggregate(Number~BOTDEPTH,subset(dat1,SEX=="M"),sum)
          names(agg.Male)[2]="N.male"
          agg.Fem=aggregate(Number~BOTDEPTH,subset(dat1,SEX=="F"),sum)
          names(agg.Fem)[2]="N.female"
          agg.sex=merge(agg.Male,agg.Fem,by=c("BOTDEPTH"),all.x=T,all.y=T)
          agg.sex$N.male=with(agg.sex,ifelse(is.na(N.male),0,N.male))
          agg.sex$N.female=with(agg.sex,ifelse(is.na(N.female),0,N.female))
          agg.sex$Tot=agg.sex$N.male+agg.sex$N.female
          agg.sex=subset(agg.sex,Tot>=4)
          agg.sex$prop=agg.sex$N.male/agg.sex$N.female
          points(agg.sex$BOTDEPTH,agg.sex$prop,bg=COL1[n],pch=PCH,cex=2)
        }
        
      } 
    }
    plot(1:max(dat.M$BOTDEPTH,na.rm=T),ylim=c(0,1),col='transparent',xaxt='n',yaxt='n',ylab='',xlab="")
    legend("center",SP,pch=PCH,col=1,pt.bg=COL1,bty='n',cex=2)
  }
  
  fn.show.block(rbind(Data.GLM[["WH"]],Data.GLM[["GM"]],Data.GLM[["TK"]],Data.GLM[["BW"]]))
  
  
  
  
  ####
  
  
  
  
  Data.GLM=vector('list',length(This.sp))
  names(Data.GLM)=This.sp
  for(i in 1:length(Data.GLM))Data.GLM[[i]]=fn.subset(subset(DATA,SPECIES%in%This.sp[i]))
  
  
  #2.1 Density distributions
  Dens.dist=function(datas,Gear)
  {
    datas=subset(datas,Method==Gear)
    name=unique(datas$COMMON_NAME)
    d <- density(datas$FL,adjust = 2)
    plot(d, type="n", main=name)
    polygon(d, col="red", border="gray")
  }
  
  Dens.dist(subset(DATA,SPECIES=="BW"),"GN")
  
  
  #2.2
  #2.2 GLM By GN (6 to 8 inch) and LL separately
  Size.fn=function(dat)
  {
    dat=subset(dat,!(is.na(BLOCK) | BLOCK==0 | is.na(FL) | FL<=20 | FL>550))
    dat$year=as.factor(dat$year)
    dat$Month=as.factor(dat$Month)
    dat$BLOCK=as.factor(dat$BLOCK)
    
    model<- glm(log(FL)~BLOCK+year+Month+SOI+Freo, data=dat, family=gaussian, maxit=500)
    
    #Anova table
    Signifcance1=anova(model,test="Chisq")
    
    #Deviance explained
    Dev.exp=Dsquared(model,adjust=F)$d3
    
    return(list(model=model,ANOVA=Signifcance1,Dev.exp=Dev.exp))
  }
  Size.fn(subset(DATA,SPECIES%in%Commercial.Sks[i]))
  
}