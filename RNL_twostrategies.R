# For BARCOO

# Before running:
# 1. Check if R and Netlogo versions work on comp (initial Mac OS and R config).
# 2. Change working dirs (~lines 44, 49-52, 101, 105, 356, 377-395)
# 3. This model runs DEB.R and onelump_varenv.R from source.

# Model run for breeding season (117 days)
# sc = simulation count. Change this number to run N number of simulations (100). 
# Run 100 sims for Optimising and 100 sims for Satisficing strategy 
# There should be nine output files for each (100) model run. 


#Step 1 run model sims

args <- (commandArgs(TRUE))
sc<-as.numeric(args[1])+1
strat<-as.numeric(args[2])
resource<-as.numeric(args[3])


# the rest of the code is the DEB, microclimate, and NL setup and the sim model  
setwd("/vlsci/VR0212/mrke/achaves/abm")
source('DEB.R')
source('onelump_varenv.R')

# read in microclimate data
tzone<-paste("Etc/GMT+",10,sep="")
metout<-read.csv('metout.csv')
soil<-read.csv('soil.csv')
shadmet<-read.csv('shadmet.csv')
shadsoil<-read.csv('shadsoil.csv')
micro_sun_all<-cbind(metout[,2:5],metout[,9],soil[,6],metout[,14:16])
colnames(micro_sun_all)<-c('dates','JULDAY','TIME','TALOC','VLOC','TS','ZEN','SOLR','TSKYC')
micro_shd_all<-cbind(shadmet[,2:5],shadmet[,9],shadsoil[,6],shadmet[,14:16])
colnames(micro_shd_all)<-c('dates','JULDAY','TIME','TALOC','VLOC','TS','ZEN','SOLR','TSKYC')



# choose a day(s) to simulate
daystart<-paste('09/09/05',sep="") # yy/mm/dd
dayfin<-paste('09/12/31',sep="") # yy/mm/dd
micro_sun<-subset(micro_sun_all, format(as.POSIXlt(micro_sun_all$dates), "%y/%m/%d")>=daystart & format(as.POSIXlt(micro_sun_all$dates), "%y/%m/%d")<=dayfin)
micro_shd<-subset(micro_shd_all, format(as.POSIXlt(micro_shd_all$dates), "%y/%m/%d")>=daystart & format(as.POSIXlt(micro_shd_all$dates), "%y/%m/%d")<=dayfin)
days<-as.numeric(as.POSIXlt(dayfin)-as.POSIXlt(daystart))
if(days==0){days=1}
# create time vectors
time<-seq(0,(days+1)*60*24,60) #60 minute intervals from microclimate output
time<-time[-1]
times2<-seq(0,(days+1)*60*24,2) #two minute intervals for prediction
time<-time*60 # minutes to seconds
times2<-times2*60 # minutes to seconds

# apply interpolation functions
velfun<- approxfun(time, micro_sun[,5], rule = 2)
Zenfun<- approxfun(time, micro_sun[,7], rule = 2)
Qsolfun_sun<- approxfun(time, micro_sun[,8], rule = 2)
Tradfun_sun<- approxfun(time, rowMeans(cbind(micro_sun[,6],micro_sun[,9])), rule = 2)
Tairfun_sun<- approxfun(time, micro_sun[,4], rule = 2)
Qsolfun_shd<- approxfun(time, micro_shd[,8]*.1, rule = 2)
Tradfun_shd<- approxfun(time, rowMeans(cbind(micro_shd[,6],micro_shd[,9])), rule = 2)
Tairfun_shd<- approxfun(time, micro_shd[,4], rule = 2)



VTMIN<- 26 # from Bundey field site (2-9-14)
VTMAX<- 35 # from max(results$Tb)

# ****************************************** read in DEB parameters *******************************

debpars=as.data.frame(read.csv('DEB_pars_Tiliqua_rugosa.csv',header=FALSE))$V1 # read in DEB pars

# set core parameters
z=debpars[8] # zoom factor (cm)
F_m = 13290 # max spec searching rate (l/h.cm^2)
kap_X=debpars[11] # digestion efficiency of food to reserve (-)
v=debpars[13] # energy conductance (cm/h)
kap=debpars[14] # kappa, fraction of mobilised reserve to growth/maintenance (-)
kap_R=debpars[15] # reproduction efficiency (-)
p_M=debpars[16] # specific somatic maintenance (J/cm3)
k_J=debpars[18] # maturity maint rate coefficient (1/h)
E_G=debpars[19] # specific cost for growth (J/cm3)
E_Hb=debpars[20] # maturity at birth (J)
E_Hp=debpars[21] # maturity at puberty (J)
h_a=debpars[22]*10^-1 # Weibull aging acceleration (1/h^2)
s_G=debpars[23] # Gompertz stress coefficient (-)

# set thermal respose curve paramters
T_REF = debpars[1]-273
TA = debpars[2] # Arrhenius temperature (K)
TAL = debpars[5] # low Arrhenius temperature (K)
TAH = debpars[6] # high Arrhenius temperature (K)
TL = debpars[3] # low temp boundary (K)
TH = debpars[4] # hight temp boundary (K)

# set auxiliary parameters
del_M=debpars[9] # shape coefficient (-) 
E_0=debpars[24] # energy of an egg (J)
mh = 1 # survivorship of hatchling in first year
mu_E = 585000 # molar Gibbs energy (chemical potential) of reserve (J/mol)
E_sm=186.03*6
gutfull <- 1
# set initial state
E_pres_init = (debpars[16]*debpars[8]/debpars[14])/(debpars[13]) # initial reserve
E_m <- E_pres_init
E_H_init = debpars[21] + 5

#### change inital size here by multiplying by < 0.85 ####
V_pres_init = (debpars[26] ^ 3) * 0.85 
d_V<-0.3
mass <- V_pres_init + V_pres_init*E_pres_init/mu_E/d_V*23.9

# ***************** end TRANSIENT MODEL SETUP ***************
#************************************************************

# initial setup complete, now run the model for a particular 2 min interval 
# **************************************************************************************************************
# ***************************************** start NETLOGO SIMULATION  ******************************************

# ****************************************** open NETLOGO ***************************************** 
Sys.setenv(NOAWT=1)
library(RNetLogo);library(adehabitatHR); library(sp)



nl.path<-"/vlsci/VR0212/mrke/achaves/abm/netlogo-5.1.0"
gui<-FALSE


NLStart(nl.path, gui=gui)
#model.path<-""/Applications/Programs/NetLogo 5.0.5/Soft foraging model/Soft foraging submodel of Tiliqua rugosa_v.5.nlogo"
model.path<-"/vlsci/VR0212/mrke/achaves/abm/Sleepy IBM_v.6.1.1_two strategies.nlogo"
NLLoadModel(model.path) 

if(strat == 0){
 NLCommand("set strategy \"Optimising\" ") 
}else{
 NLCommand("set strategy \"Satisficing\" ") 
}
if(resource == 0){
 NL_shade<-1000L         # Shade patches
 NL_food<-as.numeric(100)      # Food patches (*10, so 10x less than shade)
}else{
 NL_shade<-100000L         # Shade patches
 NL_food<-as.numeric(10000)      # Food patches (*10, so 10x less than shade)
}
# ****************************************** setup NETLOGO MODEL **********************************

# 1. update animal and env traits
month<-"sep"
#mass<-800             # Weight in grams
NL_shade<-100000L          # Shade patches
NL_food<-as.numeric(10000)      # Food patches (*10)
NL_days<-1#17         # No. of days simulated
NL_gutthresh<-0.75

# 2. update initial conditions for DEB model 
Es_pres_init<-(E_sm*gutfull)*V_pres_init
acthr<-1
Tb_init<-20
step = 1/24
debout<-DEB(step = step, z = z, del_M = del_M, F_m = F_m * 
    step, kap_X = kap_X, v = v * step, kap = kap, p_M = p_M * 
    step, E_G = E_G, kap_R = kap_R, k_J = k_J * step, E_Hb = E_Hb, 
    E_Hj = E_Hb, E_Hp = E_Hp, h_a = h_a/(step^2), s_G = s_G, 
    T_REF = T_REF, TA = TA, TAL = TAL, TAH = TAH, TL = TL, 
    TH = TH, E_0 = E_0, E_pres=E_pres_init, V_pres=V_pres_init, E_H_pres=E_H_init, acthr = acthr, breeding = 1, Es_pres = Es_pres_init, E_sm = E_sm)


# 3. calc direct movement cost
# ------- New loco cost 16-10-16 -------
#V_pres = 3.9752^3 #structure 
V_pres<-debout[2]
step<-1/24 #hourly
#step<-2/1440 #2min

p_M<-32*step #J/h
p_M<-p_M*V_pres # loco cost * structure
names(p_M)<-NULL # remove V_pres name attribute from p_M

# movement cost for time period
VO2<-0.45 # O2/g/h JohnAdler etal 1986

# multiple p_M by structure = movement cost (diff between p_M with loco cost and structure for movement period)
# p_M with loco cost 
loco<-VO2*mass*20.1 # convert ml O2 to J = J/h 
loco<-loco+p_M # add to p_M = J/h
loco<-loco/30 ; loco #J/2min


Es_pres_init<-(E_sm*gutfull)*V_pres_init
X_food<-3000
V_pres<-debout[2]
wetgonad<-debout[19]
wetstorage<-debout[20]
wetfood<-debout[21]
ctminthresh<-120000
Tairfun<-Tairfun_shd
Tc_init<-Tairfun(1)+0.1 # Initial core temperature

NL_T_b<-Tc_init       # Initial T_b
NL_T_b_min<-VTMIN         # Min foraging T_b
NL_T_b_max<-VTMAX        # Max foraging T_b
NL_ctminthresh<-ctminthresh # No. of consecutive hours below CTmin that leads to death
NL_reserve<-E_m        # Initial reserve density
NL_max_reserve<-E_m    # Maximum reserve level
NL_maint<-round(p_M, 3)               # Maintenance cost
NL_move<-round(loco, 3) 		      # Movement cost
NL_zen<-Zenfun(1*60*60)     # Zenith angle

#sc<-100 # sim count for automating writing of each sim results to file (set before NL loop)
for (i in sc:sc){ # start sc sim loop

NLCommand("set Shade-patches",NL_shade,"set Food-patches",NL_food,"set No.-of-days",NL_days,"set T_b precision",
NL_T_b, "2","set T_opt_lower precision", NL_T_b_min, "2","set T_opt_upper precision", NL_T_b_max, "2",
"setup", "set reserve-level", NL_reserve, "set Maximum-reserve", NL_max_reserve, "set Maintenance-cost", NL_maint,
"set Movement-cost precision", NL_move, "3", "set zenith", NL_zen, "set ctminthresh", NL_ctminthresh, 
"set gutthresh", NL_gutthresh, 'set gutfull', gutfull, 'set V_pres precision', V_pres, "5", 'set wetstorage precision', wetstorage, "5", 
'set wetfood precision', wetfood, "5", 'set wetgonad precision', wetgonad, "5")

#NLCommand("inspect turtle 0")

NL_ticks<-NL_days / (2 / 60 / 24) * 1 # No. of NL ticks (measurement of days)
NL_T_opt_l<-NLReport("[T_opt_lower] of turtle 0")
NL_T_opt_u<-NLReport("[T_opt_upper] of turtle 0")

# data frame setup for homerange polygon
turtles<-data.frame() # make an empty data frame
NLReport("[X] of turtle 0"); NLReport("[Y] of turtle 0")
who<-NLReport("[who] of turtle 0")

# **********************************************************
# ******************** start NETLOGO SIMULATION  ***********
debcall<-0 # check for first call to DEB
stepcount<-0 # DEB model step count

for (i in 1:NL_ticks){
stepcount<-stepcount+1
NLDoCommand(1, "go")

######### Reporting presence of shade
#if (NLReport("any? turtles")){
shade<-NLGetAgentSet("in-shade?","turtles", as.data.frame=T); shade<-as.numeric(shade) # returns an agentset of whether turtle is currently on shade patch
#food<-NLGetAgentSet("in-food?","turtles", as.data.frame=T) ; food<-as.numeric(food)
#}

# choose sun or shade
tick<-i
times3<-c(times2[tick],times2[tick+1])

if(shade==0){
  Qsolfun<-Qsolfun_sun
  Tradfun<-Tradfun_sun
  Tairfun<-Tairfun_sun
}else{
  Qsolfun<-Qsolfun_shd
  Tradfun<-Tradfun_shd
  Tairfun<-Tairfun_shd
}
if(i==1){
Tc_init<-Tairfun(1)+0.1 #initial core temperature
#Zenf<-Zenf(1*60*60)  # initial zenith angle at 1 a.m (?)
}

# ----------------------------------- 2-12-14 new one_lump_trans params
Qsol<-Qsolfun(mean(times3)); Qsol
vel<-velfun(mean(times3)) ;vel
Tair<-Tairfun(mean(times3));Tair
Trad<-Tradfun(mean(times3)); Trad
Zen<-Zenfun(mean(times3)); Zen

# calc Tb params at 2 mins interval
Tbs<-onelump_varenv(t=120,time=times3[2],Tc_init=Tc_init,thresh = 30, AMASS = mass, lometry = 3, Tairf=Tairfun,Tradf=Tradfun,velf=velfun,Qsolf=Qsolfun,Zenf=Zenfun)
Tb<-Tbs$Tc
rate<-Tbs$dTc
Tc_init<-Tb

NLCommand("set T_b precision", Tb, "2") # Updating Tb
NLCommand("set zenith", Zenfun(i*60*60)) # Updating zenith

# time spent below VTMIN
ctminhours<-NLReport("[ctmincount] of turtle 0") * 2/60 # ticks to hours
if (ctminhours == NL_ctminthresh) {NLCommand("ask turtle 0 [stop]")}

# ******************** start DEB SIMULATION  ***************

if(stepcount==1) { # run DEB loop every time step (2 mins)
stepcount<-0

# report activity state
actstate<-NLReport("[activity-state] of turtle 0")
 # Reports true if turtle is in food 
actfeed<-NLGetAgentSet("in-food?","turtles", as.data.frame=T); actfeed<-as.numeric(actfeed)
 
n<-1 # time steps
step<-2/1440 # step size (2 mins). For hourly: 1/24
# update direct movement cost
if(actstate == "S"){
	NLCommand("set Movement-cost", NL_move)
	}else{
		NLCommand("set Movement-cost", 1e-09)
		} 
# if within activity range, it's daytime, and gut below threshold 
if(Tbs$Tc>=VTMIN & Tbs$Tc<=VTMAX & Zen!=90 & gutfull<=NL_gutthresh){ 
  acthr=1 # activity state = 1 
if(actfeed==1){ # if in food patch
	X_food<-NLReport("[energy-gain] of turtle 0") # report joules intake
	}
	}else{
		X_food = 0 
		acthr=0
		}

# calculate DEB output 
if(debcall==0){
	# initialise DEB
	debout<-matrix(data = 0, nrow = n, ncol = 26)
	deb.names<-c("E_pres","V_pres","E_H_pres","q_pres","hs_pres","surviv_pres","Es_pres","cumrepro","cumbatch","p_B_past","O2FLUX","CO2FLUX","MLO2","GH2OMET","DEBQMET","DRYFOOD","FAECES","NWASTE","wetgonad","wetstorage","wetfood","wetmass","gutfreemass","gutfull","fecundity","clutches")
	colnames(debout)<-deb.names
	# initial conditions
	debout<-DEB(E_pres=E_pres_init, V_pres=V_pres_init, E_H_pres=E_H_init, acthr = acthr, Tb = Tb_init, breeding = 1, Es_pres = Es_pres_init, E_sm = E_sm, step = step, z, del_M = del_M, F_m = F_m * 
    step, kap_X = kap_X, v = v * step, kap = kap, p_M = p_M * 
    step, E_G = E_G, kap_R = kap_R, k_J = k_J * step, E_Hb = E_Hb, 
    E_Hj = E_Hb, E_Hp = E_Hp, h_a = h_a/(step^2), s_G = s_G, 
    T_REF = T_REF, TA = TA, TAL = TAL, TAH = TAH, TL = TL, 
    TH = TH, E_0 = E_0)
	debcall<-1
	}else{
		debout<-DEB(step = step, z = z, del_M = del_M, F_m = F_m * 
    step, kap_X = kap_X, v = v * step, kap = kap, p_M = p_M * 
    step, E_G = E_G, kap_R = kap_R, k_J = k_J * step, E_Hb = E_Hb, 
    E_Hj = E_Hb, E_Hp = E_Hp, h_a = h_a/(step^2), s_G = s_G, 
    T_REF = T_REF, TA = TA, TAL = TAL, TAH = TAH, TL = TL, 
    TH = TH, E_0 = E_0, 
		  X=X_food,acthr = acthr, Tb = Tbs$Tc, breeding = 1, E_sm = E_sm, E_pres=debout[1],V_pres=debout[2],E_H_pres=debout[3],q_pres=debout[4],hs_pres=debout[5],surviv_pres=debout[6],Es_pres=debout[7],cumrepro=debout[8],cumbatch=debout[9],p_B_past=debout[10])
		}
mass<-debout[22]
gutfull<-debout[24]
NL_reserve<-debout[1]
V_pres<-debout[2]
wetgonad<-debout[19]
wetstorage<-debout[20]
wetfood<-debout[21] 

#update NL wetmass properties 
NLCommand("set V_pres precision", V_pres, "5")
NLDoCommand("plot xcor ycor")
NLCommand("set wetgonad precision", wetgonad, "5")
NLDoCommand("plot xcor ycor")
NLCommand("set wetstorage precision", wetstorage, "5")
NLDoCommand("plot xcor ycor")
NLCommand("set wetfood precision", wetfood, "5")
NLDoCommand("plot xcor ycor") 
 

} #--- end DEB loop

NLCommand("set reserve-level", NL_reserve) # update reserve
NLCommand("set gutfull", debout[24])# update gut level

# ******************** end DEB SIMULATION ******************

# generate results, with V_pres, wetgonad, wetstorage, and wetfood from debout
if(i==1){
	results<-cbind(tick,Tb,rate,shade,V_pres,wetgonad,wetstorage,wetfood,NL_reserve) 
	}else{
		results<-rbind(results,c(tick,Tb,rate,shade,V_pres,wetgonad,wetstorage,wetfood,NL_reserve))
		}
results<-as.data.frame(results)

# generate data frames for homerange polygon
if (tick == NL_ticks - 1){
	X<-NLReport("[X] of turtle 0"); head(X)
	Y<-NLReport("[Y] of turtle 0"); head(Y)
	turtles<-data.frame(X,Y)
	who1<-rep(who,NL_ticks); who # who1<-rep(who,NL_ticks - 1); who 
	turtledays<-rep(1:NL_days,length.out=NL_ticks,each=720) 
	turtle<-data.frame(ID = who1,days=turtledays)
	turtles<-cbind(turtles,turtle)
	}

} # *************** end NL loop **************************

  
# get hr data
spdf<-SpatialPointsDataFrame(turtles[1:2], turtles[3]) # creates a spatial points data frame (adehabitatHR package)
homerange<-mcp(spdf,percent=95)
  
  # writing new results
  if (exists("results")){  #if results exist
    sc<-sc-1 
    nam <- paste("results", sc, sep = "") # generate new name with added sc count
    rass<-assign(nam,results) #assign new name to results. call 'results1, results2 ... resultsN'
    namh <- paste("turtles", sc, sep = "")  #generate new name with added sc count
    rassh<-assign(namh,turtles) #assign new name to results. call 'results1, results2 ... resultsN'
    nams <- paste("spdf", sc, sep = "") 
    rasss<-assign(nams,spdf) 
    namhr <- paste("homerange", sc, sep = "")  
    rasshr<-assign(namhr,homerange) 
    
# write each result for each sim to file dir getwd()
fh<-"/vlsci/VR0212/mrke/achaves/abm/results/"
	for (i in rass){
		# export all results
		write.table(results,file=paste(fh,nam,".R",sep=""))
		}
	for (i in rassh){
		# export turtle location data
		write.table(turtles,file=paste(fh,namh,".R",sep=""))
		}
#	for (i in rasss){
		# export spdf data
#		write.table(spdf,file=paste(fh,nams,".R",sep=""))  
#		}
#	for (i in rasshr){
		# export hr data
#		write.table(homerange,file=paste(fh,namhr,".R",sep=""))  
#		}#export NL plots
month<-"sep"
#spatial plot
sfh<-paste(month,NL_days,round(mass,0),NL_shade,as.integer(NL_food*10),"_",sc,"_move","",sep="");sfh
NLCommand(paste("export-plot \"Spatial coordinates of transition between activity states\" \"/vlsci/VR0212/mrke/achaves/abm/results/",sfh,".csv\"",sep=""))
#temp plot 
tfh<-paste(month,NL_days,round(mass,0),NL_shade,as.integer(NL_food*10),"_",sc,"_temp",sep="")
NLCommand(paste("export-plot \"Body temperature (T_b)\" \"/vlsci/VR0212/mrke/achaves/abm/results/",tfh,".csv\"",sep=""))
#activity budget
afh<-paste(month,NL_days,round(mass,0),NL_shade,as.integer(NL_food*10),"_",sc,"_act","",sep="");afh
NLCommand(paste("export-plot \"Global time budget\" \"/vlsci/VR0212/mrke/achaves/abm/results/",afh,".csv\"",sep=""))
#text output
xfh<-paste(month,NL_days,round(mass,0),NL_shade,as.integer(NL_food*10),"_",sc,"_txt",sep="");xfh
NLCommand(paste("export-output \"/vlsci/VR0212/mrke/achaves/abm/results/",xfh,".csv\"",sep=""))
#gut level
gfh<-paste(month,NL_days,round(mass,0),NL_shade,as.integer(NL_food*10),"_",sc,"_gut","",sep="");gfh
NLCommand(paste("export-plot \"Gutfull\" \"/vlsci/VR0212/mrke/achaves/abm/results/",gfh,".csv\"",sep=""))
#wet mass 
mfh<-paste(month,NL_days,round(mass,0),NL_shade,as.integer(NL_food*10),"_",sc,"_wetmass","",sep="");mfh
NLCommand(paste("export-plot \"Total wetmass plot\" \"/vlsci/VR0212/mrke/achaves/abm/results/",mfh,".csv\"",sep=""))
#movement cost (loco) 
lfh<-paste(month,NL_days,round(mass,0),NL_shade,as.integer(NL_food*10),"_",sc,"_loco","",sep="");lfh
NLCommand(paste("export-plot \"Movement costs\" \"/vlsci/VR0212/mrke/achaves/abm/results/",lfh,".csv\"",sep=""))
  }
} # ********************** end sc sim loop **********************
