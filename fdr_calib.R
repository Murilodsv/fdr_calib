#--- Sugarcane Water Use Data


wd = "d:/Murilo/exp_db/sugarcane/areao/fdr"
plot_swc = F
export   = F

setwd(wd)

source("fdr_calib_f.R")
#--- Read FDR data
fdr_db = read.csv("fdr_data.csv")

ret_curv = read.csv("ret_curv_all.csv")

#--- order data according to dates
fdr_db$date = as.Date(fdr_db$date, "%m/%d/%Y")
fdr_db = fdr_db[order(fdr_db$date),]

#--- Remove ignored values
fdr_data = fdr_db[fdr_db$ignore==0,]

#--- compute original scaled and absolute frequencies (sf and af)
fdr_data$sf = fdr_data$A * fdr_data$orig_meas ^ fdr_data$B + fdr_data$C
fdr_data$af = fdr_data$fa - fdr_data$sf * (fdr_data$fa - fdr_data$fw)

#--- separate equipments 
fdr_fol = fdr_data[fdr_data$equip_owner == "Folegatti",]
fdr_gep = fdr_data[fdr_data$equip_owner == "GEPEMA",]

feature = "calib_meas"
f_name  = "Scaled Frenquency (0-1)"

if(plot_swc){
  if(export){png(paste0("FDR_",feature,".png"),units="in",width=12,height=12,pointsize=15,res=300)}
  par(mfrow=c(2,2), mar = c(4.5, 4.5, 0.5, 0.5), oma = c(0, 0, 0, 0))
  
  #--- Plot FDR data of Folegatti Sensor
  boxplot(fdr_fol[,feature]~fdr_fol$depth,
          col = "lightblue",
          ylim= c(max(0,min(fdr_data[,feature])),max(fdr_data[,feature])),
          ylab= f_name,
          xlab= "Depth (cm)")

  #--- Plot FDR data of GEPEMA Sensor
  boxplot(fdr_gep[,feature]~fdr_gep$depth,
          col = "lightgreen",
          ylim= c(max(0,min(fdr_data[,feature])),max(fdr_data[,feature])),
          ylab= f_name,
          xlab= "Depth (cm)")

  fdr_fol_bp = bp_df(fdr_fol,T, y_lim = c(max(0,min(fdr_data[,feature])),max(fdr_data[,feature])), feat = feature,c("Depth (cm)",f_name))
  fdr_gep_bp = bp_df(fdr_gep,T, y_lim = c(max(0,min(fdr_data[,feature])),max(fdr_data[,feature])), feat = feature,c("Depth (cm)",f_name))
  
  if(export){dev.off()}
  
  par(mfrow=c(1,1), mar = c(4.5, 4.5, 0.5, 0.5), oma = c(0, 0, 0, 0))
}

#--- retention curve indexing with FDR depth
rc_depth = data.frame(depth_fdr   = seq(10,150,by = 10),
                      depth= c(    5, #10
                                  15, #20
                                  30, #30
                                  30, #40
                                  60, #50
                                  60, #60
                                  60, #70
                                  100,#80
                                  100,#90
                                  100,#100
                                  100,#110
                                  100,#120
                                  100,#130
                                  100,#140
                                  100))#150

#--- indexer
rc_idx = merge(ret_curv,rc_depth,by=c("depth"))

#-----------------------------
#--- GEPEMA Diviner ----------
#-----------------------------

fdr_opt = fdr_gep
fdr_opt$calib_meas_new = fdr_opt$orig_meas

#--- Assuming that the frequency on top layers are equivalent to other layers:
#--- Assuming that the pmp ocurred at least once at 10 cm depth
min_sf = aggregate(sf ~ tube,fdr_opt[fdr_opt$depth==10,],min)
pmp_sf = quantile(min_sf$sf,0.25)[1]

#--- Assuming that the saturation point is the maximun scaled frenquency measured on all series
max_sf = aggregate(sf ~ tube, fdr_opt,max)
sat_sf = quantile(max_sf$sf,0.75)[1]

#--- Original parameters
a = unique(fdr_opt$A)
b = unique(fdr_opt$B)
c = unique(fdr_opt$C)

#--- add paramters to a single vector
p = c(a,b,c)

#--- opmtimze paramters
par_abc = optim(p,opt_fdr)
par_abc$par

#--- Set of optmized parameters
a = par_abc$par[1]
b = par_abc$par[2]
c = par_abc$par[3]

#--- Compute calibrated swc
fdr_opt$calib_meas_new = 10 ^ (log10((fdr_opt$sf-c)/a)/b)

#--- write optimized parameters to original db
fdr_gep$A_calib = a
fdr_gep$B_calib = b
fdr_gep$C_calib = c
fdr_gep$calib_meas = fdr_opt$calib_meas_new


#--------------------------

fdr_opt = fdr_fol
fdr_opt$calib_meas_new = fdr_opt$orig_meas

#--- Assuming that the frequency on top layers are equivalent to other layers:
#--- Assuming that the pmp ocurred at least once at 10 cm depth
min_sf = aggregate(sf ~ tube,fdr_opt[fdr_opt$depth==10,],min)
pmp_sf = quantile(min_sf$sf,0.25)[1]

#--- Assuming that the saturation point is the maximun scaled frenquency measured on all series
max_sf = aggregate(sf ~ tube, fdr_opt,max)
sat_sf = quantile(max_sf$sf,0.75)[1]

#--- Original parameters
a = unique(fdr_opt$A)
b = unique(fdr_opt$B)
c = unique(fdr_opt$C)

#--- add paramters to a single vector
p = c(a,b,c)

#--- opmtimze paramters
par_abc = optim(p,opt_fdr)
par_abc$par

#--- Set of optmized parameters
a = par_abc$par[1]
b = par_abc$par[2]
c = par_abc$par[3]

#--- Compute calibrated swc
fdr_opt$calib_meas_new = 10 ^ (log10((fdr_opt$sf-c)/a)/b)

#--- write optimized parameters to original db
fdr_fol$A_calib = a
fdr_fol$B_calib = b
fdr_fol$C_calib = c
fdr_fol$calib_meas = fdr_opt$calib_meas_new


#--------------------------

fdr_data_calib = rbind(fdr_fol,fdr_gep)





# #--- Compute FDR data for Folegatti equipment for bottom soil (depth >80 cm)
# #--- Assume same a,b c parameters as used for 80 cm
# 
# l_tb = unique(fdr_fol$tube)
# for(tb in l_tb){
#   
#   #--- 80cm parameters
#   a   = unique(fdr_fol$A_calib[fdr_fol$tube==tb & fdr_fol$depth==80])
#   b   = unique(fdr_fol$B_calib[fdr_fol$tube==tb & fdr_fol$depth==80])
#   c   = unique(fdr_fol$C_calib[fdr_fol$tube==tb & fdr_fol$depth==80])
#   sf  = fdr_fol$sf[fdr_fol$tube==tb & fdr_fol$depth > 80] 
#   
#   fdr_fol$calib_meas[fdr_fol$tube==tb & fdr_fol$depth > 80] = 10 ^ (log10((sf-c)/a)/b)
#   
#   #--- assign abc parameters to fol data
#   fdr_fol$A_calib[fdr_fol$tube==tb & fdr_fol$depth > 80] = a
#   fdr_fol$B_calib[fdr_fol$tube==tb & fdr_fol$depth > 80] = b
#   fdr_fol$C_calib[fdr_fol$tube==tb & fdr_fol$depth > 80] = c
#   
# }



l_dp = c(10,30,50,70,90,110,130)
colors = data.frame(depth = l_dp,
                    colors= topo.colors(n = length(l_dp), alpha = 1))

colnames(fdr_data_calib)
fdr_data_calib = fdr_data_calib[fdr_data_calib$calib_meas>0,]

max(fdr_data_calib$calib_meas)

boxplot(fdr_data_calib$calib_meas)
boxplot(fdr_data_calib$sf)

for(dp in l_dp){
  ad = T
  if(dp == l_dp[1]){ad=F}
  ag_plot(fdr_data_calib,"calib_meas",dp, ad,colors$colors[colors$depth==dp],y_li = c(min(fdr_data_calib$calib_meas),max(fdr_data_calib$calib_meas)),"calib_meas")
}

legend("topright",inset = 0.03,legend = l_dp,lt = rep(1,length(colors$colors)),col = as.character(colors$colors))

#--- write optmized values
write.csv(fdr_data_calib,file = "fdr_data_calib_v3.csv", row.names = F)



















feature = "calib_meas"

for(tb in unique(fdr_opt$tube)){
  png(paste0("FDR_calib_tube_",tb,".png"),units="in",width=12,height=12,pointsize=15,res=300)
  colors = data.frame(depth = unique(fdr_opt$depth[fdr_opt$tube==tb]),
                      colors= heat.colors(n = length(unique(fdr_opt$depth[fdr_opt$tube == tb])), alpha = 1))
  
  plot(fdr_opt$calib_meas_new[fdr_opt$tube==tb & fdr_opt$depth==10]~fdr_opt$date[fdr_opt$tube==tb & fdr_opt$depth==10],
       type = "b",
       ylim = c(0.15,0.50),
       ylab = "SWC (cm3 cm-3)",
       xlab = "Date")
  
  for(d in unique(fdr_opt$depth[fdr_opt$tube==tb])){
    
    lines(fdr_opt$calib_meas_new[fdr_opt$tube==tb & fdr_opt$depth==d]~fdr_opt$date[fdr_opt$tube==tb & fdr_opt$depth==d],
          type = "b",
          ylim = c(0.15,0.50),
          col = as.character(colors$colors[colors$depth==d]))
    
  }
  
  legend("bottomright",legend = unique(fdr_opt$depth[fdr_opt$tube==tb]),col = as.character(colors$colors), lt = rep(1,length(colors$colors)))
  dev.off()
}














tb1 = 69
tb2 = 68
d = 50

t1_t2 = merge(fdr_opt[fdr_opt$tube==tb1 & fdr_opt$depth==d,c("date","calib_meas_new")],
      fdr_opt[fdr_opt$tube==tb2 & fdr_opt$depth==d,c("date","calib_meas_new")],
      by = "date")

colnames(t1_t2) = c("date","tb1","tb2")

mperf(sim = t1_t2[,2], obs = t1_t2[,3],vnam = "swc")

for(sd in l_sd){
  for(tb in l_tb){
    if(length(fdr_opt[fdr_opt$tube==tb & fdr_opt$depth==sd,"calib_meas_new"]) > 0){
      
      min_swc = min(fdr_opt[fdr_opt$tube==tb & fdr_opt$depth==sd,"calib_meas_new"])
      max_swc = max(fdr_opt[fdr_opt$tube==tb & fdr_opt$depth==sd,"calib_meas_new"])
    
    }
  }
}






write.csv(fdr_fol_bp, file = "fdr_fol_bp.csv", row.names = F)
write.csv(fdr_gep_bp, file = "fdr_gep_bp.csv", row.names = F)

fdr_sd_ag    = aggregate(freq ~ date+depth, data = fdr_data, sd)
fdr_mean_ag  = aggregate(freq ~ date+depth, data = fdr_data, mean)

dp = 10

fdr_mean  = fdr_mean_ag[fdr_mean_ag$depth==dp,]
fdr_sd    = fdr_sd_ag[fdr_sd_ag$depth==dp,]

ybot = fdr_mean[,3]-fdr_sd[,3]
yupp = fdr_mean[,3]+fdr_sd[,3]

plot(fdr_mean[,3]~fdr_mean[,1],
     ylim = c(min(ybot),max(yupp)))
     
arrows(fdr_mean[,1],
       ybot,
       fdr_mean[,1],
       yupp,
       length=0.02, angle=90, code=3)


  avg = rnorm(10,5,1)
  sdev = rnorm(10,0.1,0.01)
  plot(avg)
  
  arrows(seq(1,10), avg-sdev, seq(1,10), avg+sdev, length=0.05, angle=90, code=3)
  
  plot(fdr_data$calib_meas~fdr_data$date)

mperf(fdr_fol_bp$b1,fdr_gep_bp$b1,vnam = "swc",dchart = T, outidx = "rmse")

hist(fdr_fol$calib_meas[fdr_fol$depth==50])
hist(fdr_gep$calib_meas[fdr_gep$depth==50])








bp = boxplot(fdr_dt$calib_meas~fdr_dt$tube, plot = F)
bp_df$names

for(t in tube_l){
  

  bp_df = data.frame(tube = t,
                     b1   = bp$stats[1,])
  
}

a   = fdr_data$A
b   = fdr_data$B
c   = fdr_data$C
swc = fdr_data$Leitura.Calibrada_Tese_Murilo

sf  = a * swc ^ b + c
swc = 10 ^ (log10((sf-c)/a)/b)

plot(sf~fdr_data$date, type = "l")

unique(fdr_data$A)