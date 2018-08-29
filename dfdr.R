#--- fdr rate

wd = "d:/Murilo/exp_db/sugarcane/areao/fdr"
plot_swc = T
export   = F

setwd(wd)

source("fdr_calib_f.R")
#--- Read FDR data
fdr_db = read.csv("Diviner_data_calib_all.csv")

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

fdr_series = fdr_fol

l_tb = unique(fdr_series$tube)
for(tb in l_tb){
  l_dp = unique(fdr_series$depth[fdr_series$tube==tb])
  for(dp in l_dp){
    
    fdr = fdr_series[fdr_series$tube==tb & fdr_series$depth==dp,]
    fdr = fdr[order(fdr$date),]
    
    dfdr = fdr
    dfdr$dorig_meas = 0
    for(i in 2:length(fdr$date)){
      if(dfdr$orig_meas[i-1] > 0 & dfdr$orig_meas[i] > 0){
        dfdr$dorig_meas[i] = dfdr$orig_meas[i-1] - dfdr$orig_meas[i]
      }
    }
    
    if(dp == l_dp[1]){
      dfdr_dp = dfdr
    }else{
      dfdr_dp = rbind(dfdr_dp,dfdr)
    }
  }
  
  if(tb == l_tb[1]){
    dfdr_tb = dfdr_dp
  }else{
    dfdr_tb = rbind(dfdr_tb,dfdr_dp)
  }
}


boxplot(dfdr_tb$dorig_meas[dfdr_tb$dorig_meas<0]~dfdr_tb$depth[dfdr_tb$dorig_meas<0],
        col = "grey")
lines(c(0,0)~c(-1,20), col = "red", lty = 3)






plot(dfdr_tb$dorig_meas~dfdr_tb$date)

plot(dfdr$sf~dfdr$date,
     xaxt = "n",
     type = "b",
     ylim = c(0.7,1.0),
     ylab = "Scaled Frequency",
     xlab = "")
axis(1, dfdr$date, format(dfdr$date, "%b %d"), cex.axis = .7)
