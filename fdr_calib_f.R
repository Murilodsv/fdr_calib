

bp_df = function(fdr_dt, p_data, y_lim, feat,lab){
  
  if(missing(feat)){feat="calib_meas"}
  
  bp    = boxplot(fdr_dt[,feat]~fdr_dt$depth, plot = F)
  bp_dt = data.frame(depth = as.numeric(bp$names),
                     b1   = bp$stats[1,],
                     b2   = bp$stats[2,],
                     b3   = bp$stats[3,],
                     b4   = bp$stats[4,],
                     b5   = bp$stats[5,])
  
  
  if(missing(p_data)){p_data=F}
  if(missing(y_lim)){y_lim=c(min(bp_dt$b1),max(bp_dt$b5))}
  
  if(p_data){
    plot(bp_dt$b1~bp_dt$depth, type = "l", 
         ylim = y_lim,
         xlab = lab[1],
         ylab = lab[2],lty = 3)
    lines(bp_dt$b2~bp_dt$depth, type = "l")
    lines(bp_dt$b3~bp_dt$depth, type = "l", col = "red")
    lines(bp_dt$b4~bp_dt$depth, type = "l")
    lines(bp_dt$b5~bp_dt$depth, type = "l", lty = 3)
  }
  
  return(bp_dt)
}

mperf     = function(sim,obs,vnam,dchart,outidx){
  
  #--------------------------------------------------#
  #------------- Performance function ---------------#
  #--- Compute statistical indexes of performance ---#
  #--------------------------------------------------#
  #
  # Decription:
  # sim     - Simulated values          [Real]
  # obs     - Observed values           [Real]
  # vnam    - Name of variable          [String]
  # dchart  - Display Chart?            [T or F]
  # outidx  - Output peformance indexes [List]
  # 
  # Murilo Vianna
  # source: https://github.com/Murilodsv/R-scripts/blob/master/mperf.r
  #
  # Literature: 
  # Brun, F., Wallach, D., Makowski, D., & Jones, J. W. (2006). 
  # Working with dynamic crop models: evaluation, analysis, 
  # parameterization, and applications. Elsevier.
  #--------------------------------------------------#
  
  if(missing(sim)){stop("Missing sim argument")}
  if(missing(obs)){stop("Missing obs argument")}
  if(missing(dchart)){dchart = T}
  if(missing(outidx)){outidx = c("bias","mse","rmse","mae","rrmse","rmae","ef","r","r2","d")}
  if(missing(vnam)){
    warning("Missing vnam argument: vnam set to none.")
    vnam = ""
  }
  
  
  #--- Check Input data
  sim = as.numeric(sim)
  obs = as.numeric(obs)
  
  if(length(sim) != length(obs)){stop("Vector length of Simulated and Observed do not match.")}
  if(length(sim[is.na(sim)]) > 0 |  length(obs[is.na(obs)]) > 0){
    
    #--- remove NA values
    df_so = na.omit(data.frame(sim,obs))
    sim = df_so$sim
    obs = df_so$obs
    warning("NA values ignored in Simulated or Observed data")
  }
  
  #--- Statistical indexes
  fit   = lm(sim~obs)
  bias  = (1/length(obs)) * sum(sim-obs)
  mse   = (1/length(obs)) * sum((sim-obs)^2)
  rmse  = sqrt(mse)
  mae   = (1/length(obs)) * sum(abs(sim-obs))
  rrmse = rmse / mean(obs)
  rmae  = (1/length(obs[obs>0])) * sum(abs(sim[obs>0]-obs[obs>0])/abs(obs[obs>0]))
  ef    = 1 - (sum((sim-obs)^2) / sum((obs-mean(obs))^2))
  r     = sum((obs-mean(obs))*(sim-mean(sim)))/sqrt(sum((obs-mean(obs))^2)*sum((sim-mean(sim))^2))
  r2    = r^2
  d     = 1 - (sum((sim-obs)^2) / sum((abs(sim-mean(obs))+abs(obs-mean(obs)))^2))
  if(length(unique(sim)) > 1){
    a     = summary(fit)$coefficients["(Intercept)","Estimate"]
    b     = summary(fit)$coefficients["obs","Estimate"]
  }
  
  if(dchart){
    #--- Chart Sim ~ Obs
    varlab = vnam 
    
    mindt = min(obs,sim)
    maxdt = max(obs,sim)
    #--- Ploting limits 
    pllim = c(mindt-0.1*(maxdt-mindt),maxdt+0.1*(maxdt-mindt))
    xx = seq(min(obs),max(obs),length = (max(obs)-min(obs))*1000)
    #z = summary(fit)
    
    plot(sim~obs,
         ylab = paste("Sim - ",varlab,sep = ""),
         xlab = paste("Obs - ",varlab,sep = ""),
         ylim = pllim,
         xlim = pllim)
    
    lines(xx, predict(fit, data.frame(obs=xx)),
          col = "black",
          lty = 1,
          lwd = 1.5)
    
    l11 = seq(pllim[1]-0.5*(maxdt-mindt), pllim[2] + 0.5 * (maxdt-mindt),length = 1000)
    
    lines(l11*1~l11,
          col = "red",
          lty = 2,
          lwd = 1.5)
  }
  
  if(outidx[1] == "all"){outidx = c("bias","mse","rmse","mae","rrmse","rmae","ef","r","r2","d")}
  perf = data.frame(model = vnam,
                    bias,
                    mse,
                    rmse,
                    mae,
                    rrmse,
                    rmae,
                    ef,
                    r,
                    r2,
                    d)
  
  return(perf[,c("model",outidx)])
  
}

#--- Function to opmize parameters a,b,c
opt_fdr = function(p){
  
  #--- extract parameters from vector
  a = p[1]
  b = p[2]
  c = p[3]
  
  #--- compute soil water content based on the minimum (pmp_sf) and maximum frequency (sat_sf) 
  pmp = 10 ^ (log10((pmp_sf-c)/a)/b)
  sat = 10 ^ (log10((sat_sf-c)/a)/b)
  
  #--- compute the rmse between the pmp and sat estimated and measured (ret_curv)
  perf = mperf(sim = c(pmp,sat),obs = c(ret_curv$swc[ret_curv$depth==5 & ret_curv$h==15000], max(ret_curv$swc))*100,vnam = "SWC",dchart = F,outidx = "rmse")
  return(perf$rmse)
}


ag_plot = function(fdr_data_calib,feature,dp,ad,c,y_li,y_la){
  
  df = fdr_data_calib[,c("date","depth","tube",feature)]
  colnames(df)[4] = "feature"
  df = df[df$feature>0,]
  ag = aggregate(feature ~ date + depth, df, mean)
  ag = ag[order(ag$date),]
  if(ad){
    lines(ag$feature[ag$depth==dp]~ag$date[ag$depth==dp],
          type = "b",
          xaxt = "n",
          col = as.character(c))
  }else{
    plot(ag$feature[ag$depth==dp]~ag$date[ag$depth==dp], type = "b",xaxt = "n",
         ylim = y_li,
         ylab = y_la,
         xlab = "",
         col  = as.character(c))
    axis(1, ag$date[ag$depth==dp], format(ag$date[ag$depth==dp], "%b %d"), cex.axis = .7)  
  }
  
}
