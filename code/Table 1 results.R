########################
#print results in Table 1 in Section 3.2

library(resample)
setwd("~/Desktop/reproducibility_materials") #the path of the folder


base_method =  c("logistic", "forest", "svm")
for(i_md in 1:3 ){
print(base_method[i_md])
load(paste("output/",base_method[i_md],"_m1.Rdata", sep = ""))
error_summary = NULL
for(i in 2:12){
  error_i = NULL
  for(j in 1:50){
    min_j = max(which(list_error[[13]][j,] == min(list_error[[13]][j,] )))
    error_i = c(error_i,list_error[[i]][j,min_j ])
    
  }
  error_summary = cbind( error_summary,c(mean(error_i),sd(error_i)))
}
error_summary =  rbind(colMeans(list_error$non_np),
                       colStdevs(list_error$non_np),error_summary)
colnames(error_summary) = c("overall","error1","error12","error13","error2","error21","error23","error3","error31","error32", "error*")
rownames(error_summary) =c("MEAN","SD","MEAN(NP)","SD(NP)")

set_col = NULL
for(name in c("error1","error23","error21","error31","error32","overall")){
  set_col = c(set_col, which(colnames(error_summary) == name))
}

print("M.1")
print(round(error_summary[,set_col], digits = 3)[c(1,3),])



load(paste("output/",base_method[i_md],"_m2.Rdata", sep = ""))

error_summary = NULL
for(i in 2:12){
  error_i = NULL
  for(j in 1:50){
    min_j = max(which(list_error[[13]][j,] == min(list_error[[13]][j,] )))
    error_i = c(error_i,list_error[[i]][j,min_j ])
    
  }
  error_summary = cbind( error_summary,c(mean(error_i),sd(error_i)))
}
error_summary =  rbind(colMeans(list_error$non_np),
                       colStdevs(list_error$non_np),error_summary)
colnames(error_summary) = c("overall","error1","error12","error13","error2","error21","error23","error3","error31","error32", "error*")
rownames(error_summary) =c("MEAN","SD","MEAN(NP)","SD(NP)")

set_col = NULL
for(name in c("error1","error23","error21","error31","error32","overall")){
  set_col = c(set_col, which(colnames(error_summary) == name))
}


print("M.2")
print(round(error_summary[,set_col], digits = 3)[c(1,3),])



load(paste("output/",base_method[i_md],"_m3.Rdata", sep = ""))

library(resample)

error_summary = NULL
for(i in 2:12){
  error_i = NULL
  for(j in 1:50){
    min_j = max(which(list_error[[13]][j,] == min(list_error[[13]][j,] )))
    error_i = c(error_i,list_error[[i]][j,min_j ])
    
  }
  error_summary = cbind( error_summary,c(mean(error_i),sd(error_i)))
}
error_summary =  rbind(colMeans(list_error$non_np),
                       colStdevs(list_error$non_np),error_summary)
colnames(error_summary) = c("overall","error1","error12","error13","error2","error21","error23","error3","error31","error32", "error*")
rownames(error_summary) =c("MEAN","SD","MEAN(NP)","SD(NP)")

set_col = NULL
for(name in c("error1","error23","error21","error31","error32","overall")){
  set_col = c(set_col, which(colnames(error_summary) == name))
}

print("M.3")
print(round(error_summary[,set_col], digits = 3)[c(1,3),])




load(paste("output/",base_method[i_md],"_m4.Rdata", sep = ""))

error_summary = NULL
for(i in 2:12){
  error_i = NULL
  for(j in 1:50){
    min_j = max(which(list_error[[13]][j,] == min(list_error[[13]][j,] )))
    error_i = c(error_i,list_error[[i]][j,min_j ])
    
  }
  error_summary = cbind( error_summary,c(mean(error_i),sd(error_i)))
}
error_summary =  rbind(colMeans(list_error$non_np),
                       colStdevs(list_error$non_np),error_summary)
colnames(error_summary) = c("overall","error1","error12","error13","error2","error21","error23","error3","error31","error32", "error*")
rownames(error_summary) =c("MEAN","SD","MEAN(NP)","SD(NP)")

set_col = NULL
for(name in c("error1","error23","error21","error31","error32","overall")){
  set_col = c(set_col, which(colnames(error_summary) == name))
}


print("M.4")
print(round(error_summary[,set_col], digits = 3)[c(1,3),])}

