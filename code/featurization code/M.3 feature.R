library(nnet)
library(readr)
library(resample)
library(randomForest)
library(e1071)


#get path to the folder "reproducibility_materials"
path1 = getwd()
check_path = unlist(strsplit(path1, split = "/"))
while (check_path[length(check_path)] != "reproducibility_materials"){
  path1 = dirname(path1)  
  check_path = unlist(strsplit(path1, split = "/"))
}

setwd(path1) #the path of the folder

#################
#h-np threshold selection function
np_decision <- function(t1,t2, prob, class_set ){
  n = dim(prob)[1]
  local = colnames(prob)
  index1 = which(local == class_set[1])
  index2 = which(local == class_set[2])
  index3 = which(local == class_set[3])
  d1 = which(prob[,index1] >= t1)
  d2 = which(prob[,index2]/prob[,index3] >= t2)
  decision = rep(class_set[3],n)
  decision[d2] = class_set[2]
  decision[d1] = class_set[1]
  return(decision)
}






############
#load data and construct features

aggExprs_meta <- readRDS("data/aggExprs_meta.rds")
patient_meta = readRDS("data/patient_meta.rds")
scMerge <- readRDS( "data/aggExprs_scMerge.rds")




class_set = unique(aggExprs_meta$meta_severity )

set_cellselect = unique(aggExprs_meta$celltype)[-18]
set_patient = unique(aggExprs_meta$rp_id)




patient_pca_1 = NULL

index_cell = which(aggExprs_meta$celltype %in% set_cellselect )
for(i in set_patient){
  index_patient = which(aggExprs_meta$rp_id == i)
  index_i = intersect( index_patient, index_cell )
  patient_matrix = scMerge[,index_i]
  pca = prcomp(patient_matrix,rank. = 1)
  PC1 =  patient_matrix%*%as.matrix(abs(as.numeric(pca$rotation)), ncol = 1)
  patient_pca_1 = cbind(  patient_pca_1,PC1)
}

patient_pca_1 = t(patient_pca_1)
rownames(patient_pca_1) = set_patient
colnames(patient_pca_1) = rownames(scMerge )

patient_pca_1 = as.matrix(patient_pca_1)
head(patient_pca_1[,1:10])


var_pca = colVars(patient_pca_1)
mean_pca = colMeans(patient_pca_1)
rank_pca = rank(-var_pca/abs(mean_pca))



########################


#age and condition
aggExprs_meta$meta_age[which(substr(aggExprs_meta$meta_age,1,1) == "(" )] = substr(aggExprs_meta$meta_age,2,3)[which(substr(aggExprs_meta$meta_age,1,1) == "(" )]
aggExprs_meta$meta_age[which(substr(aggExprs_meta$meta_age,3,3) == "-")] = substr(aggExprs_meta$meta_age,1,2)[which(substr(aggExprs_meta$meta_age,3,3) == "-")]
aggExprs_meta$meta_age[which(substr(aggExprs_meta$meta_age,1,2) == ">=")] = substr(aggExprs_meta$meta_age,3,4)[which(substr(aggExprs_meta$meta_age,1,2) == ">=")]
aggExprs_meta$meta_age[which(substr(aggExprs_meta$meta_age,3,3) == "+")] = substr(aggExprs_meta$meta_age,1,2)[which(substr(aggExprs_meta$meta_age,3,3) == "+")]
aggExprs_meta$age2 = aggExprs_meta$meta_age
aggExprs_meta$age2[is.na(aggExprs_meta$age2)] = 0
aggExprs_meta$age2[which(aggExprs_meta$age2 == "unknown" )] = 0
aggExprs_meta$age2[which(aggExprs_meta$age2 == "Unknown" )] = 0
aggExprs_meta$age2 = as.integer(aggExprs_meta$age2)
aggExprs_meta$age2[which(  aggExprs_meta$age2 == 0)] =  mean(aggExprs_meta$age2[which(  aggExprs_meta$age2 != 0)])


age_set = NULL
for(i in 1:length(set_patient)){
  index_patient = which( aggExprs_meta$rp_id == set_patient[i])[1]
  age_set  = c(age_set , aggExprs_meta$age2[index_patient])
}
condition = patient_meta$meta_severity
######################




factor_index = which(rank_pca<=3000)

diagnosis_matix = data.frame(patient_pca_1)
diagnosis_matix$condition = as.factor(condition)
diagnosis_matix$age= age_set

condition_set = unique(diagnosis_matix$condition)



diagnosis_matix_1 =   diagnosis_matix[which(  diagnosis_matix$condition== "Severe/Critical"),]
diagnosis_matix_2 = diagnosis_matix[which(  diagnosis_matix$condition== "Mild/Moderate"),]
diagnosis_matix_3 = diagnosis_matix[which(  diagnosis_matix$condition== "Healthy" ),]


alpha_1 = 0.2
alpha_2 = 0.2
delta_1 = 0.2
delta_2 = 0.2


n1 = sum(condition ==  "Severe/Critical")
n2 = sum(condition ==  "Mild/Moderate")
n3 = sum(condition ==  "Healthy" )
ns1 = as.integer(n1*0.35)
nt1 = as.integer(n1*0.35)
nc1 = n1 - ns1 - nt1
ns2 = as.integer(n2*0.35)
nt2 = as.integer(n2*0.25)
ne2 = as.integer(n2*0.1)
nc2 = n2 - ns2 - nt2 -ne2
ns3 = as.integer(n3*0.35)
ne3 = as.integer(n3*0.35)
nc3 = n3 - ns3 - ne3

n_all = n1 + n2 + n3

pi_1 = n1 /n_all
pi_2 = n2 /n_all
pi_3 = n3 /n_all

c = 2
c_n = c/sqrt(nt2)
#c_n = alpha_2
delta_22 = exp(-2* nt2 * c_n^2)
delta_21 = max(delta_2 - delta_22, 0)
k1_u = qbinom(delta_1, nt1, alpha_1)
k2_l = qbinom(delta_2, nt2, alpha_2)



class_set = c("Severe/Critical","Mild/Moderate","Healthy" )








##############################################
# base method: logistic

errors_train = NULL  
errors_test = NULL
errors_test_1 = NULL
errors_test_2 = NULL
errors_test_3 = NULL


errors_np_train = NULL  
errors_np_test = NULL
errors_np_test_1 = NULL
errors_np_test_2 = NULL
errors_np_test_3 = NULL
errors_np_test_12 = NULL
errors_np_test_13 = NULL
errors_np_test_21 = NULL
errors_np_test_23 = NULL
errors_np_test_31 = NULL
errors_np_test_32 = NULL
errors_np_test_star = NULL
errors_np_test_em = NULL ############


set.seed(2023)
for(i in 1:50){
  ##############
  diagnosis_matix_1 =  diagnosis_matix_1[sample(1:n1),]
  diagnosis_matix_2 =  diagnosis_matix_2[sample(1:n2),]
  diagnosis_matix_3 =  diagnosis_matix_3[sample(1:n3),]
  
  diagnosis_train  = rbind(diagnosis_matix_1[1:ns1,],diagnosis_matix_2[1:ns2,], diagnosis_matix_3[1:ns3,] )
  diagnosis_threshold = rbind(diagnosis_matix_1[(ns1 + 1) : (ns1 + nt1),],diagnosis_matix_2[(ns2 + 1) : (ns2 + nt2),])
  
  
  diagnosis_empricial_2 = diagnosis_matix_2[(ns2 + nt2 + 1) : (ns2 + nt2 + ne2),]
  diagnosis_empricial_3 = diagnosis_matix_3[(ns3  + 1) : (ns3 + ne3),]
  diagnosis_em = rbind(diagnosis_empricial_2,diagnosis_empricial_3 )
  
  
  
  diagnosis_test_1 = diagnosis_matix_1[(ns1 + nt1 + 1) : n1,]
  diagnosis_test_2 = diagnosis_matix_2[(ns2 + nt2+ ne2 + 1) : n2,]
  diagnosis_test_3 = diagnosis_matix_3[(ns3 + ne3  + 1) : n3,]
  
  diagnosis_test = rbind(diagnosis_test_1,diagnosis_test_2, diagnosis_test_3 )
  ################
  
  
  mymodel <-  multinom(condition~. ,data = diagnosis_train,MaxNWts = 200000)
  
  errors_train = c(errors_train, sum(predict(mymodel
  ) != diagnosis_train$condition)/dim(diagnosis_train)[1])
  
  error1 = sum(predict(mymodel,newdata = diagnosis_test_1
  ) != diagnosis_test_1$condition)/dim(diagnosis_test_1)[1]
  
  error12 =  sum(predict(mymodel,newdata = diagnosis_test_1
  ) == class_set[2])/dim(diagnosis_test_1)[1]
  
  error13 =  sum(predict(mymodel,newdata = diagnosis_test_1
  ) == class_set[3])/dim(diagnosis_test_1)[1]
  
  error2 = sum(predict(mymodel,newdata = diagnosis_test_2
  ) != diagnosis_test_2$condition)/dim(diagnosis_test_2)[1]
  
  error21 = sum(predict(mymodel,newdata = diagnosis_test_2
  ) == class_set[1])/dim(diagnosis_test_2)[1]
  
  error23 = sum(predict(mymodel,newdata = diagnosis_test_2
  ) == class_set[3])/dim(diagnosis_test_2)[1]
  
  error3 = sum(predict(mymodel,newdata = diagnosis_test_3
  ) != diagnosis_test_3$condition)/dim(diagnosis_test_3)[1]
  
  error31 = sum(predict(mymodel,newdata = diagnosis_test_3
  ) == class_set[1])/dim(diagnosis_test_3)[1]
  
  error32 = sum(predict(mymodel,newdata = diagnosis_test_3
  ) == class_set[2])/dim(diagnosis_test_3)[1]
  
  
  error = sum(predict(mymodel,newdata = diagnosis_test
  ) != diagnosis_test$condition)/dim(diagnosis_test)[1]
  
  
  
  class = predict(mymodel,newdata = diagnosis_test)
  
  index_set_1 = which(diagnosis_test$condition == class_set[1])
  index_set_2 = which(diagnosis_test$condition == class_set[2] & class == class_set[3] )
  index_error = c(index_set_1,index_set_2)
  error_star = sum(class[-index_error] != diagnosis_test$condition[-index_error])/length(class)
  
  
  error_gen = c(error,error1,error12,error13,error2,error21,error23,error3,error31,error32,error_star)
  errors_test = rbind(errors_test, error_gen)
  
  
  ############################################################
  
  
  
  prob_set = predict(  mymodel, newdata = diagnosis_threshold, "probs")
  local = colnames(prob_set)
  
  
  
  
  index1 = which(local == class_set[1])
  index2 = which(local == class_set[2])
  index3 = which(local == class_set[3])
  
  prob_1 =  prob_set[1:nt1, ]
  prob_2 =  prob_set[(nt1 + 1) :(nt1 + nt2), ]
  prob_em  = predict(  mymodel, newdata = diagnosis_em, "probs") #######################
  prob_e = predict(  mymodel, newdata = diagnosis_test, "probs")
  prob_e1 = predict( mymodel, newdata = diagnosis_test_1, "probs")
  prob_e2 = predict(  mymodel, newdata = diagnosis_test_2, "probs")
  prob_e3 = predict(  mymodel, newdata = diagnosis_test_3, "probs")
  set_T_1 = sort(prob_1[,index1])
  t2_l = sort(prob_2[,index2]/prob_2[,index3])[k2_l]
  error_set = NULL
  error_set_1 = NULL
  error_set_2 = NULL
  error_set_3 = NULL
  error_set_12 = NULL
  error_set_13 = NULL
  error_set_21 = NULL
  error_set_23 = NULL
  error_set_31 = NULL
  error_set_32 = NULL
  error_set_star = NULL
  error_set_em = NULL ##########
  for (k1 in 1:k1_u){
    k2 = 0
    t1 = set_T_1[k1]
    set_T_2 = prob_2[which(prob_2[,index1] < t1),index2]/prob_2[which(prob_2[,index1] < t1),index3]
    nt22 = length(set_T_2)
    p_hat =nt22 /nt2 
    p_2 = p_hat + c_n
    alpha_21 = alpha_2/p_2
    if(alpha_21 <= 1 ){
      k2 = qbinom(delta_21, nt22, alpha_21)
      t2 = sort(set_T_2)[k2]}
    if(k2 == 0){t2 = t2_l}
    
    
    ##############
    class = np_decision(t1,t2,prob_em,class_set)
    i_class = which(diagnosis_em$condition == class_set[2])
    error2em = mean(class[i_class]== class_set[1])
    i_class = which(diagnosis_em$condition == class_set[3])
    error3em = mean(class[i_class]!= class_set[3] )
    error_em = pi_2*error2em + pi_3*error3em
    error_set_em = c(error_set_em,error_em)
    #################
    
    
    
    class = np_decision(t1,t2,prob_e,class_set)
    error =  mean(class != diagnosis_test$condition)
    error_set = c(error_set,error)
    
    index_set_1 = which(diagnosis_test$condition == class_set[1])
    index_set_2 = which(diagnosis_test$condition == class_set[2] & class == class_set[3] )
    index_error = c(index_set_1,index_set_2)
    error_star = sum(class[-index_error] != diagnosis_test$condition[-index_error])/length(class)
    error_set_star = c(error_set_star,error_star)
    
    class = np_decision(t1,t2,prob_e1,class_set)
    error = mean(class != class_set[1])
    error_set_1 = c(error_set_1,error)
    error = mean(class == class_set[2])
    error_set_12 = c(error_set_12,error)
    error = mean(class == class_set[3])
    error_set_13 = c(error_set_13,error)
    
    class = np_decision(t1,t2,prob_e2,class_set)
    error = mean(class != class_set[2])
    error_set_2 = c(error_set_2,error)
    error = mean(class == class_set[1])
    error_set_21 = c(error_set_21,error)
    error = mean(class == class_set[3])
    error_set_23 = c(error_set_23,error)
    
    
    class = np_decision(t1,t2,prob_e3,class_set)
    error = mean(class != class_set[3])
    error_set_3 = c(error_set_3,error)
    error = mean(class == class_set[1])
    error_set_31 = c(error_set_31,error)
    error = mean(class == class_set[2])
    error_set_32 = c(error_set_32,error)
  }
  
  
  
  errors_np_test = rbind(errors_np_test,error_set)
  errors_np_test_1 = rbind(errors_np_test_1,error_set_1)
  errors_np_test_2 = rbind(errors_np_test_2,error_set_2)
  errors_np_test_3 = rbind(errors_np_test_3,error_set_3)
  errors_np_test_12 = rbind(errors_np_test_12,error_set_12)
  errors_np_test_13 = rbind(errors_np_test_13,error_set_13)
  errors_np_test_21 = rbind(errors_np_test_21,error_set_21)
  errors_np_test_23 = rbind(errors_np_test_23,error_set_23)
  errors_np_test_31 = rbind(errors_np_test_31,error_set_31)
  errors_np_test_32 = rbind(errors_np_test_32,error_set_32)
  errors_np_test_star = rbind(errors_np_test_star,error_set_star)
  errors_np_test_em = rbind(errors_np_test_em,error_set_em)
}





list_error = list(non_np = errors_test,
                  overall = errors_np_test,
                  error1 = errors_np_test_1,
                  error12 = errors_np_test_12,
                  error13 = errors_np_test_13,
                  error2 = errors_np_test_2,
                  error21 = errors_np_test_21,
                  error23 = errors_np_test_23,
                  error3 = errors_np_test_3,
                  error31 = errors_np_test_31,
                  error32 = errors_np_test_32,
                  errorstar = errors_np_test_star,
                  errorem = errors_np_test_em)

save(list_error, file = "output/logistic_m3.Rdata")





##############################################
# base method: random forest


factor_index = which(rank_pca<=3000)
diagnosis_matix = data.frame(patient_pca_1[,factor_index ])
diagnosis_matix$condition = as.factor(condition)
diagnosis_matix$age= age_set


condition_set = unique(diagnosis_matix$condition)




diagnosis_matix_1 = diagnosis_matix[which(  diagnosis_matix$condition== "Severe/Critical"),]
diagnosis_matix_2 = diagnosis_matix[which(  diagnosis_matix$condition== "Mild/Moderate"),]
diagnosis_matix_3 = diagnosis_matix[which(  diagnosis_matix$condition== "Healthy" ),]
errors_train = NULL  
errors_test = NULL
errors_test_1 = NULL
errors_test_2 = NULL
errors_test_3 = NULL

errors_np_train = NULL  
errors_np_test = NULL
errors_np_test_1 = NULL
errors_np_test_2 = NULL
errors_np_test_3 = NULL
errors_np_test_12 = NULL
errors_np_test_13 = NULL
errors_np_test_21 = NULL
errors_np_test_23 = NULL
errors_np_test_31 = NULL
errors_np_test_32 = NULL
errors_np_test_star = NULL
errors_np_test_em = NULL

set.seed(2022)
for(i in 1:50){
  ##############
  diagnosis_matix_1 =  diagnosis_matix_1[sample(1:n1),]
  diagnosis_matix_2 =  diagnosis_matix_2[sample(1:n2),]
  diagnosis_matix_3 =  diagnosis_matix_3[sample(1:n3),]
  
  diagnosis_train  = rbind(diagnosis_matix_1[1:ns1,],diagnosis_matix_2[1:ns2,], diagnosis_matix_3[1:ns3,] )
  diagnosis_threshold = rbind(diagnosis_matix_1[(ns1 + 1) : (ns1 + nt1),],diagnosis_matix_2[(ns2 + 1) : (ns2 + nt2),])
  
  
  diagnosis_empricial_2 = diagnosis_matix_2[(ns2 + nt2 + 1) : (ns2 + nt2 + ne2),]
  diagnosis_empricial_3 = diagnosis_matix_3[(ns3  + 1) : (ns3 + ne3),]
  diagnosis_em = rbind(diagnosis_empricial_2,diagnosis_empricial_3 )
  
  
  
  diagnosis_test_1 = diagnosis_matix_1[(ns1 + nt1 + 1) : n1,]
  diagnosis_test_2 = diagnosis_matix_2[(ns2 + nt2+ ne2 + 1) : n2,]
  diagnosis_test_3 = diagnosis_matix_3[(ns3 + ne3  + 1) : n3,]
  
  diagnosis_test = rbind(diagnosis_test_1,diagnosis_test_2, diagnosis_test_3 )
  ################
  
  
  mymodel <-   randomForest(as.factor(condition)~., data=diagnosis_train , proximity=TRUE,ntree = 500)
  
  
  errors_train = c(errors_train, sum(predict(mymodel
  ) != diagnosis_train$condition)/dim(diagnosis_train)[1])
  
  error1 = sum(predict(mymodel,newdata = diagnosis_test_1
  ) != diagnosis_test_1$condition)/dim(diagnosis_test_1)[1]
  
  error12 =  sum(predict(mymodel,newdata = diagnosis_test_1
  ) == class_set[2])/dim(diagnosis_test_1)[1]
  
  error13 =  sum(predict(mymodel,newdata = diagnosis_test_1
  ) == class_set[3])/dim(diagnosis_test_1)[1]
  
  error2 = sum(predict(mymodel,newdata = diagnosis_test_2
  ) != diagnosis_test_2$condition)/dim(diagnosis_test_2)[1]
  
  error21 = sum(predict(mymodel,newdata = diagnosis_test_2
  ) == class_set[1])/dim(diagnosis_test_2)[1]
  
  error23 = sum(predict(mymodel,newdata = diagnosis_test_2
  ) == class_set[3])/dim(diagnosis_test_2)[1]
  
  error3 = sum(predict(mymodel,newdata = diagnosis_test_3
  ) != diagnosis_test_3$condition)/dim(diagnosis_test_3)[1]
  
  error31 = sum(predict(mymodel,newdata = diagnosis_test_3
  ) == class_set[1])/dim(diagnosis_test_3)[1]
  
  error32 = sum(predict(mymodel,newdata = diagnosis_test_3
  ) == class_set[2])/dim(diagnosis_test_3)[1]
  
  
  error = sum(predict(mymodel,newdata = diagnosis_test
  ) != diagnosis_test$condition)/dim(diagnosis_test)[1]
  
  
  
  class = predict(mymodel,newdata = diagnosis_test)
  
  index_set_1 = which(diagnosis_test$condition == class_set[1])
  index_set_2 = which(diagnosis_test$condition == class_set[2] & class == class_set[3] )
  index_error = c(index_set_1,index_set_2)
  error_star = sum(class[-index_error] != diagnosis_test$condition[-index_error])/length(class)
  
  
  error_gen = c(error,error1,error12,error13,error2,error21,error23,error3,error31,error32,error_star)
  errors_test = rbind(errors_test, error_gen)
  
  
  ############################################################
  
  
  
  prob_set = predict(  mymodel, newdata = diagnosis_threshold, "prob")
  local = colnames(prob_set)
  
  
  
  
  index1 = which(local == class_set[1])
  index2 = which(local == class_set[2])
  index3 = which(local == class_set[3])
  
  prob_1 =  prob_set[1:nt1, ]
  prob_2 =  prob_set[(nt1 + 1) :(nt1 + nt2), ]
  prob_em  = predict(  mymodel, newdata = diagnosis_em, "prob") #######################
  prob_e = predict(  mymodel, newdata = diagnosis_test, "prob")
  prob_e1 = predict( mymodel, newdata = diagnosis_test_1, "prob")
  prob_e2 = predict(  mymodel, newdata = diagnosis_test_2, "prob")
  prob_e3 = predict(  mymodel, newdata = diagnosis_test_3, "prob")
  set_T_1 = sort(prob_1[,index1])
  t2_l = sort(prob_2[,index2]/prob_2[,index3])[k2_l]
  error_set = NULL
  error_set_1 = NULL
  error_set_2 = NULL
  error_set_3 = NULL
  error_set_12 = NULL
  error_set_13 = NULL
  error_set_21 = NULL
  error_set_23 = NULL
  error_set_31 = NULL
  error_set_32 = NULL
  error_set_star = NULL
  error_set_em = NULL ##########
  for (k1 in 1:k1_u){
    k2 = 0
    t1 = set_T_1[k1]
    set_T_2 = prob_2[which(prob_2[,index1] < t1),index2]/prob_2[which(prob_2[,index1] < t1),index3]
    nt22 = length(set_T_2)
    p_hat =nt22 /nt2 
    p_2 = p_hat + c_n
    alpha_21 = alpha_2/p_2
    if(alpha_21 <= 1 ){
      k2 = qbinom(delta_21, nt22, alpha_21)
      t2 = sort(set_T_2)[k2]}
    if(k2 == 0){t2 = t2_l}
    
    
    ##############
    class = np_decision(t1,t2,prob_em,class_set)
    i_class = which(diagnosis_em$condition == class_set[2])
    error2em = mean(class[i_class]== class_set[1])
    i_class = which(diagnosis_em$condition == class_set[3])
    error3em = mean(class[i_class]!= class_set[3] )
    error_em = pi_2*error2em + pi_3*error3em
    error_set_em = c(error_set_em,error_em)
    #################
    
    
    
    class = np_decision(t1,t2,prob_e,class_set)
    error =  mean(class != diagnosis_test$condition)
    error_set = c(error_set,error)
    
    index_set_1 = which(diagnosis_test$condition == class_set[1])
    index_set_2 = which(diagnosis_test$condition == class_set[2] & class == class_set[3] )
    index_error = c(index_set_1,index_set_2)
    error_star = sum(class[-index_error] != diagnosis_test$condition[-index_error])/length(class)
    error_set_star = c(error_set_star,error_star)
    
    class = np_decision(t1,t2,prob_e1,class_set)
    error = mean(class != class_set[1])
    error_set_1 = c(error_set_1,error)
    error = mean(class == class_set[2])
    error_set_12 = c(error_set_12,error)
    error = mean(class == class_set[3])
    error_set_13 = c(error_set_13,error)
    
    class = np_decision(t1,t2,prob_e2,class_set)
    error = mean(class != class_set[2])
    error_set_2 = c(error_set_2,error)
    error = mean(class == class_set[1])
    error_set_21 = c(error_set_21,error)
    error = mean(class == class_set[3])
    error_set_23 = c(error_set_23,error)
    
    
    class = np_decision(t1,t2,prob_e3,class_set)
    error = mean(class != class_set[3])
    error_set_3 = c(error_set_3,error)
    error = mean(class == class_set[1])
    error_set_31 = c(error_set_31,error)
    error = mean(class == class_set[2])
    error_set_32 = c(error_set_32,error)
  }
  
  
  
  errors_np_test = rbind(errors_np_test,error_set)
  errors_np_test_1 = rbind(errors_np_test_1,error_set_1)
  errors_np_test_2 = rbind(errors_np_test_2,error_set_2)
  errors_np_test_3 = rbind(errors_np_test_3,error_set_3)
  errors_np_test_12 = rbind(errors_np_test_12,error_set_12)
  errors_np_test_13 = rbind(errors_np_test_13,error_set_13)
  errors_np_test_21 = rbind(errors_np_test_21,error_set_21)
  errors_np_test_23 = rbind(errors_np_test_23,error_set_23)
  errors_np_test_31 = rbind(errors_np_test_31,error_set_31)
  errors_np_test_32 = rbind(errors_np_test_32,error_set_32)
  errors_np_test_star = rbind(errors_np_test_star,error_set_star)
  errors_np_test_em = rbind(errors_np_test_em,error_set_em)
}


list_error = list(non_np = errors_test,
                  overall = errors_np_test,
                  error1 = errors_np_test_1,
                  error12 = errors_np_test_12,
                  error13 = errors_np_test_13,
                  error2 = errors_np_test_2,
                  error21 = errors_np_test_21,
                  error23 = errors_np_test_23,
                  error3 = errors_np_test_3,
                  error31 = errors_np_test_31,
                  error32 = errors_np_test_32,
                  errorstar = errors_np_test_star,
                  errorem = errors_np_test_em)

save(list_error, file = "output/forest_m3.Rdata")


##############################################
# base method: svm



factor_index = which(rank_pca<=3000)
diagnosis_matix = data.frame(patient_pca_1[,factor_index ])
diagnosis_matix$condition = as.factor(condition)
diagnosis_matix$age= age_set


condition_set = unique(diagnosis_matix$condition)





diagnosis_matix_1 =   diagnosis_matix[which(  diagnosis_matix$condition== "Severe/Critical"),]
diagnosis_matix_2 = diagnosis_matix[which(  diagnosis_matix$condition== "Mild/Moderate"),]
diagnosis_matix_3 = diagnosis_matix[which(  diagnosis_matix$condition== "Healthy" ),]

errors_train = NULL  
errors_test = NULL
errors_test_1 = NULL
errors_test_2 = NULL
errors_test_3 = NULL

errors_np_train = NULL  
errors_np_test = NULL
errors_np_test_1 = NULL
errors_np_test_2 = NULL
errors_np_test_3 = NULL
errors_np_test_12 = NULL
errors_np_test_13 = NULL
errors_np_test_21 = NULL
errors_np_test_23 = NULL
errors_np_test_31 = NULL
errors_np_test_32 = NULL
errors_np_test_star = NULL
errors_np_test_em = NULL

set.seed(2022)
for(i in 1:50){
  ##############
  diagnosis_matix_1 =  diagnosis_matix_1[sample(1:n1),]
  diagnosis_matix_2 =  diagnosis_matix_2[sample(1:n2),]
  diagnosis_matix_3 =  diagnosis_matix_3[sample(1:n3),]
  
  diagnosis_train  = rbind(diagnosis_matix_1[1:ns1,],diagnosis_matix_2[1:ns2,], diagnosis_matix_3[1:ns3,] )
  diagnosis_threshold = rbind(diagnosis_matix_1[(ns1 + 1) : (ns1 + nt1),],diagnosis_matix_2[(ns2 + 1) : (ns2 + nt2),])
  
  
  diagnosis_empricial_2 = diagnosis_matix_2[(ns2 + nt2 + 1) : (ns2 + nt2 + ne2),]
  diagnosis_empricial_3 = diagnosis_matix_3[(ns3  + 1) : (ns3 + ne3),]
  diagnosis_em = rbind(diagnosis_empricial_2,diagnosis_empricial_3 )
  
  
  
  diagnosis_test_1 = diagnosis_matix_1[(ns1 + nt1 + 1) : n1,]
  diagnosis_test_2 = diagnosis_matix_2[(ns2 + nt2+ ne2 + 1) : n2,]
  diagnosis_test_3 = diagnosis_matix_3[(ns3 + ne3  + 1) : n3,]
  
  diagnosis_test = rbind(diagnosis_test_1,diagnosis_test_2, diagnosis_test_3 )
  ################
  
  
  mymodel <-   svm(as.factor(condition)~., data=diagnosis_train, kernel =
                     "linear", probability = TRUE)
  
  
  errors_train = c(errors_train, sum(predict(mymodel
  ) != diagnosis_train$condition)/dim(diagnosis_train)[1])
  
  error1 = sum(predict(mymodel,newdata = diagnosis_test_1
  ) != diagnosis_test_1$condition)/dim(diagnosis_test_1)[1]
  
  error12 =  sum(predict(mymodel,newdata = diagnosis_test_1
  ) == class_set[2])/dim(diagnosis_test_1)[1]
  
  error13 =  sum(predict(mymodel,newdata = diagnosis_test_1
  ) == class_set[3])/dim(diagnosis_test_1)[1]
  
  error2 = sum(predict(mymodel,newdata = diagnosis_test_2
  ) != diagnosis_test_2$condition)/dim(diagnosis_test_2)[1]
  
  error21 = sum(predict(mymodel,newdata = diagnosis_test_2
  ) == class_set[1])/dim(diagnosis_test_2)[1]
  
  error23 = sum(predict(mymodel,newdata = diagnosis_test_2
  ) == class_set[3])/dim(diagnosis_test_2)[1]
  
  error3 = sum(predict(mymodel,newdata = diagnosis_test_3
  ) != diagnosis_test_3$condition)/dim(diagnosis_test_3)[1]
  
  error31 = sum(predict(mymodel,newdata = diagnosis_test_3
  ) == class_set[1])/dim(diagnosis_test_3)[1]
  
  error32 = sum(predict(mymodel,newdata = diagnosis_test_3
  ) == class_set[2])/dim(diagnosis_test_3)[1]
  
  
  error = sum(predict(mymodel,newdata = diagnosis_test
  ) != diagnosis_test$condition)/dim(diagnosis_test)[1]
  
  
  
  class = predict(mymodel,newdata = diagnosis_test)
  
  index_set_1 = which(diagnosis_test$condition == class_set[1])
  index_set_2 = which(diagnosis_test$condition == class_set[2] & class == class_set[3] )
  index_error = c(index_set_1,index_set_2)
  error_star = sum(class[-index_error] != diagnosis_test$condition[-index_error])/length(class)
  
  
  error_gen = c(error,error1,error12,error13,error2,error21,error23,error3,error31,error32,error_star)
  errors_test = rbind(errors_test, error_gen)
  
  
  ############################################################
  
  
  
  prob_set = attr(predict(  mymodel, newdata = diagnosis_threshold, probability = TRUE),"probabilities")
  local = colnames(prob_set)
  
  
  
  
  index1 = which(local == class_set[1])
  index2 = which(local == class_set[2])
  index3 = which(local == class_set[3])
  
  prob_1 =  prob_set[1:nt1, ]
  prob_2 =  prob_set[(nt1 + 1) :(nt1 + nt2), ]
  
  
  
  prob_set =attr(predict(  mymodel, newdata = diagnosis_threshold, probability = TRUE),"probabilities")
  local = colnames(prob_set)
  
  prob_em = attr(predict(  mymodel, newdata = diagnosis_em, probability = TRUE),"probabilities")
  prob_e = attr(predict(  mymodel, newdata = diagnosis_test, probability = TRUE),"probabilities")
  prob_e1 = attr(predict( mymodel, newdata = diagnosis_test_1, probability = TRUE),"probabilities")
  prob_e2 = attr(predict(  mymodel, newdata = diagnosis_test_2, probability = TRUE),"probabilities")
  prob_e3 = attr(predict(  mymodel, newdata = diagnosis_test_3, probability = TRUE),"probabilities")
  set_T_1 = sort(prob_1[,index1])
  t2_l = sort(prob_2[,index2]/prob_2[,index3])[k2_l]
  error_set = NULL
  error_set_1 = NULL
  error_set_2 = NULL
  error_set_3 = NULL
  error_set_12 = NULL
  error_set_13 = NULL
  error_set_21 = NULL
  error_set_23 = NULL
  error_set_31 = NULL
  error_set_32 = NULL
  error_set_star = NULL
  error_set_em = NULL ##########
  for (k1 in 1:k1_u){
    k2 = 0
    t1 = set_T_1[k1]
    set_T_2 = prob_2[which(prob_2[,index1] < t1),index2]/prob_2[which(prob_2[,index1] < t1),index3]
    nt22 = length(set_T_2)
    p_hat =nt22 /nt2 
    p_2 = p_hat + c_n
    alpha_21 = alpha_2/p_2
    if(alpha_21 <= 1 ){
      k2 = qbinom(delta_21, nt22, alpha_21)
      t2 = sort(set_T_2)[k2]}
    if(k2 == 0){t2 = t2_l}
    
    
    ##############
    class = np_decision(t1,t2,prob_em,class_set)
    i_class = which(diagnosis_em$condition == class_set[2])
    error2em = mean(class[i_class]== class_set[1])
    i_class = which(diagnosis_em$condition == class_set[3])
    error3em = mean(class[i_class]!= class_set[3] )
    error_em = pi_2*error2em + pi_3*error3em
    error_set_em = c(error_set_em,error_em)
    #################
    
    
    
    class = np_decision(t1,t2,prob_e,class_set)
    error =  mean(class != diagnosis_test$condition)
    error_set = c(error_set,error)
    
    index_set_1 = which(diagnosis_test$condition == class_set[1])
    index_set_2 = which(diagnosis_test$condition == class_set[2] & class == class_set[3] )
    index_error = c(index_set_1,index_set_2)
    error_star = sum(class[-index_error] != diagnosis_test$condition[-index_error])/length(class)
    error_set_star = c(error_set_star,error_star)
    
    class = np_decision(t1,t2,prob_e1,class_set)
    error = mean(class != class_set[1])
    error_set_1 = c(error_set_1,error)
    error = mean(class == class_set[2])
    error_set_12 = c(error_set_12,error)
    error = mean(class == class_set[3])
    error_set_13 = c(error_set_13,error)
    
    class = np_decision(t1,t2,prob_e2,class_set)
    error = mean(class != class_set[2])
    error_set_2 = c(error_set_2,error)
    error = mean(class == class_set[1])
    error_set_21 = c(error_set_21,error)
    error = mean(class == class_set[3])
    error_set_23 = c(error_set_23,error)
    
    
    class = np_decision(t1,t2,prob_e3,class_set)
    error = mean(class != class_set[3])
    error_set_3 = c(error_set_3,error)
    error = mean(class == class_set[1])
    error_set_31 = c(error_set_31,error)
    error = mean(class == class_set[2])
    error_set_32 = c(error_set_32,error)
  }

  
  
  errors_np_test = rbind(errors_np_test,error_set)
  errors_np_test_1 = rbind(errors_np_test_1,error_set_1)
  errors_np_test_2 = rbind(errors_np_test_2,error_set_2)
  errors_np_test_3 = rbind(errors_np_test_3,error_set_3)
  errors_np_test_12 = rbind(errors_np_test_12,error_set_12)
  errors_np_test_13 = rbind(errors_np_test_13,error_set_13)
  errors_np_test_21 = rbind(errors_np_test_21,error_set_21)
  errors_np_test_23 = rbind(errors_np_test_23,error_set_23)
  errors_np_test_31 = rbind(errors_np_test_31,error_set_31)
  errors_np_test_32 = rbind(errors_np_test_32,error_set_32)
  errors_np_test_star = rbind(errors_np_test_star,error_set_star)
  errors_np_test_em = rbind(errors_np_test_em,error_set_em)
}



list_error = list(non_np = errors_test,
                  overall = errors_np_test,
                  error1 = errors_np_test_1,
                  error12 = errors_np_test_12,
                  error13 = errors_np_test_13,
                  error2 = errors_np_test_2,
                  error21 = errors_np_test_21,
                  error23 = errors_np_test_23,
                  error3 = errors_np_test_3,
                  error31 = errors_np_test_31,
                  error32 = errors_np_test_32,
                  errorstar = errors_np_test_star,
                  errorem = errors_np_test_em)


save(list_error, file = "output/svm_m3.Rdata")



