############################
# generating Figure 2
library(nnet)

#############################
# Classifer with respect to t1 and t2
np_decision <- function(t1,t2, prob){
  n = dim(prob)[1]
  d1 = which(prob[,1] >= t1)
  d2 = which(prob[,2]/prob[,3] >= t2)
  decision = rep(3,n)
  decision[d2] = 2
  decision[d1] = 1
  return(decision)
}
#############################

#error controls and tolerances
alpha_1 = 0.05
alpha_2 = 0.05
delta_1 = 0.05
delta_2 = 0.05

###############
#the size of samplesize
n1 = 500 # class 1
n2 = 500# class 2
n3 = 500# class 3
n_all = n1 + n2 + n3

pi_1 = n1 /n_all
pi_2 = n2 /n_all
pi_3 = n3 /n_all

######################
#the size of validation set
n_check = 60000
nc1 = n_check*pi_1
nc2 = n_check*pi_2
nc3 = n_check*pi_3

ns1 = n1*0.5#the size of sample used to train scoring function for class 1
nt1 = n1*0.5 #the size of sample used to choose threshold t1
ns2 = n2*0.45
nt2 = n2*0.5
ne2 = n2*0.05
ns3 = n3*0.95
ne3 = n3*0.05

c = 2
c_n = c/sqrt(nt2)
delta_22 = exp(-2* nt2 * c_n^2)
delta_21 = max(delta_2 - delta_22, 0)
k1_u = qbinom(delta_1, nt1, alpha_1)


error_non = NULL
error_sets = NULL

error_sets_1 = NULL
error_sets_2 = NULL
error_sets_3 = NULL
error_sets_12 = NULL
error_sets_13 = NULL
error_sets_21 = NULL
error_sets_23 = NULL
error_sets_31 = NULL
error_sets_32 = NULL
error_sets_star = NULL
error_sets_e = NULL
error_sets_32_0  = NULL
error_sets_23_0 = NULL
error_sets_star_0 = NULL
error_sets_e_0  = NULL
set.seed(202201)

for(ii in 1:1000){
  
  set_1 = data.frame(Y = rep(1,n1), X = rnorm(n1, 0, 1), X2 = rnorm(n1, -1, 1))
  set_2 = data.frame(Y = rep(2,n2), X = rnorm(n2,-1, 1), X2 = rnorm(n2, 1, 1))
  set_3 = data.frame(Y = rep(3,n3), X = rnorm(n3, 1, 1), X2 = rnorm(n3, 0, 1))
  
  
  set_s = rbind(set_1[1:ns1,],set_2[1:ns2,], set_3[1:ns3,] )
  set_t = rbind(set_1[(ns1 + 1) : (ns1 + nt1),],set_2[(ns2 + 1) : (ns2 + nt2),] )
  
  set_2e = set_2[(ns2 + nt2 + 1) : n2,] 
  set_3e = set_3[(ns3  + 1) : n3,] 
  set_e = rbind(set_2e, set_3e)
  
  
  
  set_1c = data.frame(Y = rep(1,nc1), X = rnorm(nc1, 0, 1), X2 = rnorm(nc1, -1, 1))
  set_2c = data.frame(Y = rep(2,nc2), X = rnorm(nc2,-1, 1), X2 = rnorm(nc2, 1, 1))
  set_3c = data.frame(Y = rep(3,nc3), X = rnorm(nc3, 1, 1), X2 = rnorm(nc3, 0, 1))
  set_c = rbind(set_1c,set_2c, set_3c)
  
  
  model <-  multinom(as.factor(Y)~.,data = set_s)
  
  
  class = predict(model,newdata = set_c)
  error = 1 - sum(class == set_c$Y)/(nc1 + nc2 + nc3)
  i_class = which(set_c$Y == 1)
  error1 = mean(class[i_class]!= 1)
  error12 = mean(class[i_class]== 2)
  error13 = mean(class[i_class]== 3)
  
  i_class = which(set_c$Y == 2)
  error2 = mean(class[i_class]!= 2)
  error21 = mean(class[i_class]== 1)
  error23 = mean(class[i_class]== 3)
  
  
  
  i_class = which(set_c$Y == 3)
  error3 = mean(class[i_class]!= 3)
  
  error31 = mean(class[i_class]== 1)
  error32 = mean(class[i_class]== 2)
  
  error_non = rbind(error_non, c(error, error1,error12,error13,error2,error21,error23, error3,error31,error32))
  
  prob_set = predict(model, newdata = set_t, "probs")
  prob_1 =  prob_set[which( set_t$Y == 1), ]
  prob_2 =  prob_set[which( set_t$Y == 2), ]
  prob_e = predict(model, newdata = set_e, "probs")
  prob_c = predict(model, newdata = set_c, "probs")
  set_T_1 = sort(prob_1[,1])
  k2_l = qbinom(delta_2, nt2, alpha_2)
  t2_l = sort(prob_2[,2]/prob_2[,3])[k2_l]
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
  error_set_e  = NULL
  for (k1 in 1:k1_u){
    t1 = set_T_1[k1]
    set_T_2 = prob_2[which(prob_2[,1] < t1),2]/prob_2[which(prob_2[,1] < t1),3]
    nt22 = length(set_T_2)
    p_hat =nt22 /nt2 
    p_2 = p_hat + c_n
    alpha_21 = alpha_2/p_2
    if(alpha_21 <= 1) {k2 = qbinom(delta_21, nt22, alpha_21)}else{k2 = 0}
    t2 = sort(set_T_2)[k2]
    if(k2 == 0){t2 = t2_l}
    
    
    class = np_decision(t1,t2,prob_e)
    i_class = which(set_e$Y == 2)
    error2e = mean(class[i_class]== 1)
    i_class = which(set_e$Y == 3)
    error3e = mean(class[i_class]!= 3 )
    error_e = pi_2*error2e + pi_3*error3e
    error_set_e = c(error_set_e,error_e)
    
    class = np_decision(t1,t2,prob_c)
    error = 1 - sum(class == set_c$Y)/(nc1 + nc2 + nc3)
    error_set = c(error_set,error)
    i_class = which(set_c$Y == 1)
    error = mean(class[i_class]!= 1)
    error_set_1 = c(error_set_1,error)
    error = mean(class[i_class]== 2)
    error_set_12 = c(error_set_12,error)
    error = mean(class[i_class]== 3)
    error_set_13 = c(error_set_13,error)
    
    i_class = which(set_c$Y == 2)
    error = mean(class[i_class]!= 2)
    error_set_2 = c(error_set_2,error)
    error = mean(class[i_class]== 1)
    error_set_21 = c(error_set_21,error)
    error = mean(class[i_class]== 3)
    error_set_23 = c(error_set_23,error)
    
    
    i_class = which(set_c$Y == 3)
    error = mean(class[i_class]!= 3)
    error_set_3 = c(error_set_3,error)
    error = mean(class[i_class]== 1)
    error_set_31 = c(error_set_31,error)
    error = mean(class[i_class]== 2)
    error_set_32 = c(error_set_32,error)
    
    n_class =c( which(set_c$Y == 1), which(set_c$Y == 2 & class == 3))
    
    error = sum(class[-n_class]!= set_c$Y[-n_class])/(nc1 + nc2 + nc3)
    error_set_star = c(error_set_star,error)
  }
  error_sets = rbind(error_sets,error_set)
  error_sets_1 = rbind(error_sets_1,error_set_1)
  error_sets_12 = rbind(error_sets_12,error_set_12)
  error_sets_13= rbind(error_sets_13,error_set_13)
  error_sets_2 = rbind(error_sets_2,error_set_2)
  error_sets_21 = rbind(error_sets_21,error_set_21)
  error_sets_23 = rbind(error_sets_23,error_set_23)
  error_sets_3 = rbind(error_sets_3,error_set_3)
  error_sets_31 = rbind(error_sets_31,error_set_31)
  error_sets_32 = rbind(error_sets_32,error_set_32)
  error_sets_star = rbind(error_sets_star,error_set_star)
  error_sets_e = rbind(error_sets_e,error_set_e)
  
  
  error_set_32_0  = NULL
  error_set_23_0 = NULL
  error_set_star_0 = NULL
  error_set_e_0 = NULL
  for (k1 in 1:k1_u){
    t1 = set_T_1[k1]
    set_T_2 = prob_2[which(prob_2[,1] < t1),2]/prob_2[which(prob_2[,1] < t1),3]
    nt22 = length(set_T_2)
    p_hat =nt22 /nt2 
    p_2 = p_hat + c_n
    alpha_21 = alpha_2/p_2
    if(alpha_21 <= 1) {k2 = qbinom(delta_21, nt22, alpha_21)}else{k2 = 0}
    t2 = sort(set_T_2)[k2]
    k2 = 0
    if(k2 == 0){t2 = t2_l}
    
    class = np_decision(t1,t2,prob_e)
    i_class = which(set_e$Y == 2)
    error2e = mean(class[i_class]== 1)
    i_class = which(set_e$Y == 3)
    error3e = mean(class[i_class]!= 3 )
    error_e = pi_2*error2e + pi_3*error3e
    error_set_e_0 = c(error_set_e_0,error_e)
    
    class = np_decision(t1,t2,prob_c)
    
    
    i_class = which(set_c$Y == 2)
    error = mean(class[i_class]== 3)
    error_set_23_0 = c(error_set_23_0,error)
    n_class =c( which(set_c$Y == 1), which(set_c$Y == 2 & class == 3))
    
    error = sum(class[-n_class]!= set_c$Y[-n_class])/(nc1 + nc2 + nc3)
    error_set_star_0 = c(error_set_star_0,error)
    
    i_class = which(set_c$Y == 3)
    error = mean(class[i_class]== 2)
    error_set_32_0 = c(error_set_32_0,error)}
  error_sets_32_0 = rbind(error_sets_32_0,error_set_32_0)
  error_sets_23_0 = rbind(error_sets_23_0,error_set_23_0)
  error_sets_star_0 = rbind(error_sets_star_0,error_set_star_0)
  error_sets_e_0  = rbind(error_sets_e_0 ,error_set_e_0 )
  
}




#################################
# Figure 2(a)
par(mfrow =c(1,1))
quants = NULL
for(i in 1:k1_u){
  quants = c(quants, quantile(error_sets_1[,i],probs = 0.95))
}
boxplot( error_sets_1,boxwex = 0.1, staplewex = 4, ylab = "errors", xlab = "k",cex.axis=1.5,cex.lab =  1.5)
points(1:  k1_u, quants,col = "blue",cex = 1.5,pch = 18)
abline(h = 0.05, col = "red", lty = 2, lwd = 2)

#################################
# Figure 2(b)
quants = NULL
for(i in 1:k1_u){
  quants = c(quants, quantile(error_sets_23[,i],probs = 0.95))
}


boxplot( error_sets_23,boxwex = 0.1, staplewex = 4, ylab = "errors",xlab ="k" ,cex.axis=1.5,cex.lab =  1.5)
points(1:  k1_u, quants,col = "blue",cex = 1.5,pch = 18)
abline(h = 0.05, col = "red", lty = 2, lwd = 2)



#################################
# Figure 2(c)
boxplot( error_sets_star,boxwex = 0.1, staplewex = 4 , ylab = "errors",xlab ="k",cex.axis=1.5,cex.lab =  1.5)






