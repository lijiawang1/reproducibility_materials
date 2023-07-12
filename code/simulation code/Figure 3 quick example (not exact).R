########################
#generating SIMPLE EXAMPLE similar to the results in Figure 3
library(nnet)
library(randomForest)
library(e1071)
library(ggplot2)
library(gridExtra)

##################
#functions for producing half violin plots
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           })

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}
###############################



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

alpha_1 = 0.05
alpha_2 = 0.05
delta_1 = 0.05
delta_2 = 0.05
n1 = 500 # class 1
n2 = 500# class 2
n3 = 500# class 3
n_all = n1 + n2 + n3

pi_1 = n1 /n_all
pi_2 = n2 /n_all
pi_3 = n3 /n_all
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

for(ii in 1:20){
  
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





par(mfrow =c(1,1))
table_record =   matrix( round(c(colMeans(error_sets_23_0)[k1_u],colMeans(error_sets_23)[k1_u],
                                 colMeans(error_sets_32_0)[k1_u],colMeans(error_sets_32)[k1_u]), digits = 3),ncol = 2,byrow = FALSE)

rownames(table_record) = c("Prop 1","Thm 1")
colnames(table_record) = c("Error23","Error32")

table_1 = table_record



plot_error= rbind(error_sets_23,error_sets_32,error_sets_23_0,error_sets_32_0)[,k1_u]
errornames = rep(c("error23","error32"),each = dim(error_non)[1])
errornames = rep(errornames,2)
method = rep(c("Thm 1","Prop 1"),each = length(plot_error)/2)
plot_data = data.frame(errors =  plot_error , type =  errornames,method = method)



color_value = c("pink2", "#56B4E9")

plot_data$method = factor(plot_data$method, level = c("Prop 1","Thm 1"))

p4_1 <-ggplot(data = plot_data, aes( type , errors, fill = method)) + geom_split_violin(trim = TRUE)   +
  labs(title = "Logistic Regression", x = "")+
  geom_vline(xintercept = 1.5, size=1.5 , linetype="dotted")+
  scale_colour_manual(values = color_value) +  scale_fill_manual(values = color_value) + ylim(0,0.2)+
  theme(text = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 15)) + theme(legend.position="none")

p4_1
###################################



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

for(ii in 1:20){
  
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
  
  model <-  randomForest(as.factor(Y)~., data=set_s , proximity=TRUE)
  
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
  
  prob_set = predict(model, newdata = set_t, "prob")
  prob_1 =  prob_set[which( set_t$Y == 1), ]
  prob_2 =  prob_set[which( set_t$Y == 2), ]
  prob_e = predict(model, newdata = set_e, "prob")
  prob_c = predict(model, newdata = set_c, "prob")
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



table_record =   matrix( round(c(colMeans(error_sets_23_0)[k1_u],colMeans(error_sets_23)[k1_u],
                                 colMeans(error_sets_32_0)[k1_u],colMeans(error_sets_32)[k1_u]), digits = 3),ncol = 2,byrow = FALSE)

rownames(table_record) = c("Prop 1","Thm 1")
colnames(table_record) = c("Error23","Error32")

table_2 = table_record




plot_error= rbind(error_sets_23,error_sets_32,error_sets_23_0,error_sets_32_0)[,k1_u]
errornames = rep(c("error23","error32"),each = dim(error_non)[1])
errornames = rep(errornames,2)
method = rep(c("Thm 1","Prop 1"),each = length(plot_error)/2)
plot_data = data.frame(errors =  plot_error , type =  errornames,method = method)



color_value = c("pink2", "#56B4E9")

plot_data$method = factor(plot_data$method, level = c("Prop 1","Thm 1"))

p4_2 <-ggplot(data = plot_data, aes( type , errors, fill = method)) + geom_split_violin(trim = TRUE)   +
  geom_vline(xintercept = 1.5, size=1.5 , linetype="dotted")+
  labs(title = "Random Forest", x = "")+
  scale_colour_manual(values = color_value) +  scale_fill_manual(values = color_value) + ylim(0,0.2)+
  theme(text = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 15)) + theme(legend.position="none")

p4_2


###################################


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



for(ii in 1:20){
  
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
  
  
  model <-   svm(as.factor(Y)~., data=set_s, kernel =
                   "linear", probability = TRUE)
  
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
  
  
  prob_set = attr(predict(  model, newdata = set_t, probability = TRUE),"probabilities")
  
  prob_1 =  prob_set[which( set_t$Y == 1), ]
  prob_2 =  prob_set[which( set_t$Y == 2), ]
  prob_e =  attr(predict(  model, newdata = set_e, probability = TRUE),"probabilities")
  prob_c = attr(predict(  model, newdata = set_c, probability = TRUE),"probabilities")
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



table_record =   matrix( round(c(colMeans(error_sets_23_0)[k1_u],colMeans(error_sets_23)[k1_u],
                                 colMeans(error_sets_32_0)[k1_u],colMeans(error_sets_32)[k1_u]), digits = 3),ncol = 2,byrow = FALSE)

rownames(table_record) = c("Prop 1","Thm 1")
colnames(table_record) = c("Error23","Error32")

table_3 = table_record



plot_error= rbind(error_sets_23,error_sets_32,error_sets_23_0,error_sets_32_0)[,k1_u]
errornames = rep(c("error23","error32"),each = dim(error_non)[1])
errornames = rep(errornames,2)
method = rep(c("Thm 1","Prop 1"),each = length(plot_error)/2)
plot_data = data.frame(errors =  plot_error , type =  errornames,method = method)



color_value = c("pink2", "#56B4E9")

plot_data$method = factor(plot_data$method, level = c("Prop 1","Thm 1"))

p4_3 <-ggplot(data = plot_data, aes( type , errors, fill = method)) + geom_split_violin(trim = TRUE)   +
  geom_vline(xintercept = 1.5, size=1.5 , linetype="dotted")+
  scale_colour_manual(values = color_value) +  scale_fill_manual(values = color_value) + ylim(0,0.2)+
  labs(title = "SVM", x = "")+
  theme(text = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 15))  + theme(legend.position="none")

p4_3




############################
# SIMPLE EXAMPLE half violin plots similar to that in Figure 3
grid.arrange(p4_1,p4_2,p4_3)


#################
#SIMPLE EXAMPLE Table in similar to that in Figure 3
#logistic
table_1
#random forest
table_2
#svm
table_3


