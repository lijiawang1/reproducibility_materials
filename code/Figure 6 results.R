##########################
#generate the half violin plots in Figure 6 in Section 3.2
library(ggplot2)
library(gridExtra)


#get path to the root directory
path1 = getwd()
check_path = unlist(strsplit(path1, split = "/"))
while (check_path[length(check_path)] != "reproducibility_materials"){
  path1 = dirname(path1)  
  check_path = unlist(strsplit(path1, split = "/"))
}

setwd(path1) #the path of the folder

#################################
#set up the featurization method that you want to investigate
feature_method = "m1"   # "mi" stands for Featurization method M.i proposed in section 3.1. e.g., m1,m2,m3,m4.


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


############################################




load(paste("output/logistic_",feature_method,".Rdata", sep = ""))
lss = dim(list_error$non_np)[1]
plot_error = c(list_error$non_np)
quant = NULL
for ( i in 2:12 ){
  error_i = NULL
  for(j in 1:lss){
    min_j = max(which(list_error[[13]][j,] == min(list_error[[13]][j,] )))
    error_i = c(error_i,list_error[[i]][j,min_j ])
  }
  plot_error = c( plot_error,error_i)
  quant = c(quant,quantile(error_i ,probs = c(0.8)))
}

errornames = rep(c("overall","error1","error12","error13","error2","error21","error23","error3","error31","error32","error*"),each = dim(list_error$non_np)[1])
errornames = rep(errornames,2)
method = rep(c("non-NP","NP"),each = length(c(list_error$non_np)))

plot_data = data.frame(errors =  plot_error , type =  errornames,method = method)


q21 = quant[2]
q71 = quant[7]



color_value = c("#E69F00", "#56B4E9")
plot_data_4 = plot_data[which(plot_data$type %in% c("error1","error23","error21","error31","error32","overall")),]
plot_data_4$type = factor(plot_data_4$type, level = c("error1","error23","error21","error31","error32","overall"))

p1 <-ggplot(data = plot_data_4, aes( type , errors, fill = method)) + geom_split_violin(trim = TRUE)   + ylim(0,0.8)+
  geom_segment(aes(x = 0.5, y = 0.2, xend = 1.5, yend = 0.2),linetype="dashed",col = c("red"))+
  geom_segment(aes(x = 1.5, y = 0.2, xend = 2.5, yend = 0.2),linetype="dashed",col = c("red")) +
  geom_segment(aes(x = 1, y = q21, xend = 1.2, yend = q21)) +
  geom_segment(aes(x = 2, y = q71, xend = 2.2, yend = q71))+
  geom_vline(xintercept = 2.5, size=1.5 , linetype="dotted")+
  labs(x = "",title = "Logistic Regression")+
  scale_colour_manual(values = color_value) +  scale_fill_manual(values = color_value) +
  theme(text = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 15)) + theme(legend.position="none")


load(paste("output/forest_",feature_method,".Rdata", sep = ""))
lss = dim(list_error$non_np)[1]
plot_error = c(list_error$non_np)
quant = NULL
for ( i in 2:12 ){
  error_i = NULL
  for(j in 1:50){
    min_j = max(which(list_error[[13]][j,] == min(list_error[[13]][j,] )))
    error_i = c(error_i,list_error[[i]][j,min_j ])
    
  }
  plot_error = c( plot_error,error_i)
  quant = c(quant,quantile(error_i ,probs = c(0.8)))
}


errornames = rep(c("overall","error1","error12","error13","error2","error21","error23","error3","error31","error32","error*"),each = dim(list_error$non_np)[1])
errornames = rep(errornames,2)
method = rep(c("original","NP-adjusted"),each = length(c(list_error$non_np)))

plot_data = data.frame(errors =  plot_error , type =  errornames,method = method)



q22 = quant[2]
q72 = quant[7]


color_value = c("#E69F00", "#56B4E9")
plot_data_4 = plot_data[which(plot_data$type %in% c("error1","error23","error21","error31","error32","overall")),]
plot_data_4$type = factor(plot_data_4$type, level = c("error1","error23","error21","error31","error32","overall"))
plot_data_4$method = factor(plot_data_4$method, level =c("original","NP-adjusted"))
p2 <-ggplot(data = plot_data_4, aes( type , errors, fill = method)) + geom_split_violin(trim = TRUE)  + ylim(0,0.8)+
  geom_segment(aes(x = 0.5, y = 0.2, xend = 1.5, yend = 0.2),linetype="dashed",col = c("red"))+
  geom_segment(aes(x = 1.5, y = 0.2, xend = 2.5, yend = 0.2),linetype="dashed",col = c("red")) +
  geom_segment(aes(x = 1, y = q22, xend = 1.2, yend = q22)) +
  geom_segment(aes(x = 2, y = q72, xend = 2.2, yend = q72))+
  geom_vline(xintercept = 2.5, size=1.5 , linetype="dotted")+
  labs(x = "",title = "Random Forest")+
  scale_colour_manual(values = color_value) +  scale_fill_manual(values = color_value) +
  theme(text = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 15)) + theme(legend.position="none")

p2 



load(paste("output/svm_",feature_method,".Rdata", sep = ""))
lss = dim(list_error$non_np)[1]
plot_error = c(list_error$non_np)
quant = NULL
for ( i in 2:12 ){
  error_i = NULL
  for(j in 1:50){
    min_j = max(which(list_error[[13]][j,] == min(list_error[[13]][j,] )))
    error_i = c(error_i,list_error[[i]][j,min_j ])
    
  }
  plot_error = c( plot_error,error_i)
  quant = c(quant,quantile(error_i ,probs = c(0.8)))
}

errornames = rep(c("overall","error1","error12","error13","error2","error21","error23","error3","error31","error32","error*"),each = dim(list_error$non_np)[1])
errornames = rep(errornames,2)
method = rep(c("non-NP","NP"),each = length(c(list_error$non_np)))

plot_data = data.frame(errors =  plot_error , type =  errornames,method = method)





q23 = quant[2]
q73 = quant[7]
color_value = c("#E69F00", "#56B4E9")
plot_data_4 = plot_data[which(plot_data$type %in% c("error1","error23","error21","error31","error32","overall")),]
plot_data_4$type = factor(plot_data_4$type, level = c("error1","error23","error21","error31","error32","overall"))

p3 <-ggplot(data = plot_data_4, aes( type , errors, fill = method)) + geom_split_violin(trim = TRUE)  + ylim(0,0.8)+
  geom_segment(aes(x = 0.5, y = 0.2, xend = 1.5, yend = 0.2),linetype="dashed",col = c("red"))+
  geom_segment(aes(x = 1.5, y = 0.2, xend = 2.5, yend = 0.2),linetype="dashed",col = c("red")) +
  geom_segment(aes(x = 1, y = q23, xend = 1.2, yend = q23)) +
  geom_segment(aes(x = 2, y = q73, xend = 2.2, yend =  q73))+
  labs(x = "",title = "SVM")+
  geom_vline(xintercept = 2.5, size=1.5 , linetype="dotted")+
  scale_colour_manual(values = color_value) +  scale_fill_manual(values = color_value) +
  theme(text = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 15)) + theme(legend.position="none")



grid.arrange(p1,p2,p3) #generate the half violin plots for the corresponding featurization method in Figure 6 in section 3.2


