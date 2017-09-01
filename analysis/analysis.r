# Analyze the out put files
#v2.0
library("ggplot2")
library("ggthemes")
library("EbayesThresh")
library(data.table)
library(dplyr)
library(MESS)
library(sfsmisc)

# todo KDE
# todo -sense
setwd("/Users/roujia/Documents/02_dev/07_prob_model/output/test_l2_without_chromatin/")
output_dir <- "/Users/roujia/Documents/02_dev/07_prob_model/output/test_l2_without_chromatin/plots/02_all-ROC/"
#setwd("/home/ryogali/output/test_l2_without_chromatin/")
# todo move the time point file

time_point_file <- "/Users/roujia/Documents/02_dev/07_prob_model/data/vaf.txt"
# open timepoint file
time_point = read.delim2(time_point_file,sep="\t",header=FALSE,fill=TRUE,stringsAsFactors = F)
# get list of test and random file
test_prob_files <- list.files(path="./test_prob/", pattern="*.txt", full.names=T, recursive=FALSE)
train_prob_files <- list.files(path="./train_prob/", pattern="*.txt", full.names=T, recursive=FALSE)
low_prob_files <- list.files(path="./lowsup_prob/", pattern="*.txt", full.names=T, recursive=FALSE)
random_prob_files <- list.files(path="./random_prob/", pattern="*.txt", full.names=T, recursive=FALSE)


# sort the files by name
test_prob_files <- test_prob_files[order(test_prob_files)]
train_prob_files <- train_prob_files[order(train_prob_files)]
low_prob_files <- low_prob_files[order(low_prob_files)]
random_prob_files <- random_prob_files[order(random_prob_files)]


calculate_ratio <- function(vector, thresh){
    return(length(vector[vector>thresh])/length(vector))
}

k_means <- function(data){
    k <- 1
    means <- kmeans(data$prob, k)
    assign_cluster <- data.frame(data$VAF, data$prob, means$cluster)
    class_one <- assign_cluster[which(assign_cluster$means.cluster==1),]
    class_two <- assign_cluster[which(assign_cluster$means.cluster==2),]
    threshold <- (min(class_two$data.prob)+max(class_one$data.prob)) / 2

    return(threshold)
}


plots <- function(x_col, y_col, x_lim, y_lim, title, x_lab, y_lab){

    plot(random[[x_col]],random[[y_col]],xlim=x_lim,pch=18, axes = TRUE, xlab = x_lab, ylab = y_lab, lwd=2,ylim=y_lim, col=alpha("#3366CC", .8))
    # plot(random[[x_col]],random[[y_col]],xlim=x_lim,pch=18, axes = TRUE, xlab = x_lab, ylab = y_lab, xaxt = "n", col=alpha("blue", .8),lwd=2,ylim=y_lim)
    par(new= TRUE)
    plot(test[[x_col]],test[[y_col]],xlim=x_lim,main=title,axes = FALSE,pch=18, xlab="", ylab="", col="red", lwd=2, ylim=y_lim)
    par(new= TRUE)
    plot(train[[x_col]],train[[y_col]],xlim=x_lim,pch=18, axes = FALSE, xlab = "", ylab = "", col=alpha("darkgrey", .8),lwd=2,ylim=y_lim)
    par(new= TRUE)
    plot(low_sup[[x_col]],low_sup[[y_col]],xlim=x_lim,pch=18, axes = FALSE, xlab = "", ylab = "", col=alpha("darkgreen", .8),lwd=2,ylim=y_lim)
    par(new=TRUE)
    abline(v=time_points, col="darkgrey",lty=2)

}

plot_kde <- function(test_d, random_d, train_d, low_sup_d, x_lim, y_lim, title, x_lab, y_lab){
    plot(random_d,xlim=x_lim, main="",axes = TRUE, xlab = x_lab, ylab = y_lab, lwd=2,ylim=y_lim, col=alpha("#3366CC", .8))
    par(new= TRUE)
    plot(test_d,xlim=x_lim,main=title,axes = FALSE,pch=18, xlab="", ylab="", col="red", lwd=2, ylim=y_lim)
    par(new= TRUE)
    plot(train_d,xlim=x_lim,pch=18,main="", axes = FALSE, xlab = "", ylab = "", col=alpha("darkgrey", .8),lwd=2,ylim=y_lim)
    par(new= TRUE)
    plot(low_sup_d,xlim=x_lim,pch=18,main="", axes = FALSE, xlab = "", ylab = "", col=alpha("darkgreen", .8),lwd=2,ylim=y_lim)
}

calculate_roc <- function(neg_data, pos_data){

    p <- nrow(pos_data) # total positive
    n <- nrow(neg_data) # total negative

    total <- seq(0,max(c(neg_data$prob, pos_data$prob)), 0.00001)

    tp_vector <- numeric()
    fp_vector <- numeric()
    sum <- numeric()
    for (prob in total) {
        tp <- length(pos_data$prob[which(pos_data$prob >= prob)])
        tpr <- tp/p
        tp_vector <- c(tp_vector, tpr)

        fp <- length(neg_data$prob[which(neg_data$prob >= prob)])
        fpr <- fp/n
        fp_vector <- c(fp_vector, fpr)

        #find Youden-Index
        sum <- c(sum, tpr+(1-fpr))
    }
    output <- data.frame(cbind(tp_vector, fp_vector, total, sum))
    return (output)
}


train_ratio <- vector()
test_ratio <- vector()
lowsup_ratio <- vector()
random_ratio <- vector()
tumor_ids <- vector()

# number of files
total <- length(test_prob_files)
# total <- 30
count <- 0
# thresh <- 0.005

i <- 1
while (i <= total){
    tumour_name <- unlist(strsplit(basename(test_prob_files[i]), "\\."))[1]

    time_points <- time_point[time_point$V1==tumour_name,]
	filename <- paste(output_dir, tumour_name, ".pdf", sep = "")
    print(tumour_name)
	pdf(filename, width = 10, height = 15)

    par(mfrow = c(3,2))
    test_avg <- 0
    train_avg <- 0
    random_avg <- 0
    lowsup_avg <- 0
    auc <- 0
	for (j in i:(i+2)){

	    test <- read.table(test_prob_files[j], header=TRUE)
	    train <- read.table(train_prob_files[j], header=TRUE)
        low_sup <- read.table(low_prob_files[j], header=TRUE)
        random <- read.table(random_prob_files[j], header=TRUE)

        test <- test[order(test$VAF),]
        train <- train[order(train$VAF),]
        low_sup <- low_sup[order(low_sup$VAF),]
        random <- random[order(random$VAF),]


        # get density
        test_d <- density(test$prob)
        train_d <- density(train$prob)
        low_sup_d <- density(low_sup$prob)
        random_d <- density(random$prob)
#        plot_kde(test_d, random_d, train_d, low_sup_d, c(0,0.03), c(0,2000), "KDE", "probability", "density")
#        # use k-means to get threshold
#        # threshold <- (k_means(train) + k_means(test) + k_means(low_sup) + k_means(random))/4
#

#
#        ######## remove sense strand #########
#        ######## update new prob #########
#        test$prob <- log(test$p_ai) * log(test$p_ti) * log(test$p_ce) * log(test$p_si)
#        print(min(test$prob))
#        random$prob <- log(random$p_ai) * log(random$p_ti) * log(random$p_ce) * log(random$p_si)
#        train$prob <- log(train$p_ai) * log(train$p_ti) * log(train$p_ce) * log(train$p_si)
#        low_sup$prob <- log(low_sup$p_ai) * log(low_sup$p_ti) * log(low_sup$p_ce) * log(low_sup$p_si)
#   #
#        test$prob_log <- log(test$prob)
#        #print(max(test$prob))
#        random$prob_log <- log(random$prob)
#        train$prob_log <- log(train$prob)
#        low_sup$prob_log <- log(low_sup$prob)
        x_lim = c(0,2)
        y_lim = c(0,0.3)
#        # mean and std
        plots("VAF", "prob", x_lim, y_lim, "Predicted Probabilities","Avg number of mutant alleles per cancer cell", "Probabilities")
        out <- calculate_roc(random, test)
        max <- max(out$sum)
        thresh <- out$total[which(out$sum==max)]
        x <- out$fp_vector[which(out$sum == max)]
        y <- out$tp_vector[which(out$sum == max)]
        auc_fold <- integrate.xy(out$fp_vector,out$tp_vector)

        abline(h=log(thresh), col = "black", lty=2.5)

        # test_pos_ratio <- paste(toString(nrow(test)), toString(calculate_ratio(test$prob, thresh)))
        # train_pos_ratio <- paste(toString(nrow(train)), toString(calculate_ratio(train$prob, thresh)))
        # lowsup_pos_ratio <- paste(toString(nrow(low_sup)), toString(calculate_ratio(low_sup$prob, thresh)))
        # random_pos_ratio <- paste(toString(nrow(random)), toString(calculate_ratio(random$prob, thresh)))
        legend("topleft", # places a legend at the appropriate place
			c(paste("Test data: ", toString(nrow(test))),
			  paste("Low Support data: ",toString(nrow(low_sup))),
			  paste("Train data: ", toString(nrow(train))),
			  paste("Random data: ",toString(nrow(random)))), # puts text in the legend

			pch=c(18,18,18,18), # gives the legend appropriate symbols (lines)
            #
			# lwd=c(2.5,2.5,2.5,2.5),
			pt.cex=c(2,2,2,2),
			col=c("red", "darkgreen","grey", "#3366CC"))
        legend("topright", # places a legend at the appropriate place
             sprintf("%s",c(toString(calculate_ratio(test$prob, thresh)),
                            toString(calculate_ratio(train$prob, thresh)),
                            toString(calculate_ratio(low_sup$prob, thresh)),
                            toString(calculate_ratio(random$prob, thresh)))), # puts text in the legend
             pch=c(18,18,18,18), # gives the legend appropriate symbols (lines)
             pt.cex=c(2,2,2,2),
             col=c("red", "darkgreen","grey", "#3366CC"))
##
###
        plot(out$fp_vector, out$tp_vector, type="l",xlim=c(0,1), ylim=c(0,1),lwd=2.5, main="ROC - test(pos); random(neg)", xlab="FPR", ylab="TPR", col="orange")
        par(new=TRUE)
        plot(x, y, pch=19, xlim=c(0,1), ylim=c(0,1),col="black", axes=FALSE, xlab="", ylab="")
        par(new=TRUE)
        plot(c(0,1), c(0,1), type="l", lty=5, col="black", axes=FALSE, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="")


        legend("bottomright",
                sprintf("%s", c(paste("AUC: ", toString(auc_fold)))),
                col=c("Black"))

#         plots("VAF", "p_ce", x_lim, y_lim, "p_ce", "Avg number of mutant alleles per cancer cell", "Probabilities")
#         plots("VAF", "p_ai", x_lim, c(0,1), "p_ai", "Avg number of mutant alleles per cancer cell", "Probabilities")
#         plots("VAF", "p_ti", x_lim, c(0,1), "p_ti", "Avg number of mutant alleles per cancer cell", "Probabilities")
#         plots("VAF", "p_si", x_lim, c(0,1), "p_si", "Avg number of mutant alleles per cancer cell", "Probabilities")
#         plots("mut", "prob", c(0,96), y_lim, "prob")
#         axis(1, at=0:95, labels=tri_file$tri, las = 2, cex.axis = 0.8)

        test_avg <- test_avg + calculate_ratio(test$prob, thresh)
        train_avg <- train_avg + calculate_ratio(train$prob, thresh)
        lowsup_avg <- lowsup_avg + calculate_ratio(low_sup$prob, thresh)
        random_avg <- random_avg + calculate_ratio(random$prob, thresh)
#        auc <- auc + auc_fold

    }
    tumor_ids <- append(tumor_ids, tumour_name)
    test_ratio <- append(test_ratio, test_avg/3)
    train_ratio <- append(train_ratio, train_avg/3)
    lowsup_ratio <- append(lowsup_ratio, lowsup_avg/3)
    random_ratio <- append(random_ratio, random_avg/3)
    dev.off()
#    print(auc/3)

    i <- i + 3

}
print(test_ratio)
print(train_ratio)
print(lowsup_ratio)
print(random_ratio)
#
#filename <- paste(output_dir, "ratios", ".pdf", sep = "")
#pdf(filename, width = 22, height = 15)
#total <- total/3
#par(mar=c(16, 8, 7, 7) + 0.1)
#
#plot(1:total, test_ratio ,cex=2,pch=18, xaxt="n", xlab="", xaxt = "n",axes = TRUE,main="Ratios", ylab = "positive ratio", col=alpha("red", 1),lwd=2,ylim=c(0,1))
#par(new=TRUE)
#plot(1:total, train_ratio ,cex=2,pch=18, axes = FALSE, xlab = "", ylab = "", col=alpha("darkgrey", .8),lwd=2,ylim=c(0,1))
#par(new=TRUE)
#plot(1:total, lowsup_ratio ,cex=2,pch=18, axes = FALSE, xlab = "", ylab = "", col=alpha("darkgreen", .8),lwd=2,ylim=c(0,1))
#par(new=TRUE)
#plot(1:total, random_ratio ,cex=2,pch=18, axes = FALSE, xlab = "", ylab = "", col=alpha("#3366CC", .8),lwd=2,ylim=c(0,1))
#par(new=TRUE)
#abline(v=1:total, col = "lightgrey", lty=2)
#legend("topleft", # places a legend at the appropriate place
#        c("Test data",
#          "Low Support data",
#          "Train data",
#          "Random data"), # puts text in the legend
#
#        pch=c(18,18,18,18), # gives the legend appropriate symbols (lines)
#        #
#        # lwd=c(2.5,2.5,2.5,2.5),
#        pt.cex=c(2,2,2,2),
#        col=c("red", "darkgreen","grey", "#3366CC"))
#
#axis(1, at=1:total, labels=1:total)
#dev.off()
#
## plot as histogram
#filename <- paste(output_dir,"ratios_his", ".pdf", sep = "")
#pdf(filename, width = 10, height = 20)
## total <- total/3
#par(mar=c(8, 4, 3.5, 3.5) + 0.1, mfrow = c(4, 1))
#
#plot(1:total, test_ratio ,cex=2,pch=18,axes = TRUE,main="Ratios above threshold - Test data", ylab = "positive ratio", col=alpha("red", .8),lwd=2,ylim=c(0,1))
#par(new=TRUE)
#abline(h=mean(test_ratio), col = "lightgrey", lty=3)
#
#plot(1:total, train_ratio ,cex=2,pch=18, axes = TRUE, main="Ratios above threshold - Train data", ylab = "positive ratio", col=alpha("darkgrey", .8),lwd=2,ylim=c(0,1))
#par(new=TRUE)
#abline(h=mean(train_ratio), col = "lightgrey", lty=3)
#
#plot(1:total, lowsup_ratio ,cex=2,pch=18, axes = TRUE, main="Ratios above threshold - Low support data", ylab = "positive ratio", col=alpha("darkgreen", .8),lwd=2,ylim=c(0,1))
#par(new=TRUE)
#abline(h=mean(lowsup_ratio), col = "lightgrey", lty=3)
#
#plot(1:total, random_ratio ,cex=2,pch=18, axes = TRUE, main="Ratios above threshold - Random data", ylab = "positive ratio", col=alpha("#3366CC", .8),lwd=2,ylim=c(0,1))
#par(new=TRUE)
#abline(h=mean(random_ratio), col = "lightgrey", lty=3)
#dev.off()
#
## plot a heatmap
#filename <- paste(output_dir,"ratios_heat", ".pdf", sep = "")
#pdf(filename, width = 30, height = 10)
#df <- matrix(c(train_ratio, test_ratio, lowsup_ratio, random_ratio),nrow=length(train_ratio))
#df.m <- data.frame(tumor_ids, "Train" = train_ratio, "Test" = test_ratio, "Low support" = lowsup_ratio, "Random" = random_ratio)
#df.m <- melt(df.m)
#head(df.m)
#ggplot(df.m, aes(x=tumor_ids, y=variable, fill=value)) +
#        geom_tile(color = "white") +
#        scale_fill_gradient(low = "white", high = "darkblue", limit = c(0,1), name="Ratios") +
#        xlab("Tumors") +
#        theme(axis.text.x = element_blank())
#
## heatmap(df, Rowv=F, Colv=F, main="ratios above threshold")
dev.off()
