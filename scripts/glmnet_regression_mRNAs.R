# Glmnet regression
library(pls)
library(glm2)
library(glmnet)
library(ROCR)
library(Rcpp)
library(ggplot2)
library(corrplot)

# regression model for mRNAs
file_mrna = 'mRNA_allfeatures'

t <- read.delim (file_mrna, sep = '\t')
t0 = t
# eliminate columns you are not interested in
t = t[, -c (1,5,6,7,8,25)]
# build feature matrix
z= as.matrix (t[, -1])
xst = stdize (z, center = TRUE, scale =TRUE) #i try not to standardize
xst = data.frame(xst) # this has changed to non-standardized matrix
xm = cbind (xst, t$halftime)
names(xm)[19 ]= 'halftime'

#Randomly shuffle the data
data<-xm[sample(nrow(xm)),]
#Create 10 equally size folds
folds <- cut(seq(1,nrow(data)),breaks=10,labels=FALSE)
#Perform 10 fold cross validation

#create empty data frame to accomodate all 10 models: 
# one for coefficients and one vector for mean square error
df_coeff <- data.frame(matrix(ncol = 18, nrow = 0))
n = names(xm[,-19])
colnames(df_coeff)=n
R=c() # initialization vector of correlations
l=list()

for(i in 1:10){
  #segement your data by fold using the which() function 
  testIndexes <- which(folds==i,arr.ind=TRUE)
  testData <- data[testIndexes, ]
  trainData <- data[-testIndexes, ]
  #Use test and train data partitions however you desire...
  
  y=as.integer(trainData$halftime) # regression
  x=trainData[,-19]
  #x=trainData[,-1]
  cvfit=cv.glmnet(as.matrix(x),y, alpha=0.5, type.measure="mse")
  plot(cvfit)
  lam=0.5
  #lam=cvfit$lambda.min
  #cvfit$lambda.1se; #print(cvfit$lambda.1se); #lam=cvfit$lambda.1se
  coef(cvfit,s=lam)
  fit=glmnet(x,y,alpha=0.5, lambda=lam)
  print(coef(fit))
  #ii <- which(cvfit$lambda == lam); #mse.min <- cvfit$cvm[ii]
  
  # prediction and correlation calculation
  c = as.vector(coef(fit))[2:19]
  df_coeff=rbind(df_coeff,c)
  names(df_coeff)=n
  # predict on test set
  y_pred=predict(fit, newx=as.matrix(testData[,-19]), type="response", s=lam)
  #print(y_pred)
  y_test=testData$halftime
  #print(y_test)
  # compute mis-classification error
  cc = data.frame(cbind(y_pred,y_test))
  l[[i]] = cc
  corr= cor(y_pred,y_test)
  print(corr)
  R = append(R, corr)
  
}  

print(mean(R))
print(max(R))
print(min(R))

# plotting coefficients
# produce plots of coefficients
# rank coefficients y average value
df_coeff=df_coeff[,order(apply(df_coeff,2,mean))]
# now nice boxplot on the columns
pdf(file="boxplot_coeff_glmnet_regression_mRNAs.pdf", width=7, height = 7)
par(mar=c(11,5,5,2))
boxplot(df_coeff[,1], df_coeff[,2], df_coeff[,3], df_coeff[,4], df_coeff[,5], df_coeff[,6], df_coeff[,7], df_coeff[,8], df_coeff[,9], df_coeff[,10], df_coeff[,11], df_coeff[,12], df_coeff[,13], df_coeff[,14], df_coeff[,15], df_coeff[,16], df_coeff[,17], df_coeff[,18], df_coeff[,19], df_coeff[,20], df_coeff[,21], df_coeff[,22], main="Coefficients mRNA dynamics", names = names(df_coeff), las=2, col="red", border="black", notch=TRUE)
dev.off()

# plot coefficients with ggplot
h = rep(c("Beta"), dim(df_coeff)[2]*10)
i = rep(c(names(df_coeff)),10)
j=c(); j= append(j,as.vector(t(df_coeff[1,])));j= append(j,as.vector(t(df_coeff[2,])))
j= append(j,as.vector(t(df_coeff[3,]))); j= append(j,as.vector(t(df_coeff[4,])))
j= append(j,as.vector(t(df_coeff[5,]))); j= append(j,as.vector(t(df_coeff[6,])))
j= append(j,as.vector(t(df_coeff[7,]))); j= append(j,as.vector(t(df_coeff[8,])))
j= append(j,as.vector(t(df_coeff[9,]))); j= append(j,as.vector(t(df_coeff[10,])))
df_coeff_new = data.frame(h,i,j)

colnames(df_coeff_new)=c("b", "name","value")
pdf(file="boxplot_coeff_glmnet_regression_mRNAs_ggplot2.pdf", width=9, height = 8)
par(mar=c(11,5,5,2))
ggplot(df_coeff_new, aes(x=reorder(name, value, FUN=median),y=value))+
  geom_boxplot(fill="turquoise2", colour="black", notch = TRUE, lwd=1.0)+
  scale_y_continuous(name = "Coefficients' values")+
  scale_x_discrete(name="Features")+theme_bw()+
  theme(axis.text.y = element_text(size = 17, face="bold"), axis.text.x = element_text(angle = 90, size =17, vjust = 0.5, hjust=1, face="bold"), plot.title=element_text(size=20, face="bold"), axis.title=element_text(size=17, face="bold"))+
  labs(title="Regularized regression mRNAs")
dev.off()

# Scatter plot measured versus predicted y
ii <- which(R == max(R))
best_prediction_cc <- l[[ii]]

pdf(file="scatterplot_glmnet_regression_mRNAs_ggplot2.pdf", width=9, height = 7)
ggplot(best_prediction_cc, aes(x=best_prediction_cc[,1],y=best_prediction_cc[,2])) +
  geom_point(size=3, col='turquoise2') +
  labs(x= "predicted halftime", y="computed halftime")+
  geom_smooth(method=lm, col="grey30", lwd=1.1)+
  theme_bw()+
  theme(axis.text.y = element_text(size = 20), axis.text.x = element_text(size =20), plot.title=element_text(size=22, face="bold"), axis.title=element_text(size=17, face="bold"))+
  labs(title="Regularized linear regression mRNAs")
dev.off()

print(cor(l[[ii]][,1],l[[ii]][,2]))