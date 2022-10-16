# running glmnet
# load necessary libraries
library(corrplot)
library(pls)
library(glm2)
library(glmnet)
library(ROCR)
library(Rcpp)
library(ggplot2)

file_mrna = 'mRNA_allfeatures'
t <- read.delim (file_mrna, sep = '\t')
t0 = t
t = t[, -c (1,2,5,6,7,8)]
M = M = cor(t[,-c(19)])
par(mar=c(5,5,18,5))
pdf(file="correlation_matrix_mRNAs.pdf", width=6, height = 7)
# save to a pdf file with proper margins
corrplot::corrplot(M, main="Feature correlation mRNAs", order = 'hclust', tl.cex=0.9, number.cex=1.4, tl.col='black', mar = c(2, 2, 4, 2)) # default circles
dev.off()

# write data to temporary file without column 1 and 2 and for only fast and slow
t = t[t$id =='fast' | t$id == 'slow', ]
write.table(t, 'tmp', sep = '\t', quote=F, row.names=F)

# without termination index feature
z= as.matrix (t[, -19])
xst = stdize (z, center = TRUE, scale =TRUE) #i try not to standardize
xst = data.frame(xst) # this has changed to non-standardized matrix
xm = cbind (xst, t$id)
names(xm)[19 ]= 'id'
xm[xm=="slow"]<-1
xm[xm=="fast"]<-0
#xm$id <- as.factor(xm$id)  # new line of code, necessary!!

set.seed(1234)
#xm is the standardized matrix

#Randomly shuffle the data
data<-xm[sample(nrow(xm)),]
#Create 10 equally size folds
folds <- cut(seq(1,nrow(data)),breaks=10,labels=FALSE)
#Perform 10 fold cross validation

#create empty data frame to accomodate all 10 models: 
# one for coefficients and one vector for misclassification error

df_coeff = data.frame(matrix(ncol = 18, nrow = 0))
n = names(xm[,-19])
colnames(df_coeff)=n
mce=c() # initialization vector miscl. error

for(i in 1:10){
  #segement your data by fold using the which() function 
  testIndexes <- which(folds==i,arr.ind=TRUE)
  testData <- data[testIndexes, ]
  trainData <- data[-testIndexes, ]
  #Use test and train data partitions however you desire...
  
  y=as.integer(trainData$id)
  x=trainData[,-19]
  # nested corss-validation on training data to determine optimal regularization
  # parameter lambda
  cvfit=cv.glmnet(as.matrix(x),y, alpha=0.5, family="binomial", type.measure="class")
  plot(cvfit)
  cvfit$lambda.1se
  print(cvfit$lambda.1se)
  coef(cvfit,s="lambda.1se")
  fit=glmnet(x,y,alpha=0.5, family="binomial", lambda=cvfit$lambda.1se)
  print(coef(fit))
  # update coefficients in data frame
  c = as.vector(coef(fit))[2:19]
  df_coeff=rbind(df_coeff,c)
  names(df_coeff)=n
  # predict on test set
  y_pred=as.integer(predict(fit, newx=as.matrix(testData[,-19]), type="class", s="lambda.sd1"))
  print(y_pred)
  y_test=as.integer(testData$id)
  print(y_test)
  # compute mis-classification error
  cc = cbind(y_pred,y_test)
  err= (dim(cc[cc[,1]!=cc[,2],])[1] / dim(cc)[1])*100
  print(err)
  mce = append(err, mce)
}


print(mean(mce))
print(max(mce))
print(min(mce))

# produce plots of coefficients
# rank coefficients y average value
df_coeff=df_coeff[,order(apply(df_coeff,2,mean))]
# now nice boxplot on the columns
pdf(file="boxplot_coeff_glmnet_mrnas.pdf", width=7, height = 7)
par(mar=c(11,5,5,2))
boxplot(df_coeff[,1], df_coeff[,2], df_coeff[,3], df_coeff[,4], df_coeff[,5], df_coeff[,6], df_coeff[,7], df_coeff[,8], df_coeff[,9], df_coeff[,10], df_coeff[,11], df_coeff[,12], df_coeff[,13], df_coeff[,14], df_coeff[,15], df_coeff[,16], df_coeff[,17], df_coeff[,18], main="Coefficients mRNA dynamics", names = names(df_coeff), las=2, col="red", border="black", notch=TRUE)
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
pdf(file="boxplot_coeff_glmnet_mRNAs_ggplot2.pdf", width=9, height = 8)
par(mar=c(11,5,5,2))
ggplot(df_coeff_new, aes(x=reorder(name, value, FUN=median),y=value))+
  geom_boxplot(fill="turquoise2", colour="black", notch = TRUE, lwd=1.0)+
  scale_y_continuous(name = "Coefficients' values")+
  scale_x_discrete(name="Features")+theme_bw()+
  theme(axis.text.y = element_text(size = 17, face="bold"), axis.text.x = element_text(angle = 90, size =17, vjust = 0.5, hjust=1, face="bold"), plot.title=element_text(size=20, face="bold"), axis.title=element_text(size=17, face="bold"))+
  labs(title="Regularized logistic regression mRNAs")
dev.off()