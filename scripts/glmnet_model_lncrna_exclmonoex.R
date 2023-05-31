# running glmnet
# load necessary libraries
library(pls)
library(glm2)
library(glmnet)
library(ROCR)
library(Rcpp)
library(ggplot2)
library(corrplot)
file_lnc = 'lnc2077_allfeatures'

t <- read.delim (file_lnc, sep = '\t')
t = t[, -c (1,2,5,6,7,8)]
M = cor(t[,-c(19)])
par(mar=c(5,5,18,5))
pdf(file="correlation_matrix_lncRNAs.pdf", width=6, height = 7)
# save to a pdf file with proper margins
corrplot::corrplot(M, main="Feature correlation lncRNAs", order = 'hclust', tl.cex=0.9, number.cex=1.4, tl.col='black', mar = c(2, 2, 4, 2)) # default circles
dev.off()

t = t[t$id =='fast' | t$id == 'slow', ]

# build feature matrix
z= as.matrix (t[, -19])
xst = stdize (z, center = TRUE, scale =TRUE) #i try not to standardize
xst = data.frame(xst) # this has changed to non-standardized matrix
xm = cbind (xst, t$id)
names(xm)[19 ]= 'id'
xm[xm=="slow"]<-1
xm[xm=="fast"]<-0

# Model training
set.seed(1234)
#Randomly shuffle the data
data<-xm[sample(nrow(xm)),]
#Create 10 equally size folds
folds <- cut(seq(1,nrow(data)),breaks=10,labels=FALSE)
#Perform 10 fold cross validation
df_coeff <- data.frame(matrix(ncol = 18, nrow = 0))
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
  y_test=as.integer(testData$id)
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

# plot coefficients with ggplot
h = rep(c("Beta"), dim(df_coeff)[2]*10)
i = rep(c(names(df_coeff)),10)
j=c(); j= append(j,as.vector(t(df_coeff[1,])));j= append(j,as.vector(t(df_coeff[2,])))
j= append(j,as.vector(t(df_coeff[3,]))); j= append(j,as.vector(t(df_coeff[4,])))
j= append(j,as.vector(t(df_coeff[5,]))); j= append(j,as.vector(t(df_coeff[6,])))
j= append(j,as.vector(t(df_coeff[7,]))); j= append(j,as.vector(t(df_coeff[8,])))
j= append(j,as.vector(t(df_coeff[9,]))); j= append(j,as.vector(t(df_coeff[10,])))
df_coeff_new = data.frame(h,i,j)
print(df_coeff_new)
save(df_coeff_new,file="lncrna_coeff_glmnet_logistic.Rda")

colnames(df_coeff_new)=c("b", "name","value")
pdf(file="boxplot_coeff_glmnet_lncRNAs_ggplot2.pdf", width=9, height = 8)
par(mar=c(11,5,5,2))
ggplot(df_coeff_new, aes(x=reorder(name, value, FUN=median),y=value))+
  geom_boxplot(fill="indianred1", colour="black", notch = FALSE, lwd=0.5)+
scale_y_continuous(name = "Coefficients' values")+
scale_x_discrete(name="Features")+theme_bw()+
theme(axis.text.y = element_text(size = 17, face="bold"), axis.text.x = element_text(angle = 90, size =17, vjust = 0.5, hjust=1, face="bold"), plot.title=element_text(size=20, face="bold"), axis.title=element_text(size=17, face="bold"))+
labs(title="Regularized logisitc regression lncRNAs")
dev.off()

###########################################################
# Let's apply the trained lncRNA model to predict mRNA dynamics
# Let's prepare the matrix for mRNAs
file_mrna = 'mRNA_allfeatures'
t <- read.delim (file_mrna, sep = '\t')
t0 = t
t = t[, -c (1,2,5,6,7,8)]
t = t[t$id =='fast' | t$id == 'slow', ]
write.table(t, 'tmp', sep = '\t', quote=F, row.names=F)
# build feature matrix
z= as.matrix (t[, -19])
xst = stdize (z, center = TRUE, scale =TRUE) #i try not to standardize
xst = data.frame(xst) # this has changed to non-standardized matrix
xm = cbind (xst, t$id)
names(xm)[19 ]= 'id'
xm[xm=="slow"]<-1
xm[xm=="fast"]<-0

# create lncRNA test set from these data
# random shuffling
data<-xm[sample(nrow(xm),100),]
testData <- data

y_pred=as.integer(predict(fit, newx=as.matrix(testData[,-19]), type="class", s="lambda.sd1"))
y_test=as.integer(testData$id)
cc = cbind(y_pred,y_test)
err= (dim(cc[cc[,1]!=cc[,2],])[1] / dim(cc)[1])*100
cc_fast = cc[cc[,1]==0,]
cc_slow = cc[cc[,1]==1,]
err_fast = (dim(cc_fast[cc_fast[,1]!=cc_fast[,2],])[1] / dim(cc_fast)[1])*100
print(err_fast)
err_slow = (dim(cc_slow[cc_slow[,1]!=cc_slow[,2],])[1] / dim(cc_slow)[1])*100
print(err_slow)
print(err)

#############################################################################
# regularized logistic regression model for lRNAs, exlcuding monoexonic transcripts
f = 'all_txt_monoexonicexcluded_oldANDnew_id'
t <- read.delim (f, sep = '\t')
t = t[t$type =='lncRNA',] # select lncRNAs only from all transcript file
t = t[, -c (1,2,5,6,7,8,26,27,28)]
t = t[t$old_id_monoexonincluded =='fast' | t$old_id_monoexonincluded == 'slow', ]

# Check SE min and SE mean distribution after exlcuding monoexonic lncRNAs
xx = t$SE_mean[round(t$exon_density*t$txt_length/1000.0) > 1 & t$old_id_monoexonincluded=="slow"]
yy = t$SE_mean[round(t$exon_density*t$txt_length/1000.0) > 1 & t$old_id_monoexonincluded=="fast"]
boxplot(xx,yy,names=c('slow', 'fast'), col=c('red', 'blue'))

xx = t$SE_min[round(t$exon_density*t$txt_length/1000.0) > 1 & t$old_id_monoexonincluded=="slow"]
yy = t$SE_min[round(t$exon_density*t$txt_length/1000.0) > 1 & t$old_id_monoexonincluded=="fast"]
boxplot(xx,yy,names=c('slow', 'fast'), col=c('red', 'blue'))

# check after setting number of exon > 2
xx = t$SE_min[round(t$exon_density*t$txt_length/1000.0) > 2 & t$old_id_monoexonincluded=="slow"]
yy = t$SE_min[round(t$exon_density*t$txt_length/1000.0) > 2 & t$old_id_monoexonincluded=="fast"]
boxplot(xx,yy,names=c('slow', 'fast'), col=c('red', 'blue'))

# build feature matrix
z= as.matrix (t[, -19])
xst = stdize (z, center = TRUE, scale =TRUE)
xst = data.frame(xst)
xm = cbind (xst, t$old_id_monoexonincluded)
names(xm)[19 ]= 'id'
xm[xm=="slow"]<-1; xm[xm=="fast"]<-0

# Model training
#set.seed(1234)
#Randomly shuffle the data
data<-xm[sample(nrow(xm)),]
#Create 10 equally size folds
folds <- cut(seq(1,nrow(data)),breaks=10,labels=FALSE)
#Perform 10 fold cross validation
#create empty data frame to accomodate all 10 models: 
# one for coefficients and one vector for misclassification error

df_coeff <- data.frame(matrix(ncol = 18, nrow = 0))
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
  # nested cross-validation on training data to determine optimal regularization
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
  y_test=as.integer(testData$id)
  # compute mis-classification error
  cc = cbind(y_pred,y_test)
  err= (dim(cc[cc[,1]!=cc[,2],])[1] / dim(cc)[1])*100
  print(err)
  mce = append(err, mce)
}

print(mean(mce)); print(max(mce)); print(min(mce))

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
pdf(file="boxplot_coeff_glmnet_lncRNAs_exclude_monoexonic.pdf", width=9, height = 8)
par(mar=c(11,5,5,2))
ggplot(df_coeff_new, aes(x=reorder(name, value, FUN=median),y=value))+
  geom_boxplot(fill="indianred1", colour="black", notch = FALSE, lwd=0.5)+
  scale_y_continuous(name = "Coefficients' values")+
  scale_x_discrete(name="Features")+theme_bw()+
  theme(axis.text.y = element_text(size = 17, face="bold"), axis.text.x = element_text(angle = 90, size =17, vjust = 0.5, hjust=1, face="bold"), plot.title=element_text(size=20, face="bold"), axis.title=element_text(size=17, face="bold"))+
  labs(title="Regularized logisitc regression lncRNAs (excl mono)")
dev.off()
