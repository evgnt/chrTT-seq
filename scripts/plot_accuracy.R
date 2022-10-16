h = rep(c("transcript"), 20)
i = rep(c("mRNA","lncRNA"),each=10)
j=c(); j= append(j,c(24.06577,24.36472,30.23952,25.56054,27.09581,26.15845,27.35426,26.64671,25.26158,28.25112))
j= append(j,c(22.15190, 22.78481, 15.82278, 20.38217, 16.45570, 13.92405, 12.10191, 15.82278, 18.98734, 10.12658))
df_mce = data.frame(h,i,j)
colnames(df_mce)=c("t", "type","mce")
pdf(file="boxplot_mce_glmnet.pdf", width=5, height = 7)
par(mar=c(9,5,5,2))
ggplot(df_mce, aes(x=reorder(type, mce, FUN=median),y=mce))+
  geom_boxplot(fill=c("indianred1","turquoise2"), colour="black", notch=TRUE, lwd=1.4)+
  scale_y_continuous(name = "MSE (%)", limits=c(0, 50))+
  scale_x_discrete(name="")+theme_bw()+
  theme(axis.text.y = element_text(size = 17, face="bold"), axis.text.x = element_text(angle = 90, size =17, vjust = 0.5, hjust=1, face="bold"), plot.title=element_text(size=20, face="bold"), axis.title=element_text(size=17, face="bold"))+
  labs(title="Misclassification error")
dev.off()