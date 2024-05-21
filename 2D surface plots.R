library(visreg)

names(SEM_df)
names(SEM_df)[12]="Tmean"
names(SEM_df)[13]="Tseas"

#Model SR 
##-------------------------- Include interactions 
m1_I = lmer(S_perc ~ Latitude * Tmean + 
              Latitude * Tseas + 
              Tmean * SUD_PC1 + 
              Tmean * LCc_PC2 +
              Tseas * SUD_PC1 + 
              Tseas * LCc_PC2 +
              Flow_PC1 + Flow_PC2 + (1|Island/Ecoregions), REML = F, data = SEM_df) 

library(MuMIn)
options(na.action="na.fail")

d1=dredge(m1_I, rank="AICc", beta = "sd", evaluate = T)
subset(d1, delta < 4)

#models with delta.aicc < 4
ma1=model.avg(d1, subset = delta < 4)
summary(ma1)

avgmod.95p <- get.models(d1, cumsum(weight) <= .95)

#or as a 95% confidence set:
avgmod.95p <- model.avg(d1, cumsum(weight) <= .95)
conf = as.data.frame(confint(avgmod.95p))
conf$variable=rownames(conf)
conf$effect=(conf$`2.5 %` + conf$`97.5 %`)/2
varImp=sw(ma1)

#Model FRic 
##-------------------------- Include interactions 
m2_I = lmer(FRic_perc ~ Latitude * Tmean + 
              Latitude * Tseas + 
              Tmean * SUD_PC1 + 
              Tmean * LCc_PC2 +
              Tseas * SUD_PC1 + 
              Tseas * LCc_PC2 +
              Flow_PC1 + Flow_PC2 + (1|Island/Ecoregions), REML = F, data = SEM_df) 

library(MuMIn)
options(na.action="na.fail")

d2=dredge(m2_I, rank="AICc", beta = "sd", evaluate = T)
subset(d2, delta < 4)

#models with delta.aicc < 4
ma2=model.avg(d2, subset = delta < 4)
summary(ma2)

#Functional redundancy 
##-------------------------- Include interactions 
m3_I = lmer(FRed_perc ~ Latitude * Tmean + 
              Latitude * Tseas + 
              Tmean * SUD_PC1 + 
              Tmean * LCc_PC2 +
              Tseas * SUD_PC1 + 
              Tseas * LCc_PC2 +
              Flow_PC1 + Flow_PC2 + (1|Island/Ecoregions), REML = F, data = SEM_df) 

library(MuMIn)
options(na.action="na.fail")

d3=dredge(m3_I, rank="AICc", beta = "sd", evaluate = T, extra = "R^2")
subset(d3, delta < 4)
range(d3$'R^2')
mean(d3$'R^2')
#models with delta.aicc < 4
ma3=model.avg(d3, subset = delta < 4)
summary(ma3)
r.squaredGLMM(ma3)


##-------------------------------- Plots of marginal effects of Latitude and the interactions ----------------------------##
#Relationship with latitude
ylab_SR=substitute(SR~("%.decade"^-1))
p1_1=visreg(m1_I, "Latitude", type="conditional", gg = TRUE, fill=list(col="grey"), line=list(col="midnight blue")) + 
  theme_bw() +
  geom_point(size = 2, col="darkgrey") +
  ylab(ylab_SR)
p1_1

p1_2 = interactions::interact_plot(m1_I, pred = SUD_PC1, modx = Tmean, plot.points=T) + 
  theme_bw() +
  xlab("Land-use 2") +
  ylab(ylab_SR) 
p1_2  

interactions::interact_plot(m1_I, pred = SUD_PC1, modx = Tseas, plot.points=T) + 
  theme_bw() +
  xlab("Land-use 2") +
  ylab(ylab_SR)

p1_2D = visreg2d(m1_I, "SUD_PC1", "Tmean", type="conditional", plot.type="gg") +
  theme_bw() +
 # viridis::scale_fill_viridis(name = ylab_SR, option="rocket") +
  scale_fill_gradient2(name = ylab_SR) + xlab("Land-use (2)")

visreg2d(m1_I, "SUD_PC1", "Tseas", type="conditional", plot.type="gg") +
  theme_bw()

ylab_FRic=substitute(FRic~("%.decade"^-1))
p2_1=visreg(m2_I, "Latitude", type="conditional", gg = TRUE, fill=list(col="grey"), line=list(col="midnight blue")) + 
  theme_bw() +
  geom_point(size = 2, col="darkgrey") +
  ylab(ylab_FRic)
p2_1

p2_2 = interactions::interact_plot(m2_I, pred = LCc_PC2, modx = Tmean, plot.points=T) + 
  theme_bw() +
  xlab("Land-use 1") +
  ylab(ylab_FRic) 
p2_2

p2_2D = visreg2d(m2_I, "LCc_PC2", "Tmean", type="conditional", plot.type="gg") +
  theme_bw() +
#  viridis::scale_fill_viridis(name = ylab_SR, option="cividis") +
  scale_fill_gradient2(name = ylab_FRic) + xlab("Land-use (1)")

ylab_Fred=substitute(FRed~("%.decade"^-1))

visreg(m3_I, "Latitude", type="conditional", gg = TRUE, fill=list(col="grey"), line=list(col="midnight blue")) +
  theme_bw() +
  geom_point(size = 2, col="darkgrey") +
  ylab(ylab_Fred)

p3_1 = ggplot(SEM_df, aes(y=FRed_perc, x=Latitude)) + 
  geom_point(size = 2, col="darkgrey") + 
  theme_bw() + ylab(ylab_Fred)
p3_1

p3_2 = interactions::interact_plot(m3_I, pred = LCc_PC2, modx = Tseas, plot.points=T) +
  theme_bw() +
  xlab("Land-use 1") +
  ylab(ylab_FRic) 
p3_2

p3_2D = visreg2d(m3_I, "LCc_PC2", "Tseas", type="conditional", plot.type="gg") +
  theme_bw() +
  scale_fill_gradient2(name = ylab_Fred) + xlab("Land-use (1)")

grid_most=ggpubr::ggarrange(p1_1, p2_1, p3_1,
               p1_2D, p2_2D, p3_2D, nrow = 2, ncol= 3, align="hv", labels = c("A", "B", "C", "","",""))
ggsave(file="Grid of relationships with Latitude and most important interactions 2D surface.jpeg", grid_most, width=600*0.75, height=250*0.75, units="mm")

##------------------------------------- Plots of Akaike Weights 

#SR 
varImp=sw(ma1)
varImp=as.data.frame(varImp)
varImp$variable=rownames(varImp)

varImp$variable=c("Land-use 1", "Latitude", "Land-use 2", "Tmean", "Tmean x Land-use 2", "Tseas", 
                  "Tseas x Land-use 2", "Flow 1", "Latitude x Tmean", "Flow 2", "Tmean x Land-use 1")

varImp_SR_int=ggplot(varImp, aes(x=varImp, y=reorder(variable, varImp))) + 
  geom_bar(stat = "identity") +
  theme_bw() +
  xlab("Variable importance (%)") + ylab(NULL)
varImp_SR_int

varImp$Group=c("Land-use", "Latitude", "Land-use", "Climate", "Climate x Land-use", "Climate", "Climate x Land-use", "Flow", "Latitude x Climate", "Flow", "Climate x Land-use")

varImp_Group=varImp %>% group_by(Group) %>% summarise_at(vars(varImp), list(name = mean))
varImp_Group$Group=factor(varImp_Group$Group, levels=c("Climate x Land-use", "Latitude x Climate", "Flow", "Land-use", "Climate", "Latitude"))

g1 = ggplot(varImp_Group, aes(x=name, y=Group)) + 
  geom_bar(stat = "identity") +
  theme_bw() +
  xlab("Akaike weight") + ylab(NULL)

#FRic 
varImp=sw(ma2)
varImp=as.data.frame(varImp)
varImp$variable=rownames(varImp)

varImp$variable=c("Latitude", "Flow 2", "Land-use 1", "Tmean", "Flow 1", 
                  "Tmean x Land-use 1", "Tseas", "Land-use 2", "Latitude x Tmean",   
                  "Tseas x Land-use 1")

varImp_SR_int=ggplot(varImp, aes(x=varImp, y=reorder(variable, varImp))) + 
  geom_bar(stat = "identity") +
  theme_bw() +
  xlab("Akaike weight") + ylab(NULL)
varImp_SR_int

varImp$Group=c("Latitude", "Flow", "Land-use", "Climate", "Flow", "Climate x Land-use", "Climate", "Land-use", "Latitude x Climate", "Climate x Land-use")

varImp_Group=varImp %>% group_by(Group) %>% summarise_at(vars(varImp), list(name = mean))
varImp_Group$Group=factor(varImp_Group$Group, levels=c("Climate x Land-use", "Latitude x Climate", "Flow", "Land-use", "Climate", "Latitude"))

g2 = ggplot(varImp_Group, aes(x=name, y=Group)) + 
  geom_bar(stat = "identity") +
  theme_bw() +
  xlab("Akaike weight") + ylab(NULL) +
  theme(axis.text.y = element_blank())

#FRed
varImp=sw(ma3)
varImp=as.data.frame(varImp)
varImp$variable=rownames(varImp)

varImp$variable=c("Land-use 1", "Flow 2", "Tmean", "Tseas", "Latitude", "Flow 1", 
                  "Tseas x Land-use 1", "Land-use 2", "Tmean x Land-use 1", "Latitude x Tseas")

varImp_FRed_int=ggplot(varImp, aes(x=varImp, y=reorder(variable, varImp))) + 
  geom_bar(stat = "identity") +
  theme_bw() +
  xlab("Variable importance (%)") + ylab(NULL)
varImp_FRed_int

varImp$Group=c("Land-use", "Flow", "Climate", "Climate", "Latitude", "Flow", "Climate x Land-use", "Land-use", "Climate x Land-use", "Latitude x Climate")

varImp_Group=varImp %>% group_by(Group) %>% summarise_at(vars(varImp), list(name = mean))
varImp_Group$Group=factor(varImp_Group$Group, levels=c("Climate x Land-use", "Latitude x Climate", "Flow", "Land-use", "Climate", "Latitude"))

g3 = ggplot(varImp_Group, aes(x=name, y=Group)) + 
  geom_bar(stat = "identity") +
  theme_bw() +
  xlab("Akaike weight") + ylab(NULL) +
  theme(axis.text.y = element_blank())

egg::ggarrange(g1, g1_b, g2, g2_b, g3, g3_b, nrow = 1)

require(ggpubr)
grid_group = ggpubr::ggarrange(
egg::ggarrange(g1, g1_b, nrow = 1),
egg::ggarrange(g2, g2_b, nrow = 1), 
egg::ggarrange(g3, g3_b, nrow = 1), 
align = "v", nrow = 1, labels = c("A", "B", "C"), label.x = 0.1
)

ggsave(file="Grid of Akaike Weights and variable importance by group.jpeg", grid_group, width=600, height=100, units="mm")

