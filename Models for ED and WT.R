setwd("")

require(broom)
require(car)
require(curl)
require(gridExtra)
require(MCMCglmm)
require(lme4)
require(tidyverse)
require(knitr)
require(ordinal)
require(mmmgee)
require(geesmv)
require(gtsummary)
require(gt)
require(ICC)
install.packages("geepack")
citation(package = "geepack")
#Loading data
dat=read.csv("puzzlepersonality1.csv",na.strings=c("NA",""))
summary(dat)

dat <- dat %>% group_by(id) %>% mutate(edavg = mean(dat$ed))
tapply(dat$ed, dat$id, mean)
dat.scr <- dat %>% filter(id == "SCR")
mean(dat.scr$ed)
dat.nonscr <- dat %>% filter(id != "SCR")
mean(dat.nonscr$ed)
dat <- dat %>% group_by(id) %>% mutate(total1 = sum(result))
dat %>% group_by(id, total1) %>% summarize(count = n())
####test of exploration, wt, and latency to approach on first attempt####
require(lmerTest)
dat1 <- dat%>% filter(trial == 1)
dat1 <- dat1 %>% filter(!is.na(latency))

####some really strange results, reducing model to see if this helps with convergence
m1 <- glm(result ~ ed + scale(lat1),data = dat1, family = binomial)
summary(m1)
qqnorm(residuals(m1))
plot(residuals(m1))

plot(result ~ latency, data = dat1)

###can't do wt in with the analysis so doing it separately
wtm1 <- glm(result ~  scale(wt),data = dat1, family = binomial)
summary(wtm1)


wtm1bic <- BIC(wtm1)
### taking out wt
wtm1.1 <- glm(result ~1,data = dat1)
anova(wtm1,wtm1.1, test = "LRT")

wtm1.1bic <- BIC(wtm1.1)

###BF for wt

e^((wtm1.1bic - wtm1bic)/2)


###taking out ed
m1.2 <- glm(result ~ latency,data = dat1)
anova(m1,m1.2, test = "LRT")
pchisq(deviance(m1.2)-deviance(m1), 1, lower.tail=FALSE)
m1.2bic <- BIC(m1.2)
dbic <- m1bic - m1.2bic
dbic

###BF for ed
e^0.9935

e^((dbic)/2)

###taking out latency
m1.3 <- glm(result ~ ed, data = dat1, family = "binomial")

anova(m1, m1.3, test = "LRT")
pchisq(deviance(m1.3)-deviance(m1), 1, lower.tail=FALSE)

m1.3bic <- BIC(m1.3)

e^((m1bic - m1.3bic)/2)



###all trials
install.packages("mmmgee")
require(mmmgee)

install.packages("geesmv")
require(geesmv)
m2 <- glmer(result ~ ed + trial + scale(wt) + (1|id), data = dat, family = binomial (link = "logit"))
summary(m2)
model.df <- tidy(m2)
model.df <- model.df[c(-5),]
model.df %>% 
  mutate(or = exp(estimate),  # Odds ratio/gradient
         var.diag = diag(vcov(m2)),  # Variance of each coefficient
         or.se = sqrt(or^2 * var.diag))  # Odds-ratio adjusted 
??geem2
dat$id <- as.factor(dat$id)
dat <- as.data.frame(dat)
mgee <- geem2(result ~ ed  + trial * scale(wt), id = id, data = dat, family = binomial, corstr = "exchangeable", sandwich = T)
mgee.inter <- geem2(result ~ ed * trial + trial + scale(wt), id = id, data = dat, family = binomial, corstr = "exchangeable", sandwich = T)
summary(mgee)
mgee
summary(mgee.inter)

qqnorm(residuals(mgee))
require(DHARMa)
s1 <- simulate(mgee, nsim = 10)
testResiduals(s1, plot = T)
plot(s1)

citation(package = mmmgee)
citation(package = stats)
require(wgeesel)
wgee.1 <- wgee(result ~ ed + trial + scale(wt), data = dat, id = dat$id, family = "binomial", corstr = "exchangeable")
summary(wgee.1)
QIC.gee(wgee.1)
wgee.2 <- wgee(result ~ ed + trial + scale(wt), data = dat, id = dat$id, family = "binomial", corstr = "exchangeable")
summary(wgee.2)
QIC.gee(wgee.2)

require(geepack)
###USE THESE MODELS FOR REVISION
geep <- geeglm(result ~ ed + trial + scale(wt), data = dat, id = as.factor(id), family = "binomial", corstr = "exchangeable")
summary(geep)

geep.2 <- geeglm(result ~ ed + scale(wt) + trial * scale(wt), data = dat, id = as.factor(id), family = "binomial", corstr = "exchangeable")
summary(geep.2)

geep.3 <- geeglm(result ~ ed * trial + scale(wt) + trial * scale(wt), data = dat, id = as.factor(id), family = "binomial", corstr = "exchangeable")
summary(geep.3)  
m2.1 <- glm(result ~ ed + scale(wt) * trial, data = dat, family = binomial (link = "logit"))

summary(m2.1)
wl.gee <- GEE.var.wl(result~ed+scale(wt)+trial,id="id",family=binomial,dat, corstr = "exchangeable")

s <- c(-2.3498706, 1.1564908, 0.4522724, 0.5865347)
sqrtexch <- sqrt(wl.gee$cov.beta)
sqrtexch
pv <- s/(sqrt(mean(wl.gee$cov.beta)))
c <- var(wl.gee$cov.beta)/(2*mean(wl.gee$cov.beta))
d <- (2*mean(wl.gee$cov.beta^2))/var(wl.gee$cov.beta)

t <- pv/((sqrtexch/(c*d)))
t
pval <-  2*pnorm(-abs(t))
pval

var(dat$ed)
wl.gee1 <- GEE.var.wl(result~ scale(wt) +trial, id = "id", family = binomial, dat, corstr = "exchangeable")

anova(wl.gee, wl.gee1)

m2bic <- BIC(m2)

m2.1 <- glmer(result ~ scale(wt) + trial + (1|id), data = dat, family = binomial)
summary(m2.1)
anova(m2, m2.1)

m2.1bic <- BIC(m2.1)

e^((m2.1bic - m2bic)/2)

m2.2 <- glmer(result ~ ed + trial + (1|id), data = dat, family = binomial)
summary(m2.2)
anova(m2,m2.2)

m2.2bic <- BIC(m2.2)

e^((m2.2bic-m2bic)/2)

m2.3 <- glmer(result ~ ed + scale(wt) + (1|id), data = dat, family = binomial)
summary(m2.3)
anova(m2,m2.3)

m2.3bic <- BIC(m2.3)

e^((m2.3bic-m2bic)/2)

sim1 <- simulateResiduals(fittedModel = m2, n = 500)
plotSimulatedResiduals(simulationOutput = sim1)
testOverdispersion(simulationOutput = sim1)


AIC(m2)
AIC(m2.1)
AIC(m2.2)
AIC(m2.3)
###test how wt influenced changes across trial when successful###
dat2.1 <- dat %>% filter(result == 1)
unique(dat$id)
count(unique(dat2.1$id))
dat2 <- dat %>% filter(result == 1)

dat3 <- dat %>% filter(result == 0)

dat2 <- dat2 %>% group_by(id) %>% mutate(trial1 = 1:n())
dat2 <- dat2 %>% group_by(id) %>% mutate(trialsum = sum(result))

dat2 <- dat2 %>% filter(trialsum >1)

dat3 <- dat3 %>% group_by(id) %>% mutate(trial1 = 1:n())

m3 <- lmer(wtlog ~ ed + trial1 + (1|id), data = dat2)
summary(m3)

m3.1 <- lmer(scale(wt) ~ ed + trial1 + (1|id), data = dat2)
summary(m3.1)

pro <- profile(m3.1)
confint(pro, method = wald)

bicm3.1 <- BIC(m3.1)

bfm3.1 <- lmer(wtlog ~ trial1 + (1|id), data = dat2)

bicbfm3.1 <- BIC(bfm3.1)

dat2 <- as.data.frame(dat2)
length(unique(dat2$id))
hist(dat$wt)
require(mmmgee)
dat2$id <- as.numeric(dat2$id)
dat2 <- as.data.frame(dat2)
summary(dat2)
dim(dat2)
length(dat2$id)
gee.wt <- geem2(scale(wt) ~ ed + trial1, id = id, data = dat2, family = gaussian ,corstr = "exchangeable", sandwich = T)

summary(gee.wt)
###BF for ed
e^((bicbfm3.1 - bicm3.1)/2)

m4 <- lmer(wt ~ ed + trial1 + (1|id), data = dat3)
summary(m4)

hist(dat2$wtlog)
dat2$wtlog <- log(dat2$wt)

qqnorm(residuals(m3))
qqline(residuals(m3))
plot(residuals(m3))
abline(h = 0)

qqnorm(residuals(m3.1))
qqline(residuals(m3.1))
plot(residuals(m3.1))
abline(h = 0)


#####try to get at how they got into the box changes over trials####

dat.1 <- dat %>% group_by(id) %>% mutate(how.total )

######plots for glm model (m1)
newdat <- augment(m1, type.predict = "response") %>%
  select(scale.lat1., ed, .fitted, result)

newdat1 <- augment(mgee, type.predict = "response") %>%
  select(scale.latency., ed, .fitted, result)

dat$pred <- predict(geep.2, newdata = dat)
g1 <- ggplot(newdat, aes(x = scale.lat1., y = .fitted)) +
  stat_smooth(data = newdat,
              aes(y = .fitted),
              method = "glm", method.args = list(family = "binomial"), se = T, alpha = 0.25, size = 0.4, color = "blue") +
  geom_point(data = newdat, aes(y = result)) +
  labs(y = "Solved Puzzle Box", x = "Latency to Approach") +
  scale_y_continuous(limits = c(0,1), breaks = c(0,1)) +
  ggtitle('a') +
  theme_classic() + 
  theme(
    axis.line.x = element_line(colour = "grey50"),
    axis.line.y = element_line(colour = "grey50"),
    axis.title = element_text(size = 16),
    legend.position = "None")
g1

g2 <- ggplot(newdat, aes(x = scale.wt., y = result)) +
  stat_smooth(data = newdat,
              aes(y = .fitted),
              method = "glm", method.args = list(family = "binomial"), se = T, alpha = 0.25, size = 0.4, color = "blue") +
  geom_point(data = newdat, aes(y = result)) +
  labs(y = "", x = "Persistenc (wt)") +
  scale_y_continuous(limits = c(0,1), breaks = c(0,1)) +
  ggtitle('b') +
  theme_classic() + 
  theme(
    axis.line.x = element_line(colour = "grey50"),
    axis.line.y = element_line(colour = "grey50"),
    axis.title = element_text(size = 16),
    legend.position = "None")

g2

g3 <- ggplot(newdat, aes(x = ed, y = result)) +
  stat_smooth(data = newdat,
              aes(y = .fitted),
              method = "glm", method.args = list(family = "binomial"), se = T, alpha = 0.25, size = 0.4, color = "blue") +
  geom_jitter(height = 0, width = 0.1) +
  labs(y = "", x = "Exploration Diversity") +
  scale_y_continuous(limits = c(0,1), breaks = c(0,1)) +
  ggtitle('b') +
  theme_classic() + 
  theme(
    axis.line.x = element_line(colour = "grey50"),
    axis.line.y = element_line(colour = "grey50"),
    axis.title = element_text(size = 16),
    legend.position = "None")
g3

png('glmgraph.png',height=6,width=10,pointsize=11,family="serif",units="in",res=300)
grid.arrange(g1,g3, nrow=1) # Write the grid.arrange in the file
dev.off()


###plots for glmm
summary(m2)
newdat <- augment(m2) %>%
  select(trial, scale.wt., ed, .fitted, result)

dat$predict1 <- predict(geep.2, data = dat, type = "response")
dat$pred3 <- ifelse(dat$predict1 <0.5, 0, 1)
summary(dat$pred3)

dat <- dat %>% mutate(predict2 = exp(predict)/(1+exp(predict)))
summary(dat$predict1)
summary(dat$predict2)
g2.1 <- ggplot(dat, aes(x = trial, y = result)) +
  stat_smooth(data = dat,
              aes(y = predict1),
              method = "glm", method.args = list(family = "binomial"), se = T, alpha = 0.25, size = 0.4, color = "blue") +
  geom_jitter(data = dat, aes(y = result), height = 0, width = 0.1) +
  labs(y = "Solved Puzzle Box", x = "Trial") +
  scale_y_continuous(limits = c(0,1), breaks = c(0,1)) +
  ggtitle('a') +
  theme_classic() + 
  theme(
    axis.line.x = element_line(colour = "grey50"),
    axis.line.y = element_line(colour = "grey50"),
    axis.title = element_text(size = 16),
    legend.position = "None")
g2.1

dat$trial <- as.factor(dat$trial)
g2.2 <- ggplot(dat, aes(x = scale(wt), y = result, by = trial, color = trial)) +
  stat_smooth(data = dat,
              aes(y = predict1),
              method = "glm", method.args = list(family = "binomial"), se = T, alpha = 0.1, size = 0.75) +
  geom_point(data = dat, aes(y = result)) +
  labs(y = "Solved Puzzle Box", x = "Persistence (s)") +
  scale_y_continuous(limits = c(0,1), breaks = c(0,1)) +
  ggtitle('a') +
  theme_classic() + 
  theme(
    axis.line.x = element_line(colour = "grey50"),
    axis.line.y = element_line(colour = "grey50"),
    axis.title = element_text(size = 16))

g2.2

names(dat)
g2.3 <- ggplot(dat, aes(x = ed, y = result)) +
  stat_smooth(data = dat,
              aes(y = predict1),
              method = "glm", method.args = list(family = "binomial"), se = T, alpha = 0.25, size = 0.4, color = "blue") +
  geom_jitter(height = 0, width = 0.1) +
  labs(y = "", x = "Exploration Diversity") +
  scale_y_continuous(limits = c(0,1), breaks = c(0,1)) +
  ggtitle('b') +
  theme_classic() + 
  theme(
    axis.line.x = element_line(colour = "grey50"),
    axis.line.y = element_line(colour = "grey50"),
    axis.title = element_text(size = 16),
    legend.position = "None")
g2.3

png('glmmgraph1.3.24.png',height=6,width=10,pointsize=11,family="serif",units="in",res=300)
grid.arrange(g2.2,g2.3, nrow=1) # Write the grid.arrange in the file
dev.off()


###graph for wt
summary(m3)
require(broom)
install.packages("broom.mixed")
require(broom.mixed)
detach("package:lmerTest")
names(newdata)
newdata <- augment(m3.1) %>% 
  select(ed, trial1, .fitted)
newdata2 <- predict(m3.1)
dat2$predicted <- predict(m3.1)
dat2$predicted.1 <- predict(gee.wt, data = dat2, type = "response")


wtg1 <- ggplot(dat2, aes(x = trial1, y = predicted.1, group = trial1)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(height = 0.05, width = 0.05), shape = 21) +
  labs(y = "Work Time (s)", x = "Trial") +
  ggtitle('a') +
  theme_classic() + 
  theme(
    axis.line.x = element_line(colour = "grey50"),
    axis.line.y = element_line(colour = "grey50"),
    axis.title = element_text(size = 16),
    legend.position = "None")
wtg1

wtg2 <- ggplot(dat2, aes(x = ed, y = predicted.1, group = ed)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(height = 0.05, width = 0.05),shape = 21) +
  labs(y = "", x = "Exploration Diversity") +
  scale_x_continuous(breaks = c(1,2,3)) +
  ggtitle('b') +
  theme_classic() + 
  theme(
    axis.line.x = element_line(colour = "grey50"),
    axis.line.y = element_line(colour = "grey50"),
    axis.title = element_text(size = 16),
    legend.position = "None")

wtg2
hist(newdata$ed)


png('lmmwtgraph1.png',height=6,width=10,pointsize=11,family="serif",units="in",res=300)
grid.arrange(wtg1, wtg2, nrow=1) # Write the grid.arrange in the file
dev.off()


######graph for ed and trial



exg1 <- ggplot(dat2, aes(x = trial1, y = ed, group = trial1)) +
  geom_boxplot() +
  geom_point(aes(fill = ed),shape = 21) +
  stat_smooth(data = dat2,
              aes(group = 1),
              method = "lm",se = T, alpha = 0.25, size = 0.4, color = "blue") +
  labs(y = "Exploration Diversity", x = "Trial") +
  theme_classic() + 
  theme(
    axis.line.x = element_line(colour = "grey50"),
    axis.line.y = element_line(colour = "grey50"),
    legend.position = "None")

exg1
hist(newdata$ed)


png('lmmwtgraph1.png',height=6,width=10,pointsize=11,family="serif",units="in",res=300)
grid.arrange(wtg1, wtg2, nrow=1) # Write the grid.arrange in the file
dev.off()

####work time for successful vs. unsuccessful trials
summary(dat)
wtall <- ggplot(dat, aes(x = result, y = wt)) +
  geom_boxplot(aes(group = result)) +
  geom_point(aes(fill = wt)) +
  theme_classic() + 
  theme(
    axis.line.x = element_line(colour = "grey50"),
    axis.line.y = element_line(colour = "grey50"),
    legend.position = "None")
wtall




####trying out data from paper
data(dental)
data_alt <- reshape(dental, direction="long", timevar="Time",varying=names(dental)[3:6], v.names="response", times=c(8,10,12,14))
data_alt <- data_alt[order(data_alt$subject),]
data_alt$gender <- as.numeric(data_alt$gender)
data_alt$Time <- sqrt(data_alt$Time)
formula <- response~Time+gender

wl.exch <- GEE.var.wl(formula,id="subject",family=gaussian,data_alt,corstr="exchangeable") ##Exchangeable correlation structure;
wl.exch$cov.beta
sqrt(wl.exch$cov.beta)
mean(wl.exch$cov.beta)

s <- c(6.07712, 4.319197, 2.321023)

pv <- s/sqrt(mean(wl.exch$cov.beta))
c <- var(wl.exch$cov.beta)/2*mean(wl.exch$cov.beta)
d <- (2*mean(wl.exch$cov.beta^2))/var(wl.exch$cov.beta)

t <- pv/(sqrt(wl.exch$cov.beta/c*d))
t
pval <-  pnorm(-abs(t))
pval
require(ordinal)
data(soup)
dat <- subset(soup, as.numeric(as.character(RESP)) <=  24)
dat$RESP <- dat$RESP[drop=TRUE]
m1 <- clmm2(SURENESS ~ PROD, random = RESP, data = dat, link="logistic",  Hess = TRUE,doFit=T)
summary(m1)
str(dat)

pred <- sapply(as.character(1:6), function(x){ newdata1=data.frame(PROD=factor(c("Ref", "Test"), levels=levels(dat$PROD)), SURENESS=factor(c(x,x), levels=levels(dat$SURENESS))  );predict(m1, newdata=newdata1)})


edm1 <- clmm2(as.factor(ed) ~ trial, random = id, data = dat2, Hess = T)
summary(edm1)
dat2$ed <- as.factor(dat2$ed)
str(dat2)
edm2 <- clmm(as.factor(ed) ~ trial1 + (1|id), data = dat2)
summary(edm2)

pred <- sapply(as.character(1:3), function(x){ newData <- data.frame(trial = rep(1:5,8),
                      ed = factor(c(x,x), levels = levels(dat2$ed)));predict(edm2, newdata=newData)})


newData <- data.frame(trial = rep(1:5,8),
                      ed = factor(c("1","1"), levels = levels(dat2$ed)))
                    predict(edm2, newdata = newData)
pred <- predict(edm1, newdata = dat2, type = class)
pred
cbind(newData, round(predict(edm1, newdata = newData)$fit, 3),
      "class" = predict(edm1, newdata = newData, type = "class")$fit)

dat2$pred1 <- predict(edm1, data = dat2, type = "class")
summary(dat2$pred1)



require(mmmgee)
data(karatosis)

####try stacking tables####

summary(m1)
m1.tidy <- tidy(m1)
m1.tidy

summary(mgee)
tidy.gee <- function(x, conf.int = FALSE, conf.level = 0.95, ...) {
  
  result <- summary(x)$coefficients[,-3] %>%
    tibble::as_tibble(rownames = "term") %>%
    dplyr::rename(estimate = Estimate,
                  std.error = `Model SE`,
                  statistic = `Wald`,
                  p.value = `p`)
  
  if (conf.int) {
    ci <- confint(x, level = conf.level)
    result <- dplyr::left_join(result, ci, by = "term")
  }
  
  result
}

mgee.tidy <- tidy.gee(mgee)

mgee.tidy

gee.wt.tidy <- tidy.gee(gee.wt)

gee.wt.tidy


m1.tidy
m1.tidy_df <- as.data.frame(m1.tidy)

m1.tidy_df$M1 <- ""


m1.tidy_df <- m1.tidy_df[,c("M1", "term", "estimate", "std.error", "statistic", "p.value")]

m1.tidy_df

m1.tidy_df[2,2] = "Exploration Diversity"
m1.tidy_df[3,2] = "Latency"


m1.tidy_df <- m1.tidy_df %>% rename("Variable" = term, "Estimate" = estimate, "Std Error" = std.error, "Statistic" = statistic, "P value" = p.value)


m1.tidy_df

gt.m1 <- gt(m1.tidy_df) %>% 
  tab_stubhead(label = "M1") %>% 
  cols_align('center') %>% 
  fmt_number(columns = vars("Estimate", "Std Error","Statistic", "P value"),
             decimals = 3)

gt.m1
gtsave(gt.m1, "First trial.rtf")


mgee.tidy
mgee.tidy_df <- as.data.frame(mgee.tidy)

mgee.tidy_df$M2 <- ""
mgee.tidy_df <- mgee.tidy_df[,c("M2", "term", "estimate", "std.error", "statistic", "p.value")]

mgee.tidy_df

mgee.tidy_df[2,2] = "Exploration Diversity"
mgee.tidy_df[3,2] = "Trial"
mgee.tidy_df[4,2] = "Persistence (s)"


mgee.tidy_df <- mgee.tidy_df %>% rename("Variable" = term, "Estimate" = estimate, "Std Error" = std.error, "Statistic" = statistic, "P value" = p.value)

gt.mgee <- gt(mgee.tidy_df) %>% 
  cols_align('center') %>% 
  fmt_number(columns = vars("Estimate", "Std Error", "Statistic", "P value"),
             decimals = 3)

gt.mgee
gtsave(gt.mgee, "All trials.rtf")


gee.wt.tidy_df <- as.data.frame(gee.wt.tidy)

gee.wt.tidy_df$M3 <- ""
gee.wt.tidy_df <- gee.wt.tidy_df[,c("M3", "term", "estimate", "std.error", "statistic", "p.value")]


gee.wt.tidy_df

gee.wt.tidy_df[2,2] = "Exploration Diversity"
gee.wt.tidy_df[3,2] = "Trial"
gee.wt.tidy_df


gee.wt.tidy_df <- gee.wt.tidy_df %>% rename("Variable" = term, "Estimate" = estimate, "Std Error" = std.error, "Statistic" = statistic, "P value" = p.value)

gt.gee.wt <- gt(gee.wt.tidy_df) %>% 
  cols_align('center') %>% 
  fmt_number(columns = vars("Estimate", "Std Error", "Statistic", "P value"),
             decimals = 3)

gt.gee.wt

gtsave(gt.gee.wt, "Learning over successful trials.rtf")

###try to do ICC of how###

dat$how <- as.factor(dat$how)
summary(dat$how)

dat.how <- dat %>% filter(how != "Unsolved")

dat.how <- dat.how %>% mutate(how.2 = recode(how,
                                            Lever = 1,
                                            Lid = 2,
                                            String = 3))
dat.how$how.2 <- as.numeric(dat.how$how.2)
dat.how$id <- as.factor(dat.how$id)
summary(dat.how$how.2)
table(dat.how$id, dat.how$how.2)
ICCbare(x = id, y = how.2, data = dat.how)
install.packages("irr")
require(irr)
install.packages("rel")
require(rel)
dat.how.df <- dat.how %>% select(how.2, id, trial)

dat.how.df <- dat.how.df %>% pivot_wider(names_from = trial, values_from = how.2)

dat.how.df <- dat.how.df %>% select(2:6)
kappam.light(dat.how.df)

kendall(ratings = dat.how.df, TRUE)

