# Packages
library(tidyr)
library(dplyr)
library(survival)
library(ggplot2)
library(RColorBrewer)


#Creating the data set from multiple files

load("FinalCore2009.rdat")
data2009 <- temp

#Eliminating the columns not common across data files over the years
data2009 <- data2009[,(-150)]
data2009 <- data2009[,(-150)]

#Process repeated individually for remaining data files due to different columns in each year

#2010
load("WA_SIDC_2010_CORE.rdat")
data2010 <- temp
data2010 <- data2010[,-125]
data2010 <- data2010[,(-150)]
data2010 <- data2010[,(-150)]


#2011
load("WA_SIDC_2011_CORE.rdat")
data2011 <- temp
data2011 <- data2011[,(-103)]
data2011 <- data2011[,(-103)]
data2011 <- data2011[,(-111)]
data2011 <- data2011[,(-111)]
data2011 <- data2011[,(-119)]
data2011 <- data2011[,(-119)]
data2011 <- data2011[,(-122)]
data2011 <- data2011[,(-122)]
data2011 <- data2011[,(-125)]
data2011 <- data2011[,(-143)]
data2011 <- data2011[,(-150)]
data2011 <- data2011[,(-150)]


#2012
load("WA_SIDC_2012_CORE.rdat")
data2012 <- temp
data2012 <- data2012[,(-125)]
data2012 <- data2012[,(-143)]
data2012 <- data2012[,(-150)]
data2012 <- data2012[,(-150)]


#2013
load("WA_SIDC_2013_CORE.rdat")
data2013 <- temp


data2013 <- data2013[,(-8)]
data2013 <- data2013[,(-125)]
data2013 <- data2013[,(-143)]
data2013 <- data2013[,(-150)]

#Renaming column names
colnames(data2013)[153] <- "pl_rucc2003"
colnames(data2013)[154] <- "pl_uic2003"


rm(temp)

#Exporting cleaned up data

write.csv(data2009, file = "data2009.csv", row.names = FALSE)
write.csv(data2010, file = "data2010.csv", row.names = FALSE)
write.csv(data2011, file = "data2011.csv", row.names = FALSE)
write.csv(data2012, file = "data2012.csv", row.names = FALSE)
write.csv(data2013, file = "data2013.csv", row.names = FALSE)


#Combine dxn (where n: 1- 25) in a column "diagnosis"
data_comb_diag_2009 <- unite(data2009,diagnosis,dx1,dx2,dx3,dx4,dx5,dx6,dx7,dx8,dx9,dx10,dx11,dx12,dx13,dx14,dx15,dx16,dx17,dx18,dx19,dx20,dx21,dx22,dx23,dx24,dx25, sep=",")

# MRSA column generated to indicate presence of MRSA (MRSA = TRUE)
data2009$mrsa <- rep(FALSE,nrow=data2009)
data2009$mrsa[grep("(?=.*04112)", data_comb_diag_2009$diagnosis, ignore.case = TRUE, perl = TRUE)] <- TRUE 

# Eczema column generated to indicate presence of Eczema (Eczema = TRUE)
data2009$ezcema <- rep(FALSE,nrow=data2009)
data2009$ezcema[grep("(?=.*6929)", data_comb_diag_2009$diagnosis, ignore.case = TRUE, perl = TRUE)] <- TRUE 

#Column created to set flag for matched pair to enable case control analysis
data2009$match.pair <- rep(NA, nrow(data2009))


# for loop to find controls in one data set
ix<-which(data2009$mrsa)
for(i in 1:200){
  temp.age <- data2009[ix[i],"age"]
  temp.female <- data2009[ix[i],"female"]
  temp.amonth <- data2009[ix[i],"amonth"]
  
  temp.ix <- which(data2009$age==temp.age & data2009$female==temp.female & data2009$amonth == temp.amonth)
  
  data2009[ix[i], "match.pair"] <- i
  data2009[sample(temp.ix, 1), "match.pair"] <- i
}

temp.data <- data2009[!is.na(data2009$match.pair),]

# Load data - HCUP data for the years 2009 - 2014
data_6_yrs <- read.csv("case_control_6_yrs.csv")

# Make match pair year-specific
data_6_yrs$match.pair.unique <- as.numeric(paste0(data_6_yrs$match.pair,".",data_6_yrs$year))

# Research Question 1
# Conditional logistic regression
mod.1 <- clogit(mrsa ~ eczema + strata(match.pair.unique), data=data_6_yrs, method="approximate")

summary(mod.1)
exp(confint(mod.1))

#Research Question 2
mod.2 <- clogit(mrsa ~ orproc + strata(match.pair.unique), data=data_6_yrs, method="approximate")
summary(mod.2)
exp(confint(mod.2))

mod.3 <- clogit(mrsa ~ eczema + orproc + strata(match.pair.unique), data=data_6_yrs, method="approximate" )
summary(mod.3)
exp(confint(mod.3))

mod.4 <- clogit(mrsa ~ eczema + orproc + skin_infection + strata(match.pair.unique), data=data_6_yrs, method="approximate" )
summary(mod.4)
exp(confint(mod.4))


# Filter childern & adults from case control data for research research question 3
children_3 <- filter(data_6_yrs, age <=3)
adults <- filter(data_6_yrs, age >= 60)

mod.children <- clogit(mrsa ~ eczema + strata(match.pair.unique), data=children_3, method="approximate" )
summary(mod.children)

mod.adults <- clogit(mrsa ~ eczema + strata(match.pair.unique), data = adults, method = "approximate")
summary(mod.adults)


#  Plot data prep
age_graph <- data_6_yrs %>%
  group_by(age,mrsa) %>%
  filter(mrsa == TRUE) %>%
  summarize(total.count=n())

surgery_graph <- data_6_yrs %>%
  group_by(age,orproc) %>%
  filter(orproc== 1) %>%
  summarize(total.count=n())

skin_infection_graph <- data_6_yrs %>%
  group_by(age,skin_infection) %>%
  filter(skin_infection == TRUE) %>%
  summarize(total.count=n())

eczema_graph <- data_6_yrs %>%
  group_by(age,eczema) %>%
  filter (eczema == TRUE) %>%
  summarize(total.count=n())


# Plot age vs mrsa freq / age vs. orproc freq  
plot(age_graph$age,age_graph$total.count, type="o", xlab="Age", ylab = "Occurrence")
points(surgery_graph$age, surgery_graph$total.count, type = "o", col=2)

# Adding  age in eczema graph so that age's considered are uniform across age_graph, eczema_graph and surgery_graph. 
# The occurrence of eczema for these ages remains zero. 
temp <- age_graph$age[which(!age_graph$age%in%eczema_graph$age)]
eczema_graph_full <- rbind(data.frame(eczema_graph), cbind("age"=temp, "eczema"=rep(T, length(temp)), "total.count" = rep(0, length(temp))))
eczema_graph_full <- eczema_graph_full[order(eczema_graph_full$age),]
points(eczema_graph_full$age, eczema_graph_full$total.count, type = "o", col=4)
lines(age_graph$age, age_graph$total.count + eczema_graph_full$total.count)
legend("topright",legend=c("MRSA","Surgery"), lty=c(1,1), col = c(2,"black"))


plot(surgery_graph$total.count, age_graph$total.count)
abline(0,1)
summary(lm(age_graph$total.count ~ surgery_graph$total.count))

tab <- table("mrsa"=data_6_yrs$mrsa, "orproc"=data_6_yrs$orproc)
tab/sum(tab)


#Plot age vs. eczema freq / age vs. skin infection freq
plot(skin_infection_graph$age,skin_infection_graph$total.count, type = "o", col='red', ylim=c(0,40))
points(eczema_graph$age,eczema_graph$total.count, type = "o", col='green')

# Additional questions - Top primary diagnosis associated with MRSA

case_control <- data_6_yrs

valid_dx1_2009<- case_control[(!is.na(case_control$dx1)) & case_control$mrsa== TRUE & case_control$ayear== "2009", ]
dx1_groups_2009<- valid_dx1_2009 %>%
  group_by(dx1) %>%
  summarise(count= n()) %>%
  filter(count>150)
dx1_desc_2009<- arrange(dx1_groups_2009, desc(count))
dx1_desc_2009$pct <- (dx1_desc_2009$count/sum(dx1_desc_2009$count))*100


valid_dx1_2010<- case_control[(!is.na(case_control$dx1)) & case_control$mrsa== TRUE & case_control$ayear== "2010", ]
dx1_groups_2010<- valid_dx1_2010 %>%
  group_by(dx1) %>%
  summarise(count= n()) %>%
  filter(count>150)
dx1_desc_2010<- arrange(dx1_groups_2010, desc(count))
dx1_desc_2010$pct <- (dx1_desc_2010$count/sum(dx1_desc_2010$count))*100

valid_dx1_2011<- case_control[(!is.na(case_control$dx1)) & case_control$mrsa== TRUE & case_control$ayear== "2011", ]
dx1_groups_2011<- valid_dx1_2011 %>%
  group_by(dx1) %>%
  summarise(count= n()) %>%
  filter(count>150)
dx1_desc_2011<- arrange(dx1_groups_2011, desc(count))
dx1_desc_2011$pct <- (dx1_desc_2011$count/sum(dx1_desc_2011$count))*100

valid_dx1_2012<- case_control[(!is.na(case_control$dx1)) & case_control$mrsa== TRUE & case_control$ayear== "2012", ]
dx1_groups_2012<- valid_dx1_2012 %>%
  group_by(dx1) %>%
  summarise(count= n()) %>%
  filter(count>150)
dx1_desc_2012<- arrange(dx1_groups_2012, desc(count))
dx1_desc_2012$pct <- (dx1_desc_2012$count/sum(dx1_desc_2012$count))*100


valid_dx1_2013<- case_control[(!is.na(case_control$dx1)) & case_control$mrsa== TRUE & case_control$ayear== "2013", ]
dx1_groups_2013<- valid_dx1_2013 %>%
  group_by(dx1) %>%
  summarise(count= n()) %>%
  filter(count>150)
dx1_desc_2013<- arrange(dx1_groups_2013, desc(count))
dx1_desc_2013$pct <- (dx1_desc_2013$count/sum(dx1_desc_2013$count))*100

valid_dx1_2014<- case_control[(!is.na(case_control$dx1)) & case_control$mrsa== TRUE & case_control$ayear== "2014", ]
dx1_groups_2014<- valid_dx1_2014 %>%
  group_by(dx1) %>%
  summarise(count= n()) %>%
  filter(count>150)
dx1_desc_2014<- arrange(dx1_groups_2014, desc(count))
dx1_desc_2014$pct <- (dx1_desc_2014$count/sum(dx1_desc_2014$count))*100

# Stacked bar plot for top primary diagnosis with MRSA

# Rename col name
colnames(dx1_desc_2009)[3] <- "2009"  

colnames(dx1_desc_2010)[3] <- "2010"  

colnames(dx1_desc_2011)[3] <- "2011"  

colnames(dx1_desc_2012)[3] <- "2012"  

colnames(dx1_desc_2013)[3] <- "2013"  

colnames(dx1_desc_2014)[3] <- "2014"  

dx1_desc_2009 <- dx1_desc_2009[,(-2)]
dx1_desc_2010 <- dx1_desc_2010[,(-2)]
dx1_desc_2011 <- dx1_desc_2011[,(-2)]
dx1_desc_2012 <- dx1_desc_2012[,(-2)]
dx1_desc_2013 <- dx1_desc_2013[,(-2)]
dx1_desc_2014 <- dx1_desc_2014[,(-2)]

# Join years into a single data frame
temp <- full_join(dx1_desc_2009,dx1_desc_2010, by = "dx1")

temp <- full_join(temp,dx1_desc_2011, by = "dx1")

temp <- full_join(temp,dx1_desc_2012, by = "dx1")

temp <- full_join(temp,dx1_desc_2013, by = "dx1")

temp <- full_join(temp,dx1_desc_2014, by = "dx1")

data_mrsa_dx1 <- temp

rm(temp)


# Transpose
data_mrsa_dx1 <- t(data_mrsa_dx1)

# Fix column names based on icd code
colnames(data_mrsa_dx1) = data_mrsa_dx1[1, ] # the first row will be the header
data_mrsa_dx1 = data_mrsa_dx1[-1, ]   

# Transpose
data_mrsa_dx1 <- t(data_mrsa_dx1)

# Color Brewer
diverging <- brewer.pal(8,"Spectral")

# Assign 0 to Null values
data_mrsa_dx1[is.na(data_mrsa_dx1)] <- 0

barplot(as.matrix(data_mrsa_dx1),xlim=c(0, ncol(as.matrix(data_mrsa_dx1)) + 5), ylim = range(0,100), col=diverging)
legend("right",cex = 0.80,legend=c("Unspecified septicemia",
                                   "Cellulitis & Abscess of face",
                                   "Conductive hearing loss",
                                   "Cellulitis & Abscess of Upper Arm",
                                   "Cellulitis & Abscess of Buttocks",
                                   "Other post-op infection",
                                   "Cellulitis & Abscess of Trunk",
                                   "Cellulitis & Abscess of Leg"),pt.cex = 2, fill = diverging[8:1])
