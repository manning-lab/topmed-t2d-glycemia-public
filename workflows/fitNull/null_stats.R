## Sample_Stats. R
## Author: Jasen Jackson
## First version completed: 11/3/2017
## Goal: summarize characteristics of sample for poster presentations


## User specifies: trait file, cases, controls, and traits. 
# TODO: 1) Optimize parameter input process using command line args 
#       2) Make "traits" a parameter

## import ggplot and dataframe
library("ggplot2")
df <- read.table("Downloads/t2d.ped", header = T)

## remove rows with "NA" for traits of interest (BMI)
df = df[!is.na(df$last_exam_BMI),]
df = df[!is.na(df$last_exam_age),]

## define cases and controls .. split each by sex
df = df[df$T2D == 2 | df$T2D == 1 | df$T2D == 0,]
case = df[df$T2D == 2 | df$T2D == 1,]
control = df[df$T2D == 0,]

#males:
m = df[df$sex==1,]
m_case = case[case$sex==1,]
m_control = control[control$sex==1,]

#females:
f = df[df$sex==2,]
f_case = case[case$sex==2,]
f_control = control[control$sex==2,]

## Sample total for each group... 
print("Sample totals...")
nrow(df) # total in sample
nrow(case) # total # of cases (according to params)
nrow(control) # of controls
nrow(m) # total males
nrow(f) # total females
nrow(m_case) # total male cases
nrow(m_control) # total male controls
nrow(f_case) # total female cases
nrow(f_control) # total female controls

## In the future.. loop through each specified quantitative trait

### AGE ###

#-(Mean, median, std) age of each category

mean(df$last_exam_age) # mean total age
median(df$last_exam_age) # median total age
sd(df$last_exam_age) # sd of total age

mean(case$last_exam_age) # mean age of cases
median(case$last_exam_age) # median age of cases
sd(case$last_exam_age) # sd of age in cases

mean(control$last_exam_age) # mean age of controls
median(control$last_exam_age) # median age of controls
sd(control$last_exam_age) # sd of age in controls

mean(f$last_exam_age) # mean age of females
median(f$last_exam_age) # median age of females
sd(f$last_exam_age) # sd of age in females

mean(m$last_exam_age) # mean age of males
median(m$last_exam_age) # median age of males
sd(m$last_exam_age) # sd of age in males

mean(f_case$last_exam_age) # mean age of female cases
median(f_case$last_exam_age) # median age of female cases
sd(f_case$last_exam_age) # sd of age in female cases

mean(f_control$last_exam_age) # mean age of female controls
median(f_control$last_exam_age) # median age of female controls
sd(f_control$last_exam_age) # sd of age in female controls

mean(m_case$last_exam_age) # mean age of male cases
median(m_case$last_exam_age) # median age of male cases
sd(m_case$last_exam_age) # sd of age in male cases

mean(m_control$last_exam_age) # mean age of male controls
median(m_control$last_exam_age) # median age of male controls
sd(m_control$last_exam_age) # sd of age in male controls

#-Age distribution by study
plot <- ggplot(df, aes(factor(STUDY_ANCESTRY), last_exam_age))
plot + geom_violin()

# Case Age distribution by study
plot <- ggplot(case, aes(factor(STUDY_ANCESTRY), last_exam_age))
plot + geom_violin()

# Control Age distribution by study
plot <- ggplot(control, aes(factor(STUDY_ANCESTRY), last_exam_age))
plot + geom_violin()

#-Age distribution by study & sex
plot <- ggplot(df, aes(factor(STUDY_ANCESTRY), last_exam_age))
plot + geom_violin(aes(fill=factor(sex)))

#-Case distribution by study & sex
plot <- ggplot(case, aes(factor(STUDY_ANCESTRY), last_exam_age))
plot + geom_violin(aes(fill=factor(sex)))

#-Control distribution by study & sex
plot <- ggplot(df, aes(factor(STUDY_ANCESTRY), last_exam_age))
plot + geom_violin(aes(fill=factor(sex)))

### BMI ###

#-(Mean, median, std) BMI of each category

mean(df$last_exam_BMI) # mean total BMI
median(df$last_exam_BMI) # median total BMI
sd(df$last_exam_BMI) # sd of total BMI

mean(case$last_exam_BMI) # mean BMI of cases
median(case$last_exam_BMI) # median BMI of cases
sd(case$last_exam_BMI) # sd of BMI in cases

mean(control$last_exam_BMI) # mean BMI of controls
median(control$last_exam_BMI) # median BMI of controls
sd(control$last_exam_BMI) # sd of BMI in controls

mean(f$last_exam_BMI) # mean BMI of females
median(f$last_exam_BMI) # median BMI of females
sd(f$last_exam_BMI) # sd of BMI in females

mean(m$last_exam_BMI) # mean BMI of males
median(m$last_exam_BMI) # median BMI of males
sd(m$last_exam_BMI) # sd of BMI in males

mean(f_case$last_exam_BMI) # mean BMI of female cases
median(f_case$last_exam_BMI) # median BMI of female cases
sd(f_case$last_exam_BMI) # sd of BMI in female cases

mean(f_control$last_exam_BMI) # mean BMI of female controls
median(f_control$last_exam_BMI) # median BMI of female controls
sd(f_control$last_exam_BMI) # sd of BMI in female controls

mean(m_case$last_exam_BMI) # mean BMI of male cases
median(m_case$last_exam_BMI) # median BMI of male cases
sd(m_case$last_exam_BMI) # sd of BMI in male cases

mean(m_control$last_exam_BMI) # mean BMI of male controls
median(m_control$last_exam_BMI) # median BMI of male controls
sd(m_control$last_exam_BMI) # sd of BMI in male controls

#-BMI distribution by study
plot <- ggplot(df, aes(factor(STUDY_ANCESTRY), last_exam_BMI))
plot + geom_violin()

# Case BMI distribution by study
plot <- ggplot(case, aes(factor(STUDY_ANCESTRY), last_exam_BMI))
plot + geom_violin()

# Control BMI distribution by study
plot <- ggplot(control, aes(factor(STUDY_ANCESTRY), last_exam_BMI))
plot + geom_violin()

#-BMI distribution by study & sex
plot <- ggplot(df, aes(factor(STUDY_ANCESTRY), last_exam_BMI))
plot + geom_violin(aes(fill=factor(sex)))

#-Case distribution by study & sex
plot <- ggplot(case, aes(factor(STUDY_ANCESTRY), last_exam_BMI))
plot + geom_violin(aes(fill=factor(sex)))

#-Control distribution by study & sex
plot <- ggplot(df, aes(factor(STUDY_ANCESTRY), last_exam_BMI))
plot + geom_violin(aes(fill=factor(sex)))




    