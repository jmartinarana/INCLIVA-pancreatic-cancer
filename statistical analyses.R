### R CODE FOR ALL STATISTICAL ANALYSES, FIGURES AND TABLES


####Figure 3

mar <- read_excel("C:/Users/raistlol/Juanbox/Marisol Huerta/MH_1/Base pancreas corregida.xlsx", na = "NA")
mar$dif <- (mar$`FA_PostTreatment 6_8 weeks`-mar$FA_PreTreatment)/mar$FA_PreTreatment
mar_train <- subset(mar, cohort=="first" & FA_PreTreatment > 0)
mar_test <- subset(mar, cohort=="second" & FA_PreTreatment > 0)
mod1 <- coxph(Surv(PFS1_dias, staus_PFS)~dif, data=mar_test)
summary(mod1)

res.cut <- surv_cutpoint(mar_train, time = "PFS1_dias", event = "staus_PFS",
                         variables = c("dif"))

summary(res.cut)


plot(res.cut, "dif", palette = "npg")


mar_test$dif_cat1 <- car::recode(mar_test$dif,
                                 "lo:-0.8475='Low risk';
                              -0.8475:hi='High risk'")

mar_train$dif_cat1 <- car::recode(mar_train$dif,
                                  "lo:-0.8475='Low risk';
                              -0.8475:hi='High risk'")

test1 <- survfit(Surv(PFS1_dias, staus_PFS)~dif_cat1, data=mar_train)
survminer::ggsurvplot(test1, title="Training cohort", legend.title="", legend.labs= c("High risk", "Low Risk"),
                      pval = T, pval.method = T, xlab="Days", ylab="PFS", 
                      risk.table = "abs_pct", risk.table.col="strata", risk.table.title= "Patients at risk (%)",
                      risk.table.y.text = FALSE,
                      tables.y.text = FALSE,surv.median.line = "hv",
                      break.time.by = 60)

mar_test$PFS1_meses <- mar_test$PFS1_dias /30.417

test1 <- survfit(Surv(PFS1_meses, staus_PFS)~dif_cat1, data=mar_test)
survminer::ggsurvplot(test1, title="Cohorte de validación", legend.title="", legend.labs= c("Alto riesgo", "Bajo riesgo"),
                      pval = T, pval.method = T, xlab="Meses", ylab="SLP", 
                      risk.table = "abs_pct", risk.table.col="strata", risk.table.title= "Pacientes en riesgo (%)",
                      risk.table.y.text = FALSE,
                      tables.y.text = FALSE,surv.median.line = "hv",
                      break.time.by = 1)

####Sup Figure 7

library(dplyr)
library(ggplot2)
library(readxl)
library(readxl)

marisol <- read_excel("marisol.xlsx", 
                      col_types = c("text", "date", "text", 
                                    "date", "date", "text", "date", "date", 
                                    "text", "date", "date", "date", "text", 
                                    "date", "text", "date", "text", "date", 
                                    "text", "date", "text", "date", "text", 
                                    "date", "numeric", "skip"))

marisol$linea1_start <- (marisol$StartTr1 - marisol$Diagnosis)/30.417
marisol$linea2_start <- (marisol$StartTr2 - marisol$Diagnosis)/30.417
marisol$linea3_start <- (marisol$StartTr3 - marisol$Diagnosis)/30.417


marisol$linea1_fin <- (marisol$StopTr1 - marisol$Diagnosis)/30.417
marisol$linea2_fin <- (marisol$StopTr2 - marisol$Diagnosis)/30.417
marisol$linea3_fin <- (marisol$StopTr3 - marisol$Diagnosis)/30.417

marisol$L1_status <- as.factor(marisol$`Status LB1`)
marisol$L2_status <- as.factor(marisol$`Status LB2`)
marisol$L3_status <- as.factor(marisol$`Status LB3`)
marisol$L4_status <- as.factor(marisol$`Status L4`)
marisol$L5_status <- as.factor(marisol$`Status L5`)
marisol$L6_status <- as.factor(marisol$`Status L6`)

marisol$L1_date <- difftime(marisol$LB1, marisol$Diagnosis, units = "days")/30.417
marisol$L2_date <- (marisol$LB2 - marisol$Diagnosis)/30.417
marisol$L3_date <- (marisol$L3 - marisol$Diagnosis)/30.417
marisol$L4_date <- (marisol$L4 - marisol$Diagnosis)/30.417
marisol$L5_date <- (marisol$L5 - marisol$Diagnosis)/30.417
marisol$L6_date <- (marisol$L6 - marisol$Diagnosis)/30.417

marisol$T_LFU <- (marisol$`last FUP` -  marisol$Diagnosis)/30.417

marisol$status <- as.factor(marisol$`status last FUP`)


marisol %>% ggplot(aes(as.factor(Patient), T_LFU))+
  coord_flip()+
  labs(x="Patient", y="Months since diagnose", colour="", shape="")+
  geom_segment(aes(x=as.factor(Patient), xend=as.factor(Patient), y=0, yend=T_LFU), size=1, alpha=0.7)+
  geom_segment(mapping=aes(x=as.factor(Patient), xend=as.factor(Patient), y=linea1_start, yend=linea1_fin,  colour = "palegreen3"), alpha = 0.3, 
               size=5)+
  geom_segment(mapping=aes(x=as.factor(Patient), xend=as.factor(Patient), y=linea2_start, yend=linea2_fin,  colour = "orange3"), alpha = 0.3, 
               size=5)+ 
  geom_segment(mapping=aes(x=as.factor(Patient), xend=as.factor(Patient), y=linea3_start, yend=linea3_fin,  colour = "orchid2"), alpha = 0.3, 
               size=5)+ 
  geom_point(aes(as.factor(Patient), L1_date, shape=L1_status),size=3,show.legend = T, fill="blue")+
  geom_point(aes(as.factor(Patient), L2_date, shape=L2_status),size=3,show.legend = T, fill="blue")+
  geom_point(aes(as.factor(Patient), L3_date, shape=L3_status),size=3,show.legend = T, fill="blue")+
  geom_point(aes(as.factor(Patient), L4_date, shape=L4_status),size=3,show.legend = T, fill="blue")+
  geom_point(aes(as.factor(Patient), L5_date, shape=L5_status),size=3,show.legend = T, fill="blue")+
  geom_point(aes(as.factor(Patient), L6_date, shape=L6_status),size=3,show.legend = T, fill="blue")+
  geom_point(aes(as.factor(Patient), T_LFU, shape=status),size=3,show.legend = T, fill="red")+
  theme(legend.title = element_blank(),
        legend.spacing.y = unit(0, "mm"), 
        panel.border = element_rect(colour = "black", fill=NA),
        aspect.ratio = 1, axis.text = element_text(colour = 1, size = 12),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"))+
  guides(colour=guide_legend(override.aes=list(shape=NA)))+
  scale_shape_manual(values=c(1, 21,17, 18), name="direc", labels=c("Alive", "Dead", "KRAS wildtype", "KRAS mutated"))+
  scale_colour_manual(breaks=c("palegreen3",  "orange3","orchid2"),
                      values=c("palegreen3",  "orange3","orchid2"),
                      labels=c("First line","Second line","Third line"))

#Sup Figure 8


mar <- read_excel("C:/Users/raistlol/Juanbox/Marisol Huerta/MH_1/Base pancreas corregida.xlsx", na = "NA")
mar_train <- subset(mar, cohort=="first")
mar <- subset(mar, cohort=="second")

#PFS, validation cohort
mar$PFS1_meses <- mar$PFS1_dias /30.417
test1 <- survfit(Surv(PFS1_meses, staus_PFS)~FA_PreTreatment_cat, data=mar)
survminer::ggsurvplot(test1, title="Pre treatment", legend.title="", legend.labs= c("Negative", "Positive"),
                      pval = T, pval.method = T, xlab="Months", ylab="PFS", conf.int = F,
                      risk.table = "abs_pct", risk.table.col="strata", risk.table.title= "Patients at risk (%)",
                      risk.table.y.text = FALSE,
                      tables.y.text = FALSE,surv.median.line = "hv",
                      break.time.by = 6)

mod <- coxph(Surv(PFS1_meses, staus_PFS)~FA_PreTreatment_cat, data=mar)
summary(mod)


test1 <- survfit(Surv(PFS1_meses, staus_PFS)~FA_PostTreatment_6_8_weeks_cat, data=mar)
survminer::ggsurvplot(test1, title="Post treatment", legend.title="", legend.labs= c("Negative", "Positive"),
                      pval = T, pval.method = T, xlab="Days", ylab="PFS", conf.int = F,
                      risk.table = "abs_pct", risk.table.col="strata", risk.table.title= "Patients at risk (%)",
                      risk.table.y.text = FALSE,
                      tables.y.text = FALSE,surv.median.line = "hv",
                      break.time.by = 6)



#PFS, training cohort
mar_train$PFS1_meses <- mar_train$PFS1_dias /30.417
test1 <- survfit(Surv(PFS1_meses, staus_PFS)~FA_PreTreatment_cat, data=mar_train)
survminer::ggsurvplot(test1, title="Pre treatment", legend.title="", legend.labs= c("Negative", "Positive"),
                         pval = T, pval.method = T, xlab="Months", ylab="PFS", conf.int = F,
                         risk.table = "abs_pct", risk.table.col="strata", risk.table.title= "Patients at risk (%)",
                         risk.table.y.text = FALSE,
                         tables.y.text = FALSE,surv.median.line = "hv",
                         break.time.by = 6,
                         pval.method.coord = c(0.2, 0.30),
                         pval.coord = c(0.2,0.23))
mod <- coxph(Surv(PFS1_dias, staus_PFS)~FA_PreTreatment_cat, data=mar_train)
summary(mod)

test1 <- survfit(Surv(PFS1_meses, staus_PFS)~FA_PostTreatment_6_8_weeks_cat, data=mar_train)
survminer::ggsurvplot(test1, title="Post treatment", legend.title="", legend.labs= c("Negative", "Positive"),
                         pval = T, pval.method = T, xlab="Months", ylab="PFS", conf.int = F,
                         risk.table = "abs_pct", risk.table.col="strata", risk.table.title= "Patients at risk (%)",
                         risk.table.y.text = FALSE,
                         tables.y.text = FALSE,surv.median.line = "hv",
                         break.time.by = 6)

#Sup Figure 9

mar$evol <- paste(mar$FA_PreTreatment_cat, mar$FA_PostTreatment_6_8_weeks_cat)
mar_evol <- subset(mar, evol=="NEG NEG" | evol=="POS NEG" | evol=="POS POS")


test1 <- survfit(Surv(PFS1_dias, staus_PFS)~evol, data=mar_evol)
survminer::ggsurvplot(test1, title="KRAS stauts baseline to first line treatment", legend.title="", legend.labs= c("Remains negative", "Positive to negative", "Remains positive"),
                      pval = "Log-rank p<0.001", pval.method = T, xlab="Days", ylab="PFS", 
                      risk.table = "abs_pct", risk.table.col="strata", risk.table.title= "Patients at risk (%)",
                      risk.table.y.text = FALSE,
                      tables.y.text = FALSE,surv.median.line = "hv",
                      break.time.by = 60)
pairwise_survdiff(Surv(PFS1_dias, staus_PFS)~evol, data=mar_evol)
#Differences found between POS-POS vs NEG-NEG and vs POS-POS and POS-NEG.


###Table 1

mar <- read_excel("C:/Users/Raistlol/JuanBox/Marisol Huerta/MH_1/Base pancreas Juan multivariable_V2.xlsx", na = "NA")
mar$LIVER_METS <- mar$`LIVER METS`
mar$ECOG <- as.factor(mar$ECOG)
mar$N_METS_SITES <- car::recode(mar$`N METS SITES`,
                                "3:hi = 'More than 2'")
mar$PERITONEUM_METS <- mar$`PERITONEUM METS`
mar$firstline <- mar$`1ST LINE TREATMENT`
mar$LUNG_METS <- mar$`LUNG METS`
mar$dif <- (mar$`FA_PostTreatment 6_8 weeks`-mar$FA_PreTreatment)/mar$FA_PreTreatment



explanatory = c("AGE_CONT", "SEX", "ECOG", "LIVER_METS", "PERITONEUM_METS", "LUNG_METS",
                "OTHER", "N_METS_SITES", "firstline", "OS","CA_L1_cat")
dependent = "cohort"
mar %>%
  summary_factorlist(dependent, explanatory,
                     p=TRUE, add_dependent_label=TRUE, p_cat="fisher",  cont="median",total_col =  T) -> t1
knitr::kable(t1, row.names=FALSE, align=c("l", "l", "r", "r", "r"))


###Table 3

mar_test <- subset(mar, cohort=="second" & FA_PreTreatment > 0)
mar_test$INCLIVA_score <- car::recode(mar_test$dif,
                                      "lo:-0.8475='Low risk';
                              -0.8475:hi='High risk'")
mar_test$INCLIVA_score <- relevel(as.factor(mar_test$INCLIVA_score),ref="Low risk")

dependent_pfs = "Surv(PFS1_dias, staus_PFS)"
explanatory = c("AGE", "ECOG", "LIVER_METS", "INCLIVA_score", "CA_L1_cat")

mar_test_pfs <- na.omit(mar_test[,c(explanatory, "PFS1_dias", "staus_PFS")])
# 
mod1 <- coxph(Surv(PFS1_dias, staus_PFS)~ AGE + ECOG + LIVER_METS + INCLIVA_score + CA_L1_cat, data=mar_test_pfs)
summary(mod1)
step(mod1)
#  
best<-coxph(formula = Surv(PFS1_dias, staus_PFS) ~ ECOG + INCLIVA_score, 
            data = mar_test_pfs)
summary(best)
explanatory_multi = c("ECOG", "INCLIVA_score")

mar_test_pfs %>% 
  finalfit(dependent_pfs, explanatory, explanatory_multi) %>% 
  kable(row.names = F, caption= "PFS (Validation set). INCLIVA Score") # for vignette only