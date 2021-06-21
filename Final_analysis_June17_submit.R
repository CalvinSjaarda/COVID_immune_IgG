library("dplyr")
library("ggplot2")
library("ggpubr")
library("lme4")
library("gridExtra")

#INPUT DATA
Input <- read.csv("~/Projects/Collaborations/Sheth_immunology/IgG_NT_dataset_Jun11.csv")
Input$Sex <- as.factor(as.character(Input$Sex))
Input$Donor <- as.factor(as.character(Input$Donor))
Input$NTassay <- as.factor(as.character(Input$NTassay))
Input$Age.bin <- as.factor(as.character(Input$Age.bin))

#############################
#Figure 1: Age and Sex
#SEX
SEX_data <- Input %>% group_by(Donor) %>% slice_head(n = 1) %>% filter(!is.na(Sex)) %>% group_by(Donor,Sex) %>% summarise(IgG_YHLO)

colnames(SEX_data) <- c("Donor","Sex", "iFlash")
x_labels <- c("Female", "Male")
my_comparisons <- list( c("m", "f"))

FIG1A <- ggplot(SEX_data, aes(x=Sex, y=iFlash)) +
  geom_dotplot(binwidth=2,binaxis="y",stackdir="center") +
  theme_classic() +
  ylab("Anti-SARS-CoV-2 IgG level\n(AU/mL)") +
  geom_hline(yintercept = c(10), color = "black", linetype="dashed")+
  stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", width = 0.7, colour = "red") +
  scale_y_continuous(breaks = round(seq(min(0), max(140), by = 20),1)) +
  stat_compare_means(comparisons = my_comparisons, size = 0) +
  stat_compare_means(label = "p.signif", label.y = 152, label.x = 1.5) +
  scale_x_discrete(labels= x_labels) +
  theme(axis.title.x=element_blank())

#AGE
AGE_data <- Input %>% group_by(Donor) %>% slice_head(n = 1) %>% filter(!is.na(Age.bin)) %>% filter(!is.na(IgG_YHLO)) %>% group_by(Donor,Age.bin) %>% summarise(IgG_YHLO)
colnames(AGE_data) <- c("Donor","Age", "iFlash")
x_labels <- c("<20","20-29","30-39","40-49","50-59","60+")

FIG1B <- ggplot(AGE_data, aes(x=Age, y=iFlash)) +
  geom_dotplot(binwidth=2,binaxis="y",stackdir="center") +
  theme_classic() +
  geom_hline(yintercept = c(10), color = "black", linetype="dashed")+
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.7, colour = "red") +
  scale_y_continuous(breaks = round(seq(min(0), max(140), by = 20),1)) +
  stat_compare_means(label.y = 150, label.x=2) +
  scale_x_discrete(labels= x_labels) +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank())

FIG1final <- ggarrange(FIG1A, FIG1B, ncol = 2, nrow = 1, labels = c("A","B"))
FIG1final
ggsave("Figure1final.tiff", plot = FIG1final, units="in",width=5, height=3, dpi=300, compression = 'lzw')

compare_means(iFlash ~ Age,  data = AGE_data)

#Data for Table 1
summary(SEX_data) #male:female
Input %>% group_by(Donor) %>% slice_head(n = 1) %>% group_by(Sex,Age.bin) %>% summarise(n = n()) #Age.bin
Input %>% group_by(Donor) %>% slice_head(n = 1) %>% group_by(Sex,Age.bin,NTassay) %>% count(n()) %>% filter(!is.na(NTassay)) #NTassay
Input %>% group_by(Donor) %>% slice_head(n = 1) %>% filter(!is.na(IgG_YHLO)) %>% group_by(Sex,Age.bin) %>% summarise(n = n(), median = median(IgG_YHLO), min = min(IgG_YHLO), max = max(IgG_YHLO)) #YHLO IgG levels
Input %>% group_by(Donor) %>% slice_head(n = 1) %>% filter(!is.na(IgG_Abbott)) %>% group_by(Sex,Age.bin) %>% summarise(median = median(IgG_Abbott), min = min(IgG_Abbott), max = max(IgG_Abbott)) #Abbott IgG levels
Input %>% group_by(Donor) %>% filter(IgG_YHLO == max(IgG_YHLO))%>% filter(!is.na(Time)) %>% group_by(Sex,Age.bin) %>% summarise(min = min(Time), max = max(Time)) #Time to reach max IgG levels


####################################
#LINEAR MIXED MODEL
STATS_data <- tidyr::gather(Input,"Assay","AB",5:6)
STATS_input <- STATS_data %>% filter(!is.na(AB)) %>% add_count(Donor, Assay) %>% group_by(Assay, Donor) %>% filter(n>2)
STATS_input$Assay <- as.factor(as.character(STATS_input$Assay))

#MLL8 was the model that best fit the data
MLL8 = lmer(log(AB) ~ 1 + log(Time) + Age + Assay + log(Time):Age + log(Time):Assay + (1 | Donor) + (0 + log(Time) | Donor),data = STATS_input, REML=TRUE)

#Random donor intercept
MLL8_donor = lmer(log(AB) ~ 1 + log(Time) + Age + Assay + log(Time):Age + log(Time):Assay + (0 + log(Time) | Donor),data = STATS_input, REML=TRUE)
anova(MLL8,MLL8_donor)
#Random donor slope
MLL8_slope = lmer(log(AB) ~ 1 + log(Time) + Age + Assay + log(Time):Age + log(Time):Assay + (1 | Donor),data = STATS_input, REML=TRUE)
anova(MLL8,MLL8_slope)
#Main effect of Time
MLL8_time = lmer(log(AB) ~ 1 + Age + Assay + log(Time):Age + log(Time):Assay + (1 | Donor) + (0 + log(Time) | Donor),data = STATS_input, REML=TRUE)
anova(MLL8,MLL8_time)
#Main effect of Assay
MLL8_assay = lmer(log(AB) ~ 1 + log(Time) + Age + log(Time):Age + log(Time):Assay + (1 | Donor) + (0 + log(Time) | Donor),data = STATS_input, REML=TRUE)
anova(MLL8,MLL8_assay)
#Interaction between time and age
MLL8_time_age = lmer(log(AB) ~ 1 + log(Time) + Age + Assay + log(Time):Assay + (1 | Donor) + (0 + log(Time) | Donor),data = STATS_input, REML=TRUE)
anova(MLL8,MLL8_time_age)
#Interaction between time and assay
MLL8_time_assay = lmer(log(AB) ~ 1 + log(Time) + Age + Assay + log(Time):Age + (1 | Donor) + (0 + log(Time) | Donor),data = STATS_input, REML=TRUE)
anova(MLL8,MLL8_time_assay)


################################
#Interpolation of decay below LOD
DECAY_iFlash <- STATS_input %>% filter(Assay == "IgG_YHLO") %>% group_by(Donor) %>% filter(!is.na(Time)) %>% mutate(Decay = (exp(predict(lm(log(Time)~AB), newdata = data.frame(AB = 10))))) %>% slice_head(n = 1) %>% filter(Decay < 400) %>% select(Donor,Sex,Age,Assay,Decay)
DECAY_Architect <- STATS_input %>% filter(Assay == "IgG_Abbott") %>% group_by(Donor) %>% filter(!is.na(Time)) %>% mutate(Decay = (exp(predict(lm(log(Time)~AB), newdata = data.frame(AB = 1.4))))) %>% slice_head(n = 1) %>% filter(Decay < 400) %>% select(Donor,Sex,Age,Assay,Decay)

#Use t.text to mean # of days for IgG levels to decay below LOD
t.test(as.numeric(DECAY_iFlash$Decay))
t.test(as.numeric(DECAY_Architect$Decay))

#FIGURE S5: Density figure
DECAY_data <- rbind(DECAY_iFlash, DECAY_Architect)
FIGS5 <- ggplot(DECAY_data, aes(Decay)) +
  geom_density(aes(fill=factor(Assay)), alpha=0.6) +
  xlim(0,400) +
  theme_classic() +
  xlab("Time since infection (# of days)") +
  ylab("Density\n(Proportion of samples whose IgG\nlevels dropped below the LOD)") +
  scale_fill_discrete(name = "Assay", labels = c("Achitect", "iFlash")) +
  theme(
    legend.position = c(1, 1),
    legend.justification = c("right", "top"),
    legend.key = element_rect(colour = "black"),
  )
FIGS5
ggsave("FigureS5final.tiff", plot = FIGS5, units="in",width=3, height=3, dpi=300, compression = 'lzw')

##################################
#FIGURES PLOTTING IgG DECAY OVER TIME
TIME_data <- STATS_input %>% group_by(Assay, Donor) %>%mutate(ABlog=log(AB)-mean(log(AB))) #log transform AB
scaleFUN <- function(x) sprintf("%.0f", x)

#FIGURE S3: IgG decay rates over time
TIME_data_YHLO <- TIME_data %>% filter(Assay == "IgG_YHLO")
TIME_data_Abbott <- TIME_data %>% filter(Assay == "IgG_Abbott")

FIGS3A <- ggplot(TIME_data_YHLO, aes(x=Time, y=ABlog)) +
  geom_point(shape=21, colour = "dark grey",size = 1, alpha = 1) +
  geom_smooth(method = "lm", colour = "black") +
  stat_regline_equation(label.y = 1.2, label.x = 5, size = 4) +
  theme_classic() +
  xlab("Time since SARS-CoV-2 infection (# of days)") +
  ylab("Anti-SARS-CoV-2 IgG antibody levels \niFlash (Normalized AU/mL)") +
  ylim(-1.5,1.5) +
  scale_x_continuous(trans = 'log', labels=scaleFUN)
FIGS3B <- ggplot(TIME_data_Abbott, aes(x=Time, y=ABlog)) +
  geom_point(shape=21, colour = "dark grey",size = 1, alpha = 1) +
  geom_smooth(method = "lm", colour = "black") +
  stat_regline_equation(label.y = 1.2, label.x = 5, size = 4) +
  theme_classic() +
  xlab("Time since SARS-CoV-2 infection (# of days)") +
  ylab("Anti-SARS-CoV-2 IgG antibody levels \nArchitect (Normalized Index S/C)") +
  ylim(-1.5,1.5) +
  scale_x_continuous(trans = 'log', labels=scaleFUN)
FIGS3final <- ggarrange(FIGS3A, FIGS3B, ncol = 2, nrow = 1)
FIGS3final
ggsave("FigureS3final.tiff",plot = FIGS3final, units="in",width=8, height=4, dpi=300, compression = 'lzw')

#Sample size for each assay
TIME_data_YHLO %>% group_by(Donor) %>% count(n())
TIME_data_Abbott %>% group_by(Donor) %>% count(n())


#FIGURE S4: Sex
#Sample size for each assay
TIME_data_YHLO %>% group_by(Donor,Sex) %>% count(n()) %>% group_by(Sex) %>% count(n())
TIME_data_Abbott %>% group_by(Donor,Sex) %>% count(n()) %>% group_by(Sex) %>% count(n())

TIME_data_YHLO$Sex <- factor(TIME_data_YHLO$Sex, levels = c("f","m"), labels = c("Female (n=142)","Male (n=170)"))
TIME_data_Abbott$Sex <- factor(TIME_data_Abbott$Sex, levels = c("f","m"), labels = c("Female (n=75)","Male (n=95)"))

FIGS4A <- ggplot(TIME_data_YHLO, aes(x=Time, y=ABlog, colour=Sex, fill=Sex)) +
  geom_point(shape=21, colour = "black",size = 1, alpha = 1) +
  geom_smooth(method = "lm") +
  theme_classic() +
  xlab("Time since SARS-CoV-2 infection\n(# of days)") +
  ylab("Anti-SARS-CoV-2 IgG antibody levels\niFlash (Normalized AU/mL)") +
  ylim(-1.5,1.5) +
  scale_x_continuous(trans = 'log', labels=scaleFUN) +
    theme(
    legend.position = c(0.95, 0.95),
    legend.justification = c("right", "top"),
    legend.key = element_rect(colour = "black"),
    legend.title = element_blank())
FIGS4B <- ggplot(TIME_data_Abbott, aes(x=Time, y=ABlog, colour=Sex, fill=Sex)) +
    geom_point(shape=21, colour = "black",size = 1, alpha = 1) +
    geom_smooth(method = "lm") +
    theme_classic() +
    xlab("Time since SARS-CoV-2 infection\n(# of days)") +
    ylab("Anti-SARS-CoV-2 IgG antibody levels\nArchitect (Normalized Index S/C)") +
    ylim(-1.5,1.5) +
    scale_x_continuous(trans = 'log', labels=scaleFUN) +
    theme(
      legend.position = c(0.95, 0.95),
      legend.justification = c("right", "top"),
      legend.key = element_rect(colour = "black"),
      legend.title = element_blank())
FIGS4final <- ggarrange(FIGS4A, FIGS4B, ncol = 2, nrow = 1)
FIGS4final
ggsave("FigureS4final.tiff",plot = FIGS4final, units="in",width=8, height=4, dpi=300, compression = 'lzw')

###############################
#FIGURE 2: Age binned IgG decay 
#Sample size for each assay
TIME_data_YHLO %>% group_by(Donor,Age.bin) %>% count(n()) %>% group_by(Age.bin) %>% count(n())
TIME_data_Abbott %>% group_by(Donor,Age.bin) %>% count(n()) %>% group_by(Age.bin) %>% count(n())

TIME_data_YHLO$Age.bin <- factor(TIME_data_YHLO$Age.bin, levels = c("<20","20-29","30-39","40-49","50-59","60+"), labels = c("<20 (n=15)","20-29 (n=119)","30-39 (n=72)","40-49 (n=56)","50-59 (n=40)","60+ (n=10)"))
TIME_data_Abbott$Age.bin <- factor(TIME_data_Abbott$Age.bin, levels = c("<20","20-29","30-39","40-49","50-59","60+"), labels = c("<20 (n=11)","20-29 (n=61)","30-39 (n=44)","40-49 (n=29)","50-59 (n=22)","60+ (n=3)"))
FIG2A <- ggplot(TIME_data_YHLO, aes(x=Time, y=ABlog, colour=Age.bin, fill=Age.bin)) +
  geom_point(shape=21, colour = "black",size = 1, alpha = 1) +
  geom_smooth(method = "lm") +
  theme_classic() +
  stat_regline_equation(label.y = 1.7, label.x = 4.7, size = 2.5) +
  xlab("Time since SARS-CoV-2 infection\n(# of days)") +
  ylab("Anti-SARS-CoV-2 IgG level\niFlash (Normalized AU/mL)") +
  scale_x_continuous(trans = 'log', labels=scaleFUN) +
  theme(legend.position = "none") +
  ylim(-1.5,1.75) +
  facet_wrap(~Age.bin, ncol = 3)
FIG2B <- ggplot(TIME_data_Abbott, aes(x=Time, y=ABlog, colour=Age.bin, fill=Age.bin)) +
  geom_point(shape=21, colour = "black",size = 1, alpha = 1) +
  geom_smooth(method = "lm") +
  theme_classic() +
  stat_regline_equation(label.y = 1.7, label.x = 4.7, size = 2.5) +
  xlab("Time since SARS-CoV-2 infection\n(# of days)") +
  ylab("Anti-SARS-CoV-2 IgG level\nArchitect (Normalized Index S/C)") +
  ylim(-1.5,1.75) +
  scale_x_continuous(trans = 'log', labels=scaleFUN) +
  theme(legend.position = "none") +
  facet_wrap(~Age.bin, ncol = 3)

FIGURE2final <- ggarrange(FIG2A, FIG2B, ncol = 1, nrow = 2, labels = c("A","B"))
FIGURE2final
ggsave("Figure2final.tiff",plot = FIGURE2final, units="in",width=6, height=8, dpi=300, compression = 'lzw')

#######################
#FIGURE 3: NT assay
NT_data <- Input %>% group_by(Donor) %>% slice_head(n = 1) %>% filter(!is.na(NTassay)) %>% filter(Time < 101) %>% select('Donor', 'NTassay','IgG_YHLO', 'Time')
summary(NT_data)
#NUMBER OF SAMPLES WITH IgG BELOW LOD
NT_data %>% filter(IgG_YHLO <10) %>% group_by(NTassay) %>% summarize(n=n())

#FIGURE 3 image
x_labels <- c("NT negative", "NT postive")
my_comparisons <- list( c("NEG", "POS"))

FIG3final <- ggplot(NT_data, aes(x=NTassay, y=IgG_YHLO)) +
  geom_dotplot(binwidth=2,binaxis="y",stackdir="center") +
  theme_classic() +
  ylab("Anti-SARS-CoV-2 IgG antibody levels\n(AU/mL)") +
  geom_hline(yintercept = c(10), color = "black", linetype="dashed")+
  stat_summary(fun = median, fun.min = median, fun.max = median, colour = "red", geom = "crossbar", width = 0.7) +
  scale_y_continuous(breaks = round(seq(min(0), max(155), by = 20),1)) +
  stat_compare_means(comparisons = my_comparisons, size = 3, label = "p.format") +
  scale_x_discrete(labels= x_labels) +
  theme(axis.title.x=element_blank())
FIG3final
ggsave("Figure3final.tiff", plot = FIG3final, units="in",width=3, height=3, dpi=300, compression = 'lzw')

wilcox.test(IgG_YHLO ~ NTassay, data=NT_data)

#SUPPLEMENTARY FIGURE 1
donorY <- TIME_data_YHLO %>%
  filter(Donor %in% c("13921","13964","14124","13916","14094","14153","14028",
                      "14090","14192","14275","14171","13994","14288","13923","14084"))
donorA <- TIME_data_Abbott %>%
  filter(Donor %in% c("13921","13964","14124","13916","14094","14153","14028",
                      "14090","14192","14275","14171","13994","14288","13923","14084"))

#FIGURE S1: Donor plot for iFlash
FIGS1A <- ggplot(donorY, aes(x=Time, y=AB, color = Donor)) + 
  geom_point(shape=21,size = 2) +
  geom_smooth(fill = "transparent") +
  theme_classic() +
  xlab("Time since SARS-CoV-2 infection\n(# of days)") +
  ylab("Anti SARS-CoV-2 IgG level\niFlash (AU/mL)") +
  facet_wrap(~Donor, nrow = 3) +
  ylim(0,120) +
  geom_hline(yintercept = c(10), color = "grey", linetype="dashed")+
  theme(legend.position = "none")
FIGS1B <- ggplot(donorY, aes(x=Time, y=ABlog, color = Donor)) +
  geom_point(shape=21, size = 2) +
  geom_smooth(fill = "transparent") +
  geom_smooth(fill = "transparent", colour = "black",size = 1, method = "lm", linetype="dashed") +
  theme_classic() +
  xlab("Time since SARS-CoV-2 infection\n(# of days)") +
  ylab("Anti SARS-CoV-2 IgG level\niFlash (Normalized AU/mL)") +
  facet_wrap(~Donor, nrow = 3) +
  ylim(-2,2) +
  scale_x_continuous(trans = 'log', labels=scaleFUN, ) +
  theme(legend.position = "none")

FIGURES1final <- ggarrange(FIGS1A, FIGS1B, ncol = 1, nrow = 2, labels = c("A","B"))
FIGURES1final
ggsave("FigureS1final.tiff",plot = FIGURES1final, units="in",width=6, height=6, dpi=300, compression = 'lzw')

#FIGURE S2: Donor plot for Architct
FIGS2A <- ggplot(donorA, aes(x=Time, y=AB, color = Donor)) + 
  geom_point(shape=21,size = 2) +
  geom_smooth(fill = "transparent") +
  theme_classic() +
  xlab("Time since SARS-CoV-2 infection\n(# of days)") +
  ylab("Anti SARS-CoV-2 IgG level\nArchitect (Index S/C)") +
  facet_wrap(~Donor, nrow = 3) +
  ylim(0,10) +
  geom_hline(yintercept = c(1.4), color = "grey", linetype="dashed")+
  theme(legend.position = "none")
FIGS2B <- ggplot(donorA, aes(x=Time, y=ABlog, color = Donor)) +
  geom_point(shape=21, size = 2) +
  geom_smooth(fill = "transparent") +
  geom_smooth(fill = "transparent", colour = "black",size = 1, method = "lm", linetype="dashed") +
  theme_classic() +
  xlab("Time since SARS-CoV-2 infection\n(# of days)") +
  ylab("Anti SARS-CoV-2 IgG level\nArchitect (Normalized Index S/C)") +
  facet_wrap(~Donor, nrow = 3) +
  ylim(-2,2) +
  scale_x_continuous(trans = 'log', labels=scaleFUN, ) +
  theme(legend.position = "none")

FIGURES2final <- ggarrange(FIGS2A, FIGS2B, ncol = 1, nrow = 2, labels = c("A","B"))
FIGURES2final
ggsave("FigureS2final.tiff",plot = FIGURES2final, units="in",width=6, height=6, dpi=300, compression = 'lzw')
