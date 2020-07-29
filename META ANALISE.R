#https://bookdown.org/MathiasHarrer/Doing_Meta_Analysis_in_R/correlations.html




### Para estimar media e SD por meio da mediana e IRQ
#http://www.math.hkbu.edu.hk/~tongt/papers/median2mean.html 
#https://home.ubalt.edu/ntsbarsh/business-stat/otherapplets/Pooled.htm
#FORMULA TO COMBINE MEAN
library(meta)

rvmeta <- metacont(N experimental, media experimental, 
                   SD experimental, N controle, 
                   media controle, SD controle, 
                   nome do estudo, data= DATASET, comb.fixed = FALSE,
                   comb.random = TRUE,
                   method.tau = "SJ",
                   hakn = TRUE,
                   prediction = TRUE,
                   sm = "SMD") #calcula fixed-effects-model and random-effects-model
##Acima, usei The Hartung-Knapp-Sidik-Jonkman method para random-fixed-model. É o mais conservador
library(readxl)
DIAS <- read_excel("COVID/DIAS.xlsx",  sheet = "dias")

dias_meta <- metacont(N1, Dias_depois_CP_M,Dias_depois_CP_SD,
                      N2,
                    Dias_antes_CP_M,
                    Dias_antes_CP_SD, 
                   Autor, data= DIAS, comb.fixed = TRUE,
                   comb.random = TRUE,
                   method.tau = "DL",
                   hakn = TRUE,
                   prediction = TRUE,
                   sm = "ROM") 
dias_meta
##forest plot###
forest (dias_meta)
forest(dias_meta,
       sortvar=TE,
       xlim = c(-1.5,0.5),
       rightlabs = c("g","95% CI","weight"),
       leftlabs = c("Author", "N","Mean","SD","N","Mean","SD"),
       lab.e = "Intervention",
       pooled.totals = FALSE,
       smlab = "",
       text.random = "Overall effect",
       print.tau2 = FALSE,
       col.diamond = "blue",
       col.diamond.lines = "black",
       col.predict = "black",
       print.I2.ci = TRUE,
       digits.sd = 2
)

forest(dias_meta,
       layout = "RevMan5",
       digits.sd = 2) #Lsyout da chocrane

forest(dias_meta,
       layout = "JAMA",
       text.predict = "95% PI",
       col.predict = "black",
       backtransf=TRUE) ## layout do JAMA


##meta regressão com volume de plasma
library(metafor)
output.metareg<-metareg(dias_meta,CP_Vol_media, method.tau="SJ")
output.metareg

#com idade
library(metafor)
idade_metareg<-metareg(dias_meta,IdadeCasos, method.tau="SJ")
idade_metareg

#com grupo de gravidade de sintomas
sintomas_metareg<-metareg(dias_meta, N_sintomas_continuos, method.tau="DL")
sintomas_metareg

bubble(sintomas_metareg, xlab = "Number of symptoms")


interaction.model <- rma(n1i=N1, m1i= Dias_depois_CP_M,sd1i=Dias_depois_CP_SD,
                         n2i=N2,
                        m2i= Dias_antes_CP_M,
                        sd2i=Dias_antes_CP_SD,
                         data=DIAS, 
                         method = "SJ",
                        mods = ~ N_sintomas_continuos*IdadeCasos, 
                         test="z")
interaction.model
#ANALISE subgrupo em relação aos dias de internação >10 e <10 dias
DIAS$Group <- as.factor(DIAS$Group)
dias_meta_GRUPO <- metacont(N1, Dias_depois_CP_M, Dias_depois_CP_SD, 
                            N2,
                      Dias_antes_CP_M,
                      Dias_antes_CP_SD, 
                      Autor, data= DIAS, comb.fixed = TRUE,
                      comb.random = TRUE,
                      method.tau = "SJ",
                      hakn = TRUE,
                      prediction = TRUE,
                      sm = "ROM", byvar = GRUPO)

summary(dias_meta_GRUPO)
forest(dias_meta_GRUPO)

##Outra forma de análise de subgrupo
asas <- update.meta(dias_meta, 
            byvar=GRUPO, 
            comb.random = TRUE, 
            comb.fixed = TRUE)

forest(asas)
############ANALISE subgrupo em relação ao numero de sintomas USEI ESSE#####
dias_meta_sintomas <- metacont(N1, Dias_depois_CP_M, Dias_depois_CP_SD, 
                            N2,
                            Dias_antes_CP_M,
                            Dias_antes_CP_SD, 
                            Autor, data= DIAS, comb.fixed = TRUE,
                            comb.random = TRUE,
                            method.tau = "SJ",
                            hakn = TRUE,
                            prediction = TRUE,
                            sm = "ROM", byvar = Clinical_symptomatology)

summary(dias_meta_sintomas)
forest(dias_meta_sintomas, layout = "JAMA",
       text.predict = "95% PI",
       col.predict = "black",
       scientific.pval = gs("scientific.pval"),
       col.diamond = "black",
)


#ANALISE subgrupo em relação ao numero de vezes que tomou soro 
dias_meta_soro <- metacont(N1, Dias_depois_CP_M, Dias_depois_CP_SD, 
                               N2,
                               Dias_antes_CP_M,
                               Dias_antes_CP_SD, 
                               Autor, data= DIAS, comb.fixed = TRUE,
                               comb.random = TRUE,
                               method.tau = "SJ",
                               hakn = TRUE,
                               prediction = TRUE,
                               sm = "SMD", byvar = Numero_soro)

summary(dias_meta_soro)
forest(dias_meta_soro)


#####Meta regressão#
#calculando tamanho de efeito
library(MOTE)

d.ind.t(m1=16.25, m2=27, sd1= 3.4, sd2=12.14, n1=4, n2=4, a=.05)
d.ind.t(m1=8, m2=22, sd1= 2.82, sd2=5.65, n1=2, n2=2, a=.05)

d.ind.t(m1=17.83, m2=10, sd1= 8.99, sd2= 3.84, n1=6, n2=6, a=.05)
d.ind.t(m1=18.2, m2=28.4, sd1= 4.71, sd2=8.04, n1=5,n2=5, a=.05)
d.ind.t(m1=3.28, m2=11.2, sd1= 2.6, sd2=6.63, n1=25,n2=25, a=.05)
d.ind.t(m1=5.6, m2=15.7, sd1= 3.65, sd2=3.94, n1=10, n2=10, a=.05)


SE_dias <- metagen(ES_dias,
        SE_dias,
       studlab=paste(Autor),
       data=DIAS,
        comb.fixed = TRUE,
        comb.random = TRUE,
        prediction=TRUE,
        method.tau = "SJ",
        sm="ROM")
SE_dias
output.metareg <- metareg(SE_dias, CP_Vol_media)
output.metareg
##Binary outcome#####  Mantel-Haenszel method
m.bin <- metabin(Ee,
                 Ne,
                 Ec,
                 Nc,
                 data = binarydata,
                 studlab = paste(Author),
                 comb.fixed = TRUE,
                 comb.random = TRUE,
                 method.tau = "SJ",
                 hakn = TRUE,
                 prediction = TRUE,
                 incr = 0.1,
                 sm = c("RR", "OR"))
m.bin
labbe.metabin(x = m.bin,
              bg = "blue",
              studlab = TRUE,
              col.random = "red")



###### braço uncio####

####morte ####
library(readxl)
obitos <- read_excel("COVID/DIAS.xlsx", sheet = "obitos")
morte  <- metaprop(
  obito_CP,
  N,
  Autor,
  data = obitos, sm= "PFT")
morte
forest(morte)
mtprop <- metaprop(event=obito_CP, n=N, studlab=Autor, data=obitos)




####### PCR ####
library(readxl)
PCR <- read_excel("C:/Users/diego/Downloads/DIAS.xlsx", sheet = "PCR")
PCR_meta <- metacont(N1, C_reactive_POS_M,
                     C_reactive_POS_SD, 
                     N2,C_reactive_PRE_M,C_reactive_PRE_SD,
                     Autor, data= PCR, comb.fixed = TRUE,
                      comb.random = TRUE,
                      method.tau = "DL",
                      hakn = TRUE,
                      prediction = TRUE,
                      sm = "ROM") 
PCR_meta
##forest plot###
forest (PCR_meta)
forest.meta (PCR_meta,         text.predict =NULL,
             col.predict = "black",
             scientific.pval = gs("scientific.pval"),
             prediction = FALSE, col.diamond = "gray", col.study = "black", col.square="black")

##meta regressão com volume de plasma
library(metafor)
PCR_metareg<-metareg(PCR_meta,CP_ml, method.tau="DL")
PCR_metareg

#com idade
library(metafor)
PCR_idade_metareg<-metareg(PCR_meta,IdadeCasos, method.tau="DL")
PCR_idade_metareg

#com grupo de gravidade de sintomas
PCR_sintomas_metareg<-metareg(PCR_meta, Sintomas_CP_M)
PCR_sintomas_metareg
bubble(sintomas_metareg, xlab = "Number of symptoms")
##Gravidade dos sintomas, Idade, volume de CP NEM MEDICAÇÕES não estão associados com os tamanhos de efeito de PCR

#medicaçao
PCR_medicacao_metareg<-metareg(PCR_meta, IdadeCasos + Antiviral
,  method.tau="DL")
PCR_medicacao_metareg

PCR_medicacao_metareg2 <-metareg(PCR_meta, Cloroquina + antibiotico,  method.tau="DL")
PCR_medicacao_metareg2
bubble(PCR_medicacao_metareg2, regline=TRUE, studlab=TRUE,  col = "black", bg = "dodgerblue3")
bubble(PCR_medicacao_metareg2, regline=TRUE, studlab=FALSE,  col = "black", bg = "dodgerblue3", ylim = c(-3.6, -0.3))


####### CARGA VIRAL ####
library(readxl)
CARGA <- read_excel("C:/Users/diego/Downloads/DIAS.xlsx", sheet = "Carga")
CARGA_meta <- metacont(N1, CargaViral_posCP_M,
                     CargaViral_posCP_SD, 
                     N2,CargaViral_preCP_M,CargaViral_preCP_SD,
                     Autor, data= CARGA, comb.fixed = TRUE,
                     comb.random = TRUE,
                     method.tau = "DL",
                     hakn = TRUE,
                     prediction = TRUE,
                     sm = "ROM") 
CARGA_meta
CARGA_meta$I2
CARGA_meta$lower.I2
CARGA_meta$upper.I2
baujat(CARGA_meta)


##forest plot###
forest.meta (CARGA_meta,         text.predict =NULL,
             col.predict = "black",
             scientific.pval = gs("scientific.pval"),
             prediction = FALSE, col.diamond = "gray", col.study = "black", col.square="black")
CARGA_meta$upper.I2
##meta regressão com volume de plasma
library(metafor)
carga_metareg<-metareg(CARGA_meta,IdadeCasos, method.tau="DL")
carga_metareg

#com idade
library(metafor)
carga_idade_metareg<-metareg(CARGA_meta,IdadeCasos, method.tau="DL")
carga_idade_metareg
bubble(carga_idade_metareg, xlab = "Idade")

#com grupo de gravidade de sintomas
carga_sintomas_metareg<-metareg(CARGA_meta, Sintomas)
carga_sintomas_metareg
bubble(carga_sintomas_metareg, xlab = "Number of symptoms")

#medicaçao NAO DEU SIG
CARGA_medicacao_metareg<-metareg(CARGA_meta, Antiviral + antibiotico,  method.tau="DL")
CARGA_medicacao_metareg


##forest plot###
forest (CARGA_meta, layout = "JAMA",
        text.predict = "95% PI",
        col.predict = "black",
        scientific.pval = gs("scientific.pval"),
        col.diamond = "black",
)

##meta regressão com volume de plasma
library(metafor)
carga_metareg<-metareg(CARGA_meta,IdadeCasos, method.tau="DL")
carga_metareg

#com idade
library(metafor)
carga_idade_metareg<-metareg(CARGA_meta,IdadeCasos, method.tau="DL")
carga_idade_metareg
bubble(carga_idade_metareg, xlab = "Idade")

#com grupo de gravidade de sintomas
carga_sintomas_metareg<-metareg(CARGA_meta, Sintomas)
carga_sintomas_metareg
bubble(carga_sintomas_metareg, xlab = "Number of symptoms")

#medicaçao NAO DEU SIG
CARGA_medicacao_metareg<-metareg(CARGA_meta, Antiviral + antibiotico + Cloroquina,  method.tau="DL")
CARGA_medicacao_metareg
bubble(carga_sintomas_metareg, xlab = "Number of symptoms")

####### CARGA VIRAL com Li 2020####
library(readxl)
CARGA2 <- read_excel("DIAS.xlsx", sheet = "carga_nova")
CARGA_meta2 <- metacont(N1, CargaViral_posCP_M,
                       CargaViral_posCP_SD, 
                       N2,CargaViral_preCP_M,CargaViral_preCP_SD,
                       Autor, data= CARGA2, comb.fixed = TRUE,
                       comb.random = TRUE,
                       method.tau = "DL",
                       hakn = TRUE,
                       prediction = TRUE,
                       sm = "ROM") 
CARGA_meta2
CARGA_meta2$I2
CARGA_meta2$lower.I2
CARGA_meta2$upper.I2
baujat(CARGA_meta2)


##forest plot###
forest.meta (CARGA_meta2,         text.predict =NULL,
             col.predict = "black",
             scientific.pval = gs("scientific.pval"),
             prediction = FALSE, col.diamond = "gray", col.study = "black", col.square="black")

library(metafor)
carga_metareg2 <-metareg(CARGA_meta, Antiviral + Corticoide + antibiotico +	Cloroquina, method.tau="DL")
carga_metareg2 
bubble(carga_metareg2, xlab = "Idade")
bubble(carga_sintomas_metareg, xlab = "Number of symptoms")


##### Clínica #####
library(meta)
library(readxl)
Gravidade <- read_excel("C:/Users/diego/Downloads/DIAS.xlsx", sheet = "Gravidade")
Gravidade_meta <- metacont(N1, M_Score_depois,	DP_Score_depois,
                       N2,M_Score_antes,	DP_Score_antes,
                       Autor, data= Gravidade, comb.fixed = TRUE,
                       comb.random = TRUE,
                       method.tau = "DL",
                       hakn = TRUE,
                       prediction = TRUE,
                       sm = "ROM") 
Gravidade_meta
##forest plot###
forest.meta (Gravidade_meta,         text.predict =NULL,
        col.predict = "black",
        scientific.pval = gs("scientific.pval"),
        prediction = FALSE, col.diamond = "gray", col.study = "black", col.square="black")
Gravidade_meta$I2
Gravidade_meta$lower.I2
Gravidade_meta$upper.I2

#medicaçao NAO DEU SIG
Gravidademetareg <-metareg(Gravidade_meta, antibiotico + Antiviral + Corticoide + Cloroquina + IdadeCasos, method.tau="DL")
Gravidademetareg
library(RColorBrewer)
bubble(Gravidademetareg, regline=TRUE, studlab=FALSE,  col = "black", bg = "dodgerblue3", ylim = c(-1.5, 0.5))

##sEM OS ETSUDOS HETEROGENEOS
library(meta)
library(readxl)
Gravidade2 <- read_excel("DIAS.xlsx", sheet = "Gravidade2")
Gravidade_meta2 <- metacont(N1, M_Score_depois,	DP_Score_depois,
                           N2,M_Score_antes,	DP_Score_antes,
                           Autor, data= Gravidade2, comb.fixed = TRUE,
                           comb.random = TRUE,
                           method.tau = "DL",
                           hakn = TRUE,
                           prediction = TRUE,
                           sm = "ROM") 
Gravidade_meta2
forest.meta (Gravidade_meta2,         text.predict =NULL,
             col.predict = "black",
             scientific.pval = gs("scientific.pval"),
             prediction = FALSE, col.diamond = "gray", col.study = "black", col.square="black")


#### Grafico de heterogeneidade
baujat(Gravidade_meta, studlab=TRUE,  col = "black", bg = "dodgerblue3",   col.grid =NULL)

x <- c(121.1,
       5.6,
       5.5,
       11.9,
       3.4,
       1.1,
       20,
       2)

sd(x)
mean(x)
sqrt(mean(x ^ 2) - mean(x)^2)


# Carga viral dic #####

library(readxl)
carga_dic <- read_excel("C:/Users/diego/Downloads/DIAS.xlsx",  sheet = "carga_dic")
library(meta)
m.bin <- metabin(CargaViral_posCP,
                 N2,
                 CargaViral_preCP,
                 N1,
                 data = carga_dic,
                 studlab = paste(Autor),
                 comb.fixed = TRUE,
                 comb.random = TRUE,
                 method.tau = "DL",
                 hakn = TRUE,
                 prediction = TRUE,
                sm = "RR")
m.bin
summary(m.bin)
forest(m.bin,  text.predict =NULL,
       col.predict = "black",
       scientific.pval = gs("scientific.pval"),
       prediction = FALSE, col.diamond = "gray", col.study = "black", col.square="black")
