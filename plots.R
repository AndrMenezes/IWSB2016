#------------------------------------------------------------------------
#------------------------------------------------------------------------
rm(list=ls())
library(tidyverse)
library(ggplot2)
library(gridExtra)
#------------------------------------------------------------------------
setwd("D:/BETA_TRABALHO/Banco")
dados = read.table('montecarlo10000_FINAL.txt', header = T, sep=",")
head(dados)
names(dados) <- c("teste", "alpha", "n", "mu", "phi")
dados$teste <- factor(dados$teste,
                      levels = c("Table trv", "Table wald",
                                 "Table escore"),
                      labels = c("TRV", "Wald", "Escore"))  
dados$mu <- factor(dados$mu, levels = c(0.05, 0.25, 0.5, 0.75,0.95),
                   labels = c(expression(paste(mu, " = 0.05")),
                              expression(paste(mu, " = 0.25")),
                              expression(paste(mu, " = 0.5")),
                              expression(paste(mu, " = 0.75")),
                              expression(paste(mu, " = 0.95"))))
dados %>% arrange(teste)
#------------------------------------------------------------------------
legenda <- theme(legend.position = c(0.7, 0.15),
        plot.title = element_text(hjust = 0.5))
#------------------------------------------------------------------------
tema <- theme(legend.position="none") +
        theme(axis.text.x = element_text(hjust = 1,color="black"),
        axis.text.y = element_text(color="black"),
        strip.background = element_rect(fill="black"),
        strip.text.x = element_text(size=13, colour="white",face="bold"))
#------------------------------------------------------------------------
p2<-dados%>% filter(phi==15)%>%
  ggplot(aes(n, alpha, col = teste))+
  facet_wrap(~mu,labeller = label_parsed, nrow = 3) +
  geom_line()+geom_point()+ 
  scale_colour_brewer(palette = "Set1",name="Testes de\n Hipótese")+
  scale_y_continuous(breaks=seq(0,0.13,0.01))+
  geom_hline(yintercept = 0.05,linetype="dashed")+
  labs(x="\nTamanho da amostra", y="Taxa do Erro Tipo I\n", 
  title=expression(paste(phi, " = 15"))) +
  theme(legend.position="none") + tema
phi15_simul<-p2 + legenda
dir<-"D:/BETA_TRABALHO/Banco/phi15_simul.pdf"
ggsave(dir, plot = phi15_simul, width = 9, height = 6, dpi = 350)
#------------------------------------------------------------------------
p3<-dados%>% filter(phi==50)%>%
  ggplot(aes(n, alpha, col = teste))+
  facet_wrap(~mu,labeller = label_parsed, nrow = 3) +
  geom_line()+geom_point()+ 
  scale_colour_brewer(palette = "Set1",name="Testes de\n Hipótese")+
  scale_y_continuous(breaks=seq(0,0.13,0.01))+
  geom_hline(yintercept = 0.05,linetype="dashed")+
  labs(x="\nTamanho da amostra", y="Taxa do Erro Tipo I\n", 
       title=expression(paste(phi, " = 50"))) + tema
phi50_simul<-p3 + legenda
phi50_simul
dir<-"D:/BETA_TRABALHO/Banco/phi50_simul.pdf"
ggsave(dir, plot = phi50_simul, width = 9, height = 6, dpi = 350)
#------------------------------------------------------------------------
#------------------------------------------------------------------------

