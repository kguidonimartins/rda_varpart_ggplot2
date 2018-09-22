############################################################
#                                                          #
#          Autor: Karlo Gregório Guidoni Martins           #
#            E-mail: kguidonimartins@gmail.com             # 
#                                                          #
############################################################

# Baseado em: 
# http://renatabrandt.github.io/EBC2015/PCA.html

############################################################
#                                                          #
#                         Pacotes                          #
#                                                          #
############################################################

# ipak function: install and load multiple R packages.
# Check to see if packages are installed. 
# Install them if they are not, then load them into the R session.
# Forked from: https://gist.github.com/stevenworthington/3178163

ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
  {
    install.packages(new.pkg, dependencies = TRUE)
  }
  suppressPackageStartupMessages(sapply(pkg, require, character.only = TRUE))
}

install.packages("packfor", repos = "http://R-Forge.R-project.org")

ipak(packages <- c("vegan", "packfor", "ggplot2", "ggrepel"))

############################################################
#                                                          #
#                          Dados                           #
#                                                          #
############################################################

data("varespec")
data("varechem")

# Transformação e padronização dos dados
Y <- decostand(varespec, "hellinger")
X <- decostand(varechem, "standardize")

############################################################
#                                                          #
#                     RDA com seleção                      #
#                                                          #
############################################################

rda1 <- rda(Y ~ ., X) # modelo saturado
R2a.all <- RsquareAdj(rda1)$adj.r.squared
env.sel <- forward.sel(Y, X, adjR2thresh = R2a.all)
env.sign <-  sort(env.sel$order)
env.red <- X[, c(env.sign)]
rda2 <- rda(Y, env.red)

# Plot para uma rápida visualização
plot(rda2, type = "t")

############################################################
#                                                          #
#                  Plotando os resultados                  #
#                                                          #
############################################################

# Atribua o objeto da ordenação a cca.results
# A partir daqui, o script seguirá um esqueleto baseado no objeto cca.results
cca.results <- rda2

# Definindo objetos-padrão para o plot
cca.summary <- summary(cca.results)
(imp.axis.1 <- cca.summary$cont$importance[2,1]) # explicação do eixo 1
(imp.axis.2 <- cca.summary$cont$importance[2,2]) # explicação do eixo 2
(rowScores <- as.data.frame(cca.results$CCA$u))  # scores das parcelas
(colScores <- as.data.frame(cca.results$CCA$v))  # scores das espécies 
(bi.var <- as.data.frame(cca.results$CCA$biplot)) # scores das variáveis ambientais

# Para salvar a figura em alta resolução.
# Sempre adeque a resolução à altura e à largura da figura.
# tiff(filename="figura.tiff", res=600, height=600/72*600, width=600/72*600, compression= "lzw")


############################################################
#                                                          #
#              PLOTANDO AS UNIDADES AMOSTRAIS              #
#                                                          #
############################################################

ggrda1 <- ggplot() +
  # geom_point(data=rowScores, aes(x=RDA1, y=RDA2),
  #            size=4,
  #            fill = "white",
  #            shape=22) +
  geom_vline(xintercept = 0, alpha=0.9, linetype = "dashed") +
  geom_hline(yintercept = 0, alpha=0.9, linetype = "dashed") +
  geom_point(data = rowScores, aes(x = RDA1, y = RDA2),
             size = 5,
             color = "black", fill="grey50", shape=21) +
  theme(panel.background = element_rect(fill = "white", colour = NA),
        panel.grid.minor = element_line(colour = NULL),
        panel.grid.major = element_line(colour = NULL),
        plot.background = element_rect(fill = "white", colour = NA),
        panel.border = element_rect(fill = NA, colour = "black"),
        axis.text.x = element_text(color ="black", size = 12, angle = 0),
        axis.text.y = element_text(color ="black", size = 12, angle = 0))

# http://www.sthda.com/english/wiki/ggplot2-axis-ticks-a-guide-to-customize-tick-marks-and-labels
# http://www.sthda.com/english/wiki/ggplot2-line-types-how-to-change-line-types-of-a-graph-in-r-software

# Visualizar
ggrda1

############################################################
#                                                          #
#           DEFININDO EIXOS E LEGENDAS DOS EIXOS           #
#                                                          #
############################################################

ggrda2 <- ggrda1 +
  xlab(paste("RDA 1 ", "(", round(imp.axis.1*100, digits = 1), "%)")) +
  theme(axis.title.x=element_text(angle = 0, size = 15)) + #, face = "bold"
  ylab(paste("RDA 2 ", "(", round(imp.axis.2*100, digits = 1), "%)")) +
  theme(axis.title.y=element_text(angle = 90, size = 15)) #, face = "bold"

# Visualizar
ggrda2


############################################################
#                                                          #
#                   PLOTANDO AS ESPÉCIES                   #
#                                                          #
############################################################

# Geralmente as espécies estão sobrepostas. Daí tem uns truques.
# Se as espécies NÂO estiverem sobrepostas, só rode o comando abaixo.
# Se estiverem sobrepostas, bloqueie o comando geom_text e
#libere o geom_point e geom_text_repel.
# Isso vai plotar os pontos das espécies e espalhar os nomes.

ggrda3 <- ggrda2 +
  # SE AS ESPÉCIES ESTIVEREM SOBREPOSTAS, BLOQUEIE O COMANDO ABAIXO:
  # geom_text(data=colScores, aes(RDA1, RDA2),
  #           size=4,
  #           label = rownames(colScores),
  #           fontface="italic")
  
  # SE AS ESPÉCIES NÃO ESTIVEREM SOBREPOSTAS, BLOQUEIE O COMANDO ABAIXO
  geom_point(data = colScores, aes(RDA1, RDA2),
             size = 5, shape = 4)+
  geom_text_repel(aes(x = colScores$RDA1, y = colScores$RDA2,
                      label = rownames(colScores)),
                  size = 5,
                  fontface ="italic",
                  segment.color = "grey",
                  alpha = 0.8,
                  segment.size = 0.1,
                  box.padding = unit(0.3, 'lines'))#

# Visualizar
ggrda3


############################################################
#                                                          #
#            PLOTANDO AS VARIÁVEIS AMBIENTAIS              #
#                                                          #
############################################################

ggrda4 <- ggrda3 +
  geom_segment(data = bi.var, aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               size = 1,
               arrow = arrow(length = unit(0.2, "cm"),
                             type = "closed", angle = 15),
               color = "black",
               alpha = 0.5) +
  geom_text(data = bi.var, aes(RDA1*1.05, RDA2*1.05,
                               label = rownames(bi.var)),
            color = "black",
            size = 5,
            fontface = "bold")
#geom_text_repel(aes(x = colScores$PC1, y = colScores$PC2, label = rownames(colScores)))

# Visualizar
ggrda4

ggrda5 <- ggrda4 + ggtitle("a) Minha RDA no ggplot2") +
  theme(plot.title = element_text(lineheight =.8, face = "bold", size = 16))

ggrda5

# para fechar a figura
# dev.off() 
