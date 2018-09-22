############################################################
#                                                          #
#          Autor: Karlo Gregório Guidoni Martins           #
#            E-mail: kguidonimartins@gmail.com             #
#                                                          #
############################################################

# Baseado em:
# https://scriptsandstatistics.wordpress.com/2018/04/26/how-to-plot-venn-diagrams-using-r-ggplot2-and-ggforce/


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

ipak(packages <- c("vegan",
                   "packfor",
                   "tidyverse",
                   "ggforce"))

############################################################
#                                                          #
#                          Dados                           #
#                                                          #
############################################################

data("mite")
data("mite.env")
data("mite.xy")
data("mite.pcnm")

# Matriz resposta
Y  <- decostand(mite, "hellinger") # https://doi.org/10.1007/s004420100716
# Matriz ambiental
X1 <- decostand(mite.env[, 1:2], "standardize")
# Matriz de coordenadas geográficas
X2 <- mite.xy
# Matriz de filtros espaciais
X3 <- mite.pcnm

############################################################
#                                                          #
#         RDA com seleção de variáveis ambientais          #
#                                                          #
############################################################

rda_01 <- rda(Y ~ ., X1) # modelo saturado
(R2a.all <- RsquareAdj(rda_01)$adj.r.squared)
env.sel <- forward.sel(Y, X1, adjR2thresh = R2a.all)
env.sign <- sort(env.sel$order)
env.red <- X1[, c(env.sign)]
rda_02 <- rda(Y, env.red)

anova.cca(rda_02, step = 1000)
RsquareAdj(rda_02)
vif.cca(rda_02)

############################################################
#                                                          #
#     RDA com seleção de coordenadas geográficas x e y     #
#                                                          #
############################################################

rda_03 <- rda(Y ~ ., X2) # modelo saturado
(R2a.all <- RsquareAdj(rda_03)$adj.r.squared)
xy.sel <- forward.sel(Y, X2, adjR2thresh = R2a.all)
xy.sign <- sort(xy.sel$order)
xy.red <- X2[, c(xy.sign)]
rda_04 <- rda(Y, xy.red)

anova.cca(rda_04, step = 1000)
RsquareAdj(rda_04)
vif.cca(rda_04)

############################################################
#                                                          #
#           RDA com seleção de filtros espaciais           #
#                                                          #
############################################################

rda_05 <- rda(Y ~ ., X3) # modelo saturado
(R2a.all <- RsquareAdj(rda_05)$adj.r.squared)
filt.sel <- forward.sel(Y, X3, adjR2thresh = R2a.all)
filt.sign <- sort(filt.sel$order)
filt.red <- X3[, c(filt.sign)]
rda_06 <- rda(Y, filt.red)

anova.cca(rda_06, step = 1000)
RsquareAdj(rda_06)
vif.cca(rda_06)

############################################################
#                                                          #
#           Partição da variância entre matrizes           #
#                                                          #
############################################################

ab.varpart <- varpart(Y, env.red, xy.red, filt.red)
COM <- Y
X1 <- env.red
X1 %>% names()
X2 <- xy.red
X2 %>% names()
X3 <- filt.red
X3 %>% names()

plot(ab.varpart, digits = 1)

ab.varpart

############################################################
#                                                          #
#        Obtenção do p-valor para frações testáveis        #
#                                                          #
############################################################

fractions <- as.data.frame(rbind(ab.varpart$part$fract,
                                 ab.varpart$part$indfract,
                                 ab.varpart$part$contr1),
                           row.names = NULL)

fractions$fractions <- rownames(fractions)

rownames(fractions) <- NULL

fractions <- fractions[ , c(5, 1, 2, 3, 4)]

p.value <- as.character(rep("NA", nrow(fractions)))

fractions <- as.data.frame(cbind(fractions, p.value))

fractions <- transform(fractions, p.value = as.character(p.value))

# ab.varpart$part$fract[1,] : [a+d+f+g] = X1
fractions$p.value[1] <-
  anova.cca(rda(COM, X1),
            permutations = 9999)$`Pr(>F)`[1]
# ab.varpart$part$fract[2,] : [b+d+e+g] = X2
fractions$p.value[2] <-
  anova.cca(rda(COM, X2),
            permutations = 9999)$`Pr(>F)`[1]
# ab.varpart$part$fract[3,] : [c+e+f+g] = X3
fractions$p.value[3] <-
  anova.cca(rda(COM, X3),
            permutations = 9999)$`Pr(>F)`[1]
# ab.varpart$part$fract[4,] : [a+b+d+e+f+g] = X1+X2
fractions$p.value[4] <-
  anova.cca(rda(COM, cbind(X1, X2),
                permutations = 9999))$`Pr(>F)`[1]
# ab.varpart$part$fract[5,] : [a+c+d+e+f+g] = X1+X3
fractions$p.value[5] <-
  anova.cca(rda(COM, cbind(X1, X3),
                permutations = 9999))$`Pr(>F)`[1]
# ab.varpart$part$fract[6,] : [b+c+d+e+f+g] = X2+X3
fractions$p.value[6] <-
  anova.cca(rda(COM, cbind(X2, X3),
                permutations = 9999))$`Pr(>F)`[1]
# ab.varpart$part$fract[7,] : [a+b+c+d+e+f+g] = All
fractions$p.value[7] <-
  anova.cca(rda(COM, cbind(X1, X2, X3),
                permutations = 9999))$`Pr(>F)`[1]
# ab.varpart$part$indfract[1,] : [a] = X1 | X2+X3
fractions$p.value[8] <-
  anova.cca(rda(COM, X1, cbind(X2, X3)),
            permutations = 9999)$`Pr(>F)`[1]
# ab.varpart$part$indfract[2,] : [b] = X2 | X1+X3
fractions$p.value[9] <-
  anova.cca(rda(COM, X2, cbind(X1, X3)),
            permutations = 9999)$`Pr(>F)`[1]
# ab.varpart$part$indfract[3,] : [c] = X3 | X1+X2
fractions$p.value[10] <-
  anova.cca(rda(COM, X3, cbind(X1, X2)),
            permutations = 9999)$`Pr(>F)`[1]
# ab.varpart$part$contr1[1,] : [a+d] = X1 | X3
fractions$p.value[16] <-
  anova.cca(rda(COM, X1, X3,
                permutations = 9999))$`Pr(>F)`[1]
# ab.varpart$part$contr1[2,] : [a+f] = X1 | X2
fractions$p.value[17] <-
  anova.cca(rda(COM, X1, X2,
                permutations = 9999))$`Pr(>F)`[1]
# ab.varpart$part$contr1[3,] : [b+d] = X2 | X3
fractions$p.value[18] <-
  anova.cca(rda(COM, X2, X3,
                permutations = 9999))$`Pr(>F)`[1]
# ab.varpart$part$contr1[4,] : [b+e] = X2 | X1
fractions$p.value[19] <-
  anova.cca(rda(COM, X2, X1,
                permutations = 9999))$`Pr(>F)`[1]
# ab.varpart$part$contr1[5,] : [c+e] = X3 | X1
fractions$p.value[20] <-
  anova.cca(rda(COM, X3, X1,
                permutations = 9999))$`Pr(>F)`[1]
# ab.varpart$part$contr1[6,] : [c+f] = X3 | X2
fractions$p.value[21] <-
  anova.cca(rda(COM, X3, X2,
                permutations = 9999))$`Pr(>F)`[1]

# Tabela de partição formatada
knitr::kable(x = fractions,
             format = "markdown",
             digits = 3,
             caption = "Partition table")

############################################################
#                                                          #
#                  Plotando os resultados                  #
#                                                          #
############################################################

# Criando o modelo de gráfico

# Relembrando as posições das letras
showvarparts(parts = 3)

# Create a Venn Diagram Model

# begin model;

A <- c(1, 0, 0, 1, 0, 1, 1)
B <- c(0, 1, 0, 1, 1, 0, 1)
C <- c(0, 0, 1, 0, 1, 1, 1)

Counts <- paste0("[", letters[1:7], "]")

# coordinate of the letters
#       [a]   [b]   [c]   [d]   [e]   [f]   [g]
x <- c(-1.4,  1.4,  0.0,  0.0,  0.8, -0.8,  0.0)
y <- c( 1.2,  1.2, -1.0,  1.5,  0.0,  0.0,  0.5)

(df.vdc <- data.frame(A, B, C, Counts, x, y))

(df.venn <- data.frame(x = c(0, 0.866, -0.866),
                       y = c(-0.5, 1, 1)))

ggplot(df.venn) +
  geom_circle(aes(x0 = x,
                  y0 = y,
                  r = 1.5),
              alpha = .01,
              size = 0.5,
              colour = 'black') +
  coord_fixed() +
  theme_void() +
  theme(legend.position = "none") +
  annotate("text",
           x = df.vdc$x,
           y = df.vdc$y,
           label = df.vdc$Counts,
           size = 5) +
  annotate(geom = "text",
           x = c(-2, 2, 0),
           y = c(2.6, 2.6, -2.2),
           label = c("Variable 1",
                     "Variable 2",
                     "Variable 3"),
           size = 6) +
  annotate(geom = "text",
           x = 2,
           y = -2.4,
           label = "Residuals = ")
# end model;

#############################################################
#                                                           #
#     Incluindo os valores das frações no novo gráfico      #
#                                                           #
#############################################################

# Frações
df.vdc$Counts <- round(fractions[8:14, 4], digits = 2)
# Resíduos
var.resid <- round(fractions[15, 4], digits = 2)

ggplot(df.venn) +
  geom_circle(
    aes(x0 = x, y0 = y, r = 1.5),
    alpha = .01,
    size = 0.5,
    colour = 'black'
  ) +
  coord_fixed() +
  theme_void() +
  annotate(
    "text",
    x = df.vdc$x,
    y = df.vdc$y,
    label = df.vdc$Counts,
    size = 5
  ) +
  annotate(
    geom = "text",
    x = c(-2, 2, 0),
    y = c(2.7, 2.7, -2.2),
    label = c("Environment",
              "Linear spatial trend",
              "Spatial Filters"),
    size = 6
  ) +
  annotate(
    geom = "text",
    x = 2,
    y = -2.6,
    label = paste("Residuals =", var.resid)
  ) +
  xlim(c(-4, 4)) +
  ylim(c(-3, 3))

############################################################
#                                                          #
#                    Salvando a figura                     #
#                                                          #
############################################################

ggsave(filename = "varpart_3_expl_matrices.tiff",
       plot = last_plot(),
       dpi = 300)
