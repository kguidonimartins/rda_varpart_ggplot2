#' Create partition table with associated p values for each fraction.
#'
#' @param y  Response matrix
#' @param x1 First exploratory matrix
#' @param x2 Second exploratory matrix
#' @param x3 Third exploratory matrix; can be NULL
#' @param output Table output format; can "markdown" or simple "dataframe" (default) 
#'
#' @return None
#' @export
#'
#' @examples
#' # get data
#' data("mite")
#' data("mite.env")
#' data("mite.pcnm")
#' # process data
#' mite.hell <- decostand(x = mite, method = "hellinger")
#' mite.env.std <- decostand(x = mite.env[, 1:2], method = "standardize")
#' # get parition table and p values of each fraction
#' table_varpart(y = mite.hell, x1 = mite.env.std, x2 = mite.pcnm)
table_varpart <- function(y, x1, x2, x3 = NULL, output = "dataframe") {
  
  require(vegan)
  require(knitr)
  
  if (!is.null(x3)) {
  
  ab.varpart <- varpart(y, x1, x2, x3)
  
  COM <- y
  X1 <- x1
  X2 <- x2
  X3 <- x3
  
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
  
  } else {
  
  ab.varpart <- varpart(y, x1, x2)
  
  COM <- y
  X1 <- x1
  X2 <- x2
  
  fractions <- as.data.frame(rbind(ab.varpart$part$fract,
                                   ab.varpart$part$indfract),
                             row.names = NULL)
  
  fractions$fractions <- rownames(fractions)
  
  rownames(fractions) <- NULL
  
  fractions <- fractions[, c(5, 1, 2, 3, 4)]
  
  p.value <- as.character(rep("NA", nrow(fractions)))
  
  fractions <- as.data.frame(cbind(fractions, p.value))
  
  fractions <- transform(fractions, p.value = as.character(p.value))
  
  # ab.varpart$part$fract[1, ] : [a+b] = X1
  fractions$p.value[1] <-
    anova.cca(rda(COM, X1),
              permutations = 9999)$`Pr(>F)`[1]
  # ab.varpart$part$fract[2, ] : [b+c] = X2
  fractions$p.value[2] <-
    anova.cca(rda(COM, X2),
              permutations = 9999)$`Pr(>F)`[1]
  # ab.varpart$part$fract[3, ] : [a+b+c] = X1+X2
  fractions$p.value[3] <-
    anova.cca(rda(COM, cbind(X1, X2),
                  permutations = 9999))$`Pr(>F)`[1]
  # ab.varpart$part$indfract[1, ] : [a] = X1|X2
  fractions$p.value[4] <-
    anova.cca(rda(COM, X1, X2),
              permutations = 9999)$`Pr(>F)`[1]
  # ab.varpart$part$indfract[3, ] : [c] = X2|X1
  fractions$p.value[6] <-
    anova.cca(rda(COM, X2, X1),
              permutations = 9999)$`Pr(>F)`[1]
  
  }
  
  if (output == "markdown") {
    # Tabela de partição formatada
    markdown <- 
    knitr::kable(
      x = fractions,
      format = "markdown",
      digits = 3,
      caption = "Partition table"
    )  
    return(markdown)
  }
  if (output == "dataframe") {
    return(fractions)
  }
  
}
