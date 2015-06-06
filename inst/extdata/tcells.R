library(bninfo)
data(tcells)
tcell <- tcells$raw

library(signalgraph)
paste("[PKC][PKA|PKC][Raf|PKC:PKA][Mek|PKC:PKA:Raf]",
                      "[Erk|Mek:PKA][Akt|Erk:PKA][P38|PKC:PKA]",
                      "[Jnk|PKC:PKA][Plcg][PIP3|Plcg][PIP2|Plcg:PIP3]") %>%
  {bnlearn::model2network(.)} %>%
  {bnlearn::as.graphNEL(.)} %>%
  igraph.from.graphNEL 
