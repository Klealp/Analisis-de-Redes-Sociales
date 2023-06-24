library(igraph)
library(ggraph)
library(xtable)
# grilla de ggplot
library(gridExtra) 
library(egg)
library(gtable)
library(grid)
library(readr)
#modelos
library(blockmodels)
library(ergm)
library(network)
library(latentnet)

# Tema para los gráficos ------------------------------------------------------

tema <- theme(panel.background = element_blank(), 
              legend.title = element_text(color = "#FFA54F", size = 10, face = 2),
              plot.title = element_text(face = "bold", 
                                        colour = "tomato1",
                                        size = 15,                     # Tamaño de la fuente
                                        hjust = 0.5,                   # Ajuste horizontal
                                        vjust = 3,                     # Ajuste vertical
                                        lineheight = 1,                # Espacio entre líneas
                                        margin = margin(20, 0, 0, 0)), # Márgenes (t, r, b, l)
              plot.title.position = "plot")

### Primera parte -------------
### Punto 1 ----

work <- read_csv("work-edges.csv")
alliance <- read_csv("alliance-edges.csv")
nodes <- read_csv("nodes.csv")

# Creación grafos
g_work <- graph_from_data_frame(d =work, directed = F, vertices = nodes$Id)
g_alliance <- graph_from_data_frame(d =alliance, directed = F, vertices = nodes$Id)

# grafo inducido por la componente conexa work
V(g_work)$cluster <- clusters(graph = g_work)$membership
gcc_work <- induced_subgraph(graph = g_work, vids = which(V(g_work)$cluster == which.max(clusters(graph = g_work)$csize))) 

set.seed(21)
layout_1 <- create_layout(gcc_work, layout = "igraph", algorithm="dh") 

ca<-rep("turquoise",vcount(gcc_work))
ca[c(1,53,145,151,154)]<-"#CA1429"
ca[c(89,94,98,116,205,329)]<-"#5D478B"

lab<-rep(NA,vcount(gcc_work))
lab[c(1,53,89,94,98,116,145,151,154,205,329)]<-c("Uribe","Santos","A. Ordóñez","Peñalosa",
                                                 "Duque","V. Lleras","Pastrana","J.C. Restrepo",
                                                 "C. Gaviria","Petro","Fajardo")
Grado <- log(degree(gcc_work))

windows(16,12)
ggraph(layout_1)+
  geom_edge_link(edge_colour = "#CDC5BF", alpha=0.4)+
  geom_node_point(aes(size = Grado/3), color= ca, alpha=0.7)+
  geom_node_label(aes(label = lab),label.size = 0,label.padding = unit(0.2, "lines"),nudge_y = 7.5,size=4)+
  ggtitle("Red de relaciones de trabajo")+labs(size = "Grado")+
  tema +
  theme(legend.position = "none")

# grafo inducido por la componente conexa alliance
V(g_alliance)$cluster <- clusters(graph = g_alliance)$membership
gcc_alliance <- induced_subgraph(graph = g_alliance, vids = which(V(g_alliance)$cluster == which.max(clusters(graph = g_alliance)$csize))) 

#which(names(V(gcc_alliance))=="andres-pastrana-arango")

set.seed(21)
layout_a <- create_layout(gcc_alliance, layout = "igraph", algorithm="dh")

ca_a<-rep("turquoise",vcount(gcc_alliance))
ca_a[c(1,25,81,127,176)]<-"#CA1429"
ca_a[c(39,54,60,93,102,175)]<-"#5D478B"

lab_a<-rep(NA,vcount(gcc_alliance))
lab_a[c(1,25,39,54,60,81,93,102,127,175,176)]<-c("Uribe","Santos","A. Ordóñez","Peñalosa",
                                           "Duque","Vargas Lleras","Pastrana","C. Gaviria",
                                           "Petro","J.C. Restrepo","Fajardo")

Grado_al <- log(degree(gcc_alliance))

windows(16,12)
ggraph(layout_a)+
  geom_edge_link(edge_colour ="#CDC5BF", alpha=0.4)+
  geom_node_point(aes(size = Grado_al), color= ca_a, alpha=0.7)+
  geom_node_label(aes(label = lab_a),label.size = 0,label.padding = unit(0.2, "lines"),nudge_y = 5.5,size=4.2)+
  ggtitle("Red de relaciones de alianza")+labs(size = "Grado")+
  tema +
  theme(legend.position = "none")

#Análisis descriptivo -------------------------------

A_DescriptivoP1 <- function(Lista){
  
  #Calculo de estadísticas estructurales
  Densidad      <- sapply(X = Lista, FUN = edge_density)
  Coef_Agrupa   <- sapply(X = Lista, FUN = transitivity)
  Asorta <- sapply(X = Lista, FUN = assortativity_degree)
  Dist_Geod     <- sapply(X = Lista, FUN = mean_distance)
  Comp_Gigante  <- sapply(X = Lista, FUN = function(X) (max(sapply(X = decompose(X), FUN = vcount))))
  Grado    <- sapply(X = Lista, FUN = function(X) mean(degree(X)))
  Desv.est <- sapply(X = Lista, FUN = function(X) sd(degree(X)))
  Tabla         <- as.data.frame(cbind(Densidad,Coef_Agrupa,Asorta,
                                       Dist_Geod,Comp_Gigante,Grado,Desv.est))
  return(Tabla)
}

#Simplificación de las redes
gcc_work <- simplify(gcc_work)
gcc_alliance <- simplify(gcc_alliance)

tab<-A_DescriptivoP1(list(gcc_work,gcc_alliance))
rownames(tab)<-c("Trabajo","Alianza")

xtable::xtable(tab,digits = 3)

#Nodos en común para ambas componentes 

CW<-as.data.frame(V(gcc_work)$name)
CA<-as.data.frame(V(gcc_alliance)$name)
colnames(CW)<-colnames(CA)<-"NAME"
nodes<-as.data.frame(nodes)

com<-sqldf::sqldf("select Id
      from nodes
      where Id in (select NAME from CW) and Id in (select NAME from CA) ")

#Medidas de centralidad para Santos y Uribe

dc <- round(degree          (graph = gcc_work, normalized = T),2)
cc <- round(closeness       (graph = gcc_work, normalized = T),2)
bc <- round(betweenness     (graph = gcc_work, normalized = T),2)
ec <- round(eigen_centrality(graph = gcc_work, scale = T)$vector,2)

dc2 <- round(degree          (graph = gcc_alliance, normalized = T),2)
cc2 <- round(closeness       (graph = gcc_alliance, normalized = T),2)
bc2 <- round(betweenness     (graph = gcc_alliance, normalized = T),2)
ec2 <- round(eigen_centrality(graph = gcc_alliance, scale = T)$vector,2)


xtable::xtable(rbind(
  c("Grado","Cercanía","Intermediación","C. Propia","Grado","Cercanía","Intermediación","C. Propia"),
  c(dc[1],cc[1],bc[1],ec[1],dc2[1],cc2[1],bc2[1],ec2[1]),
  c(dc[53],cc[53],bc[53],ec[53],dc2[25],cc2[25],bc2[25],ec2[25])), digits = 3)

### Punto 2 -----
#### Ajuste del modelo de bloques aleatorios----

# matriz de adyacencia
MW <- as.matrix(igraph::as_adjacency_matrix(gcc_work))
dim(MW)
MA <- as.matrix(igraph::as_adjacency_matrix(gcc_alliance))
dim(MA)

### Formulación de los modelos ----
set.seed(73)
work.sbm <- blockmodels::BM_bernoulli(membership_type = "SBM_sym", adj = MW, verbosity = 0, plotting = "")
alli.sbm <- blockmodels::BM_bernoulli(membership_type = "SBM_sym", adj = MA, verbosity = 0, plotting = "")

# Estimación 
work.sbm$estimate()
alli.sbm$estimate()

# num. de grupos optimo
Q <- which.max(work.sbm$ICL)
Q2 <- which.max(alli.sbm$ICL)

#Gráficos número de grupos

windows(10,6)
par(mfrow = c(1,2), mgp = c(3,0.7,0))
plot(work.sbm$ICL, xlab = "Q", ylab = "ICL", type = "b", pch = 16,main = "Trabajo",col.main="#00868B")
lines(x = c(Q,Q), y = c(min(work.sbm$ICL),max(work.sbm$ICL)), col = "red", lty = 2)
plot(alli.sbm$ICL, xlab = "Q", ylab = "ICL", type = "b", pch = 16,main="Alianza",col.main="#00868B")
lines(x = c(Q2,Q2), y = c(min(alli.sbm$ICL),max(alli.sbm$ICL)), col = "red", lty = 2)

# probabilidades estimadas de pertenencia a las comunidades
Z <- work.sbm$memberships[[Q]]$Z
Z2 <- alli.sbm$memberships[[Q2]]$Z

# asignaciones
labs <- apply(X = Z, MARGIN = 1, FUN = which.max)
labs2 <- apply(X = Z2, MARGIN = 1, FUN = which.max)

# Tamaño y probabilidades de los grupos 
##Work
alpha <- table(labs)/vcount(gcc_work)
##Alliance
alpha2 <- table(labs2)/vcount(gcc_alliance)

xtable::xtable(cbind(rbind(table(labs),alpha),rbind(table(labs2),alpha2)))

#Matriz de probabilidades de interaccion
Pi <- work.sbm$model_parameters[[Q]]$pi
Pi2 <- alli.sbm$model_parameters[[Q2]]$pi
Pi <- 0.5*(t(Pi) + Pi)
Pi2 <- 0.5*(t(Pi2) + Pi2)

# Gráficos de las probabilidades de interacción
windows(width = 18, height = 10)
par(mfrow = c(1,2), mar = c(3,1.3,3,3))
corrplot::corrplot(corr = Pi, type = "full", col.lim = c(0,1),  method = "shade", addgrid.col = "gray90", tl.col = "black")
title("Trabajo",col.main="#00868B")
corrplot::corrplot(corr = Pi2, type = "full", col.lim = c(0,1),  method = "shade", addgrid.col = "gray90", tl.col = "black")
title("Alianza", col.main="#00868B")

### Bondad de ajuste ----
##Work ----

# Simulaciones
B<-10000
set.seed(42)
bs <- stats::rmultinom(n = B, size = vcount(gcc_work), prob = alpha)
GSIM<-NULL
for (i in 1:B){
  GSIM[[i]]<-igraph::sample_sbm(n = vcount(gcc_work), pref.matrix = Pi, block.sizes = bs[,i], directed = F)
}

### Gráficos PPP
PPP<-A_DescriptivoP1(GSIM)
P <- 0.05

windows(17,10)
ppp <- c(mean(edge_density(gcc_work) < PPP$Densidad),
         mean(transitivity(gcc_work) < PPP$Coef_Agrupa),
         mean(assortativity.degree(gcc_work) < PPP$Asorta),
         mean(mean_distance(gcc_work) < PPP$Dist_Geod),
         mean(max(sapply(X = decompose(gcc_work), FUN = vcount)) < PPP$Comp_Gigante),
         mean(mean(degree(gcc_work)) < PPP$Grado) )
ppp <- round(ppp,3)

grid.arrange(
  #Densidad
  ggplot(PPP, aes(x = Densidad)) + 
    geom_histogram(aes(y = ..density..),
                   colour = "gray86", fill = "gray86") +
    geom_vline(xintercept = edge_density(gcc_work), colour="#FF6347", lwd=1)+
    geom_vline(xintercept = quantile(PPP$Densidad, probs = c(P/2, 1-P/2)), color="turquoise", lwd=1, lty=2)+
    geom_vline(xintercept = mean(PPP$Densidad), color="dodgerblue4" ,lwd=1, lty = 4)+
    labs(y="", x="Densidad")+
    geom_label(aes(x=0.0175, y=210, label=ppp[1]), stat = "unique") +
    theme(panel.background = element_blank()),
  #Transitividad
  ggplot(PPP, aes(x = Coef_Agrupa)) + 
    geom_histogram(aes(y = ..density..),
                   colour = "gray86", fill = "gray86") +
    geom_vline(xintercept = transitivity(gcc_work) , colour="#FF6347", lwd=1)+
    geom_vline(xintercept = quantile(PPP$Coef_Agrupa, probs = c(P/2, 1-P/2)), color="turquoise", lwd=1, lty=2)+
    geom_vline(xintercept = mean(PPP$Coef_Agrupa), color="dodgerblue4" ,lwd=1, lty = 4)+
    labs(y="", x="Transitividad")+
    geom_label(aes(x=0.133, y=36.25, label=ppp[2]), stat = "unique") +
    theme(panel.background = element_blank()),
  #Asortividad
  ggplot(PPP, aes(x = Asorta)) + 
    geom_histogram(aes(y = ..density..),
                   colour = "gray86", fill = "gray86") +
    geom_vline(xintercept = assortativity_degree(gcc_work) , colour="#FF6347", lwd=1)+
    geom_vline(xintercept = quantile(PPP$Asorta, probs = c(P/2, 1-P/2)), color="turquoise", lwd=1, lty=2)+
    geom_vline(xintercept = mean(PPP$Asorta), color="dodgerblue4", lwd=1, lty = 4)+
    labs(y="", x="Asortatividad")+
    geom_label(aes(x=0, y=9.4, label=ppp[3]), stat = "unique") +
    theme(panel.background = element_blank()),
  #Dist_Geod
  ggplot(PPP, aes(x = Dist_Geod)) + 
    geom_histogram(aes(y = ..density..),
                   colour = "gray86", fill = "gray86") +
    geom_vline(xintercept = mean_distance(gcc_work) , colour="#FF6347", lwd=1)+
    geom_vline(xintercept = quantile(PPP$Dist_Geod, probs = c(P/2, 1-P/2)), color="turquoise", lwd=1, lty=2)+
    geom_vline(xintercept = mean(PPP$Dist_Geod), color="dodgerblue4", lwd=1, lty = 4)+
    labs(y="", x="Distancia Geodésica")+
    geom_label(aes(x=4.6, y=1.32, label=ppp[4]), stat = "unique") +
    theme(panel.background = element_blank()),
  #Comp_Gigante
  ggplot(PPP, aes(x = Comp_Gigante)) + 
    geom_histogram(aes(y = ..density..),
                   colour = "gray86", fill = "gray86") +
    geom_vline(xintercept = max(sapply(X = decompose(gcc_work), FUN = vcount)) , colour="#FF6347", lwd=1)+
    geom_vline(xintercept = quantile(PPP$Comp_Gigante, probs = c(P/2, 1-P/2)), color="turquoise", lwd=1, lty=2)+
    geom_vline(xintercept = mean(PPP$Comp_Gigante), color="dodgerblue4", lwd=1, lty = 4)+
    labs(y="", x="Tamaño Componente Gigante")+
    geom_label(aes(x=475, y=0.0405, label="0.000"), stat = "unique") +
    theme(panel.background = element_blank()),
  #Grado
  ggplot(PPP, aes(x = Grado)) + 
    geom_histogram(aes(y = ..density..),
                   colour = "gray86", fill = "gray86") +
    geom_vline(xintercept = mean(degree(gcc_work)) , colour="#FF6347", lwd=1)+
    geom_vline(xintercept = quantile(PPP$Grado, probs = c(P/2, 1-P/2)), color="turquoise", lwd=1, lty=2)+
    geom_vline(xintercept = mean(PPP$Grado), color="dodgerblue4", lwd=1, lty = 4)+
    labs(y="", x="Grado")+
    geom_label(aes(x=8, y=0.4305, label=ppp[6]), stat = "unique") +
    theme(panel.background = element_blank()),
  
  ncol=3,
  right = legendGrob(c("Observado", "Cuantiles", "Media"),
                     gp=gpar(lty = c(1,2,4), col= c("#FF6347", "turquoise", "dodgerblue4")))
  )

##Alliance----

# Simulaciones
set.seed(42)
bs2 <- stats::rmultinom(n = B, size = vcount(gcc_alliance), prob = alpha2)
GSIM2<-NULL
for (i in 1:B){
  GSIM2[[i]]<-igraph::sample_sbm(n = vcount(gcc_alliance), pref.matrix = Pi2, block.sizes = bs2[,i], directed = F)
  
}

### Gráficos PPP
PPP2<-A_DescriptivoP1(GSIM2)

ppp2 <- c(mean(edge_density(gcc_alliance) < PPP2$Densidad),
         mean(transitivity(gcc_alliance) < PPP2$Coef_Agrupa),
         mean(assortativity.degree(gcc_alliance) < PPP2$Asorta),
         mean(mean_distance(gcc_alliance) < PPP2$Dist_Geod),
         mean(max(sapply(X = decompose(gcc_alliance), FUN = vcount)) < PPP2$Comp_Gigante),
         mean(mean(degree(gcc_alliance)) < PPP2$Grado) )
ppp2 <- round(ppp2,3)

windows(17,10)
grid.arrange(
  #Densidad
  ggplot(PPP2, aes(x = Densidad)) + 
    geom_histogram(aes(y = ..density..),
                   colour = "gray86", fill = "gray86") +
    geom_vline(xintercept = edge_density(gcc_alliance), colour="#FF6347", lwd=1)+
    geom_vline(xintercept = quantile(PPP2$Densidad, probs = c(P/2, 1-P/2)), color="turquoise", lwd=1, lty=2)+
    geom_vline(xintercept = mean(PPP2$Densidad), color="dodgerblue4" ,lwd=1, lty = 4)+
    labs(y="", x="Densidad")+
    geom_label(aes(x=0.0165, y=240, label=ppp2[1]), stat = "unique") +
    theme(panel.background = element_blank()),
  #Transitividad
  ggplot(PPP2, aes(x = Coef_Agrupa)) + 
    geom_histogram(aes(y = ..density..),
                   colour = "gray86", fill = "gray86") +
    geom_vline(xintercept = transitivity(gcc_alliance) , colour="#FF6347", lwd=1)+
    geom_vline(xintercept = quantile(PPP2$Coef_Agrupa, probs = c(P/2, 1-P/2)), color="turquoise", lwd=1, lty=2)+
    geom_vline(xintercept = mean(PPP2$Coef_Agrupa), color="dodgerblue4" ,lwd=1, lty = 4)+
    labs(y="", x="Transitividad")+
    geom_label(aes(x=0.22, y=17.2, label=ppp2[2]), stat = "unique") +
    theme(panel.background = element_blank()),
  #Asortividad
  ggplot(PPP2, aes(x = Asorta)) + 
    geom_histogram(aes(y = ..density..),
                   colour = "gray86", fill = "gray86") +
    geom_vline(xintercept = assortativity_degree(gcc_alliance) , colour="#FF6347", lwd=1)+
    geom_vline(xintercept = quantile(PPP2$Asorta, probs = c(P/2, 1-P/2)), color="turquoise", lwd=1, lty=2)+
    geom_vline(xintercept = mean(PPP2$Asorta), color="dodgerblue4", lwd=1, lty = 4)+
    labs(y="", x="Asortatividad")+
    geom_label(aes(x=0.34, y=4.348, label=ppp2[3]), stat = "unique") +
    theme(panel.background = element_blank()),
  #Dist_Geod
  ggplot(PPP2, aes(x = Dist_Geod)) + 
    geom_histogram(aes(y = ..density..),
                   colour = "gray86", fill = "gray86") +
    geom_vline(xintercept = mean_distance(gcc_alliance) , colour="#FF6347", lwd=1)+
    geom_vline(xintercept = quantile(PPP2$Dist_Geod, probs = c(P/2, 1-P/2)), color="turquoise", lwd=1, lty=2)+
    geom_vline(xintercept = mean(PPP2$Dist_Geod), color="dodgerblue4", lwd=1, lty = 4)+
    labs(y="", x="Distancia Geodésica")+
    geom_label(aes(x=5.22, y=1.2, label=ppp2[4]), stat = "unique") +
    theme(panel.background = element_blank()),
  #Comp_Gigante
  ggplot(PPP2, aes(x = Comp_Gigante)) + 
    geom_histogram(aes(y = ..density..),
                   colour = "gray86", fill = "gray86") +
    geom_vline(xintercept = max(sapply(X = decompose(gcc_alliance), FUN = vcount)) , colour="#FF6347", lwd=1)+
    geom_vline(xintercept = quantile(PPP2$Comp_Gigante, probs = c(P/2, 1-P/2)), color="turquoise", lwd=1, lty=2)+
    geom_vline(xintercept = mean(PPP2$Comp_Gigante), color="dodgerblue4", lwd=1, lty = 4)+
    labs(y="", x="Tamaño Componente Gigante")+
    geom_label(aes(x=323, y=0.036, label="0.000"), stat = "unique") +
    theme(panel.background = element_blank()),
  #Grado
  ggplot(PPP2, aes(x = Grado)) + 
    geom_histogram(aes(y = ..density..),
                   colour = "gray86", fill = "gray86") +
    geom_vline(xintercept = mean(degree(gcc_alliance)) , colour="#FF6347", lwd=1)+
    geom_vline(xintercept = quantile(PPP2$Grado, probs = c(P/2, 1-P/2)), color="turquoise", lwd=1, lty=2)+
    geom_vline(xintercept = mean(PPP2$Grado), color="dodgerblue4", lwd=1, lty = 4)+
    labs(y="", x="Grado")+
    geom_label(aes(x=5.3, y=0.76, label=ppp2[6]), stat = "unique") +
    theme(panel.background = element_blank()),
  
  ncol=3,
  right = legendGrob(c("Observado", "Cuantiles", "Media"),
                     gp=gpar(lty = c(1,2,4), col= c("#FF6347", "turquoise", "dodgerblue4")))
)

### Gráficos de las agrupaciones -----

colores <- c("#8EE5EE","yellow1","mediumorchid1","#00CD66","red","#8B4513")

lab<-rep(NA,vcount(gcc_work))
lab[c(1,53)]<-c("Uribe","Santos")
windows(15,11)
ggraph(layout_1)+
  geom_edge_link(edge_colour = "#CDC5BF", alpha=0.15)+
  geom_node_point(aes(size = Grado/3, shape=as.factor(labs), color=as.factor(labs)))+
  scale_color_manual(values = colores)+
  geom_node_label(aes(label = lab),label.size = 0,label.padding = unit(0.2, "lines"),nudge_y = 7.5,size=4)+
  ggtitle("Partición inducida por el modelo para 'Trabajo' ")+
  tema +
  guides(size = "none", shape=guide_legend("Grupo"), color=guide_legend("Grupo"))

lab_a<-rep(NA,vcount(gcc_alliance))
lab_a[c(1,25,81)]<-c("Uribe","Santos","Vargas Lleras")
windows(15,11)
ggraph(layout_a)+
  geom_edge_link(edge_colour ="#CDC5BF", alpha=0.4)+
  geom_node_point(aes(size = Grado_al, shape=as.factor(labs2), color=as.factor(labs2)))+
  scale_color_manual(values = colores)+
  geom_node_label(aes(label = lab_a),label.size = 0,label.padding = unit(0.2, "lines"),nudge_y = 5.5,size=4.2)+
  ggtitle("Partición inducida por el modelo para 'Alianzas'")+
  tema +
  guides(size = "none", shape=guide_legend("Grupo"), color=guide_legend("Grupo"))
