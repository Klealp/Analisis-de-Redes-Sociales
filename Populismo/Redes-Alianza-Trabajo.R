
setwd("C:/Users/kevin/Desktop/UNIVERSIDAD/8 semestre/Analisis de Redes Sociales/Parcial II")
#setwd("C:/Users/Isabella/Documents/Redes/Parcial 2")
library(igraph)
library(ggraph)
library(xtable)
require(gridExtra) # grilla de ggplot
library(readr)
#modelos
library(blockmodels)
library(ergm)
library(network)
library(ROCR)
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
g_work <- graph_from_data_frame(d =work, directed = "F", vertices = nodes$Id)
g_alliance <- graph_from_data_frame(d =alliance, directed = "F", vertices = nodes$Id)

# grafo inducido por la componente conexa work
V(g_work)$cluster <- clusters(graph = g_work)$membership
gcc_work <- induced_subgraph(graph = g_work, vids = which(V(g_work)$cluster == which.max(clusters(graph = g_work)$csize)))

set.seed(21)
layout_1 <- create_layout(gcc_work, layout = "igraph", algorithm="dh") 

ca<-rep("turquoise",vcount(gcc_work))
ca[c(1,53,145,151,154)]<-"#CA1429"
lab<-rep(NA,vcount(gcc_work))
lab[c(1,53,145,151,154)]<-c("Uribe","Santos","Pastrana","J.C. Restrepo","Gaviria")
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

set.seed(21)
layout_a <- create_layout(gcc_alliance, layout = "igraph", algorithm="dh")

ca_a<-rep("turquoise",vcount(gcc_alliance))
ca_a[c(1,25,81,127,176)]<-"#CA1429"
lab_a<-rep(NA,vcount(gcc_alliance))
lab_a[c(1,25,81,127,176)]<-c("Uribe","Santos","Vargas Lleras","Petro","Fajardo")

Grado_al <- log(degree(gcc_alliance))

windows(16,12)
ggraph(layout_a)+
  geom_edge_link(edge_colour ="#CDC5BF", alpha=0.4)+
  geom_node_point(aes(size = Grado_al), color= ca_a, alpha=0.7)+
  geom_node_label(aes(label = lab_a),label.size = 0,label.padding = unit(0.2, "lines"),nudge_y = 5.5,size=4.2)+
  ggtitle("Red de relaciones de alianza")+labs(size = "Grado")+
  tema +
  theme(legend.position = "none")

#Análisis descriptivo

A_DescriptivoP1 <- function(Lista){
  
  #Calculo de estadísticas estructurales
  Densidad      <- sapply(X = Lista, FUN = edge_density)
  Coef_Agrupa   <- sapply(X = Lista, FUN = transitivity)
  Asorta <- sapply(X = Lista, FUN = assortativity_degree)
  Dist_Geod     <- sapply(X = Lista, FUN = mean_distance)
  Comp_Gigante  <- sapply(X = Lista, FUN = function(X) (max(sapply(X = decompose(X), FUN = vcount))))
  Grado    <- sapply(X = Lista, FUN = function(X) mean(degree(X)))
  Tabla         <- as.data.frame(cbind(Densidad,Coef_Agrupa,Asorta,
                                       Dist_Geod,Comp_Gigante,Grado))
  return(Tabla)
}

tab<-A_DescriptivoP1(list(gcc_work,gcc_alliance))
rownames(tab)<-c("Trabajo","Alianza")

xtable::xtable(tab,digits = 3)

#Nodos en común para ambas componentes 

CW<-as.data.frame(V(gcc_work)$name)
colnames(CW)<-colnames(CA)<-"NAME"
CA<-as.data.frame(V(gcc_alliance)$name)
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

#Simplificación de las redes
gcc_work <- simplify(gcc_work)
gcc_alliance <- simplify(gcc_alliance)

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

windows(17,10)
grid.arrange(
  #Densidad
  ggplot(PPP, aes(x = Densidad)) + 
    geom_histogram(aes(y = ..density..),
                   colour = 1, fill = "white") +
    geom_density(lwd = 1, colour = "#00C5CD",
                 fill = "#00C5CD", alpha = 0.25)+
    geom_vline(xintercept = edge_density(gcc_work) , colour="#FF6347", lwd=1)+
    labs(y="", x="Densidad")+
    theme(panel.background = element_blank()),
  #Transitividad
  ggplot(PPP, aes(x = Coef_Agrupa)) + 
    geom_histogram(aes(y = ..density..),
                   colour = 1, fill = "white") +
    geom_density(lwd = 1, colour = "#00C5CD",
                 fill = "#00C5CD", alpha = 0.25)+
    geom_vline(xintercept = transitivity(gcc_work) , colour="#FF6347", lwd=1)+
    labs(y="", x="Transitividad")+
    theme(panel.background = element_blank()),
  #Asortividad
  ggplot(PPP, aes(x = Asorta)) + 
    geom_histogram(aes(y = ..density..),
                   colour = 1, fill = "white") +
    geom_density(lwd = 1, colour = "#00C5CD",
                 fill = "#00C5CD", alpha = 0.25)+
    geom_vline(xintercept = assortativity_degree(gcc_work) , colour="#FF6347", lwd=1)+
    labs(y="", x="Asortividad")+
    theme(panel.background = element_blank()),
  #Dist_Geod
  ggplot(PPP, aes(x = Dist_Geod)) + 
    geom_histogram(aes(y = ..density..),
                   colour = 1, fill = "white") +
    geom_density(lwd = 1, colour = "#00C5CD",
                 fill = "#00C5CD", alpha = 0.25)+
    geom_vline(xintercept = mean_distance(gcc_work) , colour="#FF6347", lwd=1)+
    labs(y="", x="Dist_Geod")+
    theme(panel.background = element_blank()),
  #Comp_Gigante
  ggplot(PPP, aes(x = Comp_Gigante)) + 
    geom_histogram(aes(y = ..density..),
                   colour = 1, fill = "white") +
    geom_density(lwd = 1, colour = "#00C5CD",
                 fill = "#00C5CD", alpha = 0.25)+
    geom_vline(xintercept = max(sapply(X = decompose(gcc_work), FUN = vcount)) , colour="#FF6347", lwd=1)+
    labs(y="", x="Comp_Gigante")+
    theme(panel.background = element_blank()),
  #Grado
  ggplot(PPP, aes(x = Grado)) + 
    geom_histogram(aes(y = ..density..),
                   colour = 1, fill = "white") +
    geom_density(lwd = 1, colour = "#00C5CD",
                 fill = "#00C5CD", alpha = 0.25)+
    geom_vline(xintercept = mean(degree(gcc_work)) , colour="#FF6347", lwd=1)+
    labs(y="", x="Grado")+
    theme(panel.background = element_blank()),
  ncol=3)

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

windows(17,10)
grid.arrange(
  #Densidad
  ggplot(PPP2, aes(x = Densidad)) + 
    geom_histogram(aes(y = ..density..),
                   colour = 1, fill = "white") +
    geom_density(lwd = 1, colour = "#00C5CD",
                 fill = "#00C5CD", alpha = 0.25)+
    geom_vline(xintercept = edge_density(gcc_alliance) , colour="#FF6347", lwd=1)+
    labs(y="", x="Densidad")+
    theme(panel.background = element_blank()),
  #Transitividad
  ggplot(PPP2, aes(x = Coef_Agrupa)) + 
    geom_histogram(aes(y = ..density..),
                   colour = 1, fill = "white") +
    geom_density(lwd = 1, colour = "#00C5CD",
                 fill = "#00C5CD", alpha = 0.25)+
    geom_vline(xintercept = transitivity(gcc_alliance) , colour="#FF6347", lwd=1)+
    labs(y="", x="Transitividad")+
    theme(panel.background = element_blank()),
  #Asortividad
  ggplot(PPP2, aes(x = Asorta)) + 
    geom_histogram(aes(y = ..density..),
                   colour = 1, fill = "white") +
    geom_density(lwd = 1, colour = "#00C5CD",
                 fill = "#00C5CD", alpha = 0.25)+
    geom_vline(xintercept = assortativity_degree(gcc_alliance) , colour="#FF6347", lwd=1)+
    labs(y="", x="Asortividad")+
    theme(panel.background = element_blank()),
  #Dist_Geod
  ggplot(PPP2, aes(x = Dist_Geod)) + 
    geom_histogram(aes(y = ..density..),
                   colour = 1, fill = "white") +
    geom_density(lwd = 1, colour = "#00C5CD",
                 fill = "#00C5CD", alpha = 0.25)+
    geom_vline(xintercept = mean_distance(gcc_alliance) , colour="#FF6347", lwd=1)+
    labs(y="", x="Dist_Geod")+
    theme(panel.background = element_blank()),
  #Comp_Gigante
  ggplot(PPP2, aes(x = Comp_Gigante)) + 
    geom_histogram(aes(y = ..density..),
                   colour = 1, fill = "white") +
    geom_density(lwd = 1, colour = "#00C5CD",
                 fill = "#00C5CD", alpha = 0.25)+
    geom_vline(xintercept = max(sapply(X = decompose(gcc_alliance), FUN = vcount)) , colour="#FF6347", lwd=1)+
    labs(y="", x="Comp_Gigante")+
    theme(panel.background = element_blank()),
  #Grado
  ggplot(PPP2, aes(x = Grado)) + 
    geom_histogram(aes(y = ..density..),
                   colour = 1, fill = "white") +
    geom_density(lwd = 1, colour = "#00C5CD",
                 fill = "#00C5CD", alpha = 0.25)+
    geom_vline(xintercept = mean(degree(gcc_alliance)) , colour="#FF6347", lwd=1)+
    labs(y="", x="Grado")+
    theme(panel.background = element_blank()),
  ncol=3)

### Agrupación no probabilística ------
# Aplicación algoritmos
Ag_louvain <- cluster_louvain(gcc_work) 
Ag_eigen <- cluster_leading_eigen(gcc_alliance)

#Dist. agrupaciones
sizes(AgW2)
sizes(kc_louvain)

#### Funciones para los mapas de calor ----
# para ordenar la matriz de adyacencia respecto a una partición
get_adjacency_ordered <- function(xi, A) 
{
  xi2 <- xi[order(xi)]
  indices <- order(xi)
  d <- NULL
  for (i in 1:(length(xi)-1)) if (xi2[i] != xi2[i+1]) d <- c(d, i)
  list(A = A[indices,indices], d = d)
}
# para graficar la matriz de adyacencia
heat.plot0 <- function (mat, show.grid = FALSE, cex.axis, tick, labs, col.axis, ...)
{ 
  JJ <- dim(mat)[1]
  colorscale <- c("white", rev(heat.colors(100)))
  if(missing(labs))     labs <- 1:JJ
  if(missing(col.axis)) col.axis <- rep("black", JJ)
  if(missing(cex.axis)) cex.axis <- 1
  if(missing(tick))     tick <- TRUE
  ## adjacency matrix
  image(seq(1, JJ), seq(1, JJ), mat, axes = FALSE, xlab = "", ylab = "", col = colorscale[seq(floor(100*min(mat)), floor(100*max(mat)))], ...)
  for(j in 1:JJ){
    axis(1, at = j, labels = labs[j], las = 2, cex.axis = cex.axis, tick, col.axis = col.axis[j], col.ticks = col.axis[j])
    axis(2, at = j, labels = labs[j], las = 2, cex.axis = cex.axis, tick, col.axis = col.axis[j], col.ticks = col.axis[j])
  }
  box()
  if(show.grid) grid(nx = JJ, ny = JJ)
}

Mapa_modelo <- function(Z,Matriz){
  # asignaciones de grupos
  xi <- apply(X = Z, MARGIN = 1, FUN = which.max)
  # matriz de adyacencia ordenada y lineas divisorias de acuerdo con las comunidades
  tmp <- get_adjacency_ordered(xi = xi, A = Matriz)
  # A
  heat.plot0(mat = tmp$A, tick = F, labs = NA, main="Bloques estocásticos", col.main="#008B8B")
  abline(v = tmp$d+.5, h = tmp$d+.5)
}

Mapa_agru <- function(agru,Matriz){
  #Mapa de calor de la jerárquica
  # asignaciones 
  xi <- agru$membership
  # asignaciones ordenadas 
  xi2 <- xi[order(xi)]
  # matriz de adyacencia ordenada y lineas divisorias de acuerdo con las comunidades
  tmp <- get_adjacency_ordered(xi = xi, A = Matriz)
  # grafico
  heat.plot0(mat = tmp$A, labs=NA, tick = F, main="Agrup. no probabilística", col.main="#008B8B")
  abline(v = tmp$d+.5, h = tmp$d+.5)}

#### Mapas de calor----

#Work
windows(12,6.5)
par(mfrow=c(1,2))
Mapa_modelo(Z,MW)
title("\n \n Trabajo", col.main="#00868B")
Mapa_agru(Ag_louvain,MW)
title("\n \n Trabajo", col.main="#00868B")

#Alliance
windows(12,6.5)
par(mfrow=c(1,2))
Mapa_modelo(Z2,MA)
title("\n \n Alianza", col.main="#00868B")
Mapa_agru(Ag_eigen,MA)
title("\n \n Alianza", col.main="#00868B")

### Comparación de agrupamientos ----

#Work
med<-c(round(igraph::compare(comm1 = Ag_louvain$membership, comm2 = labs, method = "rand"), 3),
       round(igraph::compare(comm1 = Ag_louvain$membership, comm2 = labs, method = "adjusted.rand"), 3),
       round(igraph::compare(comm1 = Ag_louvain$membership, comm2 = labs, method = "nmi"), 3))
#Alliance
med2<-c(round(igraph::compare(comm1 = Ag_eigen$membership, comm2 = labs2, method = "rand"), 3),
        round(igraph::compare(comm1 = Ag_eigen$membership, comm2 = labs2, method = "adjusted.rand"), 3),
        round(igraph::compare(comm1 = Ag_eigen$membership, comm2 = labs2, method = "nmi"), 3)
)
xtable::xtable(cbind(c("AR","ARI","NMI"),med,med2))

#Tablas de contingencia
xtable::xtable(addmargins(table(labs,Ag_louvain$membership)),digits = 0)
xtable::xtable(addmargins(table(labs2,Ag_eigen$membership)),digits = 0)

### Gráficos -----

colores <- c("#8EE5EE","yellow1","mediumorchid1","#00CD66","red","#8B4513")
col_ver <- colores[labs]

windows(15,11)
ggraph(layout_1)+
  geom_edge_link(edge_colour = "#CDC5BF", alpha=0.15)+
  geom_node_point(aes(size = Grado/3, shape=as.factor(labs)), color= col_ver)+
  ggtitle("Partición inducida por el modelo para 'Trabajo' ")+
  labs(shape="Grupo",color="Grupo")+
  tema +
  guides(size = "none")

col_ver2 <- colores[labs2]
windows(15,11)
ggraph(layout_a)+
  geom_edge_link(edge_colour ="#CDC5BF", alpha=0.4)+
  geom_node_point(aes(size = Grado_al, shape=as.factor(labs2)),color= col_ver2)+
  ggtitle("Partición inducida por el modelo para 'Alianzas'")+labs(size = "Grado")+
  labs(shape="Grupo",color="Grupo")+
  tema +
  guides(size = "none")

#### Gráfico final ----
# Work
g.cl <- graph_from_adjacency_matrix(adjmatrix = Pi, mode = "undirected", weighted = TRUE)
# parametros del grafico
vsize    <- 100*sqrt(alpha)
ewidth   <- 10*E(g.cl)$weight
PolP     <- Ag_louvain$membership
class.by.PolP <- as.matrix(table(labs,PolP))
pie.vals <- lapply(1:Q, function(i) as.vector(class.by.PolP[i,]))
my.cols  <- paletteer::paletteer_c("ggthemes::Sunset-Sunrise Diverging", length(unique(PolP)))

#Alliance
g.cl2 <- graph_from_adjacency_matrix(adjmatrix = Pi2, mode = "undirected", weighted = TRUE)
# parametros del grafico
vsize2    <- 100*sqrt(alpha2)
ewidth2   <- 10*E(g.cl2)$weight
PolP2     <- Ag_eigen$membership
class.by.PolP2 <- as.matrix(table(labs2,PolP2))
pie.vals2 <- lapply(1:Q2, function(i) as.vector(class.by.PolP2[i,]))
my.cols2  <- paletteer::paletteer_c("ggthemes::Sunset-Sunrise Diverging", length(unique(PolP2)))
# grafico 
windows(11,6)
par(mfrow=c(1,2), mai= c(0, 0.3,0.5,0))
set.seed(21)
plot(g.cl, edge.width = ewidth, vertex.shape = "pie", vertex.pie = pie.vals,
     vertex.pie.color = list(my.cols), vertex.size = vsize,
     vertex.label.dist = 0.1*vsize, vertex.label.degree = pi, vertex.label.color="#00868B")
title("Red de trabajo", col.main="#00868B")
plot(g.cl2, edge.width = ewidth2, vertex.shape = "pie", vertex.pie = pie.vals2,
     vertex.pie.color = list(my.cols2), vertex.size = vsize2,
     vertex.label.dist = 0.1*vsize2, vertex.label.degree = pi, vertex.label.color="#00868B")
title("Red de alianzas", col.main="#00868B")

