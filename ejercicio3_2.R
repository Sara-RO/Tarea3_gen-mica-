#RED BOOLEANA# 
library(igraph)
g<- make_empty_graph(n=6, directed=TRUE)
V(g)$color = "aquamarine"
V(g)$shape = "sphere"
plot(g)
g<-add.edges( g, c(1,2, 1,6, 2,3, 2,4, 3,1, 5,4,4,3,6,5))
plot(g)
#SE HIZO LA RED DE REGULACIÓN TRANSCRIPCIONAL, 
#INHIBICIONES 2--->3 , 3-->1, 4-->3 

##TABLAS DE VERDAD 
### A NOT C      B, A      C (B, D)| NOT B OR D    D (B, E )| B OR E   F (A)   E (F)
### T | F        T, T      F  T  F    F            T  T  F     T       F  F    F  F 
### F | T        F, F      T  F  F    T            F  T F      F       T  T    T  T 
###                        F  F  T    F            T  F T      T 
###                        F  F  F    F            F F F       F 





#CARGAR LIBRERÍAS ADECUADAS 
library(BoolNet)
r<-loadNetwork("REDES/red.txt")
r
plotNetworkWiring(r)
atractores<-getAttractors(r)
atractores
plotAttractors(atractores)
plotStateGraph(atractores) #MUESTRA LA GRÁFICA DE ESTADOS A LOS QUE SE PUEDE IR 
#SON DOS ATRACTORES PROBABLES 
#EL ATRACTOR UNO ES (A,B,C,D,E,F)--> DONDE (0,0,1,0,0,0), SE VA EL 78.12% 
#EL SEGUNDO ATRACTOR CONSISTE EN TRES ESTADOS (COMO SE MUESTRA EN LA FIGURA ), SE VA EL 21.88% 
# CON LO ANTERIOR SE PUEDE DEFINIR QUE EL ESTADO MÁS PROBABLE ES EL DEL ATRACTOR 1 QUE ES (0,0,1,0,0,0)