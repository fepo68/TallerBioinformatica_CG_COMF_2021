
## SKAT PARA 1 (O POCOS GENES)

library(SKAT)
setwd("/cloud/project/taller SKAT secuencia")

RH <- read.delim("RH_seq_ALL.txt",stringsAsFactors=F,na.strings=c("","."))
head(RH)
dim(RH)
names(RH)
table(RH$GENE)# vemos los recuentos de variantes para cada gen 
table(RH$GENE, RH$Polyphen2_HDIV_pred)# distribucion de variantes de cada gen en funcion de polyphen

# no hace falta (lo arregla SKAT) pero si queremos la MAF real en nuestro estudio:
RH$cHET<-apply(RH[,16:43],1,function(x) length(which(x==1)))
RH$cHOM<-apply(RH[,16:43],1,function(x) length(which(x==2)))
RH$MAF<-(RH$cHET+(2*RH$cHOM))/(2*28)

# RH[RH$MAF>=0.5,c(16:25,43:48)]# visualizamos los casos que pueden ser problematicos (m?s abajo)
# si quisieramos dar la vuelta a esto
# RH[RH$MAF>0.5,16:43]<-lapply(RH[RH$MAF>0.5,16:43],function(x) ifelse(x==2,0,ifelse(x==0,2,1)))

# PARA DEFINIR PESOS EN FUNCION DE POLYPHEN POR EJ...

RH$W_POLYPHEN<-ifelse(RH$Polyphen2_HDIV_pred=='B', 0.25,
                      ifelse(RH$Polyphen2_HDIV_pred=='NA', 0.55,
                      ifelse(RH$Polyphen2_HDIV_pred=='P', 0.75,
                      ifelse(RH$Polyphen2_HDIV_pred=='D', 0.85,-99))))
table(RH$W_POLYPHEN)
head(RH)

# ANALISIS PARA LOS 51 MARCADORES DE RP1L1 - TAMBIEN SE PODR?A ANALIZAR TODO JUNTO PORQUE HAY UNA HIP?TESIS/MOTIVO 
# (genes implicados en la homeostasis de retina)
# Seleccionamos las filas correspondientes a ese gen y generamos la matriz de genotipos (Z) TRASPUESTA

names(RH)
R_RP1L1<-RH[RH$GENE=='RP1L1',]# seleccionamos las filas de ese gen en el archivo original
R_RP1L1[,c(16:44,48,49)]
RP1L1<-R_RP1L1[,16:43]# seleccionamos las columnas con los genotipos en las filas (variantes) de ese gen
head(RP1L1)
RP1L1[1:15,]

## IMPORTANTE:la funcion SKAT necesita que la matriz de genotipos est? al reves que como la hemos generado
# trasponemos 
Z1<-t(RP1L1) 
dim(Z1)
Z1# ahora los 28 individuos est?n en filas y las 51 variantes de RP1L1 son las columnas

# si quisiesemos otro..:
# NPHP4<-RH[RH$GENE=='NPHP4',16:43]
# NPHP4[1:11,1:8]
# Z2<-t(NPHP4)

# ABRIMOS EL ARCHIVO FENOTIPICO
feno_RH <- read.delim("FENOTIPO_RH.txt",stringsAsFactors=F,na.strings=c("","."))
head(feno_RH)

# si especificamos aqui las variables que vamos a usar en el modelo ya ponemos luego la f?rmula directamente:
# si no hay que indicar feno_RH$EDAD etc...

y.b<-feno_RH$GRAVEDAD # gravedad definida en dos categor?as - situaci?n an?loga a CASO/CONTROL
y.q<-feno_RH$GRAV_CUANT # medida de gravedad num?rica, variable cuantitativa
edad<-feno_RH$EDAD
sex<-feno_RH$SEX

# ajustamo el modelo nulo solo con las covariables o sin nada...

modelo1<-SKAT_Null_Model(y.b ~ edad+sex, out_type="D")
modelo0<-SKAT_Null_Model(y.b ~ 1, out_type="D")

# modelo1<-SKAT_Null_Model(feno_RH$GRAVEDAD ~ feno_RH$SEX+feno_RH$EDAD, out_type="D")# especificando el origen de las variables...

########################################################################################################################################
# SKAT "normal". Opciones por defecto: weights.beta=c(1,25), r.corr=0...

?SKAT

# sin covariables:
out_SKAT<-SKAT(Z1,modelo0)# avisan de la incongruencia de tener marcadores en los que la codificacion no encaja
                          # con lo esperado (0/1/2 en funcion de la maf propia - PERO skat ya le da la vuelta
out_SKAT # se ve todo el output en la consola
out_SKAT$p.value # lo importante: la probabilidad de SKAT para el gen analizado 
out_SKAT$param$n.marker # numero de marcadores iniciales
out_SKAT$param$n.marker.test # numero de marcadores con los que se hace test (marcadores NO MONOMORFICOS)

# con COVARIABLES (lo m?s normal):
out_SKATc<-SKAT(Z1,modelo1)
out_SKATc$p.value
out_SKATc$param$n.marker 
out_SKATc$param$n.marker.test 

# para hacer un test m?s b?sico, BURDEN: 

out_BURDENc<-SKAT(Z1, modelo1, r.corr=1)# para hacer BURDEN
out_BURDENc$p.value

#########################################################################################################################################
# JUGANDO CON LOS PESOS - por defecto SKAT da mucho mayor peso a las variantes raras:

out_SKATw1<-SKAT(Z1, modelo1, weights.beta=c(1,1))# con peso=1 a todos
out_SKATw1$p.value

r_weights<-R_RP1L1$W_POLYPHEN
out_SKATpolyphen<-SKAT(Z1, modelo1, weights=r_weights)
out_SKATpolyphen$p.value

###########################################################################################################################################
# PARA HACER SKAT solo con las variantes raras (opci?n max_maf)

out_RARES05<-SKAT(Z1, modelo1, max_maf=0.05)# para hacer SKAT SOLO RARAS DEFINIDAS COMO MAF<0.05 (se puede cambiar)
out_RARES05$p.value
out_RARES05$param$n.marker.test

## OJO! estos son procedimientos muucho m?s lentos, aqui sin problema (cuidado en GENOME):
# SKAT-O

out_SKATO<-SKAT(Z1, modelo1, method="SKATO")
out_SKATO
out_SKATO$p.value

# CON RESAMPLING BOOTSTRAP - ojo! tarda aun mas que SKATO
# EN ESTO HACEMOS PERMUTACIONES (se tiene que especificar en el modelo)

modelo_RESAMPLING<-SKAT_Null_Model(y.b ~ edad+sex, out_type="D", n.Resampling=1000)
out_SKATboot<-SKAT(Z1, modelo_RESAMPLING, weights.beta=c(1,25),r.corr=0)# SKAT con pesos por defecto
names(out_SKATboot)
out_SKATboot$p.value

