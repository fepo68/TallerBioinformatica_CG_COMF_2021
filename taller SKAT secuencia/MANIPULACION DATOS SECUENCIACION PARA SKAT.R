
### ESTABLECER DIRECTORIO DE TRABAJO  

setwd(setwd("/cloud/project/taller SKAT secuencia"))

# DATOS DE EJEMPLO PARA UN INDIVIDUO (1 ARCHIVO)- solo para visualizar y entender lo que se hace dentro del bucle de m?s abajo:
# abrimos datos gen?ticos

M1 <- read.delim("Rho_ind1.txt",stringsAsFactors=F,na.strings=c("","."))
names(M1)# para ver donde esta situada cada variable y como se llama (EXACTAMENTE)
head(M1)# visualizamos las primeras filas
M1[1:10,1:8]# visualizamos un rango (aqui las 10 primeras filas y las 8 primeras columnas)
dim(M1)# en total este archivo (ahora df) tiene 36 filas y 153 variables...

### SELECCION DE FILAS 
# nos quedamos con las variantes PASS (filtro calidad, se puede a?adir lo que se quiera) y acortamos variables...

qM1<-M1[M1$GATK.Filter=='PASS',c(1:5,8,13,14,16,17,45,51,52,56,80,83)] # SELECCION DE VARIABLES DE INTERES O QUE SE VAN A USAR...
dim(qM1)# ME QUEDO SOLO CON PASS Y REDUZCO EL NUMERO DE VARIABLES 
table(M1$GATK.Filter)# confirmo que solo hab?a 1 fila que no era PASS en M1
names(qM1)

### CAMBIO EL NOMBRE DE VARIABLES 
# manias mias, cambio el de las freq en 1000G, muy largo y empieza por numero, por eso R a?ade una X...
names(qM1)[names(qM1) == "X1000g2015aug_all"] <- "MG_ALL"
names(qM1)[names(qM1) == "X1000g2015aug_eur"] <- "MG_EUR"

### RECODIFICACION VARIABLES
# IMPORTANTE: defino una variable que sea el genotipo de los individuos, como 1 (het) o 2 (homoz), el resto se rellenar? con 0
# la llamo como el individuo porque luego se van a juntar todos los archivos (cada individuo una columna inicialmente)

qM1$ind_M1<-ifelse(qM1$GATK.Zygosity=='het',1,2)

# OTRAS MANIPULACIONES Y FILTRADOS UTILES...

# PARA SELECCIONAR VARIANTES EN FUNCION DE LA FREQ:
# se selecciona si la freq es menor al umbral que sea (0.01) PERO tambien hay que considerar el otro extremo (0.99) si el individuo es het

rM1 = subset(qM1,qM1$MG_ALL <= 0.01 | (qM1$MG_ALL >= 0.99 & qM1$GATK.Zygosity=="het"))
dim(rM1)# 6 variantes raras en este individuo en esos genes
names(qM1)

# PROBLEMA/DILEMA para variantes comunes.... 
# SKAT se puede hacer realmente con todos PERO los genotipos tienen que ir codificados en aditivo (0/1/2)
# en SNPs de freq intermedias no es f?cil estar seguro de si un homoz es 0 o 2. Los resultados podr?an variar...


# PARA SELECCIONAR VARIANTES EXONICAS O DE SPLICING
table(qM1$Func.refGene)# confirmo las categorias que hay en Func refGene
unique(qM1$Func.refGene)# solo para confirmar cuantas categor?as hay, sin recuentos

qM1_ES = subset(qM1,grepl("exonic|splicing",Func.refGene)) # se eliminan variantes que en Func.refGene contengan "exonic" o (|) "splicing"
unique(qM1_ES$Func.refGene)# se eliminan las variantes intr?nicas: 31 VARIANTES EXONICAS Y/O SPLICING

# CUANDO EL NUMERO DE CLASES ES MUY SUPERIOR: el patr?n/expresi?n regular se complica... si por ej queremos eliminar los ncRNA:
# qM1_ESno = subset(qM1,grepl("exonic|splicing",Func.refGene)&!grepl("ncRNA",Func.refGene))
# no podemos hacer esto (!grepl ncRNA) pq entonces eliminar?a todas las que lo incluyeran (y algunas pueden ser ademas exonic)
# si por ej confirmamos que exonic o splicing siempre est?n al principio o 1 caso de splicing en medio: con ^ se indica que es al comienzo:
# qM1_ESsnc<- subset(qM1,grepl("^exonic|^splicing|;splicing",Func.refGene))

# tras seleccionar exonicas y/o splicing seleccionamos las NO SINONIMAS, excluimos el resto
unique(qM1_ES$ExonicFunc.refGene)

qM1_func <- subset(qM1_ES,(!qM1_ES$ExonicFunc.refGene %in% c("synonymous SNV", "unknown","synonymous SNV;unknown")))
dim(qM1_func)# 24 VARIANTES son exonic o splicing tras eliminar las synonymous o unknown
qM1_func[,1:9]

# PARA DESCOMPONER UNA COLUMNA EN VARIAS
# Gene.refGene con frecuencia tiene m?s de un gen, separado por ;. Para descomponer esta columna en varias (una por gen):
library(splitstackshape)
gM1 <- concat.split(data = qM1, split.col = "Gene.refGene", sep = ";", drop = TRUE)
head(gM1)

gM1[!is.na(gM1$Gene.refGene_2),]# hay 3 variantes que se adscriben a dos o mas genes

# si tuviesemos varios genes y diversas funciones (no hay correspondencia u orden claro) se podr?a coger el nombre del gen
# de AAChange.refGene (m?s complicado pero posible...)
# g2M1_func2 <- concat.split(data = g2M1_func, split.col = "AAChange.refGene", sep = ";", drop = TRUE)
# g2M1_func3 <- concat.split(data = g2M1_func2, split.col = c("AAChange.refGene_1","AAChange.refGene_2"), sep = ",", drop = TRUE)
# names(g2M1_func3)
# g2M1_func3$rgene<-sub(":.*","",g2M1_func3$AAChange.refGene_1_01)
# head(g2M1_func3)

## A VECES... BUCLE!

# si se va a hacer SKAT con cada gen, una variante que se adscribe a dos deber?a aparecer dos veces, una en cada gen
# peque?o bucle para hacer eso: 

# se mira en qu? posiciones est?n las dos (en este caso) columnas con nombre de genes
# para cada una de las columnas de gene nos quedamos con las filas que tengan dato 
# primero para el gene1 (Gene_refGene_1, aqu? en la columna 16): todas las filas tendr?n algo en esa columna
# y despues para Gene_refGene_2 (columna 17), donde solo hab?a dos filas que ten?an nombre (los dos casos con 2 genes)

names(gM1)

for (i in 17:18) {
     GENO<-as.data.frame(gM1)
     GENO$GENE<-ifelse(!is.na(GENO[,i]), as.character(GENO[,i]),'NA')
     RGENO<-GENO[GENO$GENE!='NA',c(1:16,19)]# nos quedamos con las columnas originales y la ultima (GENE)     
     assign (paste("F_", i, sep = ""), RGENO)# creamos un dataframe al final de cada ejecuci?n
}

names(F_17)
F_17[,c(1:8,17)]# visualizamos las primeras columnas y la ?ltima con "GENE" (ojo! hay que indicar su posicion tras el bucle)
F_18[,c(1:8,17)]

# AHORA FUSIONAMOS LAS FILAS UNAS DEBAJO DE OTRAS (asi un marcador asignado a dos genes tendr? dos filas, una para cada gen)
genes_M1<-rbind(F_17, F_18)
dim(genes_M1)# antes teniamos 35 variantes PASS, ahora hay 38 filas pq 3 marcadores est?n adscritos a dos genes

genes_M1[,c(1:8,17)]


#####################################################################################################################
### MANIPULACION GENERAL PARA COMBINAR LOS ARCHIVOS DE CADA INDIVIUO:

for (i in c(1:28))
 {
     IND<-read.delim(paste("Rho_ind",i,".txt", sep = ""), stringsAsFactors=F,na.strings=c("","."))
     qIND<-IND[IND$GATK.Filter=='PASS',c(1:5,8,13,14,16,17, 45,51,52,56,80,83)] # SELECCION MIA, A MODIFICAR...

# defino genotipo como 1 o 2 aqui para tener ya la variable antes de seleccionar o no
     qIND$GIND<-ifelse(qIND$GATK.Zygosity=='het',1,2)

# le cambio el nombre para que sea igual que el del individuo:
     names(qIND)[names(qIND) == "GIND"] <- paste("IND",i, sep="")
     names(qIND)[names(qIND) == "X1000g2015aug_all"] <- "MG_ALL"
     names(qIND)[names(qIND) == "X1000g2015aug_eur"] <- "MG_EUR"
     rIND<-qIND[,c(1:5,7:17)]

# asigno un nombre a cada dataframe
     assign (paste("ALL_IND", i, sep = ""), rIND)
   }

head(ALL_IND1)# visualizamos las primeras filas del primer df generado en el bucle (para el individuo 1)

# MERGE de todos los archivos CON TODO:
# OJO!!! no puedo incluir en el by GATK.Zygosity porque var?a entre individuos!!!

names(ALL_IND1)

IND_ALL<-Reduce(function(x,y) merge(x, y, by=c("Chr","Start","End","Ref","Alt",
                "Func.refGene","Gene.refGene","ExonicFunc.refGene","AAChange.refGene","avsnp150",
                "MG_ALL","MG_EUR", "MaxPopFreq","SIFT_pred","Polyphen2_HDIV_pred"),all=TRUE), 
                 list(ALL_IND1,ALL_IND2,ALL_IND3,ALL_IND4,ALL_IND5,ALL_IND6,ALL_IND7,ALL_IND8,ALL_IND9,ALL_IND10,
                      ALL_IND11,ALL_IND12,ALL_IND13,ALL_IND14,ALL_IND15,ALL_IND16,ALL_IND17,ALL_IND18,ALL_IND19,ALL_IND20,
                      ALL_IND21,ALL_IND22,ALL_IND23,ALL_IND24,ALL_IND25,ALL_IND26,ALL_IND27,ALL_IND28))


dim(IND_ALL)#109 marcadores en total, 28 individuos cada uno en una columna
head(IND_ALL)
names(IND_ALL)

# MANIPULACIONES IMPORTANTES PARA SKAT
# relleno con 0 todos los missings (se supone que si tenian variante ser?n homoz para el alelo freq)
# las columnas correspondientes a los genotipos para cada individuo son de la 16 a la 43:

IND_ALL[,16:43][is.na(IND_ALL[,16:43])] = 0

# USO DE IFELSE
# si queremos una variable completa identificadora del marcador (avsnp150 no est? completa, solo cuando hay rs...)
# defino "SNP" como avsnp150 y, en los casos en los que no hay rs, lo creo combinando chr y position 

attach(IND_ALL)
IND_ALL$SNP<-ifelse(IND_ALL$avsnp150=='NA',paste("chr",Chr, Start, sep=":"),IND_ALL$avsnp150)

IND_ALL[,c(1:8,44)] # visualizo las primeras columnas y la ?ltima que acabo de crear
dim(IND_ALL)

library(splitstackshape)# no har?a falta volver a abrirla porque ya se hizo antes...
GIND_ALL <- concat.split(data = IND_ALL, split.col = "Gene.refGene", sep = ";")# HAY FILAS CON 2 GENES..
GIND_ALL[!is.na(GIND_ALL$Gene.refGene_2),]# hay 4 marcadores en los que la segunda columna para el gen tiene nombre

names(GIND_ALL)# para confirmar donde estan las columnas con el nombre del gen descompuesto

for (i in 45:46) {
     class(GIND_ALL)
     GENO<-as.data.frame(GIND_ALL)
     GENO$GENE<-ifelse(!is.na(GENO[,i]), as.character(GENO[,i]),'NA')
     RGENO<-GENO[GENO$GENE!='NA',c(1:44,47)]     
     assign (paste("G_", i, sep = ""), RGENO)
}
dim(G_45)# los 109 marcadores tienen informacion para al menos 1 gen
dim(G_46)# hay 4 que se asignan tambien a un segundo gen (Gene.refGene tiene dos genes)


# PARA REPETIR LOS MARCADORES EN CADA GENE QUE APAREZCA ANOTADO:

REUNO_ALL<-rbind(G_45, G_46)
dim(REUNO_ALL)# 113 filas (repitiendo cada variante-gen en una fila, 109+4)

# SI SE QUIERE FILTRAR LAS FUNCIONALES 

unique(REUNO_ALL$Func.refGene)
REUNO_ES<- subset(REUNO_ALL,grepl("^exonic|^splicing|;splicing",Func.refGene))
dim(REUNO_ES)# BAJA A 96 VARIANTES
unique(REUNO_ES$Func.refGene)# ya no hay intronicas
unique(REUNO_ES$ExonicFunc.refGene)

REUNO_func <- subset(REUNO_ES,(!REUNO_ES$ExonicFunc.refGene %in% c("synonymous SNV","NA", "unknown","synonymous SNV;unknown")))
dim(REUNO_func)# 70 son exonic o splicing tras eliminar las synonymous o unknown
#OJO, PUEDE HABER ALGUN CASO COMPLICADO (2 genes con exonic;intronic), posibilidad usar AAChange.ref...
 
# EN SKAT SE PUEDE HACER FILTRO EN FUNCI?N DE LA MAF, PERO ES LA CALCULADA EN EL ESTUDIO, SI SE QUIERE FILTRAR ANTES:
REUNO_raras<-REUNO_func[REUNO_func$MaxPopFreq<=0.01,]
dim(REUNO_raras)# filtrando en 0.01 solo tenemos 7 VARIANTES RARAS EN GENERAL 
  
REUNO_raras10<-REUNO_func[REUNO_func$MaxPopFreq<=0.1,]
dim(REUNO_raras10)# filtrando en 0.1 tenemos 26 VARIANTES RARAS EN GENERAL 

write.table(REUNO_ALL,"RH_seq_ALL.txt", sep="\t",col.names=TRUE, quote=FALSE,row.names=FALSE)
write.table(REUNO_func,"RH_seq_FUNC.txt", sep="\t",col.names=TRUE, quote=FALSE, row.names=FALSE)


