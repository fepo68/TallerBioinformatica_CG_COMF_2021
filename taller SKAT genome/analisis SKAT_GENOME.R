library(SKAT)
library(MetaSKAT)
library(splitstackshape)

setwd("/cloud/project/taller SKAT genome")

# datos de ejemplo de 2026 individuos genotipados con el SpainBA array
# solo chr 22 (BINARIOS DE PLINK): "SPAINBA_chr22.bed"/"SPAINBA_chr22.bim"/"SPAINBA_chr22.fam"

##############################################################################################################################
# MANIPULACI?N INICIAL: es necesario crear/tener un archivo con informaci?n de los genes a los que se asignan los SNPs (fichero .SetID)

## opcion1 - desde el fichero de anotaciones de affy (cojo solo el primer gen, analogo a illumina...)

anotaciones<-read.delim("SpainBA_anotaciones_22.txt", sep="\t", header=TRUE)
head(anotaciones)
names(anotaciones)
table(anotaciones$Chromosome)

anotaciones2<-anotaciones[,c(1:3,5,6,12:13,16)]
head(anotaciones2)
# split concatenated column by `//` SEPARAMOS TODA LA INFORMACION INCLUIDA EN Associated Gene"
anotaciones3 <- concat.split(data = anotaciones2, split.col = "Associated.Gene", sep = "//", drop = TRUE)
head(anotaciones3)

# GENERO ESTE dataframe INFO por si queremos en algun momento repetir el an?lisis limitando los SNPs en funci?n de "type"
INFO<-anotaciones3[,c(1,3:7,10,16)]
head(INFO)
names(INFO)<-c("Probeset_ID","SNP","chr","position","Allele_A","Allele_B","type","Gene")

chr22_SetID<-INFO[,c(8,1)]# solo puede incluir gene y SNP
Rchr22_SetID<-chr22_SetID[!as.character(chr22_SetID$Gene)=="---",]# eliminamos los casos sin gene asignado
dim(chr22_SetID)# 13413 filas
dim(Rchr22_SetID)# 11503 marcadores con gene

# ordenamos en funcion de los genes (si no da un warning y genera un temporal):
R_chr22_SetID<-Rchr22_SetID[order(Rchr22_SetID$Gene),]
head(R_chr22_SetID)
unique(R_chr22_SetID$Gene)# hay 484 genes

write.table(R_chr22_SetID, "SetID_SPAINBA_chr22.SetID", sep="\t", row.names=FALSE, quote=FALSE, col.names=FALSE)

############################################################################################################################
## opcion2: coger el rango de posiciones en la versi?n de genoma correspondiente. Aqui hg19, descargada de recursos de Plink:

lista_hg19<-read.delim("listado_rangos_chr22.txt", sep="\t", header=TRUE)
head(lista_hg19)
dim(lista_hg19)# aqui ya estaba filtrada por mi para tener solo el 22 pero en principio se descarga toda

	SNPs_array<-INFO
dim(SNPs_array) # hay 13413 marcadores genotipados en nuestros datos
	lista19<-lista_hg19[lista_hg19$chr==22,] 
dim(lista19) # hay 583 genes cogidos de la lista hg19 (Plink)
      asignados<-outer(INFO$position, lista_hg19$pos_ini, ">=") & outer(INFO$position, lista_hg19$pos_fin, "<=")
	index<-which(asignados==TRUE, arr.ind=T)
head(index)
dim(index)# 9687 marcadores est?n dentro de esos limites de genes
	indexcol<-as.vector(index[,2])
head(indexcol)
	track<-lista19[indexcol,]# nos quedamos con los genes de la lista a los que se ha asignado algo, repetidos
head(track)
dim(lista19)
dim(track)
	indexrow<-as.vector(index[,1])
head(indexrow)
	boop<-SNPs_array[indexrow,]# nos quedamos con los marcadores asignados a algun gen
	info<-cbind(track,boop)# unimos y aqui seleccionamos las columnas necesarias para el .SetID y exportamos o usamos...
head(info)# aqui se puede comparar el gene del fichero de anotaciones y el que asignamos

##############################################################################################################################
# SKAT - EN PRIMER LUGAR SE GENERA EL OBJETO .SSD e .Info a partir de binarios generados desde plink

SpainBA_chr22<-Generate_SSD_SetID("SPAINBA_chr22.bed", "SPAINBA_chr22.bim", "SPAINBA_chr22.fam", "SetID_SPAINBA_chr22.SetID", 
                                  "EJ_SPAINBA22.SSD", "EJ_SPAINBA22.Info", Is.FlipGenotype=TRUE)

# 484 Sets / 11503 SNPs

# INCORPORACI?N INFORMACI?N FENOTIPICA
# MEJOR OPCION: abrir a la vez el fam y el archivo de covariables CODIFICADAS NUMERICAMENTE

# flag1: 0 represents the default coding of unaffected/affected (1/2) (default=0), and 1 represents 0/1 coding
FAM_cov<-Read_Plink_FAM_Cov("SPAINBA_chr22.fam", "RCOVARIABLES_EJ_SPAINBA.txt",Is.binary=TRUE, cov_header=TRUE, flag1=0)
FAM_cov[1:20,]

# si no hay covariables se abre solo el .fam:
# FAM<-Read_Plink_FAM("SPAINBA_chr22.fam", Is.binary=TRUE, flag1=0)

# si especificamos aqui las variables que vamos a usar en el modelo ya ponemos luego la f?rmula directamente:

y<-FAM_cov$Phenotype
EDAD<-FAM_cov$edad
SEX<-FAM_cov$SEX
TIPO<-FAM_cov$TIPO

SSD.INFO<-Open_SSD("EJ_SPAINBA22.SSD", "EJ_SPAINBA22.Info")
modelo1<-SKAT_Null_Model(y ~ EDAD+SEX, out_type="D")
modelo0<-SKAT_Null_Model(y ~ 1, out_type="D")

# rmodelo1<-SKAT_Null_Model(FAM_cov$Phenotype ~ FAM_cov$SEX+FAM_cov$edad, out_type="D")# especificando el origen de las variables...
# modelo2<-SKAT_Null_Model(y ~ 1+EDAD+SEX+TIPO, out_type="D")# si definimos modelos adicionales, probando diferentes covariables...
# modelo3<-SKAT_Null_Model(y ~ 1+EDAD, out_type="D")# si definimos modelos adicionales, probando diferentes covariables...

########################################################################################################################################
# SKAT "normal". Opciones por defecto: weights.beta=c(1,25), r.corr=0

out_SKAT0<-SKAT.SSD.All(SSD.INFO, modelo0)# sin covariables
SKATsc<-out_SKAT0$results 
SKATsc # visualizo los resultados para los 484 genes

out_SKAT<-SKAT.SSD.All(SSD.INFO, modelo1)# con covariables
head(out_SKAT$results)

# para hacer un test m?s b?sico, BURDEN: 

out_BURDEN<-SKAT.SSD.All(SSD.INFO, modelo1, r.corr=1)# para hacer BURDEN
head(out_BURDEN$results)

#########################################################################################################################################
# JUGANDO CON LOS PESOS - por defecto SKAT da mucho mayor peso a las variantes raras:
# out1<-SKAT.SSD.All(SSD.INFO, modelo1, weights.beta=c(1,25))# SKAT con pesos por defecto (NO HACE FALTA INDICARLO)
# out2<-SKAT.SSD.All(SSD.INFO, modelo1, weights.beta=c(0.5,0.5))# con pesos de Madsen and Browning

out_SKATw1<-SKAT.SSD.All(SSD.INFO, modelo1, weights.beta=c(1,1))# con peso=1 a todos
head(out_SKATw1$results)

#obj.SNPWeight<-Read_SNP_WeightFile("weights_cte.txt")# para usar pesos especificados en archivo aparte, AQUI CONSTANTE=1
#out3<-SKAT.SSD.All(SSD.INFO, modelo1, obj.SNPWeight=obj.SNPWeight)# 
#head(out3$results)# al ser peso=1 tiene que salir lo mismo que en out_SKATw1

# otras alternativas de SKAT (que dan pesos m?s igualitarios a variantes comunes y raras):

# Combined sumed test (SKAT-C y BURDEN-C)
out1CS<-SKAT_CommonRare.SSD.All(SSD.INFO, modelo1)# para hacer SKAT-C
out1CB<-SKAT_CommonRare.SSD.All(SSD.INFO, modelo1, r.corr.rare=1, r.corr.common=1)# para hacer BURDEN-C

# Adaptive test (SKAT-A, BURDEN-A)
# out1AS<-SKAT_CommonRare.SSD.All(SSD.INFO, modelo1, method="A")# para hacer SKAT-A
# out1AB<-SKAT_CommonRare.SSD.All(SSD.INFO, modelo1, method="A", r.corr.rare=1, r.corr.common=1)# para hacer BURDEN-A


###########################################################################################################################################
# PARA HACER SKAT solo con las variantes raras (opci?n max_maf)

out_RARES<-SKAT.SSD.All(SSD.INFO, modelo1, max_maf=0.05)# para hacer SKAT SOLO RARAS DEFINIDAS COMO MAF<0.05 (se puede cambiar)
head(out_RARES$results)

# funcion binaria: La diferencia con SKAT es el m?todo de remuestreo para calcular los pvalues y el m?todo de imputaci?n de missings
# out_SKATB<-SKATBinary.SSD.All(SSD.INFO, modelo1)
# head(out_SKATB$results)

## OJO! estos son procedimientos muucho m?s lentos:
# SKAT-O

out_SKATO<-SKAT.SSD.All(SSD.INFO, modelo1, method="SKATO")
ROUT_SKATO<-out_SKATO$results
ROUT_SKATO[ROUT_SKATO$P.value<0.000001,]
head(out_SKATO$results)

#write.table(OUT_SKATO, "OUTPUT_SKATO_EJ_SPAINBA.txt", sep="\t", col.names=TRUE, row.names=FALSE)#exporto este porque cuesta...

# CON RESAMPLING BOOTSTRAP - ojo! tarda aun mas que SKATO
# EN ESTO HACEMOS PERMUTACIONES (se tiene que especificar en el modelo)
# se muestran luego los que pasan el umbral de FWER 0.05 - FUNCIONA??

modelo_RESAMPLING<-SKAT_Null_Model(y ~ 1, out_type="D", n.Resampling=1000)
out_SKATboot<-SKAT.SSD.All(SSD.INFO, modelo_RESAMPLING, weights.beta=c(1,25),r.corr=0)# SKAT con pesos por defecto
names(out_SKATboot)
head(out_SKATboot$results)
Resampling_FWER(out_SKATboot,FWER=0.05)

# MANIPULACION DE OUTPUTS DE SKAT/RARAS/BURDEN/SKATO

output_SKAT<-out_SKAT$results
output_RARES<-out_RARES$results
output_BURDEN<-out_BURDEN$results
output_SKATO<-out_SKATO$results

# MANIPULACION DE OUTPUTS DE MODALIDADES CON PESOS M?S EQUILIBRADOS
output_SKAT_C<-out1CS$results
output_BURDEN_C<-out1CB$results
output_SKAT_W1<-out3$results

names(output_SKAT)<-c("Gene","p_SKAT","N_markers","N_SKAT")
names(output_RARES)<-c("Gene","p_SKAT_rares","N_markers","N_rares")
names(output_BURDEN)<-c("Gene","p_Burden","N_markers","N_Burden")
names(output_SKATO)<-c("Gene","p_SKATO","N_markers","N_SKATO")

names(output_SKAT_W1)<-c("Gene","p_SKATw1","N_markers","N_SKATw1")
names(output_SKAT_C)<-c("Gene","p_SKATC","Q_SKATC","N_markers","N_SKATC","NRARESS","NSNPS")
names(output_BURDEN_C)<-c("Gene","p_BurdenC","Q_BURDENC","N_markers","N_BurdenC","NRARESB","NSNPB")

RESULTS_SKAT1<-Reduce(merge, list(output_SKAT, output_RARES, output_BURDEN, output_SKATO))
RESULTS_SKAT2<-Reduce(merge, list(output_SKAT_W1, output_SKAT_C, output_BURDEN_C))

attach(RESULTS_SKAT1)
RESULTS_SKAT1$min_prob<-pmin(p_SKAT, p_SKAT_rares, p_Burden, p_SKATO, na.rm=T)
attach(RESULTS_SKAT2)
RESULTS_SKAT2$min_prob<-pmin(p_SKATw1, p_SKATC, p_BurdenC, na.rm=T)

TOP_TIPO1<-RESULTS_SKAT1[RESULTS_SKAT1$min_prob<0.000001,]
TOP_TIPO2<-RESULTS_SKAT2[RESULTS_SKAT2$min_prob<0.000001,]
TOP_TIPO1
TOP_TIPO2

write.table(RESULTS_SKAT1, "RESULTS_SKAT1.txt", sep="\t", row.names=FALSE, col.names=TRUE)
write.table(RESULTS_SKAT2, "RESULTS_SKAT2.txt", sep="\t", row.names=FALSE, col.names=TRUE)

Close_SSD()


##############################################################################################################################
## META-ANALISIS SKAT
## DATOS EJEMPLO SIMULANDO 3 ESTUDIOS : son los mismos datos de antes pero divididos en 3 "poblaciones"
## POP1 tiene 626 MUESTRAS: 300 casos y 326 controles
## POP2 tiene 750 MUESTRAS: 500 casos y 250 controles
## POP3 tiene 650 MUESTRAS: 400 casos y 250 controles
## nos vale el mismo archivo .SetID de antes, generado arriba. La informaci?n de marcadores y genes a los que pertenecen no var?a

setwd("/cloud/project/taller SKAT meta")

FAM1_cov<-Read_Plink_FAM_Cov("POP1_SPAINBA_chr22.fam", "COV_POP1_SPAINBA.txt",Is.binary=TRUE, cov_header=TRUE, flag1=0)
FAM2_cov<-Read_Plink_FAM_Cov("POP2_SPAINBA_chr22.fam", "COV_POP2_SPAINBA.txt",Is.binary=TRUE, cov_header=TRUE, flag1=0)
FAM3_cov<-Read_Plink_FAM_Cov("POP3_SPAINBA_chr22.fam", "COV_POP3_SPAINBA.txt",Is.binary=TRUE, cov_header=TRUE, flag1=0)

y1<-FAM1_cov$Phenotype
y2<-FAM2_cov$Phenotype
y3<-FAM3_cov$Phenotype
EDAD1<-FAM1_cov$edad
EDAD2<-FAM2_cov$edad
EDAD3<-FAM3_cov$edad
SEX1<-FAM1_cov$SEX
SEX2<-FAM2_cov$SEX
SEX3<-FAM3_cov$SEX

modelo01<-SKAT_Null_Model(y1 ~ SEX1+EDAD1, out_type="C")
modelo02<-SKAT_Null_Model(y2 ~ SEX2+EDAD2, out_type="C")
modelo03<-SKAT_Null_Model(y3 ~ SEX3+EDAD3, out_type="C")

META1<-Generate_Meta_Files(modelo01, "POP1_SPAINBA_chr22.Bed", "POP1_SPAINBA_chr22.Bim", "SetID_SPAINBA_chr22.SetID",
                           "pop1.MSSD", "pop1.MInfo", N.Sample=626,
                            File.Permu = NULL, data=NULL, impute.method="fixed")

META2<-Generate_Meta_Files(modelo02, "POP2_SPAINBA_chr22.Bed", "POP2_SPAINBA_chr22.Bim", "SetID_SPAINBA_chr22.SetID",
                           "pop2.MSSD", "pop2.MInfo", N.Sample=750,
                            File.Permu = NULL, data=NULL, impute.method="fixed")

META3<-Generate_Meta_Files(modelo03, "POP3_SPAINBA_chr22.Bed", "POP3_SPAINBA_chr22.Bim", "SetID_SPAINBA_chr22.SetID",
                           "pop3.MSSD", "pop3.MInfo", N.Sample=650,
                            File.Permu = NULL, data=NULL, impute.method="fixed")

cohortes<-Open_MSSD_File_2Read(c("pop1.MSSD","pop2.MSSD","pop3.MSSD"), c("pop1.MInfo","pop2.MInfo","pop3.MInfo"))

META_SKAT<-MetaSKAT_MSSD_ALL(cohortes)
META_SKAT1<-MetaSKAT_MSSD_ALL(cohortes, is.separate=TRUE)
?MetaSKAT_MSSD_ALL
names(META_SKAT)
head(META_SKAT)
head(META_SKAT1)
dim(META_SKAT)
META_SKAT[META_SKAT$p.value<0.00001,]
META_SKAT[META_SKAT$SetID=="TRABD",]

## SI EN CADA POBLACI?N/REPLICA REPETIMOS EL PROCESO BASICO PARA VER QU? DAR?A SKAT EN CADA REPLICA:
## GENERAMOS LOS ARCHIVOS .SSD e info 
## ABRIMOS EL FAM Y LAS COV 
## DEFINIMOS EL MODELO PARA DEFINIR EL MODELO

POP1_SpainBA_chr22<-Generate_SSD_SetID("POP1_SPAINBA_chr22.bed", "POP1_SPAINBA_chr22.bim", "POP1_SPAINBA_chr22.fam", "SetID_SPAINBA_chr22.SetID", 
                                       "POP1_SPAINBA22.SSD", "POP1_SPAINBA22.Info", Is.FlipGenotype=TRUE)

POP2_SpainBA_chr22<-Generate_SSD_SetID("POP2_SPAINBA_chr22.bed", "POP2_SPAINBA_chr22.bim", "POP2_SPAINBA_chr22.fam", "SetID_SPAINBA_chr22.SetID", 
                                       "POP2_SPAINBA22.SSD", "POP2_SPAINBA22.Info", Is.FlipGenotype=TRUE)

POP3_SpainBA_chr22<-Generate_SSD_SetID("POP3_SPAINBA_chr22.bed", "POP3_SPAINBA_chr22.bim", "POP3_SPAINBA_chr22.fam", "SetID_SPAINBA_chr22.SetID", 
                                       "POP3_SPAINBA22.SSD", "POP3_SPAINBA22.Info", Is.FlipGenotype=TRUE)

FAM_cov1<-Read_Plink_FAM_Cov("POP1_SPAINBA_chr22.fam", "COV_POP1_SPAINBA.txt",Is.binary=TRUE, cov_header=TRUE, flag1=0)
FAM_cov2<-Read_Plink_FAM_Cov("POP2_SPAINBA_chr22.fam", "COV_POP2_SPAINBA.txt",Is.binary=TRUE, cov_header=TRUE, flag1=0)
FAM_cov3<-Read_Plink_FAM_Cov("POP3_SPAINBA_chr22.fam", "COV_POP3_SPAINBA.txt",Is.binary=TRUE, cov_header=TRUE, flag1=0)

SSD.INFO1<-Open_SSD("POP1_SPAINBA22.SSD", "POP1_SPAINBA22.Info")
out_SKAT1<-SKAT.SSD.All(SSD.INFO1, modelo01)
SSD.INFO2<-Open_SSD("POP2_SPAINBA22.SSD", "POP2_SPAINBA22.Info")
out_SKAT2<-SKAT.SSD.All(SSD.INFO2, modelo02)
SSD.INFO3<-Open_SSD("POP3_SPAINBA22.SSD", "POP3_SPAINBA22.Info")
out_SKAT3<-SKAT.SSD.All(SSD.INFO3, modelo03)
res_POP1<-out_SKAT1$results
res_POP2<-out_SKAT2$results
res_POP3<-out_SKAT3$results
res_POP1[res_POP1$SetID%in%c("TRABD","NOL12","MTMR3","TUBGCP6"),]
res_POP2[res_POP2$SetID%in%c("TRABD","NOL12","MTMR3","TUBGCP6"),]
res_POP3[res_POP3$SetID%in%c("TRABD","NOL12","MTMR3","TUBGCP6"),]
