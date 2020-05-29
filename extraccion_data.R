library(dplyr)

# cargar datos

counts = read.csv("counts.csv",header=T,sep = ";")
target = read.csv("targets.csv",header=T,sep = ",")


# primero extraer las que son del mismo grupo a variables separadas

NIT <- subset(target,grepl("^(NIT)", target$Group))
SFI <- subset(target,grepl("^(SFI)", target$Group))
ELI <- subset(target,grepl("^(ELI)", target$Group))

# seleccionar aleatoriamente las 10 muestras de cada grupo

set.seed(12)

NIT_10 <- sample_n(NIT, 10)
SFI_10 <- sample_n(SFI, 10)
ELI_10 <- sample_n(ELI, 10)


# extraer con las variables anteriores los datos con los que se va a trabajar

Grupo1 <- counts[,NIT_10$Sample_Name]
Grupo2 <- counts[,SFI_10$Sample_Name]
Grupo3 <- counts[,ELI_10$Sample_Name]

datos <- cbind(Grupo1,Grupo2,Grupo3)

write.csv(datos,"datos.csv", row.names = FALSE)