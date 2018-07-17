library("rJava")
library("XLConnect")

filenames <- list.files(".", pattern="*.xlsx")

total_v <- vector()

for (i in 1:length(filenames)){
  
  files <- readWorksheetFromFile(filenames[i], sheet=1, startRow=7, endRow=25)
  labID <- c(files[1:8,1])
  labID.1 <- c(files[1:8,8])
  labID.2 <- c(files[1:8,15])
  
  labID.3 <- c(files[11:18,1])
  labID.4 <- c(files[11:18,8])
  labID.5 <- c(files[11:18,15])
  
  v <- c(labID, labID.1, labID.2, labID.3, labID.4, labID.5)
  
  total_v <- c(total_v, v)

    
}

for (i in 1:length(filenames)){
  
  files2 <- readWorksheetFromFile(filenames[i], sheet=2, startRow=7, endRow=25)
  labID <- c(files2[1:8,1])
  labID.1 <- c(files2[1:8,8])
  labID.2 <- c(files2[1:8,15])
  
  labID.3 <- c(files2[11:18,1])
  labID.4 <- c(files2[11:18,8])
  labID.5 <- c(files2[11:18,15])
  
  v2 <- c(labID, labID.1, labID.2, labID.3, labID.4, labID.5)
  
  total_v <- c(total_v, v2)
  
  
}

total_v <- total_v[!is.na(total_v)]
write.table(total_v, file="out.txt", quote=FALSE, col.names=F, row.names=F)

