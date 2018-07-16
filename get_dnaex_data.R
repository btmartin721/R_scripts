library("rJava")
library("XLConnect")

filenames <- list.files(".", pattern="*.xls")

df_total <- data.frame()

for (i in 1:length(filenames)){
  
  files <- readWorksheetFromFile(filenames[i], sheet=1, startRow=11, endRow=33)
  popID <- tail(files$Col3, -1)
  indID <- tail(files$Col4, -1)
   
  fieldID <- tail(files$field.., -1)

  tissue <- tail(files$tissue, -1)

  dna <- tail(files$DNA, -1)
  
  df <- data.frame(popID, indID, fieldID, tissue, dna, stringsAsFactors = FALSE)

  df_total <- rbind(df_total, df)

}

write.csv(df_total, file="out.csv", row.names=FALSE)
