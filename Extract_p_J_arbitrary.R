K = 5
file_pathK = "/home/ch1158/DatenSTRUCTURE/STRUCTURESoftware/console/example_j3_f"

read_table_data <- function(file_path, K, M) {
  table_data_list <- list()
  for(m in 1:M){
    decimal_vector = rep(1,K)
    table_data_list[[length(table_data_list) + 1]] <- matrix(decimal_vector, ncol = K, byrow = TRUE)
  }
  lines <- readLines(file_path)
  start_indices <- grep("0.0% missing data", lines)
  for (j in 1:length(start_indices)) {
    start_index <- start_indices[j]
    end_index <- start_index + which(lines[start_index:length(lines)] == "")[1] - 1
    text <- paste(lines[(start_index+1):(end_index-1)], collapse = "\n")
    if(j == 1){
      text_final = text
    }else{
      text_final = rbind(text_final, text)
    }
  }
  temp = read.table(text = text_final)
  liste = c()
  for(i in 1:nrow(temp)){
    if(temp$V1[i] == 0){
      liste =c(liste, i)
    }
  }
  temp = temp[-liste, ]
  return(temp[1:nrow(temp), 3:ncol(temp)])  
}

p = read_table_data(file_pathK, K, 55)

