K = 3 # Number of populations
M = 100 # Number of markers
file_pathK = "CEU_IBS_TSI_enhanced_corr_f" # File Name

read_table_data <- function(file_path, K, M) {
  lines <- readLines(file_path)
  start_indices <- grep("0.0% missing data", lines)
  res = matrix(nrow = 1, ncol = K)
  for (j in 1:length(start_indices)) {
    start_index <- start_indices[j] + 1
    end_index <- start_index + which(lines[start_index:length(lines)] == "")[1] - 1
    text <- paste(lines[(start_index+1):(end_index-1)], collapse = "\n")
    
    if(start_index - end_index < -2){
      modified_text <- gsub("\\s+[0-9]+\\s*\\(0\\.[0-9]+\\)", "", text)
      values <- strsplit(modified_text, " ")[[1]]
      values <- values[values != ""]
      numeric_values <- as.numeric(values)
      res1 = numeric_values[1:K] + numeric_values[(K+1):(2*K)]
      print(res1)
      res = rbind(res,res1)
    }
    if(j == 1){
      text_final = text
    }else{
      text_final = rbind(text_final, text)
    }
  }
  temp = read.table(text = text_final)
  rownames(res) = NULL
  return(list(temp[1:nrow(temp), 3:ncol(temp)], res[2:nrow(res),]))  
}

p = read_table_data(file_pathK, K, M)


output_directory <- "/home" # Change to your output directory

file_name <- paste("p_CEU_IBS_TSI_K2")

output_filepath <- file.path(output_directory, file_name)

write.table(p[[1]], file = output_filepath, sep = " ", row.names = FALSE, col.names = FALSE)


file_name <- paste("p_CEU_IBS_TSI_K2_J")

output_filepath <- file.path(output_directory, file_name)

write.table(p[[2]], file = output_filepath, sep = " ", row.names = FALSE, col.names = FALSE)
