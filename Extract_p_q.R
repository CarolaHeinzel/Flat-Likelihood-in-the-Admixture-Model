# Input
K = 3 # Number of Populations
file_pathK = "3Island_mi5_mu1_N_f" # Output of STRUCTURE


# function to get the required format of p and q
read_table_data <- function(file_path, K, M) {
  table_data_list <- list()
  for(m in 1:M){
    decimal_vector = c(1,1,1,0,0,0)
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
    }
    text_final = rbind(text_final, text)
  }
  temp = read.table(text = text_final)
  C = ncol(temp)
  liste = c()
  for(i in 1:nrow(temp)){
    if(temp$V1[i] == 0){
      liste =c(liste, i)
    }
  }
  temp = temp[-liste, (C-K+1):C]
  return(temp[2:nrow(temp), ])  
}
extrakt_q <- function(filepath, K){
  lines <- readLines(filepath)
  start_index <- grep("Inferred ancestry of individuals:", lines)
  end_index <- grep("Estimated Allele Frequencies in each cluster",
                    lines)
  extracted_lines <- lines[(start_index + 2):(end_index - 1)]
  res = read.table(text = extracted_lines)
  C = ncol(res)
  res_final = res[,(C-K+1):C]
  res = 0*res_final
  colnames(res_final) = NULL
  
  return(res_final)
}

q = extrakt_q(file_pathK,K)
p = read_table_data(file_pathK, K, 50)

# Save the results
file_name_q <- paste("q_STRUCTURE")
file_name_q <- paste("q_STRUCTURE")


write.table(q, file = file_name, sep = " ", row.names = FALSE, col.names = FALSE)
write.table(p, file = file_name_p, sep = " ", row.names = FALSE, col.names = FALSE)

# Plot of the estimated IAs
barplot(t(as.matrix(matrix_data2)), col=c("purple", "green", "orange", "blue", "red"),
        xlab="Individuals", ylab="Estimated Ancestry", border=NA)

