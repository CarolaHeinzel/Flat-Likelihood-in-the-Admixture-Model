

input_directory1 = "/home/ch1158/Ergebnisse_STRUCTURE/SAS_EAS/Pong_ohne_Marker/FIN_IBS_ohne2.txt"
input_directory2 = "/home/ch1158/Ergebnisse_STRUCTURE/SAS_EAS/Pong_ohne_Marker/FIN_IBS_ohne7.txt"
input_directory3 = "/home/ch1158/Ergebnisse_STRUCTURE/SAS_EAS/Pong_ohne_Marker/FIN_IBS_ohne15.txt"
input_directory4 = "/home/ch1158/Ergebnisse_STRUCTURE/SAS_EAS/Pong_ohne_Marker/FIN_IBS_ohne9.txt"
file_path <- file.path("/home/ch1158/Ergebnisse_STRUCTURE/IBS_GBR/r8")
matrix_data_r1 = read.table(input_directory1)
matrix_data_r8 = read.table(input_directory2)
matrix_data_r17 = read.table(input_directory3)
matrix_data_r4 = read.table(input_directory4)


m1 = matrix_data[1:99,1]
m2 = matrix_data_r8[100:206,1]
order_indices <- order(m1)
order_indices1 <- 99 + order(m2)


q1 = matrix_data_r1[,1]
res_max = rep(0, length(q1))
res_min = rep(0, length(q1))

tbl1 <- matrix(c(q1,1-q1),ncol = 2)

t = as.matrix(q1[1:99])
t
t1 = as.matrix(q1[100:206])
sort(t)
t2 = c(sort(t), sort(t1))
tbl1[,1] =t2
tbl1[,2] = (1-tbl1[,1])
a_max =  -0.05249863
a_min = 0.1027696
b_max = 0.4490661
b_min = -0.07908124

for (i in 1:length(q1)){
  res_max[i] = tbl1[i,1]*(1-a_max) + (1-tbl1[i,1])*b_max
  res_min[i] = tbl1[i,1]*(1-a_min) + (1-tbl1[i,1])*b_min
}

x_values <- seq_along(matrix_data_r1[,1])

m1 = matrix_data_r1[1:99,1]
m2 = matrix_data_r1[100:206,1]
order_indices <- order(m1)
order_indices1 <- 99 + order(m2)


# Plotten der Punkte

list1_sorted <- c(matrix_data_r1[order_indices1,1], matrix_data_r1[order_indices,1])
list2_sorted <-  c(matrix_data_r8[order_indices1,1], matrix_data_r8[order_indices,1])
list3_sorted <-  c(matrix_data_r17[order_indices1,2], matrix_data_r17[order_indices,2])
list4_sorted <-  c(matrix_data_r4[order_indices1,1], matrix_data_r4[order_indices,1])
plot(x_values, res_max, type = "l", 
     xlab = "Individuals", ylab = "Estimated Ancestry", 
     pch = 19, col = "darkred", ylim = c(0, 1), xlim = c(1, 206),xaxt = 'n')
abline(v = 99.5, col = "black", lwd = 3)
lines(x_values, res_min,col = "red")

lines(x_values, list2_sorted,col = "lightblue")
lines(x_values, list3_sorted,col = "1E90FF")
lines(x_values, list1_sorted,col = "darkblue")
lines(x_values, tbl1[,1],col = "blue")

polygon(c(x_values, rev(x_values)), c(res_max, rev(res_min)), col = rgb(1, 0, 0, 0.3), border = NA)

mtext("FIN", side = 1, line = 1, adj = 0.25)
mtext("IBS", side = 1, line = 1, adj = 0.75)
