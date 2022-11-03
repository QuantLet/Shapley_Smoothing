subsets = function(X){
d = ncol(X)
N_subs = 2^d
seq = 1:d
my_list = vector(mode = "list", length = d)
for (k in 1:d){
  my_list[[k]] = combn(seq, m=k)
}

return(my_list)
}







