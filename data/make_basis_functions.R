library(maps)
library(sf)
library(sp)
library(spdep)
library(dplyr)

##BANERJEE BOOK reference bfs
us_state <- map("state", fill=TRUE, plot=FALSE)
us_state_sf <- st_as_sf(us_state)
sf_use_s2(FALSE)
us_nb <- poly2nb(us_state_sf)
W <- nb2mat(us_nb, style="B")
#W <- W/rowSums(W)
W_eig <- eigen(W)
banerjee_bf <- W_eig$vectors[,which(W_eig$values > 0)]
rownames(banerjee_bf) <- rownames(W)

## map state names -> state abbreviations
name_to_abb <- setNames(state.abb, tolower(state.name))
name_to_abb <- append(name_to_abb, "DC", after = 8)
names(name_to_abb)[9] <- "district of columbia"
## get abbreviations for W_ref ordering
abb <- name_to_abb[rownames(W)]

## reorder rows and columns
W_new <- W[order(abb), order(abb)]
W_new <- W_new/rowSums(W_new)
W_eig_new <- eigen(W_new)
scaled_bf <- W_eig_new$vectors[,which(W_eig_new$values > 0)]

basis_funcs <- scaled_bf
save(basis_funcs, file="data/scaled_basis_functions.RData")

W_new <- W[order(abb), order(abb)]
W_eig_new <- eigen(W_new)
unscaled_bf <- W_eig_new$vectors[,which(W_eig_new$values > 0)]

basis_funcs <- unscaled_bf
save(basis_funcs, file="data/unscaled_basis_functions.RData")
