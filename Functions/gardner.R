#############
# diallels
#############

gardner <- function(data){

#preparing dataset  
  data <- data[order(data[,1], data[,2], decreasing = F),]
  traits <- names(data)[3:ncol(data)]
  data$parent <- data$female == data$male
  data$female <- as.factor(data$female)
  data$male <- as.factor(data$male)
  data$combination <- paste(data$female, data$male, sep = "")
  
  X <- matrix(rep(1, nrow(data)))
  Z1 <- model.matrix(~female -1, data = data)/2 + model.matrix(~male -1, data = data)/2 
  Z2 <- as.numeric(data$female != data$male) 
  Z3 <- model.matrix(~female -1, data = data) + model.matrix(~male -1, data = data) 
  Z3[Z3 == 2] <- 0  
  Z4 <- model.matrix(~ combination -1, data = data)[,!data$parent]
  Z4[data$parent,] <- 0

  solution <- list()

require(EMMREML)

for (i in 1:length(traits)) {

Y <- data[,traits[i]]

#multikernel
K1 = diag(ncol(Z1))
K2 = diag(ncol(Z2))
K3 = diag(ncol(Z3))
K4 = diag(ncol(Z4))

#single kernel
model1 <- emmreml(y = as.matrix(Y), X = X, Z = Z1, K = diag(ncol(Z1)), varbetahat = T,varuhat = T, PEVuhat = T, test = T)
model2 <- emmreml(y = as.matrix(Y), X = X, Z = cbind(Z1, Z2), K = diag(ncol(cbind(Z1, Z2))), varbetahat = T,varuhat = T, PEVuhat = T, test = T)
model3 <- emmreml(y = as.matrix(Y), X = X, Z = cbind(Z1, Z2, Z3), K = diag(ncol(cbind(Z1, Z2, Z3))), varbetahat = T,varuhat = T, PEVuhat = T, test = T)
model4 <- emmreml(y = as.matrix(Y), X = X, Z = cbind(Z1, Z2, Z3, Z4), K = diag(ncol(cbind(Z1, Z2, Z3, Z4))), varbetahat = T,varuhat = T, PEVuhat = T, test = T)

all.models <- list(model1, model2, model3, model4)

solution <- c(solution, list(traits[i], all.models))

}

return(solution)

}

