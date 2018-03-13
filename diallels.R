#############
# diallels
#############

griffing <- function(data){

  data <- data[order(data[,1], data[,2], decreasing = F),]
  
  
#preparing dataset  
traits <- names(data)[3:ncol(data)]
data$reciprocal <- data$female > data$male
data$parent <- data$female == data$male
data$female <- as.factor(data$female)
data$male <- as.factor(data$male)

data$femaleR <- data$male
data$maleR <- data$female
data$combination <- paste(data$female, data$male, sep = "")
data$combination2 <- paste(data$femaleR, data$maleR, sep = "") 
data$combination3 <- data$combination2 
data$combination3[which(data$reciprocal == TRUE)] <- data$combination[data$reciprocal == TRUE]

model <- NULL

solution <- list()

require(EMMREML)

for (i in 1:length(traits)) {

if (length(unique(data$reciprocal)) > 1 & length(unique(data$parent)) > 1) {
#######################################################
# Griffing Method 1 - parents, F1's and reciprocals
model <- "model I"
  
data <- data

Y <- data[,traits[i]]
X <- matrix(rep(1, nrow(data)))

# multi kernel
Z1 <- cbind((model.matrix(~ female -1, data = data) + model.matrix(~ female -1, data = data)))
Z2 <- model.matrix( ~ combination3 -1, data = data)
Z3 <- model.matrix( ~ combination -1, data = data)[,data$parent == FALSE & data$reciprocal == FALSE] - model.matrix( ~ combination -1, data = data)[,match(data$combination[!(data$reciprocal | data$parent)], data$combination2)]
#multikernel
Kp = diag(ncol(Z1))
Kc = diag(ncol(Z2))
Kr = diag(ncol(Z3))

#single kernel
Z <- cbind(Z1, Z2, Z3)

solution <- c(solution, list(traits[i], emmreml(y = as.matrix(Y), X = X, Z = Z, K = diag(ncol(Z)), varbetahat = T,varuhat = T, PEVuhat = T, test = T)))

# multi kernel
#emmremlMultiKernel(y = as.matrix(Y), X = X, Z = list(Z1, Z2, Z3), K = list(Kp, Kc, Kr), varbetahat = T, varuhat = T, PEVuhat = T, test = T)

}


  if (length(unique(data$reciprocal)) == 1 & length(unique(data$parent)) > 1) {  
#############################################
# Griffing Method 2 - parents and F1's
    
    model <- "model II"
    
    data <- data[data$reciprocal == FALSE,]

Y <- data[,traits[i]]
X <- matrix(rep(1, nrow(data)))


#multikernel
Zp <- cbind((model.matrix(~ female -1, data = data) + model.matrix(~ male -1, data = data)))
Zc <- diag(length(Y))
Kp = diag(ncol(Zp))
Kc = diag(ncol(Zc))

# single kernel
Z <- cbind(Zp, Zc)

#single kernel
solution <- c(solution, list(traits[i], emmreml(y = as.matrix(Y), X = X, Z = Z, K = diag(ncol(Z)), varbetahat = T, varuhat = T, PEVuhat = T, test = T)))

# multi kernel
#emmremlMultiKernel(y = as.matrix(Y), X = X, Z = list(Zp, Zc), K = list(Kp, Kc), varbetahat = T, varuhat = T, PEVuhat = T, test = T)

}

  if (length(unique(data$reciprocal)) == 1 & length(unique(data$parent)) == 1) {
#######################################
# Griffing Method 4 - F1's
    
    model <- "model IV"

data <- data[data$reciprocal == FALSE & data$parent == FALSE,]

Y <- data[,traits[i]]
X <- matrix(rep(1, nrow(data)))

# single kernel
Z <- cbind((model.matrix(~ female -1, data = data) + model.matrix(~ male -1, data = data)), model.matrix(~combination -1, data = data))

#multikernel
p <- length(unique(c(data$female, data$male)))
c <- ncol(Z) - p
Zp = as.matrix(Z[,1:p])
Zc = Z[,(p+1):ncol(Z)]
Kp = diag(p)
Kc = diag(c)

# single kernel
solution <- c(solution, list(traits[i], emmreml(y = as.matrix(Y), X = X, Z = Z, K = diag(p+c), varbetahat = T,varuhat = T, PEVuhat = T, test = T)))

# multi kernel
#emmremlMultiKernel(y = as.matrix(Y), X = X, Z = list(Zp, Zc), K = list(Kp, Kc), varbetahat = T, varuhat = T, PEVuhat = T, test = T)

 }

}

outcome <- list(model, solution)

return(outcome)

}





