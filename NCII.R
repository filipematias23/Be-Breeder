#########################################
# NCII or factorial (interpopulational) 
#######################################

factorial <- function(data){

data <- data[order(data[,1], data[,2], decreasing = F),]
trait <- colnames(data)[4:ncol(data)]
data$fm <- paste(data$female, data$male, sep = "")


NCII <- function(data, trait){
  
  library(breedR)
  
  soluction <- remlf90(fixed = data[,trait] ~ 1, random = ~ rep + female + male + fm, data = data, method = "em")
  
  soluction2 <- list(trait, 
                     soluction[[7]],
                     soluction[[8]],
                     soluction[[10]][1])
  
  return(soluction2)
  
}

result <- list()

for (i in 1:length(trait)) {

sol <- NCII(data, trait[i])
 
result <- c(result, list(sol))

 }

return(result)

}






