
missVg <- function(Ne1, Ne2, cycles){
  
  Va0 = 100
  
  decline1 <- function(Va0, Ne1, cycles){
    Vat1 <- Va0*(1 - 1/(2*Ne1))^cycles
    return(Vat1)
  }
  
  decline2 <- function(Va0, Ne2, cycles){
    Vat2 <- Va0*(1 - 1/(2*Ne2))^cycles
    return(Vat2)
  }
  
  miss1 <- decline1(Va0, Ne1, cycles)
  miss2 <- decline2(Va0, Ne2, cycles)
  
  print(c(miss1, miss2))
  
  par(mfrow=c(1,2))
  
  curve(decline1(Va0, Ne1, cycles = x), from = 1, to = 30, ylab = "Genetic variance maintained",
        xlab = "Cycle of selection", col = "red", main = "Group 1")
  # abline(cicles, col = "blue")
  abline(v = cycles, col = "blue")
  text(cycles, miss1, paste(round(miss1, 2), "%", sep = ""), cex = 1, col = "black")
  
  curve(decline2(Va0, Ne2, cycles = x), from = 1, to = 30, ylab = "Genetic variance maintained",
        xlab = "Cycle of selection", col = "red", main = "Group 2")
  # abline(cicles, col = "blue")
  abline(v = cycles, col = "blue")
  text(cycles, miss2, paste(round(miss2, 2), "%", sep = ""), cex = 1, col = "black")
  
}
