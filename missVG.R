missVg <- function(Ne, cycles){
  
  Va0 = 100
  
  decline <- function(Va0, Ne, cycles){
    Vat <- Va0*(1 - 1/(2*Ne))^cycles
    return(Vat)
  }
  
  miss <- decline(Va0, Ne, cycles)
  
  print(miss)
  
  curve(decline(Va0, Ne, cycles = x), from = 1, to = 30, ylab = "Genetic variance maintained", xlab = "Cycle of selection", col = "red")
  # abline(cicles, col = "blue")
  abline(v = cycles, col = "blue")
  text(cycles, miss, paste(round(miss, 2), "%", sep = ""), cex = 1, col = "black")
  
}