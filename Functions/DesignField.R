
# DIC

trt <- c("A", "B", "C")
repeticion <- c(4, 3, 4)
outdesign <- design.crd(trt,r=repeticion,seed=777,serie=0)
book1 <- outdesign$book
head(book1)

# DBC

trt <- c("A", "B", "C","D","E")
repeticion <- 4
outdesign <- design.rcbd(trt,r=repeticion, seed=-513, serie=2)
book2<- zigzag(outdesign) 
return(book2)

# Quadrado Latino

trt <- c("A", "B", "C", "D")
outdesign <- design.lsd(trt, seed=543, serie=2)
print(outdesign$sketch)

book <- zigzag(outdesign)
print(matrix(book[,1],byrow = TRUE, ncol = 4))

# Lattice designs

trt<-letters[1:9]
outdesign <-design.lattice(trt, r = 3, serie = 2, seed = 33, kinds = "Super-Duper")
book7 <- outdesign$book
outdesign$parameters
outdesign$sketch

# Balanced Incomplete Block Designs

trt <- c("A", "B", "C", "D", "E" )
k <- 4
outdesign <- design.bib(trt,k, seed=543, serie=2)

book5 <- outdesign$book
outdesign$statistics
outdesign$sketch

# Alpha designs

trt <- letters[1:15]
outdesign <- design.alpha(trt,k=3,r=2,seed=543)
book8 <- outdesign$book
outdesign$statistics
outdesign$sketch

# Augmented block designs

trt1 <- c("A", "B", "C", "D")
trt2 <- c("t","u","v","w","x","y","z")
outdesign <- design.dau(trt1, trt2, r=5, seed=543, serie=2)
book9 <- outdesign$book
with(book9,by(trt, block,as.character))
book <- zigzag(outdesign)
with(book,by(plots, block, as.character))

# Split plot designs

trt1<-c("A","B","C","D")
trt2<-c("a","b","c")
outdesign <-design.split(trt1,trt2,r=3,serie=2,seed=543)
book10 <- outdesign$book
head(book10)

# Factorial

trt<-c(3,2) # factorial 3x2
outdesign <-design.ab(trt, r=3, serie=2)
book12 <- outdesign$book
head(book12)


