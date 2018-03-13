##########################
### GWASpoly functions ###
##########################

GWASpoly<-function (data, models, traits = NULL, params = NULL, n.core = 1, 
                    quiet = F) 
{
  stopifnot(inherits(data, "GWASpoly.K"))
  if (is.null(params)) {
    params <- set.params()
  }
  if (!is.null(params$fixed)) {
    stopifnot(is.element(params$fixed, colnames(data@fixed)))
  }
  if (is.null(traits)) {
    traits <- colnames(data@pheno)[-1]
  }
  stopifnot(is.element(traits, colnames(data@pheno)[-1]))
  geno.gid <- rownames(data@geno)
  m <- nrow(data@map)
  data@K <- data@K[geno.gid, geno.gid]
  n.gid <- length(geno.gid)
  params$models <- models
  n.model <- length(models)
  dom.models <- grep("dom", models, fixed = T)
  if (length(dom.models) > 0) {
    dom.orders <- as.integer(substr(models[dom.models], 
                                    1, 1))
    if (max(dom.orders) > data@ploidy/2) {
      stop("Maximum dominance model is ploidy/2")
    }
    dom.models2 <- sort(c(paste(models[dom.models], "ref", 
                                sep = "-"), paste(models[dom.models], "alt", sep = "-")))
  }
  else {
    dom.models2 <- character(0)
  }
  other.models <- setdiff(1:n.model, dom.models)
  if (!all(is.element(models[other.models], c("additive", 
                                              "general", "diplo-general", "diplo-additive")))) {
    stop("Invalid model")
  }
  models <- c(models[other.models], dom.models2)
  n.model <- length(models)
  n.trait <- length(traits)
  if (params$n.PC > 0) {
    eig.vec <- eigen(data@K)$vectors
  }
  all.scores <- vector("list", n.trait)
  names(all.scores) <- traits
  all.effects <- all.scores
  for (j in 1:n.trait) {
    trait <- traits[j]
    if (!quiet) {
      cat(paste("Analyzing trait:", trait, "\n"))
    }
    not.miss <- which(!is.na(data@pheno[, trait]))
    y <- data@pheno[not.miss, trait]
    pheno.gid <- data@pheno[not.miss, 1]
    n <- length(y)
    Z <- matrix(0, n, n.gid)
    Z[cbind(1:n, match(pheno.gid, geno.gid))] <- 1
    X <- matrix(1, n, 1)
    if (!is.null(params$fixed)) {
      for (i in 1:length(params$fixed)) {
        if (params$fixed.type[i] == "factor") {
          xx <- factor(data@fixed[not.miss, params$fixed[i]])
          if (length(levels(xx)) > 1) {
            X <- cbind(X, model.matrix(~x, data.frame(x = xx))[, 
                                                               -1])
          }
        }
        else {
          X <- cbind(X, data@fixed[not.miss, params$fixed[i]])
        }
      }
    }
    if (params$n.PC > 0) {
      X <- cbind(X, Z %*% eig.vec[, 1:params$n.PC])
    }
    X2 <- .make.full(X)
    if (params$P3D) {
      if (!quiet) {
        cat("P3D approach: Estimating variance components...")
      }
      Hinv <- mixed.solve(y = y, X = X2, Z = Z, K = data@K, 
                          return.Hinv = TRUE)$Hinv
      if (!quiet) {
        cat("Completed \n")
      }
    }
    else {
      Hinv <- NULL
    }
    scores <- matrix(NA, m, n.model)
    colnames(scores) <- models
    rownames(scores) <- colnames(data@geno)
    betas <- scores
    for (k in 1:n.model) {
      if (!quiet) {
        cat(paste("Testing markers for model:", models[k], 
                  "\n"))
      }
      if ((n.core > 1) & requireNamespace("parallel", 
                                          quietly = TRUE)) {
        it <- split(1:m, factor(cut(1:m, n.core, labels = FALSE)))
        score.list <- parallel::mclapply(it, .score.calc, 
                                         y, Z, X2, data@K, data@geno, Hinv, data@ploidy, 
                                         models[k], params$min.MAF, params$max.geno.freq, 
                                         mc.cores = n.core)
        scores[, k] <- unlist(lapply(score.list, function(el) {
          el$score
        }))
        betas[, k] <- unlist(lapply(score.list, function(el) {
          el$beta
        }))
      }
      else {
        ans <- .score.calc(1:m, y, Z, X2, data@K, data@geno, 
                           Hinv, data@ploidy, models[k], params$min.MAF, 
                           params$max.geno.freq)
        scores[, k] <- ans$score
        betas[, k] <- ans$beta
      }
    }
    all.scores[[j]] <- data.frame(scores, check.names = F)
    all.effects[[j]] <- data.frame(betas, check.names = F)
  }
  return(new("GWASpoly.fitted", data, scores = all.scores, 
             effects = all.effects, params = params))
}

read.GWASpoly<-function (ploidy, pheno.file, geno.file, format, n.traits, delim = ",") 
{
  if (format == "ACTG") {
    format <- "ACGT"
  }
  if (!is.element(format, c("AB", "numeric", "ACGT"))) {
    stop("Invalid genotype format.")
  }
  bases <- c("A", "C", "G", "T")
  get.ref <- function(x, format) {
    if (format == "numeric") {
      ref.alt <- c(0, 1)
    }
    if (format == "AB") {
      ref.alt <- c("A", "B")
    }
    if (format == "ACGT") {
      y <- paste(na.omit(x), collapse = "")
      ans <- apply(array(bases), 1, function(z, y) {
        length(grep(z, y, fixed = T))
      }, y)
      if (sum(ans) > 2) {
        stop("Error in genotype matrix: More than 2 alleles")
      }
      if (sum(ans) == 2) {
        ref.alt <- bases[which(ans == 1)]
      }
      if (sum(ans) == 1) {
        ref.alt <- c(bases[which(ans == 1)], NA)
      }
    }
    return(ref.alt)
  }
  geno <- read.table(file = geno.file, header = T, as.is = T, 
                     check.names = F, sep = delim)
  map <- data.frame(Marker = geno[, 1], Chrom = factor(geno[, 
                                                            2], ordered = T), Position = geno[, 3], stringsAsFactors = F)
  markers <- as.matrix(geno[, -(1:3)])
  rownames(markers) <- geno[, 1]
  tmp <- apply(markers, 1, get.ref, format)
  map$Ref <- tmp[1, ]
  map$Alt <- tmp[2, ]
  if (is.element(format, c("AB", "ACGT"))) {
    M <- apply(cbind(map$Ref, markers), 1, function(x) {
      y <- gregexpr(pattern = x[1], text = x[-1], fixed = T)
      ans <- as.integer(lapply(y, function(z) {
        ifelse(z[1] < 0, ploidy, ploidy - length(z))
      }))
      return(ans)
    })
  }
  else {
    M <- t(markers)
  }
  gid.geno <- colnames(geno)[-(1:3)]
  rownames(M) <- gid.geno
  bad <- length(which(!is.element(na.omit(M), 0:ploidy)))
  if (bad > 0) {
    stop("Invalid marker calls.")
  }
  MAF <- apply(M, 2, function(x) {
    AF <- mean(x, na.rm = T)/ploidy
    MAF <- ifelse(AF > 0.5, 1 - AF, AF)
  })
  polymorphic <- which(MAF > 0)
  M <- M[, polymorphic]
  map <- map[polymorphic, ]
  map <- map[order(map$Chrom, map$Position), ]
  M <- M[, map$Marker]
  m <- nrow(map)
  cat(paste("Number of polymorphic markers:", m, "\n"))
  impute.mode <- function(x) {
    ix <- which(is.na(x))
    if (length(ix) > 0) {
      x[ix] <- as.integer(names(which.max(table(x))))
    }
    return(x)
  }
  missing <- which(is.na(M))
  if (length(missing) > 0) {
    cat("Missing marker data imputed with population mode \n")
    M <- apply(M, 2, impute.mode)
  }
  pheno <- read.table(file = pheno.file, header = T, as.is = T, 
                      check.names = F, sep = delim)
  gid.pheno <- unique(pheno[, 1])
  gid <- intersect(gid.pheno, gid.geno)
  pheno <- pheno[is.element(pheno[, 1], gid), ]
  M <- M[gid, ]
  N <- length(gid)
  cat(paste("N =", N, "individuals with phenotypic and genotypic information \n"))
  n.fixed <- ncol(pheno) - n.traits - 1
  if (n.fixed > 0) {
    fixed <- data.frame(pheno[, (n.traits + 2):ncol(pheno)], 
                        stringsAsFactors = F)
    fixed.names <- colnames(pheno)[(n.traits + 2):ncol(pheno)]
    colnames(fixed) <- fixed.names
    pheno <- data.frame(pheno[, 1:(1 + n.traits)], stringsAsFactors = F)
    cat(paste("Detected following fixed effects:\n", paste(fixed.names, 
                                                           collapse = "\n"), "\n", sep = ""))
  }
  else {
    fixed <- data.frame(NULL)
  }
  traits <- colnames(pheno)[-1]
  cat(paste("Detected following traits:\n", paste(traits, 
                                                  collapse = "\n"), "\n", sep = ""))
  return(new("GWASpoly", map = map, pheno = pheno, fixed = fixed, 
             geno = M, ploidy = ploidy))
}

set.threshold<-function (data, method, level = 0.05, n.permute = 1000, n.core = 1) 
{
  stopifnot(inherits(data, "GWASpoly.fitted"))
  traits <- names(data@scores)
  n.trait <- length(traits)
  models <- colnames(data@scores[[1]])
  n.model <- length(models)
  methods <- c("Bonferroni", "FDR", "permute")
  stopifnot(is.element(method, methods))
  threshold <- matrix(NA, n.trait, n.model)
  colnames(threshold) <- models
  rownames(threshold) <- traits
  for (i in 1:n.trait) {
    trait <- traits[i]
    if (method == "permute") {
      print(paste("Trait:", trait), quote = F)
      y <- data@pheno[, trait]
      ix <- which(!is.na(y))
      max.scores <- matrix(NA, n.permute, n.model)
      colnames(max.scores) <- models
      for (q in 1:n.permute) {
        print(paste("Permutation", q), quote = F)
        data2 <- data
        data2@pheno[ix, trait] <- sample(y[ix])
        data2 <- GWASpoly(data2, models = data@params$models, 
                          traits = trait, params = data@params, quiet = T, 
                          n.core = n.core)
        for (j in 1:n.model) {
          max.scores[q, j] <- max(data2@scores[[trait]][, 
                                                        models[j]], na.rm = T)
        }
      }
    }
    for (j in 1:n.model) {
      model <- models[j]
      scores <- as.vector(na.omit(data@scores[[trait]][, 
                                                       model]))
      m <- length(scores)
      if (method == "Bonferroni") {
        threshold[i, j] <- -log10(level/m)
      }
      if (method == "FDR") {
        tmp <- cbind(10^(-scores), .qvalue(10^(-scores)))
        tmp <- tmp[order(tmp[, 2]), ]
        if (tmp[1, 2] > level) {
          threshold[i, j] <- -log10(tmp[1, 1]) * 1.2
        }
        else {
          k <- max(which(tmp[, 2] < level))
          threshold[i, j] <- -log10(mean(tmp[k:(k + 
                                                  1), 1]))
        }
      }
      if (method == "permute") {
        threshold[i, j] <- sort(max.scores[, model], 
                                decreasing = TRUE)[floor(level * n.permute)]
      }
    }
  }
  return(new("GWASpoly.thresh", data, threshold = threshold))
}

get.QTL<-function (data, traits = NULL, models = NULL) 
{
  stopifnot(inherits(data, "GWASpoly.thresh"))
  if (is.null(traits)) {
    traits <- names(data@scores)
  }
  else {
    stopifnot(is.element(traits, names(data@scores)))
  }
  if (is.null(models)) {
    models <- colnames(data@scores[[1]])
  }
  else {
    stopifnot(is.element(models, colnames(data@scores[[1]])))
  }
  n.model <- length(models)
  n.trait <- length(traits)
  output <- data.frame(NULL)
  for (i in 1:n.trait) {
    for (j in 1:n.model) {
      ix <- which(data@scores[[traits[i]]][, models[j]] > 
                    data@threshold[traits[i], models[j]])
      n.ix <- length(ix)
      output <- rbind(output, data.frame(Trait = rep(traits[i], 
                                                     n.ix), Model = rep(models[j], n.ix), Threshold = round(rep(data@threshold[traits[i], 
                                                                                                                               models[j]], n.ix), 2), data@map[ix, ], Score = round(data@scores[[traits[i]]][ix, 
                                                                                                                                                                                                             models[j]], 2), Effect = round(data@effects[[traits[i]]][ix, 
                                                                                                                                                                                                                                                                      models[j]], 2), stringsAsFactors = F, check.names = F))
    }
  }
  return(output)
}
