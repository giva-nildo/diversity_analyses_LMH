# =========================================================
# LD in autotetraploid from 0-4 dosage matrix
# r_Delta = cor(dose_locus1, dose_locus2)
# r2_Delta = cor(dose_locus1, dose_locus2)^2
# para genótipos não faseados e loci bialélicos
# =========================================================
# Givanildo Rodrigues da Silva
# March 09th, 2026

install.packages("scam")
library(scam)


final_snp_qc_maf_missing_refalt <- readxl::read_excel("final_snp_qc_maf_missing_refalt.xlsx")
G <- metan::column_to_rownames(final_snp_qc_maf_missing_refalt, var = "marker_ID")
Gmx <- G[,c(7:455)]

X <- t(Gmx) # it should be in n (individuals) x markers 
X <- as.data.frame(X)
map_df <- final_snp_qc_maf_missing_refalt[,1:3] # it should be adequately sorted! 
map_df <- map_df %>% 
  relocate(marker_ID)
map_df <- rename(map_df, 
                 marker = marker_ID, 
                 chr = CHROM,
                 pos = POS)

LD_plot_dose <- function(dose,
                         map,
                         max.pair = 1e4,
                         dof = 8,
                         max.loci = 2000, # it was not considered
                         position = c("bp", "cM"),
                         maf.min = 0.05,
                         missing.max = 0.50,
                         sample.seed = 123) {
  
  position <- match.arg(position)
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }
  if (!requireNamespace("scam", quietly = TRUE)) {
    stop("Package 'scam' is required.")
  }
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    stop("Package 'mgcv' is required.")
  }
  
  #-----------------------------
  # 1) Prepare genotype matrix
  #-----------------------------
  X <- as.matrix(dose)
  storage.mode(X) <- "double"
  
  if (is.null(colnames(X))) {
    stop("dose must have marker names in colnames().")
  }
  
  req_cols <- c("marker", "chr", "pos")
  if (!all(req_cols %in% names(map))) {
    stop("map must contain columns: marker, chr, pos")
  }
  
  map <- as.data.frame(map)
  map <- map[match(colnames(X), map$marker), req_cols]
  
  if (any(is.na(map$marker))) {
    stop("Some markers in dose are missing in map.")
  }
  
  #-----------------------------
  # 2) SNP filtering
  #-----------------------------
  miss <- colMeans(is.na(X))
  mean_dose <- colMeans(X, na.rm = TRUE)
  p <- mean_dose / 4
  maf <- pmin(p, 1 - p)
  vv <- apply(X, 2, var, na.rm = TRUE)
  
  keep <- miss <= missing.max & maf >= maf.min & !is.na(vv) & vv > 0
  
  X <- X[, keep, drop = FALSE]
  map <- map[keep, , drop = FALSE]
  
  if (ncol(X) < 2) {
    stop("Too few markers after filtering.")
  }
  
  #-----------------------------
  # 3) Per chromosome LD sampling
  #-----------------------------
  set.seed(sample.seed)
  chroms <- unique(map$chr)
  result_list <- vector("list", length(chroms))
  names(result_list) <- chroms
  
  for (i in seq_along(chroms)) {
    chr_now <- chroms[i]
    ix <- which(map$chr == chr_now)
    
    if (length(ix) < 2) next
    
    if (!is.null(max.loci) && length(ix) > max.loci) {
      ix <- sort(sample(ix, max.loci))
    }
    
    Xc <- X[, ix, drop = FALSE]
    pos <- map$pos[ix]
    
    # remover marcadores sem variância dentro do subconjunto
    vv_chr <- apply(Xc, 2, var, na.rm = TRUE)
    keep_chr <- !is.na(vv_chr) & vv_chr > 0
    
    Xc <- Xc[, keep_chr, drop = FALSE]
    pos <- pos[keep_chr]
    
    m <- ncol(Xc)
    if (m < 2) next
    
    # correlação entre doses = r_Delta
    R2 <- suppressWarnings(stats::cor(Xc, use = "pairwise.complete.obs"))^2
    
    # upper triangle only, excluding diagonal
    ut <- which(upper.tri(R2), arr.ind = TRUE)
    if (nrow(ut) == 0) next
    
    r2_vec <- R2[ut]
    d_vec  <- abs(pos[ut[, 2]] - pos[ut[, 1]])
    
    ok <- !is.na(r2_vec) & !is.na(d_vec)
    r2_vec <- r2_vec[ok]
    d_vec  <- d_vec[ok]
    
    if (length(r2_vec) == 0) next
    
    if (position == "bp") {
      d_vec <- d_vec / 1e6
    }
    
    tmp <- data.frame(
      chr = chr_now,
      d = d_vec,
      r2 = r2_vec
    )
    
    if (!is.null(max.pair) && nrow(tmp) > max.pair) {
      tmp <- tmp[sample(seq_len(nrow(tmp)), max.pair), , drop = FALSE]
    }
    
    result_list[[i]] <- tmp
  }
  
  result <- do.call(rbind, result_list)
  
  if (is.null(result) || nrow(result) < 10) {
    stop("Too few LD pairs available for plotting.")
  }
  
  # opcional: remover d = 0
  result <- result[is.finite(result$d) & is.finite(result$r2), , drop = FALSE]
  
  if (nrow(result) < 10) {
    stop("Too few valid LD pairs after removing NA/Inf.")
  }
  
  #-----------------------------
  # 4) Monotone decreasing spline
  #-----------------------------
  fit <- scam::scam(
    r2 ~ s(d, bs = "mdcx", k = dof),
    data = result
  )
  
  d_grid <- seq(0, max(result$d, na.rm = TRUE), length.out = 500)
  pred <- stats::predict(fit, newdata = data.frame(d = d_grid))
  
  spline_data <- data.frame(d = d_grid, r2 = pred)
  
  #-----------------------------
  # 5) Plot
  #-----------------------------
  p <- ggplot2::ggplot(spline_data, ggplot2::aes(x = d, y = r2)) +
    ggplot2::geom_line() +
    ggplot2::theme_bw() +
    ggplot2::ylab(expression(r[Delta]^2))
  
  if (position == "bp") {
    p <- p + ggplot2::xlab("Distance (Mb)")
  } else {
    p <- p + ggplot2::xlab("Distance (cM)")
  }
  
  return(list(
    plot = p,
    sampled_pairs = result,
    spline_data = spline_data,
    fit = fit
  ))
}

ans <- LD_plot_dose(
  dose = X,
  map = map_df,
  max.loci = NULL, #all markers per chromosome
  position = "bp",
  maf.min = 0.01,
  missing.max = 0.50
)

ans$plot

ans$plot + coord_cartesian(xlim = c(0, 40))

r_data <- ans$sampled_pairs
