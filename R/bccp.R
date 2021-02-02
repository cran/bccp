rtype1 <- function(n, P, T, param, mle, cdf.expression = FALSE, pdf.expression = TRUE, cdf, pdf, lb = 0)
{
  d <- length(mle)
  m <- length(T)
  if( length(param) != d ) stop("The length of parameter vector and ML estimators must be the same.")
  if( length(P) != m ) stop("The length of censoring times and percent of removals must be the same.")
  xsum <- rsum <- X <- R <- D <- C <- rep( NA, m)
  for(k in 1:d) assign(param[k], mle[k])
  if (cdf.expression == TRUE)
  {
    CDF <- function(x){}
    for(k in 1:d) assign(param[k], mle[k])
    body(CDF) <- bquote(.(cdf))
    if ( is.nan( CDF(lb) ) == TRUE ) stop( "Try for another choice of lower bound T0." )
    C[1] <- CDF(T[1]) - CDF(lb)
    X[1] <- rbinom(1, n, C[1])
    R[1] <- floor(P[1]*(n-X[1]))
    xsum[1] <- X[1]
    rsum[1] <- R[1]
    for (i in 2:m)
    {
      D1 <- CDF(T[i]) - CDF(T[i-1])
      D2 <- CDF(T[i-1]) - CDF(lb)
      D[i] <- D1/(1-D2)
      X[i] <- suppressWarnings( rbinom(1, n-xsum[i-1]-rsum[i-1], D[i]) )
      R[i] <- suppressWarnings( floor( P[i]*(n-xsum[i-1]-rsum[i-1]-X[i]) ) )
      xsum[i] <- xsum[i-1] + X[i]
      rsum[i] <- rsum[i-1] + R[i]
    }
  }
  if (pdf.expression == TRUE)
  {
    integrand <- function(x){}
    body(integrand) <- bquote(.(pdf))
    if (is.nan( integrand(lb) ) == TRUE ) stop( "Try for another choice of lower bound T0." )
    C[1] <- quadinf(integrand, lb, T[1])$Q
    X[1] <- rbinom(1, n, C[1])
    R[1] <- floor(P[1]*(n-X[1]))
    xsum[1] <- X[1]
    rsum[1] <- R[1]
    for (i in 2:m)
    {
      D1 <- quadinf(integrand, T[i-1], T[i])$Q
      D2 <- quadinf(integrand, lb, T[i-1])$Q
      D[i] <- D1/(1-D2)
      X[i] <- suppressWarnings( rbinom(1, n-xsum[i-1]-rsum[i-1], D[i]) )
      R[i] <- suppressWarnings( floor( P[i]*(n-xsum[i-1]-rsum[i-1]-X[i]) ) )
      xsum[i] <- xsum[i-1] + X[i]
      rsum[i] <- rsum[i-1] + R[i]
    }
  }
  R[is.na(R)] <- 0
  X[is.na(X)] <- 0
  return( plan = data.frame(T = T, X = X, R = R, P = P) )
}
###################################################################################################
fitype1 <- function(plan, param, mle, cdf.expression = FALSE, pdf.expression = TRUE, cdf, pdf, lb = 0)
{
  d <- length(mle)
  T <- plan$T
  X <- plan$X
  R <- plan$R
  P <- plan$P
  m <- length(T)
  n <- sum(X+R)
  K <- rep(NA, m)
  D2_e <- D2_o <- matrix (NA, nrow = d, ncol = d)
  if(cdf.expression == TRUE)
  {
    d1_CDF <- d1_CDFc <- d1i_CDF <- d1i_CDFc <- matrix(NA, nrow = d, ncol = m)
    d2_CDF <- d2_CDFc <- d2i_CDF <- array(NA, dim = c(d, m, d) )
    CDF_0 <- CDF_1 <- CDF_2 <- function(x){}
    for(k in 1:d) assign(param[k], mle[k])
    body(CDF_0) <- bquote(.(cdf))
    if ( is.nan( CDF_0(lb) ) == TRUE ) stop( "Try for another choice of lower bound lb" )
    C0 <- CDF_0(T)
    d0_CDF <- diff( c( CDF_0(lb), C0) )
    d0_CDFc <- 1-C0
    for(i in 1:m) K[i] <- prod( (1-P[1:i-1]) )
    E_X <- d0_CDF*K
    E_R <- P*d0_CDFc*K
    first <- sapply(1:d, function(i) D(cdf, param[i]))
    for (i in 1:d)
    {
      body(CDF_1) <- bquote(.(first[[i]]))
      if (is.nan( CDF_1(lb) ) == TRUE ) stop( "Try for another choice of lower bound lb" )
      C1 <- CDF_1(T)
      d1_CDF[i,] <- diff( c( CDF_1(lb), C1 ) )
      d1i_CDF[i,] <-  C1
      d1_CDFc[i,] <- -C1
    }
    for (i in 1:d)
    {
      for(j in i:d)
      {
        second <- D(first[[i]], param[j])
        body(CDF_2) <- bquote(.(second))
        if (is.nan( CDF_2(lb) ) == TRUE ) stop( "Try for another choice of lower bound lb" )
        C2 <- CDF_2(T)
        d2_CDF[i,,j] <- diff( c( CDF_2(lb), C2 ) )
        d2_CDF[j,,i] <- d2_CDF[i,,j]
        d2_CDFc[i,,j] <- -C2
        d2_CDFc[j,,i] <- d2_CDFc[i,,j]
        d2i_CDF[i,,j] <- C2
        d2i_CDF[j,,i] <- d2i_CDF[i,,j]
        D2_o[i,j] <- sum( -X*d2_CDF[i,,j]/d0_CDF + X*d1_CDF[i,]*d1_CDF[j,]/d0_CDF^2 +
                           R*d2i_CDF[i,,j]/d0_CDFc + R*d1i_CDF[i,]*d1i_CDF[j,]/d0_CDFc^2, na.rm = TRUE )
        D2_o[j,i] <- D2_o[i,j]
        D2_e[i,j] <- -n*sum( na.omit( K*( d2_CDF[i,,j] - d1_CDF[i,]*d1_CDF[j,]/d0_CDF -
                                      P*d2i_CDF[i,,j] - P*d1i_CDF[i,]*d1i_CDF[j,]/d0_CDFc) ) )
        D2_e[j,i] <- D2_e[i,j]
      }
    }
  }
  if (pdf.expression==TRUE)
  {
    d1_PDF <- d1_PDFc <- matrix(NA, nrow = d, ncol = m)
    d2_PDF<- d2_PDFc <- array(NA, dim = c(d, m, d) )
    integrand0 <- integrand1 <- integrand2 <- integrand3 <- function(x){}
    for(k in 1:d) assign(param[k], mle[k])
    body(integrand0) <- bquote(.(pdf))
    if (is.nan( integrand0(lb) ) == TRUE ) stop( "Try for another choice of lower bound lb" )
    I0 <- Vectorize( function(w) integrate(integrand0, lower = lb, upper = w)$value, "w" )
    C0 <- I0(T)
    d0_PDF <- diff( c( 0, C0 ) )
    d0_PDFc <- 1-C0
    for(i in 1:m) K[i] <- prod( (1-P[1:i-1]) )
    E_X <- d0_PDF*K
    E_R <- P*d0_PDFc*K
    first <- sapply(1:d, function(i) D(pdf, param[i]))
    for (i in 1:d)
    {
      body(integrand1) <- bquote(.(first[[i]]))
      #  if (is.nan( integrand1(lb) ) == TRUE ) stop( "Try for another choice of lower bound T0." )
      I1 <- Vectorize( function(w) quadinf(integrand1, lb, w)$Q, "w" )
      C1 <- I1(T)
      d1_PDF[i,] <- diff( c( 0, C1 ) )
      d1_PDFc[i,] <- -C1
    }
    for (i in 1:d)
    {
      for(j in i:d)
      {
        second <- D(first[[i]], param[j])
        body(integrand2) <- bquote(.(second))
        #   if (is.nan( integrand2(lb) ) == TRUE ) stop( "Try for another choice of lower bound T0." )
        I2 <- Vectorize( function(w) quadinf(integrand2, lb, w)$Q, "w" )
        C2 <- I2(T)
        d2_PDF[i,,j] <- diff( c( 0, C2 ) )
        d2_PDF[j,,i] <- d2_PDF[i,,j]
        d2_PDFc[i,,j] <- -C2
        d2_PDFc[j,,i] <- d2_PDFc[i,,j]
        D2_o[i,j] <- -sum( ( X*( d2_PDF[i,,j]/d0_PDF - d1_PDF[i,]*d1_PDF[j,]/d0_PDF^2 ) +
                             R*( d2_PDFc[i,,j]/d0_PDFc - d1_PDFc[i,]*d1_PDFc[j,]/d0_PDFc^2 )  )  )
        D2_o[j,i] <- D2_o[i,j]
        D2_e[i,j] <- -n*sum( ( E_X*( d2_PDF[i,,j]/d0_PDF - d1_PDF[i,]*d1_PDF[j,]/d0_PDF^2 ) +
                               E_R*( d2_PDFc[i,,j]/d0_PDFc - d1_PDFc[i,]*d1_PDFc[j,]/d0_PDFc^2 )  )  )
        D2_e[j,i] <- D2_e[i,j]
      }
    }
  }
  #  if ( sum( is.nan(D2) > 0 ) ) stop( "Try for another censoring scheme." )
  #  if(any(eigen(D2)$values <= 10e-15)) stop("The Hessian matrix is not invertible.")
  colnames(D2_e) <- colnames(D2_o) <- param
  rownames(D2_e) <- rownames(D2_o) <- param
  return(list( "FI.expected" = D2_e, "FI.observed" = D2_o ) )
}
###################################################################################################
mletype1 <- function(plan, param, start, cdf.expression = FALSE, pdf.expression = TRUE, cdf, pdf, method = "Nelder-Mead", lb = 0, ub = Inf, level = 0.05)
{
  T  <- plan$T
  R  <- plan$R
  R0 <- c(0, R)
  X  <- plan$X
  k  <- length(param);
  T0 <- c(lb, T)
   m <- length(R)
  R1 <- c( R[1], rep(0, m) )
  if( length(start) != k ) stop("The length of parameter vector and initial values must be the same.")
  f  <- function(x, par)
  {
    for(i in 1:k) assign(param[i], par[i])
    -sum( X*log( diff(eval(cdf)) )  ) - sum( R0*log(1 - eval(cdf)) ) - sum( R1*log(1 - eval(cdf)) )
  }
  out <- suppressWarnings( optim(start, fn = f, x = T0, method = method)$par )
  if (pdf.expression == TRUE)
  {
    D2 <- fitype1(plan, param, out, cdf.expression = FALSE, pdf.expression = TRUE, cdf, pdf, lb)$FI.expected
  }
  if (cdf.expression == TRUE)
  {
    D2 <- fitype1(plan, param, out, cdf.expression = TRUE, pdf.expression = FALSE, cdf, pdf, lb)$FI.expected
  }
  out1 <- cbind( out, sqrt( diag(solve(D2)) ), out + sqrt(diag(solve(D2)))*qnorm(level/2), out + sqrt(diag(solve(D2)))*qnorm(1 - level/2) )
  colnames(out1) <- c("estimate", "std. error", "lower bound", "upper bound")
  rownames(out1) <- param
  return(out1)
}
###################################################################################################
mletype2 <- function(plan, param, start, cdf, pdf, method = "Nelder-Mead", lb = 0, ub = Inf, N = 100, level = 0.05)
{
  X <- plan$X
  R <- plan$R
  k <- length(param)
  if( length(start) != k ) stop("The length of parameter vector and initial values must be the same.")
  f <- function(par,x)
  {
    for(k in 1:k) assign(param[k], par[k])
    -sum( log( eval(pdf)*(1-eval(cdf))^(eval(R)) ) )
  }
  out <- suppressWarnings( optim(start, fn = f, x = X, method = method)$par )
  D2 <- fitype2(plan, param, out, cdf, pdf, lb = lb, ub = ub, N = N)$FI.expected
  out1 <- cbind( out, sqrt( diag(solve(D2)) ),  out + sqrt(diag(solve(D2)))*qnorm(level/2), out + sqrt(diag(solve(D2)))*qnorm(1 - level/2) )
  colnames(out1) <- c("estimate", "std. error", "lower bound", "upper bound")
  rownames(out1) <- param
  return(out1)
}
###################################################################################################
coxbctype1 <- function(plan, param, mle, cdf.expression = FALSE, pdf.expression = TRUE, cdf, pdf, lb = 0)
{
  d <- length(mle)
  T <- plan$T
  P <- plan$P
  m <- length(T)
  n <- sum(plan$X + plan$R)
  if( length(param) != d ) stop("The length of parameter vector and ML estimators must be the same.")
  if( length(P) != m ) stop("The length of censoring times and percent of removals must be the same.")
  if( cdf.expression == TRUE & pdf.expression == TRUE ) stop("Both of cdf.expression and pdf.expression cannot be TRUE.")
  if( cdf.expression == FALSE & pdf.expression == FALSE ) stop("Both of cdf.expression and pdf.expression cannot be FALSE.")
  K <- rep(NA, m)
  D2 <- matrix (NA, nrow = d, ncol = d)
  D3 <- D3_1 <- array(NA, dim = c(d, d, d) )
  if(cdf.expression == TRUE)
  {
    d1_CDF <- d1_CDFc <- d1i_CDF <- d1i_CDFc <- matrix(NA, nrow = d, ncol = m)
    d2_CDF <- d2_CDFc <- d2i_CDF <- array(NA, dim = c(d, m, d) )
    CDF_0 <- CDF_1 <- CDF_2 <- CDF_3 <- function(x){}
    for(k in 1:d) assign(param[k], mle[k])
    body(CDF_0) <- bquote(.(cdf))
    if ( is.nan( CDF_0(lb) ) == TRUE ) stop( "Try for another choice of lower bound T0." )
    C0 <- CDF_0(T)
    d0_CDF <- diff( c( CDF_0(lb), C0) )
    d0_CDFc <- 1-C0
    for(i in 1:m) K[i] <- prod( (1-P[1:i-1]) )
    E_X <- d0_CDF*K
    E_R <- P*d0_CDFc*K
    first <- sapply(1:d, function(i) D(cdf, param[i]))
    for (i in 1:d)
    {
      body(CDF_1) <- bquote(.(first[[i]]))
      if (is.nan( CDF_1(lb) ) == TRUE ) stop( "Try for another choice of lower bound lb" )
      C1 <- CDF_1(T)
      d1_CDF[i,] <- diff( c( CDF_1(lb), C1 ) )
      d1i_CDF[i,] <-  C1
      d1_CDFc[i,] <- -C1
    }
    for(k in 1:d)
    {
      for (i in 1:d)
      {
        for(j in i:d)
        {
          second <- D(first[[i]], param[j])
          body(CDF_2) <- bquote(.(second))
          if (is.nan( CDF_2(lb) ) == TRUE ) stop( "Try for another choice of lower bound lb" )
          C2 <- CDF_2(T)
          d2_CDF[i,,j] <- diff( c( CDF_2(lb), C2 ) )
          d2_CDF[j,,i] <- d2_CDF[i,,j]
          d2_CDFc[i,,j] <- -C2
          d2_CDFc[j,,i] <- d2_CDFc[i,,j]
          d2i_CDF[i,,j] <- C2
          d2i_CDF[j,,i] <- d2i_CDF[i,,j]
          D2[i,j] <- -n*sum( na.omit( K*(
            d2_CDF[i,,j]-d1_CDF[i,]*d1_CDF[j,]/d0_CDF-P*d2i_CDF[i,,j]-
              P*d1i_CDF[i,]*d1i_CDF[j,]/d0_CDFc
          ) 	  )  	)
          D2[j,i] <- D2[i,j]
          third <- D(second, param[k])
          body(CDF_3) <- bquote(.(third))
          C3 <- CDF_3(T)
          d3_CDF <- diff( c( CDF_3(lb), C3 ) )
          d3_CDFc <- -C3
          D3[i,j,k] <- n*sum( na.omit( K*(
            d3_CDF-d2_CDF[i,,j]*d1_CDF[k,]/d0_CDF-(d2_CDF[i,,k]*
                                                     d1_CDF[j,]+d2_CDF[j,,k]*d1_CDF[i,])/d0_CDF+2*d1_CDF[k,]*
              d1_CDF[i,]*d1_CDF[j,]/d0_CDF^2-P*( C3+d1i_CDF[k,]*
                                                   d2i_CDF[i,,j]/d0_CDFc+(d2i_CDF[i,,k]*d1i_CDF[j,]+
                                                                            d2i_CDF[j,,k]*d1i_CDF[i,])/d0_CDFc+2*d1i_CDF[k,]*
                                                   d1i_CDF[i,]*d1i_CDF[j,]/d0_CDFc^2)
          )        )   )
          D3_1[i,j,k] <- n*sum( na.omit( K*(
            d3_CDF-(d2_CDF[j,,k]*d1_CDF[i,]+d2_CDF[i,,k]*d1_CDF[j,])/
              d0_CDF+d1_CDF[k,]*d1_CDF[i,]*d1_CDF[j,]/d0_CDF^2-P*(C3+
                                                                    (d2i_CDF[i,,k]*d1i_CDF[j,]+d2i_CDF[j,,k]*d1i_CDF[i,])/
                                                                    d0_CDFc+d1i_CDF[i,]*d1i_CDF[j,]*d1i_CDF[k,]/d0_CDFc^2)
          )       )   )
          D3[j,i,k] <- D3[i,j,k]
          D3_1[j,i,k] <- D3_1[i,j,k]
        }
      }
    }
  }
  if (pdf.expression == TRUE)
  {
    d1_PDF <- d1_PDFc <- matrix(NA, nrow = d, ncol = m)
    d2_PDF<- d2_PDFc <- array(NA, dim = c(d, m, d) )
    integrand0 <- integrand1 <- integrand2 <- integrand3 <- function(x){}
    for(k in 1:d) assign(param[k], mle[k])
    body(integrand0) <- bquote(.(pdf))
    #   if (is.nan( integrand0(lb) ) == TRUE ) stop( "Try for another choice of lower bound lb" )
    I0 <- Vectorize( function(w) integrate(integrand0, lower = lb, upper = w)$value, "w" )
    C0 <- I0(T)
    d0_PDF <- diff( c( 0, C0 ) )
    d0_PDFc <- 1-C0
    for(i in 1:m) K[i] <- prod( (1-P[1:i-1]) )
    E_X <- d0_PDF*K
    E_R <- P*d0_PDFc*K
    first <- sapply(1:d, function(i) D(pdf, param[i]))
    for (i in 1:d)
    {
      body(integrand1) <- bquote(.(first[[i]]))
      #   if (is.nan( integrand1(lb) ) == TRUE ) stop( "Try for another choice of lower bound lb" )
      I1 <- Vectorize( function(w) quadinf(integrand1, lb, w)$Q, "w" )
      C1 <- I1(T)
      d1_PDF[i,] <- diff( c( 0, C1 ) )
      d1_PDFc[i,] <- -C1
    }
    for(k in 1:d)
    {
      for (i in 1:d)
      {
        for(j in i:d)
        {
          second <- D(first[[i]], param[j])
          body(integrand2) <- bquote(.(second))
          #       if (is.nan( integrand2(lb) ) == TRUE ) stop( "Try for another choice of lower bound lb." )
          I2 <- Vectorize( function(w) quadinf(integrand2, lb, w)$Q, "w" )
          C2 <- I2(T)
          d2_PDF[i,,j] <-diff( c( 0, C2 ) )
          d2_PDF[j,,i] <- d2_PDF[i,,j]
          d2_PDFc[i,,j] <- -C2
          d2_PDFc[j,,i] <- d2_PDFc[i,,j]
          D2[i,j] <- -n*sum( ( E_X*( d2_PDF[i,,j]/d0_PDF - d1_PDF[i,]*d1_PDF[j,]/d0_PDF^2 ) +
                               E_R*( d2_PDFc[i,,j]/d0_PDFc - d1_PDFc[i,]*d1_PDFc[j,]/d0_PDFc^2 )  )  )
          D2[j,i] <- D2[i,j]
          third <- D(second, param[k])
          body(integrand3) <- bquote(.(third))
          I3 <- Vectorize( function(w) quadinf(integrand3, lb, w)$Q, "w" )
          C3 <- I3(T)
          d3_PDF <- diff( c( 0, C3 ) )
          d3_PDFc <- -C3
          D3[i,j,k] <- n*sum(
            E_X*( d3_PDF/d0_PDF-d2_PDF[i,,j]*d1_PDF[k,]/d0_PDF^2-
                    ((d2_PDF[j,,k]*d1_PDF[i,]+d2_PDF[i,,k]*d1_PDF[j,])*
                       d0_PDF-2*d1_PDF[k,]*d1_PDF[i,]*d1_PDF[j,])/d0_PDF^3)+
              E_R*( d3_PDFc/d0_PDFc-d2_PDFc[i,,j]*d1_PDFc[k,]/d0_PDFc^2-
                      ((d2_PDFc[j,,k]*d1_PDFc[i,]+d2_PDFc[i,,k]*d1_PDFc[j,])*
                         d0_PDFc-2*d1_PDFc[k,]*d1_PDFc[i,]*d1_PDFc[j,])/d0_PDFc^3)
          )
          D3_1[i,j,k] <- n*sum(
            K*( d3_PDF-(d2_PDF[j,,k]*d1_PDF[i,]+d2_PDF[i,,k]*d1_PDF[j,])/
                  d0_PDF+d1_PDF[k,]*d1_PDF[i,]*d1_PDF[j,]/d0_PDF^2+
                  P*( d3_PDFc-(d2_PDFc[i,,k]*d1_PDFc[j,]+d2_PDFc[j,,k]*
                                 d1_PDFc[i,])/d0_PDFc+d1_PDFc[k,]*d1_PDFc[i,]*
                        d1_PDFc[j,]/d0_PDFc^2 ) )
          )
          D3[j,i,k] <- D3[i,j,k]
          D3_1[j,i,k] <- D3_1[i,j,k]
        }
      }
    }
  }
  if ( sum( is.nan(D2) > 0 ) ) stop( "Try for another censoring scheme." )
#  if(any(eigen(D2)$values <= 10e-15)) stop("The Hessian matrix is not invertible.")
  bias <- as.vector( solve(D2)%*%matrix( cbind( D3_1 - D3/2 ), nrow = d , ncol = d*d )%*%c(solve(D2)) )
  mle.corrected <- mle - bias
  colnames(D2) <- param
  rownames(D2) <- param
  out1 <- as.matrix( rbind(mle, bias, mle.corrected), ncol = d, nrow = 3, byrow = TRUE)
  colnames(out1) <- param
  if (cdf.expression == TRUE)
  {
    measure_uncorrected <- goftype1(plan, param, mle, cdf.expression = TRUE, pdf.expression = FALSE, cdf, pdf, lb)
    measure_corrected <- goftype1(plan, param, mle.corrected, cdf.expression = TRUE, pdf.expression = FALSE, cdf, pdf, lb)
    FI.corrected <- fitype1(plan, param, mle.corrected, cdf.expression = TRUE, pdf.expression = FALSE, cdf, pdf, lb)$FI.expected
  }
  else
  {
    measure_uncorrected <- goftype1(plan, param, mle, cdf.expression = FALSE, pdf.expression = TRUE, cdf, pdf, lb)
    measure_corrected <- goftype1(plan, param, mle.corrected, cdf.expression = FALSE, pdf.expression = TRUE, cdf, pdf, lb)
    FI.corrected <- fitype1(plan, param, mle.corrected, cdf.expression = FALSE, pdf.expression = TRUE, cdf, pdf, lb)$FI.expected
  }
  out2 <- rbind( c(measure_uncorrected[1], measure_uncorrected[2], measure_uncorrected[3]),
                 c(  measure_corrected[1],   measure_corrected[2],   measure_corrected[3]) )
  colnames(out2) <- c("AD", "CVM", "KS")
  rownames(out2) <- c("uncorrected", "corrected")
  return( list( "cov" = solve(D2),  "cov.corrected" = solve(FI.corrected), "estimates" = out1, "measures" = out2 ) )
}
###################################################################################################
goftype1 <- function(plan, param, mle, cdf.expression = FALSE, pdf.expression = TRUE, cdf, pdf, lb = 0)
{
  d <- length(mle)
  T <- plan$T
  X <- plan$X
  R <- plan$R
  m <- length(T)
  n <- sum(X+R)
  PDF <- CDF <- function(x){}
  for(k in 1:d) assign(param[k], mle[k])
  body(CDF) <- bquote(.(cdf))
  body(PDF) <- bquote(.(pdf))
  stat_A <- stat_C <- rep(NA, m-1)
  Alpha <- U <- rep(NA, m)
  Alpha[1] <- X[1]/n
  for (i in 2:m)
  {
    Prod <- 1
    for (j in 2:i)
    {
      K <- n-sum(R[1:(j-1)])-sum(X[1:(j-1)])
      if (X[j]==0 && K==0)
      {
        Prod <- Prod
      }else{
        Prod <- Prod*( 1 - X[j]/K )
      }
    }
    Alpha[i] <- 1 - Prod*( 1-X[1]/n )
  }
  U <- CDF(T)
  for (i in 1:(m-1))
  {
    stat_A[i] <- Alpha[i]^2*log( U[i+1]*( 1-U[i] )/( U[i]*( 1-U[i+1] ) ) ) +
      2*Alpha[i]*log( ( 1-U[i+1] )/( 1-U[i] ) )
    stat_C[i] <- Alpha[i]*( U[i+1]-U[i] )*( Alpha[i]-U[i+1]-U[i] )
  }
  Anderson <- n*sum(stat_A) - n*log( 1 - U[m] ) - n*U[m]
  Cramer <- n*sum(stat_C) + n/3*U[m]^3
  KS <- max( abs( Alpha - U ) )
  return( list(AD = Anderson, CVM = Cramer, KS = KS ))
}
###################################################################################################
goftype2 <- function(plan, param, mle, cdf, pdf)
{
  d <- length(mle)
  R <- plan$R
  X <- plan$X
  m <- length(R)
  n <- sum(R) + length(R)
  if( length(param) != d ) stop("The length of parameter vector and ML estimators must be the same.")
  PDF <- CDF <- function(x){}
  for(k in 1:d) assign(param[k], mle[k])
  body(CDF) <- bquote(.(cdf))
  body(PDF) <- bquote(.(pdf))
  stat_A <- stat_C <- Alpha <- rep(NA, m-1)
  U <- rep(NA, m)
  U <- CDF(X)
  for (i in 1:(m-1))
  {
    Prod1 <- 1
    for ( j in (m-i+1):m ) Prod1 <- Prod1*( j+sum(R[(m-j+1):(m)]) )/( j+1+sum(R[(m-j+1):(m)]) )
    Alpha[i] <- 1 - Prod1
    stat_A[i] <- Alpha[i]^2*log( U[i+1]*( 1-U[i] )/( U[i]*( 1-U[i+1] ) ) ) + 2*Alpha[i]*log( ( 1-U[i+1] )/( 1-U[i] ) )
    stat_C[i] <- Alpha[i]*( U[i+1]-U[i] )*( Alpha[i]-U[i+1]-U[i] )
  }
  Anderson <- n*sum(stat_A) - n*log( 1 - U[m] ) - n*U[m]
  Cramer <- n*sum(stat_C) + n/3*U[m]^3
  #     Prod2 <- 1
  #  for (j in 1:(m-1)) Prod2 <- Prod2*( n-sum(R[1:j])-j )
  #  C_s  <- n*Prod2
  #  log.like <- 0
  #  for (i in 1:m) log.like <- log.like + log( PDF(X[i]) ) + R[i]*log( 1-U[i] )
  Prod2 <- 1
  for ( j in 1:m ) Prod2 <- Prod2*( j+sum(R[(m-j+1):(m)]) )/( j+1+sum(R[(m-j+1):(m)]) )
  KS <- max( abs( c(Alpha, 1 - Prod2) - U ) )
  return( list(AD = Anderson, CVM = Cramer, KS = KS ))
}
###################################################################################################
rtype2 <- function(n, R, param, mle, cdf, lb = 0, ub = Inf)
{
  m <- length(R)
  d <- length(mle)
  X <- rep(NA, m)
  ub <- ifelse(ub ==  Inf,  10e+16, ub)
  lb <- ifelse(lb == -Inf, -10e+16, lb)
  if( length(R) + sum(R) != n ) stop("The plan size in not compatible with the removed items.")
  if( length(param) != d ) stop("The length of parameter vector and ML estimators must be	the same.")
  CDF <- function(x){}
  V <- rep(NA, m)
  U <- V
  W <- runif(m)
  for (i in 1:m) V[i] <- W[i]^(1/( i+sum(R[(m-i+1):m]) ) )
  for (i in 1:m)
  {
    U[i] <- 1 - prod( V[(m-i+1):m] )
    body(CDF) <- bquote( .(cdf) - U[i] )
    for(k in 1:d) assign(param[k], mle[k]);
    X[i] <- uniroot(CDF, lower = lb, upper = ub )$root
  }
  return( data.frame(X = X, R = R))
}
###################################################################################################
fitype2 <- function(plan, param, mle, cdf, pdf, lb = 0, ub = Inf, N = 100)
{
  d <- length(mle)
  R <- plan$R
  X <- plan$X
  m <- length(R)
  n <- sum(R) + length(R)
  if( length(param) != d ) stop("The length of parameter vector and ML estimators must be the same.")
  D2_e <- D2_o <- matrix (NA, nrow = d, ncol = d)
  integrand0 <- integrand1 <- integrand2 <- function(x){}
  for(k in 1:d) assign(param[k], mle[k])
  first0 <- sapply(1:d, function(i) D( bquote(log(.(pdf))), param[i]))
  first1 <- sapply(1:d, function(i) D( bquote(log(.(pdf))), param[i]))
# first2 <- sapply(1:d, function(i) D( bquote(log(1 - .(cdf))), param[i]))
  for (r in 1:d)
  {
    for (k in r:d)
    {
      I1 <- 0
      for (s in 0:(m-1))
      {
        for (i in 0:s)
        {
          if ((i == s) & (i == 0)) { C_is <- 1
          } else if ((i == s) & (i != 0)) {
            Prod <- 1
            for (j in 1:i) Prod <- Prod*sum(R[(s-i+1):(s-i+j)]+1)
            C_is <- (-1)^i/Prod
          } else if ((i != s) & (i == 0)) {
            Prod <- 1
            for (j in 1:(s-i)) Prod <- Prod*sum(R[j:(s-i)]+1)
            C_is <- (-1)^i/Prod
          } else if ((i == s) & (i >= 1)) {
            Prod <- 1
            for (j in 1:i) Prod <- Prod*sum(R[(s-i+1):(s-i+j)]+1)
            C_is <- (-1)^i/Prod
          } else {
            Prod1 <- 1; Prod2 <- 1
            for (j in 1:i) Prod1 <- Prod1*sum(R[(s-i+1):(s-i+j)]+1)
            for (j in 1:(s-i)) Prod2 <- Prod2*sum(R[j:(s-i)]+1)
            C_is <- (-1)^i/(Prod1*Prod2)
          }
          R_is <- ifelse( i == s, n - s + i, n - s + i - sum(R[1:(s-i)]) )
          if (s == 0) { C_s <- n
          } else {
            Prod <- 1
            for (j in 1:s) Prod <- Prod*(n-sum(R[1:j])-j)
            C_s  <- n*Prod
          }
          body(integrand0) <- bquote(.(first0[[r]])*.(first0[[k]])*.(pdf)*(1 - .(cdf))^(R_is-1))
          dd <- simpson(integrand0, lb, ub, N)
          I1 <- I1 + C_s*C_is*dd
        }
      }
      D2_e[r , k] <- I1
      D2_e[k , r] <- I1
  #      body(integrand1) <- bquote(.(first1[[r]])*.(first1[[k]]))
  #      body(integrand2) <- bquote(.(first2[[r]])*.(first2[[k]]))
  #      D2_o[r , k] <- sum( integrand1(X) + R*integrand2(X) )
  #      D2_o[k , r] <- D2_o[r , k]
    }
  }
  #  if ( sum( is.nan(D2) > 0 ) ) stop( "Try for another censoring scheme." )
  #  if(any(eigen(D2)$values <= 10e-15)) stop("The Hessian matrix is not invertible.")
  colnames(D2_e) <- param
  # colnames(D2_o) <- param
  rownames(D2_e) <- param
  # colnames(D2_o) <- param
  # return(list( "FI.expected" = D2_e, "FI.observed" = D2_o ) )
   return(list( "FI.expected" = D2_e ) )
}
###################################################################################################
coxbctype2 <- function(plan, param, mle, cdf, pdf, lb = 0, ub = Inf, N = 100)
{
  d <- length(mle)
  R <- plan$R
  X <- plan$X
  m <- length(R)
  n <- sum(R) + length(R)
  if( length(param) != d ) stop("The length of parameter vector and ML estimators must be the same.")
  D2 <- matrix (NA, nrow = d, ncol = d)
  D3_1 <- D3 <- array(NA, dim = c(d, d, d) )
  integrand1 <- integrand2 <- integrand3 <- function(x){}
  for(k in 1:d) assign(param[k], mle[k])
  for (w in 1:d)
  {
    for (r in 1:d)
    {
      for (k in r:d)
      {
        I1 <- 0
        I2 <- 0
        I3 <- 0
        for (s in 0:(m-1))
        {
          for (i in 0:s)
          {
            if ((i == s) & (i == 0)) { C_is <- 1
            } else if ((i == s) & (i != 0)) {
              Prod <- 1
              for (j in 1:i) Prod <- Prod*sum(R[(s-i+1):(s-i+j)]+1)
              C_is <- (-1)^i/Prod
            } else if ((i != s) & (i == 0)) {
              Prod <- 1
              for (j in 1:(s-i)) Prod <- Prod*sum(R[j:(s-i)]+1)
              C_is <- (-1)^i/Prod
            } else if ((i == s) & (i >= 1)) {
              Prod <- 1
              for (j in 1:i) Prod <- Prod*sum(R[(s-i+1):(s-i+j)]+1)
              C_is <- (-1)^i/Prod
            } else {
              Prod1 <- 1; Prod2 <- 1
              for (j in 1:i) Prod1 <- Prod1*sum(R[(s-i+1):(s-i+j)]+1)
              for (j in 1:(s-i)) Prod2 <- Prod2*sum(R[j:(s-i)]+1)
              C_is <- (-1)^i/(Prod1*Prod2)
            }
            R_is <- ifelse( i == s, n -s + i, n - s + i - sum(R[1:(s-i)]) )
            if (s == 0)
            {
              C_s <- n
            } else {
              Prod <- 1
              for (j in 1:s) Prod <- Prod*(n-sum(R[1:j])-j)
              C_s  <- n*Prod
            }
            first0 <- sapply(1:d, function(h) D( bquote(log(.(pdf)/(1-.(cdf)))), param[h]))
            first1 <- sapply(1:d, function(h) D( bquote(log(1-.(cdf))), param[h]))
            first2 <- sapply(1:d, function(h) D( bquote(log(.(pdf))), param[h]))
            second1 <- D(first1[[r]], param[k])
            third1 <- D(second1, param[w])
            second2 <- D(first2[[r]], param[k])
            third2 <- D(second2, param[w])
            body(integrand1) <- bquote(.(third1)*.(pdf)*(1- .(cdf))^(R_is-1))
            body(integrand2) <- bquote(.(third2)*.(pdf)*(1- .(cdf))^(R_is-1))
            third_1 <- D( bquote(.(first0[[r]])*.(first0[[k]])*.(pdf)*(1- .(cdf))^(R_is-1)), param[w])
            body(integrand3) <- bquote(.(third_1))
            I1 <- I1 + R[s+1]*C_s*C_is*simpson(integrand1, lb, ub, N)
            I2 <- I2 + C_s*C_is*simpson(integrand2, lb, ub, N)
            I3 <- I3 + C_s*C_is*simpson(integrand3, lb, ub, N)
          }
        }
        D3[r , k, w] <- I1  + I2
        D3[k , r, w] <- I1  + I2
        D3_1[r , k, w] <- I3
        D3_1[k , r, w] <- I3
      }
    }
  }
  D2 <- fitype2(plan, param, mle, cdf, pdf, lb = lb, ub = ub, N = N)$FI.expected
  # if ( sum( is.nan(D2) > 0 ) ) stop( "Try for another censoring scheme." )
  # if(any(eigen(D2)$values <= 10e-15)) stop("The Hessian matrix is not invertible.")
  bias <- as.vector( solve(D2)%*%matrix( cbind( -D3_1 - D3/2 ), nrow = d , ncol = d*d )%*%c(solve(D2)) )
  mle.corrected <- mle - bias
  colnames(D2) <- param
  rownames(D2) <- param
  out1 <- as.matrix( rbind(mle, bias, mle.corrected), ncol = d, nrow = 3, byrow = TRUE)
  colnames(out1) <- param
  measure_uncorrected <- goftype2(plan = plan, param = param, mle = mle, cdf = cdf, pdf = pdf)
  measure_corrected   <- goftype2(plan = plan, param = param, mle = mle.corrected, cdf = cdf, pdf = pdf)
  FI.corrected <- fitype2(plan = plan, param = param, mle = mle.corrected, cdf = cdf, pdf = pdf, lb = lb, ub = ub, N = N)$FI.expected
  out2 <- rbind( c(measure_uncorrected[1], measure_uncorrected[2], measure_uncorrected[3]),
                 c(measure_corrected[1]  , measure_corrected[2]  , measure_corrected[3]) )
  colnames(out2) <- c("AD", "CVM", "KS")
  rownames(out2) <- c("uncorrected", "corrected")
  return( list( "cov" = solve(D2),  "cov.corrected" = solve(FI.corrected), "estimates" = out1, "measures" = out2 ) )
}
###################################################################################################
simpson <- function(fun, lb, ub, N = 100) {
  if (lb == -Inf & ub == Inf) {
    f <- function(t) ( fun( (1-t)/t ) + fun( (t-1)/t ) )/t^2
    s <- simpson(f, 0, 1, N)
  } else if (lb == -Inf & ub != Inf) {
    f <- function(t) fun( ub-(1-t)/t )/t^2
    s <- simpson(f, 0, 1, N)
  } else if (lb != -Inf & ub == Inf) {
    f <- function(t) fun( lb + (1-t)/t )/t^2
    s <- simpson(f, 0, 1, N)
  } else {
    h <- (ub-lb)/N
    x <- seq(lb, ub, by = h)
    y <- fun(x)
    y[is.nan(y)] <- 0
    s <- y[1] + y[N+1] + 2*sum(y[seq(2, N, by = 2)]) + 4 *sum( y[seq(3, N-1, by = 2)] )
    s <- s*h/3
  }
  return(s)
}
###################################################################################################
bootbctype2 <- function(plan, param, mle, cdf, pdf, lb = 0, ub = Inf, nboot = 200, coverage = 0.95)
{
  d <- length(mle)
  m <- length(plan$R)
  n <- sum(plan$R) + m
  ML <- matrix(NA, nboot, d)
  if( length(param) != d ) stop("The length of parameter vector and ML estimators must be the same.")
  for (i in 1:nboot)
  {
    plan1 <- rtype2(n = n, R = plan$R, param = param, mle = mle, cdf = cdf, lb = lb, ub = ub)
    R <- plan1$R
    X <- plan1$X
    g <- function(par, x)
    {
      for (j in 1:d) assign(param[j], par[j])
      -sum( log( eval(pdf)*(1 - eval(cdf))^R ) )
    }
    ML[i,]<- suppressWarnings( optim(mle, fn = g, x = X, method = "Nelder-Mead")$par )
  }
  bias.mean <- as.matrix( apply(ML, 2, mean) - mle )
  bias.median <- as.matrix( apply(ML, 2, median) - mle )
  LPCI <- sapply(1:d, function(i)quantile( ML[,i], (1 - coverage)/2 ) )
  UPCI <- sapply(1:d, function(i)quantile(ML[,i], 1 - (1 - coverage)/2 ) )
  cov <- cov(ML)
  colnames(cov) <- param
  rownames(cov) <- param
  out1 <- cbind( bias.mean, bias.median, LPCI, UPCI )
  colnames(out1) <- c("bias from mean", "bias from median", "LPCI", "UPCI")
  rownames(out1) <- param
  return(list ( summary = out1, cov = cov) )
}
###################################################################################################
bootbctype1 <- function(plan, param, mle, cdf, lb = 0, ub = Inf, nboot = 200, coverage = 0.95)
{
  d <- length(mle)
#  R <- plan$R
#  X <- plan$X
  n <- sum(plan$X) + length(plan$R)
  ML <- matrix(NA, nboot, d)
  if( length(param) != d ) stop("The length of parameter vector and ML estimators must be the same.")
  for (i in 1:nboot)
  {
    f <- function(x, par)
    {
      for(j in 1:d) assign(param[j], par[j])
      -sum( X*log( diff(eval(cdf)) ) ) - sum( R0*log(1 - eval(cdf)) ) - sum( R1*log(1 - eval(cdf)) )
    }
    plan1 <- rtype1(n = n, P = plan$P, T = plan$T, param = param, mle = mle, TRUE, FALSE, cdf = cdf, pdf = pdf, lb = lb)
    X <- plan1$X
    T <- plan$T
    R <- plan1$R
    R0 <- c(0, R)
    T0 <- c(lb, T)
    m <- length(R)
    R1 <- c( R[1], rep(0, m) )
    ML[i,] <- suppressWarnings( optim(mle, fn = f, x = T0, method = "Nelder-Mead")$par )
  }
  bias.mean <- as.matrix( apply(ML, 2, mean) - mle )
  bias.median <- as.matrix( apply(ML, 2, median) - mle )
  LPCI <- sapply(1:d, function(i)quantile( ML[,i], (1 - coverage)/2 ) )
  UPCI <- sapply(1:d, function(i)quantile(ML[,i], 1 - (1 - coverage)/2 ) )
  cov <- cov(ML)
  colnames(cov) <- param
  rownames(cov) <- param
  out1 <- cbind( bias.mean, bias.median, LPCI, UPCI )
  colnames(out1) <- c("bias from mean", "bias from median", "LPCI", "UPCI")
  rownames(out1) <- param
  return(list ( summary = out1, cov = cov) )
}

