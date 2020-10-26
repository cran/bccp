bctypei <- function(plan, param, mle, cdf.expression=TRUE, pdf.expression=FALSE, cdf, pdf, lower = 0)
{
  d <- length(mle)
  T <- plan$T
  P <- plan$P
  m <- length(T)
  n <- sum(plan$X + plan$R)
  if( length(param) != d ) stop("The length of parameter vector and ML estimators must be the same.")
  if( length(P) != m ) stop("The length of censoring times and percent of removals must be the same.")
  K <- rep(NA, m)
  D2 <- matrix (NA, nrow = d, ncol = d)
  D3 <- D3_1 <- array(NA, dim = c(d, d, d) )
  if(cdf.expression==TRUE)
  {
    d1_CDF <- d1_CDFc <- d1i_CDF <- d1i_CDFc <- matrix(NA, nrow = d, ncol = m)
    d2_CDF <- d2_CDFc <- d2i_CDF <- array(NA, dim = c(d, m, d) )
    CDF_0 <- CDF_1 <- CDF_2 <- CDF_3 <- function(x){}
    for(k in 1:d) assign(param[k], mle[k])
    body(CDF_0) <- bquote(.(cdf))
    if ( is.nan( CDF_0(lower) ) == TRUE ) stop( "Try for another choice of lower bound T0." )
    C0 <- CDF_0(T)
    d0_CDF <- diff( c( CDF_0(lower), C0) )
    d0_CDFc <- 1-C0
    for(i in 1:m) K[i] <- prod( (1-P[1:i-1]) )
    E_X <- d0_CDF*K
    E_R <- P*d0_CDFc*K
    first <- sapply(1:d, function(i) D(cdf, param[i]))
    for (i in 1:d)
    {
      body(CDF_1) <- bquote(.(first[[i]]))
      if (is.nan( CDF_1(lower) ) == TRUE ) stop( "Try for another choice of lower bound T0." )
      C1 <- CDF_1(T)
      d1_CDF[i,] <- diff( c( CDF_1(lower), C1 ) )
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
          if (is.nan( CDF_2(lower) ) == TRUE ) stop( "Try for another choice of lower bound T0." )
          C2 <- CDF_2(T)
          d2_CDF[i,,j] <- diff( c( CDF_2(lower), C2 ) )
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
          d3_CDF <- diff( c( CDF_3(lower), C3 ) )
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
  if (pdf.expression==TRUE)
  {
    d1_PDF <- d1_PDFc <- matrix(NA, nrow = d, ncol = m)
    d2_PDF<- d2_PDFc <- array(NA, dim = c(d, m, d) )
    integrand0 <- integrand1 <- integrand2 <- integrand3 <- function(x){}
    for(k in 1:d) assign(param[k], mle[k])
    body(integrand0) <- bquote(.(pdf))
    if (is.nan( integrand0(lower) ) == TRUE ) stop( "Try for another choice of lower bound T0." )
    I0 <- Vectorize( function(w) integrate(integrand0, lower = lower, upper = w)$value, "w" )
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
      if (is.nan( integrand1(lower) ) == TRUE ) stop( "Try for another choice of lower bound T0." )
      I1 <- Vectorize( function(w) quadinf(integrand1, lower, w)$Q, "w" )
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
          if (is.nan( integrand2(lower) ) == TRUE ) stop( "Try for another choice of lower bound T0." )
          I2 <- Vectorize( function(w) quadinf(integrand2, lower, w)$Q, "w" )
          C2 <- I2(T)
          d2_PDF[i,,j] <-diff( c( 0, C2 ) )
          d2_PDF[j,,i] <- d2_PDF[i,,j]
          d2_PDFc[i,,j] <- -C2
          d2_PDFc[j,,i] <- d2_PDFc[i,,j]
          D2[i,j] <- -n*sum( ( E_X*( d2_PDF[i,,j]/d0_PDF-d1_PDF[i,]*d1_PDF[j,]/d0_PDF^2 ) +
                                 E_R*( d2_PDFc[i,,j]/d0_PDFc-d1_PDFc[i,]*d1_PDFc[j,]/d0_PDFc^2 )  )  )
          D2[j,i] <- D2[i,j]
          third <- D(second, param[k])
          body(integrand3) <- bquote(.(third))
          I3 <- Vectorize( function(w) quadinf(integrand3, lower, w)$Q, "w" )
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
  if ( sum( is.nan(D2) >0 ) ) stop( "Try for another censoring scheme." )
  if(any(eigen(D2)$values <= 10e-15)) stop("The Hessian matrix is not invertible.")
  bias <- as.vector( solve(D2)%*%matrix( cbind( D3_1-D3/2 ), nrow = d , ncol = d*d )%*%c(solve(D2)) )
  bias.corrected <- mle-bias
  colnames(D2) <- param
  rownames(D2) <- param
  out1 <- as.matrix( rbind(mle, bias, bias.corrected), ncol = d, nrow = 3, byrow = TRUE)
  colnames(out1) <- param
  if (cdf.expression==TRUE)
  {
    measure_uncorrected   <- goftypei(plan, param, mle, TRUE, FALSE, cdf, pdf, lower)
    measure_corrected <- goftypei(plan, param, bias.corrected, TRUE, FALSE, cdf, pdf, lower)
    FI.corrected <- fitypei(plan, param, bias.corrected, TRUE, FALSE, cdf, pdf, lower)
  }
  else
  {
    measure_uncorrected   <- goftypei(plan, param, mle, FALSE, TRUE, cdf, pdf, lower)
    measure_corrected <- goftypei(plan, param, bias.corrected, FALSE, TRUE, cdf, pdf, lower)
    FI.corrected <- fitypei(plan, param, bias.corrected, FALSE, TRUE, cdf, pdf, lower)
  }
  out2 <- rbind( c(measure_uncorrected[1], measure_uncorrected[2]),
                 c(measure_corrected[1], measure_corrected[2]) )
  colnames(out2) <- c("AD", "CVM")
  rownames(out2) <- c("uncorrected", "corrected")
  return( list( "Cov" = solve(D2),  "cov.corrected" = solve(FI.corrected), "estimates" = out1, "measures" = out2 ) )
}

################################################################
fitypei <- function(plan, param, mle, cdf.expression=TRUE, pdf.expression=FALSE, cdf, pdf, lower = 0)
{
  d <- length(mle)
  T <- plan$T
  P <- plan$P
  m <- length(T)
  n <- sum(plan$X + plan$R)
  K <- rep(NA, m)
  D2 <- matrix (NA, nrow = d, ncol = d)
  if(cdf.expression==TRUE)
  {
    d1_CDF <- d1_CDFc <- d1i_CDF <- d1i_CDFc <- matrix(NA, nrow = d, ncol = m)
    d2_CDF <- d2_CDFc <- d2i_CDF <- array(NA, dim = c(d, m, d) )
    CDF_0 <- CDF_1 <- CDF_2 <- function(x){}
    for(k in 1:d) assign(param[k], mle[k])
    body(CDF_0) <- bquote(.(cdf))
    if ( is.nan( CDF_0(lower) ) == TRUE ) stop( "Try for another choice of lower bound T0." )
    C0 <- CDF_0(T)
    d0_CDF <- diff( c( CDF_0(lower), C0) )
    d0_CDFc <- 1-C0
    for(i in 1:m) K[i] <- prod( (1-P[1:i-1]) )
    E_X <- d0_CDF*K
    E_R <- P*d0_CDFc*K
    first <- sapply(1:d, function(i) D(cdf, param[i]))
    for (i in 1:d)
    {
      body(CDF_1) <- bquote(.(first[[i]]))
      if (is.nan( CDF_1(lower) ) == TRUE ) stop( "Try for another choice of lower bound T0." )
      C1 <- CDF_1(T)
      d1_CDF[i,] <- diff( c( CDF_1(lower), C1 ) )
      d1i_CDF[i,] <-  C1
      d1_CDFc[i,] <- -C1
    }
    for (i in 1:d)
    {
      for(j in i:d)
      {
        second <- D(first[[i]], param[j])
        body(CDF_2) <- bquote(.(second))
        if (is.nan( CDF_2(lower) ) == TRUE ) stop( "Try for another choice of lower bound T0." )
        C2 <- CDF_2(T)
        d2_CDF[i,,j] <- diff( c( CDF_2(lower), C2 ) )
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
    if (is.nan( integrand0(lower) ) == TRUE ) stop( "Try for another choice of lower bound T0." )
    I0 <- Vectorize( function(w) integrate(integrand0, lower = lower, upper = w)$value, "w" )
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
      if (is.nan( integrand1(lower) ) == TRUE ) stop( "Try for another choice of lower bound T0." )
      I1 <- Vectorize( function(w) quadinf(integrand1, lower, w)$Q, "w" )
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
        if (is.nan( integrand2(lower) ) == TRUE ) stop( "Try for another choice of lower bound T0." )
        I2 <- Vectorize( function(w) quadinf(integrand2, lower, w)$Q, "w" )
        C2 <- I2(T)
        d2_PDF[i,,j] <-diff( c( 0, C2 ) )
        d2_PDF[j,,i] <- d2_PDF[i,,j]
        d2_PDFc[i,,j] <- -C2
        d2_PDFc[j,,i] <- d2_PDFc[i,,j]
        D2[i,j] <- -n*sum( ( E_X*( d2_PDF[i,,j]/d0_PDF-d1_PDF[i,]*d1_PDF[j,]/d0_PDF^2 ) +
                               E_R*( d2_PDFc[i,,j]/d0_PDFc-d1_PDFc[i,]*d1_PDFc[j,]/d0_PDFc^2 )  )  )
        D2[j,i] <- D2[i,j]
      }
    }
  }
  if ( sum( is.nan(D2) >0 ) ) stop( "Try for another censoring scheme." )
  if(any(eigen(D2)$values <= 10e-15)) stop("The Hessian matrix is not invertible.")
  return( D2 )
}

rtypei <- function(n, plan, param, mle, cdf.expression = TRUE, pdf.expression = FALSE, cdf, pdf, lower = 0)
{
  d <- length(mle)
  T <- plan$T
  P <- plan$P
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
    if ( is.nan( CDF(lower) ) == TRUE ) stop( "Try for another choice of lower bound T0." )
    C[1] <- CDF(T[1])-CDF(lower)
    X[1] <- rbinom(1, n, C[1])
    R[1] <- floor(P[1]*(n-X[1]))
    xsum[1] <- X[1]
    rsum[1] <- R[1]
    for (i in 2:m)
    {
      D1 <- CDF(T[i])-CDF(T[i-1])
      D2 <- CDF(T[i-1])-CDF(lower)
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
    if (is.nan( integrand(lower) ) == TRUE ) stop( "Try for another choice of lower bound T0." )
    C[1] <- quadinf(integrand, lower, T[1])$Q
    X[1] <- rbinom(1, n, C[1])
    R[1] <- floor(P[1]*(n-X[1]))
    xsum[1] <- X[1]
    rsum[1] <- R[1]
    for (i in 2:m)
    {
      D1 <- quadinf(integrand, T[i-1], T[i])$Q
      D2 <- quadinf(integrand, 0, T[i-1])$Q
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
################################################################
goftypei <- function(plan, param, mle, cdf.expression = TRUE, pdf.expression = FALSE, cdf, pdf, lower = 0)
{
  d <- length(mle)
  m <- length(plan$T)
  gamma <- C1 <- C2 <- rep(NA, m)
  A <- rep(NA, m+1)
  X <- plan$X
  R <- plan$R
  n <- sum(X+R)
  gamma[1] <- X[1]/n
  for(i in 2:m ) gamma[i] <- ( sum( X[1:i] ) + sum( R[1:i-1] ) )/n
  if (cdf.expression==TRUE)
  {
    CDF <- function(x){}
    for(k in 1:d) assign(param[k], mle[k])
    body(CDF) <- bquote(.(cdf))
    if ( is.nan( CDF(lower) ) == TRUE ) stop( "Try for another choice of lower bound T0." )
    A <- c( CDF(lower), CDF(T) )
    for (i in 1:m)
    {
      C1[i] <- gamma[i]^2*log( A[i+1]*(1-A[i])/( A[i]*(1-A[i+1])) ) +
        2*gamma[i]*log( (1-A[i+1])/(1-A[i]) )
    }
    AD <- n*( sum(C1) - (A[m+1]-A[1]) - log( (1-A[m+1])/(1-A[1]) ) + 1-A[m+1]-log(A[m+1]) )
    for (i in 1:m) C2[i] <- gamma[i]^2*(A[i+1]-A[i]) - gamma[i]*( A[i+1]^2-A[i]^2 )
    CVM <- n*( sum(C2) + (A[m+1]^3-A[1]^3)/3 + (1-A[m+1])^3/3 )
    #log_like <- sum( X*log(diff(A)) + R*(1-A[2:(m+1)]) )
  }
  if (pdf.expression==TRUE)
  {
    integrand <- function(x){}
    for(k in 1:d) assign(param[k], mle[k])
    body(integrand) <- bquote(.(pdf))
    if (is.nan( integrand(lower) ) == TRUE ) stop( "Try for another choice of lower bound T0." )
    I2 <- Vectorize( function(w) quadinf(integrand, 0, w)$Q, "w" )
    A <- I2( c(lower, T) )
    for (i in 1:m)
    {
      C1[i] <- gamma[i]^2*log( A[i+1]*(1-A[i])/( A[i]*(1-A[i+1])) ) +
        2*gamma[i]*log( (1-A[i+1])/(1-A[i]) )
    }
    AD <- n*( sum(C1) - (A[m+1]-A[1]) - log( (1-A[m+1])/(1-A[1]) ) + 1-A[m+1]-log(A[m+1]) )
    for (i in 1:m) C2[i] <- gamma[i]^2*(A[i+1]-A[i]) - gamma[i]*( A[i+1]^2-A[i]^2 )
    CVM <- n*( sum(C2) + (A[m+1]^3-A[1]^3)/3 + (1-A[m+1])^3/3 )
    #log_like <- sum( X*log(diff(A)) + R*(1-A[2:(m+1)]) )
  }
  return( c(AD, CVM) )
}
