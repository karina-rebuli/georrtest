#'
#' Internal MSE between variogs for convoluted permuted data.
#'
#' @export
#'
.mmc.variog.core <-
  function( data, coords
            , vg
            , nConv, kernelConv, h.min, h.max
            , nSimMMC, log ){

    ##----------------------------------------------------------
    ##
    ## 1. Data
    ##
    nData <- length( data )
    ##-------------


    ##----------------------------------------------------------
    ##
    ## 2. nSimMMC permutations and their restored spatial pattern
    ##
    new.data <- matrix( NA, nrow = nData, ncol = nSimMMC )
    if( log ) residuals.perm <- numeric( nSimMMC )
    for( s in 1:nSimMMC ){
      ##
      ## I - Permute
      ##
      data.perm <- data[ sample( 1:nData, replace = FALSE ) ]

      ##
      ## II - Convolution
      ##
      ## Evaluate fitted values for locfit regression object
      fit.values <- lapply( .conv.permuted( data = data.perm
                                            , coords = coords
                                            , nConv = nConv, kernelConv = kernelConv
                                            , h.min = h.min, h.max = h.max )
                            , fitted )

      ##
      ## III - Scale
      ##
      ## Variogram of permuted and smoothed data
      vg.sm <- lapply( fit.values
                       ## For each trial ( nn or h )
                       , function( x ){
                         variog( data = x, coords = coords, uvec = vg$u, max.dist = vg$max.dist, messages = FALSE )
                       } )

      ## Linear regression between smoothed-permuted and original data variograms
      linear.fit <- lapply( vg.sm
                            ## For each trial ( nn or h )
                            , function( x ){
                              lm( vg$v ~ ( 1 + x$v ) )
                            } )

      ## Scaled smoothed permuted data
      ## Transformation linear model -> y = alpha (=coef[1]) + beta (=coef[2] * x )
      data.perm.sm.scaled <- lapply( 1:length( linear.fit )
                                     ## For each trial (nn or h )
                                     , function( n ){
                                       fit.values[[ n ]] * sqrt( abs( linear.fit[[ n ]]$coef[2] ) ) + rnorm( nData ) * sqrt( abs( linear.fit[[ n ]]$coef[1] ) )
                                     } )

      ## Residuals
      residuals <- sapply( 1:length( data.perm.sm.scaled )
                           ## For each trial ( nn or h )
                           , function( n ){
                             sum( ( vg$v - variog( data = data.perm.sm.scaled[[ n ]], uvec = vg$uvec, max.dist = vg$max.dist
                                                   , coords = coords, messages = FALSE )$v )^2 )
                           } )

      ## Choosing of delta (nn) or distance (h) for local regression based on least squares of residuals
      k <- which.min( residuals )
      if( log ) residuals.perm[ s ] <- residuals[ k ]

      ##
      ## IV - Record permutation
      ##
      new.data[, s] <- data.perm.sm.scaled[[ k ]]

    }
    ##-------------


    ##----------------------------------------------------------
    ##
    ## 3. Explicit return
    ##
    res <- list( new.data = new.data )
    if( log ){
      res <- list( new.data = new.data
                   , residuals = residuals.perm ) }

    return( res )

  }
