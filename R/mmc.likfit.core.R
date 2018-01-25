#'
#' Internal spatial MLE for convoluted permuted data.
#'
#' @export
#'
.mmc.likfit.core <-
  function( data, coords
            , ml
            , nConv, kernelConv
            , h.min, h.max, nn.min
            , nSimMMC, log ){

    ##----------------------------------------------------------
    ##
    ## 1. nSimMMC permutations and their restored spatial pattern
    ##
    ml.tausd <- sqrt( ml$tausq )
    new.data <- matrix( NA, nrow = length( data ), ncol = nSimMMC )
    ml.conv.phi <- numeric( length = nSimMMC )
    ml.conv.sigmasq <- numeric( length = nSimMMC )
    ml.conv.tausq <- numeric( length = nSimMMC )
    ##------------------

    ##
    if( log ) differences.phis <- numeric( nSimMMC )
    for( s in 1:nSimMMC ){ ## s <- 1

      ## I - Permutation
      data.perm <- data[ sample( 1:length( data ), replace = FALSE ) ]

      ## II - Convolution
      diff.phis <- 0
      k <- nConv
      ml.conv <- list()
      h <- seq( h.min, h.max, l = nConv )

      for( conv.i in 1:nConv ){

        ## Evaluate fitted values for locfit regression object
        fit.values <- lapply( .conv.permuted( data = data.perm
                                              , coords = coords
                                              , nConv = 1, kernelConv = kernelConv
                                              , h.min = h.min, h.max = h.max )
                              , fitted )[[1]]

        ## III - Permuted and convoluted data model estimation
        ## Convolution parameter estimation
        ini <- variofit( variog( data = fit.values, coords = coords, messages = F ), messages = F )$cov.pars
        if( ini[1] <= 0.1 | ini[2] <= 0.03 | ini[2] > 1 ) ini <- c(1,0.15)
        ml.conv[[ conv.i ]] <- likfit( data = fit.values
                                       , coords = coords, fix.nugget = T, messages = F
                                       , ini.cov.pars = ini )

        ## Checks if the difference between target phi and this convolution phi is getting bigger
        ## If so, choose this covolution and stop for
        if( conv.i >= 2 ){

          ## If the phi difference is getting bigger, it should stop and use the previous convolution
          if( diff.phis <= abs( ml.conv[[ conv.i ]]$phi - ml$phi ) ){
            k <- ( conv.i - 1 )
            fit.values <- fit.values.prev
            break
          }else{
            ## If the phi difference is getting smaller, if should be updated and to continue
            diff.phis <- abs( ml.conv[[ conv.i ]]$phi - ml$phi )
            fit.values.prev <- fit.values
          }

        }else{
          ## First step
          diff.phis <- abs( ml.conv[[ conv.i ]]$phi - ml$phi )
          fit.values.prev <- fit.values
        }

      }

      ## IV - Adjust variance
      ##
      ## Adjust sill
      fit.varcor <- fit.values * sqrt( ml$sigmasq / ml.conv[[ k ]]$sigmasq )
      ## Add nugget
      fit.varcor <- fit.varcor + rnorm( length( data ), mean = 0, sd = ml.tausd )

      ## IV - Record permutation
      new.data[, s] <- fit.varcor

      ## V - Record model parameters for choosen permutation
      differences.phis[ s ] <- diff.phis

    }
    ##------------------


    ##----------------------------------------------------------
    ##
    ## 2. Explicit return
    ##
    res <- list( new.data = new.data )
    if( log ){
      res <- list( new.data = new.data
                   , differences.phis = differences.phis ) }
    return( res )

  }
