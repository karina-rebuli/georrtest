#'
#' Non-parametrical Monte Carlo method for testing Pearson correlation coeficient (r) significance under null
#' hypothesis (\eqn{H_0: r = 0}) between spatial autodependent data, proposed on Viladomat, et al. 2014.
#' The algorithm follow the steps:
#'
#' 1. Permute one (or both) autocorrelated data.
#' 2. Repeat 3 to 5 N iteration times
#' 3. Convolute permuted data to reinsert autodependence on it.
#' 4. Choose the convoluted data which variogram is closer (min least square errors) to the origianl data variogram.
#' 5. Evaluate correlation between data.
#'
#' Viladomat, Juia and Mazumder, Rahul and McInturff, Alex and McCauley, Douglas J and Hastie, Trevor. 2014.
#' Assessing the significance of global and local correlations under spatial autocorrelation: a nonparametric approach.
#' Biometrics, 2, 409--18, vol 70.
#'
#' @param geodata An object of geodata class. Mandatory if X or Y are not setted.
#' @param X X data vector. Mandatory.
#' @param Y Y data vector. Mandatory.
#' @param coords Data grid (locations). Mandatory.
#' @param shuffle If it should be permuted one or both datasets. Allowed values: "single", "s" or 1 and "double", "d" or 2. Default( 1 ).
#' @param nConv Number of different window size convolutions must be done. Allowed values are integers between 3 and 20. Default( 10 ).
#' @param vg.X A variog\{geoR\} X data return object. Default( NULL ).
#' @param vg.Y A variog\{geoR\} Y data return object. Default( NULL ).
#' @param vg.uvec.X A vector with distances should be used to build empirical variograms of X variable. Default( NULL ).
#' @param vg.uvec.Y A vector with distances should be used to build empirical variograms of Y variable. Default( NULL ).
#' @param vg.maxdist.X The maximum distance that should be used to build empirical variograms of X variable. Default( NULL ).
#' @param vg.maxdist.Y The maximum distance that should be used to build empirical variograms of Y variable. Default( NULL ).
#' @param h.min Smallest locfit delta (h) parameter. If it's not setted and if nn.max ins't setted
#' , it will be the value corresponding to 10th percentil of distances vector.
#' @param h.max Biggest locfit delta (h) parameter. If it's not setted and if nn.max ins't setted
#' , it will be the value corresponding to 50th percentil of distances vector.
#' @param kernelConv Weight function for convolutions; check man locfit package to see allowed values.
#' @param nSimMMC Number of MMC simulations to be done. Allowed integers from 1e2 to 1e4. Default( 1e3 ).
#' @param logVars logical, if the returned object have a vector with estimated MMC r. Default( TRUE ).
#'
#' @return r.hat Sample Pearson correlation coefficient.
#' @return p.value The empirical p-value obtained.
#' @return r.perm A vector with estimated Pearson correlation coefficient in each MMC permutation.
#' @return residuals A vector with least squared residuals from original data empirical variogram and empirical
#' variogram of choosen permuted and transformed data used as MMC permutation.
#'
#' @author Karina Brotto rebuli and Paulo Justiniano Ribeiro Jr.
#'
#' @note Warning: The function can be time consuming!
#' @note Warning: If variog._ is not provided, this function will try to estimate spatial
#' model parameters with default settings with variofit\{geoR\} function.
#' It's not a recommended use! Since the autodependence structure affects the correlation coefficient distribution
#' , a bad parameters estimation can lead to misleading results.
#'
#' @keywords mmc, correlation coeficient significance, geostatistical data, empirical variograms
#'
#' @export
#'
mmc.variog.rtest <-
function( geodata, X, Y, coords
                           , shuffle = 1, nConv = 10
                           , vg.X = NULL, vg.Y = NULL
                           , vg.uvec.X = NULL, vg.uvec.Y = NULL
                           , vg.maxdist.X = NULL, vg.maxdist.Y = NULL
                           , h.min = NULL, h.max = NULL
                           , kernelConv = "gauss"
                           , nSimMMC = 1e3
                           , logVars = TRUE, ... ){

    ##----------------------------------------------------------
    ## Settings and args checking


    ## Either geodata or X-Y must be informed
    if( missing( geodata ) & missing(X) ) stop( "\nMissing geodata object or X and Y data vectors.\n\n" )

    ## If X and Y are missing, they are in geodata object
    if( missing(X) ){
        if( class( geodata ) != "geodata" ) stop( "\ngeodata object must be an object of geodata class.\n\n" )
        if( !all( X <- geodata$data[, 1] ) ) stop( "\ngeodata object must have the data element as a 2 columns matrix.\n\n" )
        if( !all( Y <- geodata$data[, 2] ) ) stop( "\ngeodata object must have the data element as a 2 columns matrix.\n\n" )
    }else{
        if( !is.numeric(X) ) stop( "\nX must be provided as a numeric vector.\n" )
        if( missing(Y) ) stop( "\nMissing Y (it must be provided as a numeric vector).\n\n" )
        if( !is.numeric(Y) ) stop( "\nY must be provided as a numeric vector.\n\n" )
    }

    ## Data lenght
    ##
    nData <- ifelse( length( X ) == length( Y ), length( X ), stop("\nX and Y data vectors must have the same lenght.\n\n") )

    ## 'coords'
    ##
    if( missing(coords) ){
        if( length( geodata$coords ) == 0 ){
            stop( "\nMissing coords or geodata$coords argument.\n\n" )
        }else{
            coords <- geodata$coords
        }
    }
    if( !( all( dim( coords ) == c( nData, 2 ) ) ) ) stop( "\ncoords argument must have length(X) rows and 2 columns.\n\n" )
    if( !( all( apply( coords, 2, class ) %in% c("numeric", "integer") ) ) ) stop( "\ncoords argument must be a numeric or an integer matrix with 2 columns.\n\n" )
    dists <- as.vector( dist( coords ) )

    ## 'shuffle'
    ##
    if( shuffle %in% c("single", "s", 1) ) shuffle <- 1
    if( shuffle %in% c("double", "d", 2) ) shuffle <- 2
    if( !( shuffle %in% c(1,2) ) ){
        shuffle <- 1
        cat("\nWarning: Shuffling only one dataset.\n\n")
    }

    ## 'h.max', 'h.min' - max and min values for convolution window size
    ##
    ## dists <- seq(0.1, 2, by = 0.1)
    ## h.max <- 0.2
    ## h.min <- 0.1

    if( is.null( h.max ) ){
        h.max.X <- h.max.Y <- quantile( dists, probs = 0.5 )
        cat("\nWarning: h.max not setted, using 50th percentil. It could return an error on locfit function.\n Consider to set this argument.\n")
    }else{
        if( length( h.max ) >= 2 ){
            if( !is.numeric( h.max[1] ) ){
                h.max.X <- quantile( dists, probs = 0.5 )
                cat("\nWarning: h.max X wrong value, using 50th percentil. It could return an error on locfit function.\n Consider to set this argument.\n\n")
            }else{
                if( h.max[1] < 0 || h.max[1] > max( dists ) ){
                    h.max.X <- quantile( dists, probs = 0.5 )
                    cat("\nWarning: h.max X wrong value, using 50th percentil. It could return an error on locfit function.\n Consider to set this argument.\n\n")

                }else{
                    h.max.X <- h.max[1]
                }
            }

            if( !is.numeric( h.max[2] ) ){
                h.max.Y <- quantile( dists, probs = 0.5 )
                cat("\nWarning: h.max Y wrong value, using 50th percentil. It could return an error on locfit function.\n Consider to set this argument.\n\n")
            }else{
                if( h.max[2] < 0 || h.max[2] > max( dists ) ){
                    h.max.Y <- quantile( dists, probs = 0.5 )
                    cat("\nWarning: h.max Y wrong value, using 50th percentil. It could return an error on locfit function.\n Consider to set this argument.\n\n")

                }else{
                    h.max.Y <- h.max[2]
                }
            }

        }else{
            if( h.max[1] < 0 || h.max[1] > max( dists ) ){
                h.max.X <- h.max.Y <- quantile( dists, probs = 0.5 )
                cat("\nWarning: h.min X wrong value, using 50th percentil. It could return an error on locfit function.\n Consider to set this argument.\n\n")

            }else{
                h.max.X <- h.max.Y <- h.max[1]
            }
        }
    }

    ## h.min
    if( is.null( h.min ) ){
        h.min.X <- h.min.Y <- quantile( dists, probs = 0.1 )
        cat("\nWarning: h.min not setted, using 10th percentil. It could return an error on locfit function.\n Consider to set this argument.\n\n")
    }else{
        if( length( h.min ) >= 2 ){

            if( !is.numeric( h.min[1] ) ){
                h.min.X <- quantile( dists, probs = 0.1 )
                cat("\nWarning: h.min X wrong value, using 10th percentil. It could return an error on locfit function.\n Consider to set this argument.\n\n")
            }else{

                if( h.min[1] < 0 || h.min[1] > max( dists ) ){
                    h.min.X <- quantile( dists, probs = 0.1 )
                    cat("\nWarning: h.min X wrong value, using 10th percentil. It could return an error on locfit function.\n Consider to set this argument.\n\n")

                }else{

                    h.min.X <- h.min[1]
                }
            }

            if( !is.numeric( h.min[2] ) ){
                h.min.Y <- quantile( dists, probs = 0.1 )
                cat("\nWarning: h.min Y wrong value, using 10th percentil. It could return an error on locfit function.\n Consider to set this argument.\n\n")
            }else{

                if( h.min[2] < 0 || h.min[2] > max( dists ) ){
                    h.min.Y <- quantile( dists, probs = 0.1 )
                    cat("\nWarning: h.min Y wrong value, using 10th percentil. It could return an error on locfit function.\n Consider to set this argument.\n\n")

                }else{

                    h.min.Y <- h.min[2]
                }

            }
        }else{

            if( h.min[1] < 0 || h.min[1] > max( dists ) ){
                h.min.X <- h.min.Y <- quantile( dists, probs = 0.1 )
                cat("\nWarning: h.min X wrong value, using 10th percentil. It could return an error on locfit function.\n Consider to set this argument.\n\n")

            }else{

                h.min.X <- h.min.Y <- h.min[1]
            }
        }
    }
    if( h.min.X >= h.max.X ) stop( "\nh.max for X variable must be smaller than h.min." )
    if( h.min.Y >= h.max.Y ) stop( "\nh.max for Y variable must be smaller than h.min." )

    ## 'kernelConv' - kernel choices, it must be accepted by locfit.raw()
    ##
    if( shuffle == 2 ){
        if( length( kernelConv ) != 2 ) kernelConv <- rep( kernelConv, times = 2 )
        if( !( kernelConv[1] %in% c("tcub", "rect", "trwt", "tria", "epan", "bisq", "gauss") ) ){
            stop( "\nInvalid convolution kernel. See locfit.raw() man.\n\n" )
        }
        if( !( kernelConv[2] %in% c("tcub", "rect", "trwt", "tria", "epan", "bisq", "gauss") ) ){
            stop( "\nInvalid convolution kernel. See locfit.raw() man.\n\n" )
        }

    }else{
        if( !( kernelConv[1] %in% c("tcub", "rect", "trwt", "tria", "epan", "bisq", "gauss") ) ){
            stop( "\nInvalid convolution kernel. See locfit.raw() man.\n\n" )
        }
    }
    ##----------------------------------------------------------
    ## Variog settings
    ##
    ## X var
    if( !is.null( vg.X ) ){
        if( !( any(class( vg.X ) %in% c("variogram") ) ) )
            stop( "\nvg.X must be a returned object from variog{geoR} function.\n\n" )

        vg.target.X <- vg.X

    }else{
        vg.target.X <- variog( data = X, coords = coords
                            , uvec = if( !is.null( vg.uvec.X ) ){ vg.uvec.X }else{ "default" }
                            , max.dist = ifelse( !is.null( vg.maxdist.X ), vg.maxdist.X, max( dists ) )
                            , messages = F )
    }
    ## Y var
    if( shuffle == 2 ){
        if( !is.null( vg.Y ) ){
            if( !( any(class( vg.Y ) %in% c("variogram") ) ) )
                stop( "\nvg.Y must be a returned object from variog{geoR} function.\n\n" )

            vg.target.Y <- vg.Y

        }else{
            vg.target.Y <- variog( data = Y, coords = coords
                                , uvec = if( !is.null( vg.uvec.Y ) ){ vg.uvec.Y }else{ "default" }
                                , max.dist = ifelse( !is.null( vg.maxdist.Y ), vg.maxdist.Y, max( dists ) )
                                , messages = F)
        }
    }


    ## 'nConv'
    ##
    if( nConv < 3 | nConv > 20 ){
        stop( "\nnConv must be an integer between 3 and 20. Using default value 10.\n\n" )
    }

    ## 'nSimMMC' check
    ##
    if( nSimMMC < 1e2 | nSimMMC > 1e4 ){
        nSimMMC <- 1e3
        cat( "\nWarning: Specified nSimMMC isn't valid. Using default value, 1e3.\n\n" )
    }

    ## 'logVars' check
    ##
    if( !is.logical( logVars ) ){
        logVars <- TRUE
        cat( "\nWarning: logVars must be logical. Using default value, TRUE.\n\n" )
    }

    ## r.hat
    r.hat <- cor( X, Y )
    f.call <- match.call( expand.dots = TRUE )

    ##------------------

    ##----------------------------------------------------------
    ##
    ## 1. MMC r values
    ##
    res.core <- .mmc.variog.core( data = X, coords = coords
                               , vg = vg.target.X
                               , nConv = nConv, kernelConv = kernelConv[1]
                               , h.min = h.min.X, h.max = h.max.X
                               , nSimMMC = nSimMMC, log = logVars )
    permuted.X <- res.core$new.data  ## nData x nSim matrix

    if( shuffle == 1 ){
        r.perm <- apply( permuted.X, 2
                      , function( x ){
                          cor( x, Y )
                      } )

    }else{
        res.core.Y <- .mmc.variog.core( data = Y, coords = coords
                                     , vg = vg.target.Y
                                     , nConv = nConv, kernelConv = kernelConv[2]
                                     , h.min = h.min.Y, h.max = h.max.Y
                                     , nSimMMC = nSimMMC, log = logVars )
        permuted.Y <- res.core.Y$new.data

        r.perm <- sapply( 1:nSimMMC
                       , function( s ){
                           cor( permuted.X[, s], permuted.Y[, s] )
                       } )

    }
    ##-------------


    ##----------------------------------------------------------
    ##
    ## 2. Evaluate MMC p-value
    ##
    pv <- length( which(
        ## As extreme as
        abs( r.perm ) >= abs( r.hat )
    ) ) / nSimMMC
    ##-------------


    ##----------------------------------------------------------
    ##
    ## 3. Return
    ##
    res <- list()
    res$r.hat <- r.hat

    ## P-value
    res$p.value <- pv

    ## r from permutations
    res$r.perm <- r.perm

    ## residuals
    if( logVars ){
        if( shuffle == 1 ){
            res$residuals <- res.core$residuals
        }else{
            res$residuals <- res.core.Y$residuals
        }
    }

    ## Explicit return
    return( res )
    ##-------------


}
