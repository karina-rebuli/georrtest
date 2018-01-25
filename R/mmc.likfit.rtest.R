#'
#' Parametrical Monte Carlo method for testing Pearson correlation coeficient (r) significance under null
#' hypothesis (\eqn{H_0: r = 0}) between spatial autodependent data. Permuted and convoluted data is choosen
#' as MMC step by comparing its MLE spatial model parameters  with those from original data.
#'
#' The algorithm follow the steps:
#'
#' 1. Permute one (or both) autocorrelated data.
#' 2. Repeat 3 to 5 N iteration times
#' 3. Convolute permuted data to reinsert autodependence on it.
#' 4. Choose the convoluted data which MLE parameters is closer (difference) to the origianl data MLE parameters.
#' 5. Evaluate correlation between data.
#'
#'
#' @param geodata An object of geodata class. Mandatory if X or Y are not setted.
#' @param X X data vector. Mandatory.
#' @param Y Y data vector. Mandatory.
#' @param coords Data grid (locations). Mandatory.
#' @param shuffle If it should be permuted one or both datasets. Allowed values: "single", "s" or 1 and "double", "d" or 2. Default( 1 ).
#' @param nConv Number of different window size convolutions must be done. Allowed values are integers between 3 and 20. Default( 10 ).
#' @param dist.vector A vector with distances. Default( NULL ).
#' @param likfit.X likfit{geoR} X data returned object.
#' @param likfit.Y likfit{geoR} Y data returned object.
#' @param h.min Smallest locfit delta (h) parameter. If it's not setted and if nn.max ins't setted
#' , it will be the value corresponding to 10th percentil of distances vector.
#' @param h.max Biggest locfit delta (h) parameter. If it's not setted and if nn.max ins't setted
#' , it will be the value corresponding to 50th percentil of distances vector.
#' @param kernelConv Weight function for convolutions; check man locfit package to see allowed values.
#' @param nSimMMC Number of MMC simulations to be done. Allowed integers from 1e2 to 1e4. Default( 1e3 ).
#' @param logVars logical, if the returned object have a vector with estimated MMC r. Default( TRUE ).
#'
#'
#' @return r.hat Sample Pearson correlation coefficient.
#' @return p.value The empirical p-value obtained.
#' @return r.perm A vector with estimated Pearson correlation coefficient in each MMC permutation.
#'
#' @author Karina Brotto rebuli and Paulo Justiniano Ribeiro Jr.
#'
#' @note If neither likfit._ or cov.pars._ is provided, the function will try to estimate spatial
#' model parameters using default settings with likfit\{geoR\} function. It's not a recommended use! Since the
#' autodependence structure affects the correlation coefficient distribution, a bad parameters estimation can lead to misleading results.
#'
#' @keywords mmc, correlation coeficient significance, geostatistical data, empirical variograms
#'
#' @export
#'
mmc.likfit.rtest <-
function( geodata, X, Y, coords
                           , dist.vector = NULL
                           , shuffle = 1, nConv = 10
                           , likfit.X = NULL
                           , likfit.Y = NULL
                           , h.min = NULL, h.max = NULL
                           , kernelConv = "gauss"
                           , nSimMMC = 1e3
                           , logVars = TRUE, ... ){

    ##----------------------------------------------------------
    ## Settings
    ##

    ## X and Y data
    ##

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
    dist.vector <- if( is.null( dist.vector ) ) as.vector( dist( coords ) )
    max.dist <- max( dist.vector )
    min.dist <- min( dist.vector )

    ## 'shuffle'
    ##
    if( shuffle %in% c("single", "s", 1) ) shuffle <- 1
    if( shuffle %in% c("double", "d", 2) ) shuffle <- 2
    if( !( shuffle %in% c(1,2) ) ){
        shuffle <- 1
        cat("\nWarning: Shuffling only one dataset.\n\n")
    }

    ## 'nConv'
    ##
    if( nConv < 3 | nConv > 20 ){
        stop( "\nnConv must be an integer between 3 and 20. Using default value 10.\n\n" )
    }

    ## 'h.max', 'h.min' - max and min values for convolution window size
    ##
    if( is.null( h.max ) ){
        h.max.X <- h.max.Y <- quantile( dist.vector, probs = 0.5 )
        cat("\nWarning: h.max not setted, using 50th percentil. It could be return error on locfit function.\n Consider to set this argument.\n\n")
    }else{
        if( length( h.max ) >= 2 ){
            h.max.X <- as.list( h.max )[[1]]
            h.max.Y <- as.list( h.max )[[2]]
        }else{
            h.max.X <- h.max.Y <- as.list( h.max )[[1]]
        }
    }
    if( is.null( h.min ) ){
        h.min.X <- h.min.Y <- quantile( dist.vector, probs = 0.1 )
        cat("\nWarning: h.min not setted, using 10th percentil. It could be return error on locfit function.\n Consider to set this argument.\n\n")
    }else{
        if( length( h.min ) >= 2 ){
            h.min.X <- as.list( h.min )[[1]]
            h.min.Y <- as.list( h.min )[[2]]
        }else{
            h.min.X <- h.min.Y <- as.list( h.min )[[1]]
        }
    }

    ## 'kernelConv' - kernel choices, it must be accepted by locfit.raw()
    ##
    if( !( kernelConv %in% c("tcub", "rect", "trwt", "tria", "epan", "bisq", "gauss") ) ) stop( "\nInvalid convolution kernel. See locfit.raw() man.\n\n" )

    ##
    ## Models check
    ##
    ## There are two options:
    ## 1 - likfit() return object
    ## 2 - Automatic estimation. It will use likfit(), but is not recommended
    ##

    ## 1 - 'likfit.X' - likfit() return object
    ##
    if( !is.null( likfit.X ) ){
        if( !( any(class( likfit.X ) %in% c("likGRF", "variomodel") ) ) )
            stop( "\nlikfit.X must be a returned objetc from likfit{geoR} function.\n\n" )
    }else{

        ## 2 - Run automatic likfit
        ##

        ## Variogram fitting to start cov.pars on ML pars estimation
        variog.X <- variog( data = X, coords = coords, messages = FALSE )
        ini.pars.X <- variofit( variog.X, messages = FALSE )$cov.pars
        if( ini.pars.X[1] <= 0.1 | ini.pars.X[2] <= 0.03 ){
            stop( "\nAutomatic parameters estimations has failed.\nRun model parameters estimations manually.\n\n" ) }

        ## MLE parameters estimation
        likfit.X <- likfit( data = X, coords = coords
                         , ini.cov.pars = ini.pars.X, messages = FALSE )
        if( likfit.X$sigmasq < 1e-4 | likfit.X$sigmasq < 1e-4 ){
            stop( "\nAutomatic parameters estimations for X variable has failed.\nRun model parameters estimations manually.\n\n" ) }

        ## WARNING
        msg.war <- "\nUsing automatic estimated parameters, the analysis can be spoilt.\nConsider to estimate manually the model parameters.\n\n"
        cat( msg.war )
        warning( msg.war )

    }

    ## 1 - 'likfit.Y' - likfit() return object
    ##
    if( shuffle == 2 ){
        if( !is.null( likfit.Y ) ){
            if( !( any(class( likfit.Y ) %in% c("likGRF", "variomodel") ) ) )
                stop( "\nlikfit.Y must be a returned objetc from likfit{geoR} function.\n\n" )
        }else{

            ## 2 - Run automatic likfit
            ##

            ## Variogram fitting to start cov.pars on ML pars estimation
            variog.Y <- variog( data = Y, coords = coords, messages = FALSE )
            ini.pars.Y <- variofit( variog.Y, messages = FALSE )$cov.pars
            if( ini.pars.Y[1] <= 0.1 | ini.pars.Y[2] <= 0.03 ){
                stop( "\nAutomatic parameters estimations for Y variable has failed.\nRun model parameters estimations manually.\n\n" ) }

            ## MLE parameters estimation
            likfit.Y <- likfit( data = Y, coords = coords
                             , ini.cov.pars = ini.pars.Y, messages = FALSE )
            if( likfit.Y$sigmasq < 1e-4 | likfit.Y$sigmasq < 1e-4 ){
                stop( "\nAutomatic parameters estimations has failed.\nRun model parameters estimations manually.\n\n" ) }

            ## WARNING
            msg.war <- "\nUsing automatic estimated parameters, the analysis can be spoilt.\nConsider to estimate manually the model parameters.\n\n"
            cat( msg.war )
            warning( msg.war )

        }
    }

    ## 'nSimMMC' check
    if( nSimMMC < 1e2 | nSimMMC > 1e4 ){
        nSimMMC <- 1e3
        cat( "\nSpecifeid nSimMMC isn't valid. Using default value, 1e3.\n\n" )
    }

    ## 'logVars' check
    if( !is.logical( logVars ) ){
        logVars <- TRUE
        cat( "\nlogVars must be logical. Using default value, TRUE.\n\n" )
    }

    ## r.hat
    r.hat <- cor( X, Y )
    f.call <- match.call( expand.dots = TRUE )
    ##-------------


    ##----------------------------------------------------------
    ##
    ## 1. MMC r values
    ##
    res.core <- .mmc.likfit.core( data = X, coords = coords
                               , ml = likfit.X
                               , nConv = nConv, kernelConv = kernelConv
                               , h.min = h.min.X
                               , h.max = h.max.X
                               , nSimMMC = nSimMMC, log = logVars )
    permuted.X <- res.core$new.data

    if( shuffle == 1 ){
        r.perm <- apply( permuted.X, 2
                        , function( x ){
                            cor( x, Y )
                        } )
    }else{
        res.core.Y <- .mmc.likfit.core( data = Y, coords = coords
                                     , ml = likfit.Y
                                     , nConv = nConv, kernelConv = kernelConv
                                     , h.min = h.min.Y
                                     , h.max = h.max.Y
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
            res$differences.phis <- res.core$differences.phis
        }else{
            res$differences.phis <- res.core.Y$differences.phis
        }
    }

    ## Explicit return
    return( res )
    ##-------------

}
