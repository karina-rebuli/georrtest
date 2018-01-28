#'
#' Model-based MMC for testing significance of Pearson correlation coeficient on spatial
#' autodependent data
#'
#' Model-based Monte Carlo method for testing Pearson correlation coeficient (r) significance under null
#' hypothesis (\eqn{H_0: r = 0}) between spatial autodependent data. It uses grf\{geoR\} to simulate GRFs
#' with specified parameters and evaluates r between X and Y to each simulation.
#'
#'
#' @param geodata An object of geodata class. Mandatory if X or Y are not setted.
#' @param X X data vector. Mandatory.
#' @param Y Y data vector. Mandatory.
#' @param coords Data grid (locations). Mandatory.
#' @param shuffle If it should be permuted one or both datasets. Allowed values: "single", "s" or 1 and "double", "d" or 2. Default \code{1}.
#' @param lambda Box-Cox transformation parameter. A numeric with shuffle lenght. Default \code{1}.
#' @param likfit.X likfit{geoR} X data returned object. Default \code{NULL}.
#' @param likfit.Y likfit{geoR} Y data returned object. Default \code{NULL}.
#' @param cov.model A character with shuffle lenght. It's a string indicating the type of the
#' correlation function. Allowed values: see cov.spatial{geoR} manual. Default \code{"exp"}.
#' @param kappa A numeric with shuffle lenght, it's the smoothness parameter of the correlation function.
#' Only required by the following cov.model: "matern", "powered.exponential", "cauchy", "gencauchy"
#' and "gneiting.matern". Default \code{NULL}.
#' @param nugget tau^2 parameter for cov.model. A numeric with shuffle lenght. Default \code{0}.
#' @param cov.pars.X numeric(2) vector with sigma^2 and phi from X model. Default \code{NULL}.
#' @param cov.pars.Y numeric(2) vector with sigma^2 and phi from Y model. Default \code{NULL}.
#' @param nSimMMC Number of MMC simulations to be done. Allowed integers from 1e2 to 1e4.Default \code{1e3}.
#'
#'
#' @return r.hat Sample Pearson correlation coefficient.
#' @return p.value The empirical p-value obtained.
#' @return r.perm A vector with estimated Pearson correlation coefficient in each MMC permutation.
#' @return differences.phis A vector with differences between phi OLS estimated with original data and with choosen permuted and
#' transformed data used as MMC permutation.
#'
#' @author Karina Brotto rebuli and Paulo Justiniano Ribeiro Jr.
#'
#' @note Warning: The function can be time consuming!
#' @note Warning: If likfit._ is not provided, this function will try to estimate spatial model parameters with default settings
#' with likfit\{geoR\} function. It's not a recommended use! Since the autodependence structure affects the correlation coefficient
#' distribution, a bad parameters estimation can lead to misleading results.
#'
#' @keywords mmc, correlation coeficient significance, geostatistical data, empirical variograms
#'
#' @export
#'
mmc.modelbased.rtest <-
function( geodata, X, Y, coords
       , shuffle = 1
       , lambda = 1
       , likfit.X = NULL
       , likfit.Y = NULL
       , cov.model = "exp"
       , kappa = NULL
       , nugget = 0
       , cov.pars.X = NULL
       , cov.pars.Y = NULL
       , nSimMMC = 1e3, ... ){

    ##----------------------------------------------------------
    ## Settings and args checking
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


    ## 'shuffle'
    ##
    if( shuffle %in% c("single", "s", 1) ) shuffle <- 1
    if( shuffle %in% c("double", "d", 2) ) shuffle <- 2
    if( !( shuffle %in% c(1,2) ) ){
        shuffle <- 1
        cat("\nWarning: Shuffling only one dataset.\n\n")
    }

    ##
    ## Models check
    ##
    ## There are three options:
    ## 1 - likfit() return object
    ## 2 - Manually parameters setted
    ## 3 - Automatic estimation. It will use likfit(), but is not recommended
    ##

    ## 1 - 'likfit.X' - likfit() return object
    ##
    if( !is.null( likfit.X ) ){
        if( !( any(class( likfit.X ) %in% c("likGRF", "variomodel") ) ) )
            stop( "\nlikfit.X must be a returned objetc from likfit{geoR} function.\n\n" )

    }else{

        ## 2 - Manually model parameters settings
        ##
        if( !is.null( cov.pars.X ) ){

            ## Make a list with model parameters with the same struture as likfit{geoR} return.
            likfit.X <- list()

            ## 'cov.model'
            if( !( cov.model[1] %in% c( "exp", "matern", "exponential", "gaussian", "spherical", "circular", "cubic", "wave", "power", "powered.exponential", "cauchy", "gencauchy", "gneiting", "gneiting.matern", "pure.nugget" ) ) ){
                cov.model[1] <- "exp"
                cat("\nWarning: cov.model for X variable wasn't well setted. Using 'exp' cov.model. To see allowed values, see cov.spatial{geoR} manual.\n\n")
            }
            likfit.X$cov.model <- cov.model[1]

            ## 'kappa'
            if( cov.model[1] %in% c("matern", "powered.exponential", "cauchy", "gencauchy", "gneiting.matern" ) ){
                if( is.null( kappa[1] ) ){
                    if( likfit.X$cov.model == "exp" ){
                        kappa[1] <- 0.5
                    } else{
                        stop( "\nTo \"matern\", \"powered.exponential\", \"cauchy\", \"gencauchy\", \"gneiting.matern\" cov.models it's required to set kappa value.\n\n" )
                    }
                }
            }
            likfit.X$kappa <- kappa[1]

            ## 'cov.pars'
            if( length( cov.pars.X ) != 2 ){
                stop( "\nBad cov.pars.X setting. It must be a numeric(2), corresponding to sigma^2 and phi parameters.\n\n" )
            }
            likfit.X$cov.pars <- cov.pars.X

            ## 'lambda' and 'nugget'
            likfit.X$lambda <- lambda[1]
            likfit.X$nugget <- nugget[1]

        }else{

            ## 3 - Run automatic likfit
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

    }

    ## 1 - 'likfit.Y' - likfit() return object
    ##
    if( shuffle == 2 ){
        if( !is.null( likfit.Y ) ){
            if( !( any(class( likfit.Y ) %in% c("likGRF", "variomodel") ) ) )
                stop( "\nlikfit.Y must be a returned objetc from likfit{geoR} function.\n\n" )

        }else{

            ## 2 - Manually model parameters settings
            ##
            if( !is.null( cov.pars.Y ) ){

                ## Make a list with model parameters with the same struture as likfit{geoR} return.
                likfit.Y <- list()

                ## 'cov.model' check
                if( !( cov.model[2] %in% c( "exp", "matern", "exponential", "gaussian", "spherical", "circular", "cubic", "wave", "power", "powered.exponential", "cauchy", "gencauchy", "gneiting", "gneiting.matern", "pure.nugget" ) ) ){
                    cov.model[2] <- "exp"
                    cat("\nWarning: cov.model for Y variable wasn't well setted. Using 'exp' cov.model. To see allowed values, see cov.spatial{geoR} manual.\n\n")
                }
                likfit.Y$cov.model <- cov.model[2]

                ## 'kappa' check
                if( cov.model[2] %in% c("matern", "powered.exponential", "cauchy", "gencauchy", "gneiting.matern" ) ){
                    if( is.null( kappa[2] ) ){
                        if( likfit.Y$cov.model == "exp" ){
                            kappa[2] <- 0.5
                        } else{
                            stop( "\nTo \"matern\", \"powered.exponential\", \"cauchy\", \"gencauchy\", \"gneiting.matern\" cov.models it's required to set kappa value.\n\n" )
                        }
                    }
                }
                likfit.Y$kappa <- kappa[2]

                ## 'nugget' check
                if( length( nugget ) != 2 ){
                    nugget <- rep( nugget[1], 2 )
                    cat( "\nWarning: Using the same nugget value for both models.\n\n" )
                }
                likfit.Y$nugget <- nugget[2]

                ## 'cov.pars' checks
                if( length( cov.pars.Y ) != 2 ){
                    stop( "\nBad cov.pars.Y setting. It must be a numeric(2), corresponding to sigma^2 and phi parameters.\n\n" )
                }
                likfit.Y$cov.pars <- cov.pars.Y

                ## 'lambda'
                likfit.Y$lambda <- lambda[2]

            }
            else{

                ## 3 - Run automatic likfit
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
    }

    ## 'nSimMMC' check
    ##
    if( nSimMMC < 1e2 | nSimMMC > 1e4 ){
        nSimMMC <- 1e3
        cat( "\nWarning: Specifeid nSimMMC isn't valid. Using default value, 1e3.\n\n" )
    }
    ##-------------


    ##----------------------------------------------------------
    ## 1. Model based simulation
    ##
    pars.X <- as.numeric( c( likfit.X$cov.pars[1], likfit.X$cov.pars[2] ) )
    X.sim <- grf( n = nData, cov.pars = pars.X, nsim = nSimMMC, grid = coords
               , cov.model = likfit.X$cov.model, kappa = likfit.X$kappa,  nugget = likfit.X$nugget
               , lambda = likfit.X$lambda, messages = FALSE )
    if( shuffle == 2 ){
        pars.Y <- as.numeric( c( likfit.Y$cov.pars[1], likfit.X$cov.pars[2] ) )
        Y.sim <- grf( n = nData, cov.pars = pars.Y, nsim = nSimMMC, grid = coords
               , cov.model = likfit.Y$cov.model, kappa = likfit.Y$kappa,  nugget = likfit.Y$nugget
               , lambda = likfit.Y$lambda, messages = FALSE )
    }

    ## h.rat
    r.hat <- cor( X, Y )
    f.call <- match.call( expand.dots = TRUE )
    ##-------------


    ##----------------------------------------------------------
    ##
    ## 2. Record estimated r with this 's' permutation
    ##
    if( shuffle == 1 ){
        r.perm <- sapply( 1:nSimMMC
                       , function( i ){
                           cor( X.sim$data[, i], Y ) } )
    }
    if( shuffle == 2 ){
        r.perm <- sapply( 1:nSimMMC
                       , function( i ){
                           cor( X.sim$data[, i], Y.sim$data[, i] ) } )
    }
    ##-------------


    ##----------------------------------------------------------
    ##
    ## 3. Evaluate MMC p-value
    ##
    pv <- length( which(
        ## As extreme as
        abs( r.perm ) >= abs( r.hat )
    ) ) / nSimMMC
    ##-------------


    ##----------------------------------------------------------
    ##
    ## 4. Return
    ##
    res <- list()
    res$r.hat <- r.hat
    res$p.value <- pv
    res$r.perm <- r.perm

    ## Explicit return
    return( res )
    ##-------------

}
