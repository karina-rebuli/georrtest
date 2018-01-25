#'
#' Internal convolution of permuted data function.
#'
#' @export
#'
.conv.permuted <-
  function( data, coords = NULL
            , dists = if( is.null( coords ) ) stop( "\nError on conv.permuted() arguments.\n" )
            , h.min = NULL
            , h.max = NULL
            , nConv = 10, kernelConv = "gauss" ){

    ##------------------------------------------------------------
    ## Settings
    if( is.null( dists ) ) dists <- as.vector( dist( coords ) )
    nn.max <- NULL

    ##------------------------------------------------------------
    ## Check max values for convolution window size
    if( is.null( h.max ) && is.null( nn.max ) ){
      h.max <- quantile( dists, probs = 0.5 )
    }
    if( is.null( h.min ) & is.null( nn.max ) ){
      h.min <- quantile( dists, probs = 0.1 )
      warning("\nh.min not setted, using 10th percentil. It could be return error on locfit function.\n Consider to set this argument.\n")
    }

    ##----------------------------------------------------------
    ## Check kernel choices
    ## It must be accepted by locfit.raw()
    if( !( kernelConv %in% c("tcub", "rect", "trwt", "tria", "epan", "bisq", "gauss") ) ) stop( "\nInvalid convolution kernel. See locfit.raw() man.\n" )

    ##------------------------------------------------------------
    ## Ia. Delta kernel convolution
    if( !is.null( h.max ) ){
      h <- seq( h.min, h.max, l = nConv )
      fit <- lapply( h
                     , function(x){
                       locfit( data ~ lp( coords[, 1], coords[, 2], h = x, deg = 0 ), kern = kernelConv )
                     } )

      ##------------------------------------------------------------
      ## Ib. Nearest neighbors convolution
    }else{

      nn <- seq( ( 1/length(data) )*5, nn.max, l = nConv )
      fit <- lapply( nn
                     , function(x){
                       locfit( data ~ lp( coords[, 1], coords[, 2], nn = x, deg = 0 ), kern = kernelConv )
                     } )
    }

    ##------------------------------------------------------------
    ## II. Explicit return
    return( fit )

  }
