# Calculates the weighted correlation

# @name WC
# @description Function for computing a weighted Pearson's Correlation
WC = function( x, y, w = rep(1,length(x)) ) {
    stopifnot( length(x) == dim(y)[2] )
    w = w / sum(w)
    # Center x and y, using the weighted means
    x  = x - sum(x*w)
    ty = t(y - colSums(t(y) * w))
    # Compute the variance
    vy = sum( w * x * x )
    vy = colSums(w * ty * ty)
    # Compute the covariance
    vxy = colsums(ty * x * w)
    # Compute the correlation
    vxy / sqrt(vx * vy)
}
