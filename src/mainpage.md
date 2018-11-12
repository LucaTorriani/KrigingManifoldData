# What is Manifoldgstat Package #
R Package to make inference and prediction for manifold-valued data analysis. This package provides a C++ implementation of the functions to create a model, for spatial dependent	manifold valued data, in order to perform kriging.
In each location, specified by a vector of coordinates ([lat,long], [x,y] or [x,y,z]), the datum is supposed to be a symmetric and positive definite matrix (possibly a correlation matrix). The algorithm exploits a projection of these data on a tangent space, where the tangent point is either provided by the user or computed as intrinsic mean of the data in input.

### References
JD. Pigoli, A. Menafoglio & P. Secchi (2016). Kriging prediction for manifold-valued random fields. Journal of Multivariate Analysis, 145, 117-131.
