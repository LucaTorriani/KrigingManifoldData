# What is Manifoldgstat Package #
Predictive analysis for manifold-valued data.
This package provides	a C++ implementation of functions to create a model for spatial dependent manifold-valued data, in order to perform kriging. In each location, specified by a vector of coordinates ([lat,long], [x,y] or [x,y,z]), the datum is supposed to be a symmetric positive definite matrix.
The user is provided with three main functions: model_kriging, full_RDD, mixed_RDD, each designed to deliver kriging predictions following the corresponding algorithm (GlobalModel, FullRDD and MixedRDD), as presented in the reference dissertation. They exploit, to different extents, tangent space approximations, Random Domain Decomposition and advanced differential geometry concepts like parallel transport.

### References
Ilaria Sartori Luca Torriani (2019): Mixed Random Domain Decomposition: an innovative approach for kriging prediction of manifold valued data, Master Degree Thesis
