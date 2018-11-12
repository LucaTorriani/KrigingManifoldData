# Manifoldgstat Package [<img style="float: right;" src="https://www.r-project.org/logo/Rlogo.png" width=70 alt="R logo"/>](https://www.r-project.org)
### Course on Advanced programming for scientific computing
#### *Politecnico di Milano* (ITALY)
[<img style="float: right;" src="https://github.com/fpavone/pacs_spline_density/blob/ROpenMP/rePortACS/pictures/logopoli.png" width=50 alt="Polimi logo"/>](http://www.polimi.it/)

***Author*** : Ilaria Sartori, Luca Torriani

***Mailto*** : ilaria.sartori@mail.polimi.it, luca1.torriani@mail.polimi.it

***Date***   : November 2018

The output is an R Package that performs kriging prediction for manifold-valued random fields [Manifoldgstat]

The source code is written in C++ and linked to R throught the API `RcppEigen` and `.Call`

#### Subfolder structure

- `src` contains all C++ code and a special file named Makevars necessary to build and install the R package
- `R` contains the R functions that wrap the C++ calls
- `data` contains .RData files useful for testings
- `tests` contains basic R script to run tests

To install the package, please make sure that you have the package [devtools](https://cran.r-project.org/web/packages/devtools/index.html) already installed.

From the root folder then type

    R -e "library(devtools); install()" —-silent
    R -e "library(devtools); document()" —-silent

or in another way download source code from the github repository, unzip the file and run from the terminal:

    R CMD INSTALL <path name of the package to be installed> -l <path name of the R library tree>

Otherwise you can install directly from github:

    devtools::install_github("LucaTorriani/KrigingManifoldData", ref = "master")

To run an example in R:

        example(model_GLS, "Manifoldgstat")
        example(kriging, "Manifoldgstat")
        example(model_kriging, "Manifoldgstat")
