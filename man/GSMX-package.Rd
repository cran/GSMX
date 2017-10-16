\name{GSMX-package}
\alias{GSMX-package}
\alias{GSMX}
\docType{package}
\title{
Multivariate Genomic Selection
}
\description{
The package \pkg{GSMX} consists of the functions to estimate genetic variance-covariances and calculate trait heritability, and handle overfitting by calculating the variance components and the heritability through cross validation.
}
\details{
\tabular{ll}{
Package: \tab GSMX\cr
Type: \tab GSMX\cr
Version: \tab 1.3\cr
Date: \tab 2017-10-16\cr
License: GPL>=2 \cr
}
}

\section{Control overfitting heritability in genomic selection through cross validation}{
Genomic selection (GS) is a form of marker-assisted selection (MAS) where markers across the entire genome are used such that all quantitative trait loci (QTL) are in linkage disequilibrium with at least one marker. 
Increased number of markers and their density along with increased sample size improve the resolution of QTL mapping, and therefore improve estimated breeding values and genetic gain. 
However, the genetic variance-covariances were estimated from training samples using a large number of markers including many trait-irrelevant markers, and then being used for calculating trait heritability, leading to severe overfitting. 
In this package, we developed an algorithm to handle such overfitting by calculating the variance components and the heritability through cross validation. 
This method provides an accurate estimation of trait heritability (equivalent to trait predictability), and objectively reflects the level of applicability of the GS models to other breeding materials..  \cr
}

\examples{
library(GSMX)
data(pseudo.kin)
data(pseudo.data)
myfit=gsm(pseudo.data, pseudo.kin, nfold=5)
}

\references{Control overfitting heritability in genomic selection through cross validation}

\author{Zhenyu Jia
}

\keyword{package}
