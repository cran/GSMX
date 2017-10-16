\name{gsm}
\alias{gsm}
\title{
	Esimate genetic variance-covariances and calculate heritability and predictability using cross validation
}
\description{
The function estimate genetic variance-covariances and calculate heritability and predictability for multivariate genetic selection using cross validation
}
\usage{
gsm(mydata, mykin, nfold)
}
\arguments{
  \item{mydata}{
dataset with two traits
}
  \item{mykin}{
kinship matrix
}
  \item{nfold}{
number of folds for cross validation}
}

\value{
\item{res}{Results}
}

\examples{
library(GSMX)
data(pseudo.kin)
data(pseudo.data)
myfit=gsm(pseudo.data, pseudo.kin, nfold=5)
}
\keyword{Multivariate Genomic Selection}
\keyword{Overfitting}
\keyword{Cross Validation}
