\name{galts-package}
\alias{galts-package}
\alias{galts}
\docType{package}
\title{
Genetic algorithms and C-steps based LTS algorithm.
}
\description{
This package includes the ga.lts function that estimates the LTS estimator using genetic algorithms and C-steps.
}

\details{
\tabular{ll}{
Package: \tab galts\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2010-11-11\cr
License: \tab GPL\cr
LazyLoad: \tab yes\cr
}
}
\author{
Mehmet Hakan Satman

Maintainer: Mehmet Hakan Satman <mhsatman@istanbul.edu.tr>
}
\references{
Rousseeuw and Van Driessen (2006), "Computing LTS Regression for Large Data Sets", Data Mining and Knowledge Discovery, 12, 29-45
}
\keyword{ package }

\examples{
n<-300;
x1<-rnorm(n);
x2<-rnorm(n);
x3<-rnorm(n);
x4<-rnorm(n);
x5<-rnorm(n);
x<-cbind(1,x1,x2,x3,x4,x5);
e<-rnorm(n);
p<-dim(x)[2];
betas<-rep(5,p);
y<-(x\%*\%betas)+e;
h<-floor(n/2)+floor((p+1)/2);
y[(h+1):n]<-1000;


cat("Starting algorihm with configuration:\n");
cat("n=",n,"\n");
cat("p=",p,"\n");
cat("Regression parameters=",betas,"\n");
cat("Last ",n-h," observations are outliers\n");
cat("Doing...\n");

lts<-ga.lts(formula=y~x1+x2+x3+x4+x5, lower=-20, upper=20);
cat("Genetic algorithm LTS with csteps\n");
cat("Objective function: ",lts$crit,"\n");
cat("Coefficients: ",lts$coefficients,"\n");
}
