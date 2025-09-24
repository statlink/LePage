\name{lepage}
\alias{lepage}
\title{
LePage type tests
}

\description{
LePage type tests
}

\usage{
lepage(x, y)
}

\arguments{
\item{x}{
A numerical vector with data of the first sample.
}
\item{y}{
A numerical vector with data of the second sample.
}
}

\details{
The classical LePage and 5 variants of it are performed, as proposed by Hussain and Tsagris (2025).
}

\value{
A matrix with 6 columns and two rows. The first row refers to the 6 test statistics and the second row to their corresponding p-values.
}

\references{
Hussain A. and Tsagris M. (2025). Modified Lepage-type test statistics for the weak null hypothesis. 
https://arxiv.org/pdf/2509.19126
}

\author{
Abid Hussain and Michail Tsagris.

R implementation and documentation: Abid Hussain \email{abid0100@gmail.com} and Michail Tsagris \email{mtsagris@uoc.gr}.
}

\seealso{
\code{ \link{col.lepage}
}
}

\examples{
x <- rnorm(30)
y <- rnorm(40)
lepage(x, y)
}

