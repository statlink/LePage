\name{col.lepage}
\alias{col.lepage}
\title{
Column-wise LePage type tests
}

\description{
Column-wise LePage type tests
}

\usage{
col.lepage(x, y)
}

\arguments{
\item{x}{
A numerical matrix with data, where each column corresponds to the first sample.
}
\item{y}{
A numerical matrix with data, where each column corresponds to the second sample.
}
}

\details{
The classical LePage and 5 variants of it are performed, as proposed by Hussain and Tsagris (2025). 
The test is performed for each pair of columns of the x and y matrices, between the i-th column of and y.
}

\value{
 A list including:
\item{statistic}{
A matrix with 6 columns, where each column contains each test statistic and each row refers to the pair of groups.
} 
\item{pvalue}{
A matrix with 6 columns, where each column contains each p-value and each row refers to the pair of groups.
} 
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
\code{ \link{lepage}
}
}

\examples{
x <- matrix( rnorm(30 * 5), ncol = 5 )
y <- matrix( rnorm(31 * 5), ncol = 5 )
col.lepage(x, y)
}

