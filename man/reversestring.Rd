\name{reversestring}
\alias{reversestring}

\title{
Reverse given string
}
\description{
This is the function for reverse string
}
\usage{
reversestring(string, n = 1)
}

\arguments{
  \item{string}{
Given string or seq
}
  \item{n}{
By which n-plets we should reverse the given string
}
}
\details{
This function is uesd to reverse given string or seq
}
\value{
Reversed string or seq
}

\author{
Md. Bahadur Badsha <mbbadshar@gmail.com>
}

\examples{
x1 <- c("R1C1", "R2C2", "R3C3")
reversestring(x1, 2)
}
