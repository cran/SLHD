\name{maximinSLHD}
\alias{maximinSLHD}

\title{
Maximin-Distance (Sliced) Latin Hypercube Designs
}
\description{
Generate the optimal Latin hypercube designs when \eqn{t=1} and the optimal sliced Latin hypercube designs when \eqn{t>1} for computer experiments. The maximin distance criterion is adopted as the optimality criterion.
}
\usage{
maximinSLHD(t, m, k, power = 15, nstarts = 1, itermax = 100, total_iter = 1e+06)
}

\arguments{
  \item{t}{
The number of slices. If \eqn{t=1}, it leads to a standard space-filling Latin hybercube design
}
  \item{m}{
The number of design points (runs) in each slice. When \eqn{t=1}, the number of design points for the whole design is just \eqn{m}; when \eqn{t>1}, the number of design points for the whole design is \eqn{mt}
}
  \item{k}{
The number of input factors (variables)
}
  \item{power}{
Optional, default is ``15''. The power parameter \eqn{r} in the average reciprocal inter-point distance measure. When \eqn{r\rightarrow \infty}, minimizing the average reciprocal inter-point distance measure is equivalent to maximizing the minimum distance among the design points.
}
  \item{nstarts}{
Optional, default is ``1''. The number of random starts
}
  \item{itermax}{
Optional, default is ``100''. The maximum number of non-improving searches allowed under each temperature. Lower this parameter if you want the algorithm to converge faster 
}
  \item{total_iter}{
Optional, default is ``1e+06''.The maximum total number of iterations. Lower this number if the design is prohibitively large and you want to terminate the algorithm prematurely to report the best design found so far. 
}
}
\details{
This function utilizes a version of the simulated annealing algorithm and several computational shortcuts to efficiently generate the optimal Latin Hypercube Designs (LHDs) and the optimal Sliced Latin Hypercube Designs (SLHDs). The maximin distance criterion is adopted as the optimality criterion. Please refer to Ba (2013) for details of the algorithm. 


When \eqn{t=1}, the maximin-distance LHD is popularly used for designing computer experiments with quantitative factors.  

When \eqn{t>1}, the maximin-distance SLHD is a special class of LHD which can be partitioned into several slices (blocks), each of which is also a LHD of smaller size. The optimal SLHD structure ensures the uniformity (space-filling property) in each slice as well as in the whole design. The SLHD is very important in designing computer experiments with quantitative and qualitative factors, where each slice is used as a design for quantitative factors under one of the t different level combinations of qualitative factors.
}
\value{
The value returned from the function is a list containing the following components:
\item{D}{Design matrix}
\item{DS}{A list containing the \eqn{t} small designs for the \eqn{t} slices}
\item{standD}{Standardized design matrix which is scaled into the \eqn{[0,1]^k} unit region}
\item{standDS}{A list containing the \eqn{t} standardized small designs for the \eqn{t} slices}
\item{temp0}{Initial temperature}
\item{measure}{The average reciprocal inter-point distance measure}
\item{time_rec}{Time to complete the search}
}
\references{
Ba, S. (2013), "Optimal Sliced Latin Hypercube Designs," revised for \emph{Technometrics}.
}
\author{
Shan Ba<shanbatr@gmail.com>
}


%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{

#Maximin-Distance Latin Hypercube Designs
D1<-maximinSLHD(t = 1, m = 10, k = 3) 
D1$D
D1$standD

#Maximin Distance Sliced Latin Hypercube Designs
D2<-maximinSLHD(t = 3, m = 4, k = 2) 
D2$DS
D2$D
}

\keyword{ Design of Experiments }
\keyword{ Computer Experiments }
