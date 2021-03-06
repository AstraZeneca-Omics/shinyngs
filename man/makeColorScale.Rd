% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/colormaker.R
\name{makeColorScale}
\alias{makeColorScale}
\title{Make a color palette of a specified length}
\usage{
makeColorScale(ncolors, palette = "Set1")
}
\arguments{
\item{ncolors}{Integer specifying the number of colors}

\item{palette}{RColorBrewer palette name. (default: 'Set1')}
}
\value{
output Character vector of colors
}
\description{
Given an integer, make a palette with a specified number of colors using
palettes from RColorBrewer, and interpolation where necessary.
}
\examples{
makeColorScale(10)
[1] '#999999' '#EC83BA' '#B75F49' '#E1C62F' '#FFB716' '#D16948' '#7E6E85' '#48A462' '#4A72A6' '#E41A1C'
}
