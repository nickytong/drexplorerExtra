\name{plotOneExp}
\alias{plotOneExp}
\title{plot dose response curve from fitOneExp() result}
\usage{
  plotOneExp(fitRes, ind2plot = NA, cols = NA,
    type = "plot", style = "full", h = c(0.3, 0.5, 0.7),
    tag = NA, main = NA, cex.main = 1, xlab = NA,
    ylab = NA, ylim = NA, xlim = NA, show = "both",
    cexLegend = NA, showTopN = NA, lwd = 2)
}
\arguments{
  \item{fitRes}{return value from fitOneExp()}

  \item{ind2plot}{index for the models that will be
  plotted; default is NA which leads to all curves
  available; when specified as 'best', the best model is
  selected}

  \item{col}{color for the lines}

  \item{type}{either plot or line; when specified as line,
  it will only adds to an existing figure; When
  length(ind2plot)>1, type will be reset to plot which
  means the first curve will be made with plot() and
  additional ones with lines()}

  \item{h}{horizontal line added to the figure, i.e.
  indicating IC50, IC70}

  \item{tag}{tag before main}

  \item{main}{main}

  \item{cex.main}{cex.main to adjust main title size}

  \item{ylim}{ylim}

  \item{xlim}{xlim}

  \item{xlab}{xlab}

  \item{ylab}{ylab}

  \item{style}{if style == 'full', the outlier status as
  well as its legend will be shown; if style=='simple', the
  points and legend for outlier status will be removed; if
  style=='points', only points will be shown; this is
  useful if to compare multiple curves from different
  drugs.}

  \item{show}{either RSE, IC50 or both indicating if we
  need to show RSE and/or IC50 in addition to the model}

  \item{cexLegend}{legend cex}

  \item{showTopN}{if specified show best N model in figure
  to avoid busy plotting; otherwise show all successful
  models.}

  \item{lwd}{line width for the curves}
}
\description{
  plot dose response curve from fitOneExp() result
}

