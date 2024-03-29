% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cin_metrics.R
\name{countingBreakPoints}
\alias{countingBreakPoints}
\title{countingBreakPoints}
\usage{
countingBreakPoints(
  cnvData,
  sample_column = 1,
  segmentMean = 0.2,
  numProbes = NA
)
}
\arguments{
\item{cnvData}{dataframe containing following columns: Sample, Start, End, Num_Probes, Segment_Mean}

\item{sample_column}{numerical value for the sample column number, default is 1}

\item{segmentMean}{numerical value for the minimum segment_mean cutoff/ threshold. Default is 0.2}

\item{numProbes}{Number of Probes}
}
\value{
Number of Break points for each unique sample
}
\description{
The Break Point calculation takes all the CNV data for a single patient and first filters it for segmentation mean of > 0.2 and, if specified, the minimum number of probes
covering that area. Then it counts the number of rows of data and multiplies it by 2. This represents the break points at the 5' and 3' ends of each segment.
\deqn{
Number\ of \ Break\ Points = \sum^{R}_{i = 1} (n_i \cdot 2)\ where\ |\bar{y}_{S_i}| \ge 0.2
}
}
