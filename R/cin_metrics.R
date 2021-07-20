
#' Total Aberration Index
#'
#' Total Aberration Index calculation takes the sum of lengths of each segment
#' times its segmentation mean for each sample and divides it by the sum of the
#' lengths of each sample.
#'
#' The Total Aberration Index (TAI) \href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3553118/}{(Baumbusch LO, et. al.)} is ``a measure of the abundance of genomic size of copy number changes in a tumour".
#' It is defined as a weighted sum of the segment means
#' \deqn{
#' Total\ Aberration\ Index =
#'  \frac
#' {\sum^{R}_{i = 1} {d_i} \cdot |{\bar{y}_{S_i}}|}
#' {\sum^{R}_{i = 1} {d_i}}\ \
#' where |\bar{y}_{S_i}| \ge |\log_2 1.7|
#' }
#'
#' @param cnvData dataframe containing following columns: Sample, Start, End, Num_Probes, Segment_Mean
#' @param segmentMean numerical value for the minimum segment_mean cutoff/ threshold. Default is 0.2
#' @param numProbes Number of Probes
#' @return Average of lengths weighted by segmentation mean for each unique sample
#' @examples tai(cnvData = maskCNV_BRCA)
#' @export
tai <- function(cnvData, segmentMean = 0.2, numProbes = NA ){
  unique_id <- unique(cnvData$Sample)
  tai.output <- stats::setNames(data.frame(matrix(ncol = 2, nrow = length(unique_id)), stringsAsFactors = FALSE),c("sample_id","tai"))
  for (i in 1:length(unique_id)){
    id <- unique_id[i]
    subsetSample <- subset(cnvData, cnvData$Sample == id)
    subsetSample <- subset(subsetSample, abs(subsetSample$Segment_Mean) >= segmentMean)
    if (!is.na(numProbes)){
      subsetSample <- subset(subsetSample, abs(subsetSample$Num_Probes) >= numProbes)
    }
    Length <- subsetSample$End - subsetSample$Start
    num <- Length*abs(subsetSample$Segment_Mean)
    #num <- Length*abs(subsetSample$Segment_Mean)
    den <- Length
    tai.calc <- sum(sum(num)/sum(den))
    tai.output$tai[i] <- tai.calc
    tai.output$sample_id[i] <- id
  }
  return(tai.output)
}

#' Modified Total Aberration Index
#'
#' Modified Total Aberration Index calculation takes the sum of lengths of each segment
#' times its segmentation mean for each sample and divides it by the sum of the
#' lengths of each sample.
#'
#' Modified Total Aberration Index uses all sample values instead of those in aberrant copy number state, thus does not remove the directionality from the score.
#' \deqn{
#' Modified\ Total\ Aberration\ Index =
#' \frac
#' {\sum^{R}_{i = 1} {d_i} \cdot {\bar{y}_{S_i}}}
#' {\sum^{R}_{i = 1} {d_i}}
#' }
#'
#' @seealso \code{\link{tai}}
#' @param cnvData dataframe containing following columns: Sample, Start, End, Num_Probes, Segment_Mean
#' @param segmentMean numerical value for the minimum segment_mean cutoff/ threshold. Default is 0.2
#' @param numProbes Number of Probes
#' @return Average of lengths weighted by segmentation mean for each unique sample
#' @examples taiModified(cnvData = maskCNV_BRCA)
#' @export
taiModified <- function(cnvData, segmentMean = 0, numProbes = NA ){
  unique_id <- unique(cnvData$Sample)
  tai.output <- stats::setNames(data.frame(matrix(ncol = 2, nrow = length(unique_id)), stringsAsFactors = FALSE),c("sample_id","modified_tai"))
  for (i in 1:length(unique_id)){
    id <- unique_id[i]
    subsetSample <- subset(cnvData, cnvData$Sample == id )
    subsetSample <- subset(subsetSample, abs(subsetSample$Segment_Mean) >= segmentMean)
    if (!is.na(numProbes)){
      subsetSample<-subset(subsetSample, abs(subsetSample$Num_Probes) >= numProbes)
    }
    Length <- subsetSample$End - subsetSample$Start
    num <- Length*subsetSample$Segment_Mean
    den <- Length
    tai.calc <- sum(sum(num)/sum(den))
    tai.output$modified_tai[i] <- tai.calc
    tai.output$sample_id[i] <- id
  }
  return(tai.output)
}

#' Copy Number Aberration
#'
#' Calculates the number of copy number aberrations
#'
#' Copy Number Aberrations (CNA) \href{https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0079079}{(Davidson JM, et al)}, are defined as a segment with copy number outside the pre-defined range of 1.7-2.3 \deqn{(\log_2 1.7 -1) \le \bar{y}_{S_i} \le (\log_2 2.3 -1)} that is not contiguous with an adjacent independent CNA of identical copy number. For our purposes, we have adapted the range to be \deqn{|\bar{y}_{S_i}| \ge |\log_2 1.7|}, which is only slightly larger than the original.
#' It is nearly identical to countingBreakPoints, except this one calculates breaks as adjacent segments that have a difference in segment means of \eqn{\ge 0.2}.
#' \deqn{Total\ Copy\ Number\ Aberration = \sum^{R}_{i = 1} n_i \ where \
#' \bar{y}_{S_i}| \ge |\log_2{1.7}|, \
#' \bar{y}_{S_{i-1}} - \bar{y}_{S_i}| \ge 0.2, \
#' d_i \ge 10}
#'
#' @seealso \code{\link{countingBreakPoints}}
#' @param cnvData dataframe containing following columns: Sample, Start, End, Num_Probes, Segment_Mean
#' @param segmentMean numerical value for the minimum segment_mean cutoff/ threshold. Default is 0.2
#' @param numProbes Number of Probes
#' @param segmentDistance Segment distance threshold
#' @param minSegSize Minimum segment size
#' @return Number of copy number aberrations between segments
#' @examples cna(cnvData = maskCNV_BRCA)
#' @export
cna <- function(cnvData, segmentMean = (log(1.7,2)-1), numProbes = NA, segmentDistance = 0.2, minSegSize = 10){
  unique_id <- unique(cnvData$Sample)
  cna.output <- stats::setNames(data.frame(matrix(ncol = 2, nrow = length(unique_id)), stringsAsFactors = FALSE),c("sample_id","cna"))
  for (i in 1:length(unique_id)){
    id <- unique_id[i]
    exSample <- subset(cnvData, cnvData$Sample == id )
    exSample <- subset(exSample, abs(exSample$Segment_Mean) >= abs(segmentMean) & !is.na(exSample$Segment_Mean))
    exSample <- subset(exSample, abs(exSample$End - exSample$Start) >= minSegSize)
    if (!is.na(numProbes)){
      exSample<-subset(exSample, abs(exSample$Num_Probes) >= numProbes)
    }
    breakpointNumber <- 0
    curSig <- NA
    for (segment in 1:nrow(exSample)){
      if (is.na(curSig)){
        curSig <- exSample$Segment_Mean[segment]
      }
      else {
        if (!is.na(curSig) & abs(curSig - exSample$Segment_Mean[segment]) >= segmentDistance){
          breakpointNumber = breakpointNumber + 1
          curSig <- exSample$Segment_Mean[segment]
        }
      }
    }
    cna.output$cna[i] <- breakpointNumber
    cna.output$sample_id[i] <- id
  }
  return(cna.output)
}


#' countingBaseSegments
#'
#' Function for counting altered base segments
#'
#' The Altered Base Segment calculation takes all the CNV data for a single patient and first filters it for a segmentation mean of > 0.2 and, if specified, the minimum number of probes
#' covering that area. Then, it calculates the sums of the lengths of each segment for a particular patient and outputs that.
#' \deqn{
#' Number\ of\ Altered\ Bases = \sum^{R}_{i = 1} d_i\ where\ |\bar{y}_{S_i}| \ge 0.2
#' }
#'
#' @param cnvData dataframe containing following columns: Sample, Start, End, Num_Probes, Segment_Mean
#' @param segmentMean numerical value for the minimum segment_mean cutoff/ threshold. Default is 0.2
#' @param numProbes Number of Probes
#' @return Number of Base segments for each unique sample
#' @examples countingBaseSegments(cnvData = maskCNV_BRCA)
#' @export
countingBaseSegments <- function(cnvData, segmentMean = 0.2, numProbes = NA ){
  unique_id <- unique(cnvData$Sample)
  NumBases <- stats::setNames(data.frame(matrix(ncol = 2, nrow = length(unique_id)), stringsAsFactors = FALSE),c("sample_id","base_segments"))
  for (i in 1:length(unique_id)){
    id <- unique_id[i]
    exSample<-subset(cnvData, cnvData$Sample == id )
    exSample<-subset(exSample, abs(exSample$Segment_Mean) >= segmentMean)
    if (!is.na(numProbes)){
      exSample<-subset(exSample, abs(exSample$Num_Probes) >= numProbes)
    }
    segSizes<-exSample$End - exSample$Start
    count <- sum(segSizes)
    NumBases$base_segments[i] <- count
    NumBases$sample_id[i] <- id
  }
  return(NumBases)
}


#' countingBreakPoints
#'
#' The Break Point calculation takes all the CNV data for a single patient and first filters it for segmentation mean of > 0.2 and, if specified, the minimum number of probes
#' covering that area. Then it counts the number of rows of data and multiplies it by 2. This represents the break points at the 5' and 3' ends of each segment.
#' \deqn{
#' Number\ of \ Break\ Points = \sum^{R}_{i = 1} (n_i \cdot 2)\ where\ |\bar{y}_{S_i}| \ge 0.2
#' }
#'
#' @param cnvData dataframe containing following columns: Sample, Start, End, Num_Probes, Segment_Mean
#' @param segmentMean numerical value for the minimum segment_mean cutoff/ threshold. Default is 0.2
#' @param numProbes Number of Probes
#' @return Number of Break points for each unique sample
#' @example countingBreakPoints(cnvData = maskCNV_BRCA)
#' @export
countingBreakPoints <- function(cnvData, segmentMean = 0.2, numProbes = NA){
  unique_id <- unique(cnvData$Sample)
  NumBpt <- stats::setNames(data.frame(matrix(ncol = 2, nrow = length(unique_id)), stringsAsFactors = FALSE),c("sample_id","break_points"))
  for (i in 1:length(unique_id)){
    id <- unique_id[i]
    exSample<-subset(cnvData, cnvData$Sample == id )
    exSample<-subset(exSample, abs(exSample$Segment_Mean) >= segmentMean)
    if (!is.na(numProbes)){
      exSample<-subset(exSample, abs(exSample$Num_Probes) >= numProbes)
    }
    number <- nrow(exSample) * 2 # Multiplied by 2 because of counting 5' and 3' breakpoints
    NumBpt$break_points[i] <- number
    NumBpt$sample_id[i] <- id
  }
  return(NumBpt)
}


#' Fraction Genome Altered
#'
#' Fraction Genome Altered looks at the fraction of the genome that deviates from a diploid state
#' fga calculates the fraction of the genome altered (FGA; [Chin SF, et. al.](https://www.ncbi.nlm.nih.gov/pubmed/17925008)), measured by taking the sum of the number of bases altered and dividing it by the genome length covered ($G$). Genome length covered was calculated by summing the lengths of each probe on the Affeymetrix 6.0 array. This calculation **excludes** sex chromosomes.
#' \deqn{
#' Fraction\ Genome\ Altered =
#' \frac
#' {\sum^{R}_{i = 1} d_i}
#' {G}
#' \ \ where\  |\bar{y}_{S_i}| \ge 0.2
#' }
#'
#' @param cnvData dataframe containing following columns: Sample, Start, End, Num_Probes, Segment_Mean
#' @param segmentMean numerical value for the minimum segment_mean cutoff/ threshold. Default is 0.2
#' @param numProbes Number of Probes
#' @param genomeSize Size of the genome derived from Affymetrix 6.0 array probe. Default is 2873203431 calculated based on hg38 **excluding sex chromosomes**
#' @return Fraction of the genome altered
#' @examples fga(cnvData = maskCNV_BRCA)
#' @export
fga <- function(cnvData, segmentMean = 0.2, numProbes = NA, genomeSize = 2873203431){
  # genomeSize derived from Affymetrix 6.0 array probe information. The default values is 2873203431 based on hg38
  unique_id <- unique(cnvData$Sample)
  fgaOutput <- stats::setNames(data.frame(matrix(ncol = 2, nrow = length(unique_id)), stringsAsFactors = FALSE),c("sample_id","fga"))
  for (i in 1:length(unique_id)){
    id <- unique_id[i]
    subsetSample <- subset(cnvData, cnvData$Sample == id )
    subsetSample <- subset(subsetSample, subsetSample$Chromosome %in% c(1:22))
    subsetSample <- subset(subsetSample, abs(subsetSample$Segment_Mean) >= segmentMean)
    if (!is.na(numProbes)){
      subsetSample<-subset(subsetSample, abs(subsetSample$Num_Probes) >= numProbes)
    }
    Length <- subsetSample$End - subsetSample$Start
    Sum <- sum(Length)
    fga.calc <- Sum/genomeSize
    fgaOutput$fga[i] <- fga.calc
    fgaOutput$sample_id[i] <- id
  }
  return(fgaOutput)
}


#' CINmetrics
#'
#' Calculate all CINmetrics on a given dataframe
#'
#' @param cnvData dataframe containing following columns: Sample, Start, End, Num_Probes, Segment_Mean
#' @param segmentMean_tai numerical value for the minimum segment_mean cutoff/ threshold for Total Aberration Index calculation. Default is 0.2
#' @param segmentMean_cna numerical value for the minimum segment_mean cutoff/ threshold for Copy Number Aberration calculation. Default is 0.2
#' @param segmentMean_base_segments numerical value for the minimum segment_mean cutoff/ threshold for Base segments calculation. Default is 0.2
#' @param segmentMean_break_points numerical value for the minimum segment_mean cutoff/ threshold for Break points calculation. Default is 0.2
#' @param segmentMean_fga numerical value for the minimum segment_mean cutoff/ threshold for Fraction of genome altered calculation. Default is 0.2
#' @param numProbes Number of Probes
#' @param segmentDistance_cna Segment distance threshold
#' @param minSegSize_cna Minimum segment size
#' @param genomeSize_fga Size of the genome derived from Affymetrix 6.0 array probe. Default is 2873203431 calculated based on hg38 **excluding sex chromosomes**
#' @return All Chromosomal INstability metrics
#' @examples CINmetrics(cnvData = maskCNV_BRCA)
#' @export
CINmetrics <- function(cnvData, segmentMean_tai = 0.2, segmentMean_cna = (log(1.7,2)-1), segmentMean_base_segments = 0.2, segmentMean_break_points = 0.2, segmentMean_fga = 0.2, numProbes = NA, segmentDistance_cna = 0.2, minSegSize_cna = 10, genomeSize_fga = 2873203431){
  # Calculate all Chromosomal Instability metrics as a single data frame
  unique_id <- unique(cnvData$Sample)
  cinmetrics <- data.frame(matrix(ncol = 6, nrow = length(unique_id)), stringsAsFactors = FALSE)
  tai <- CINmetrics::tai(cnvData = cnvData, segmentMean = segmentMean_tai, numProbes = numProbes)
  cna <- CINmetrics::cna(cnvData = cnvData, segmentMean = segmentMean_cna, numProbes = numProbes, segmentDistance = segmentDistance_cna, minSegSize = minSegSize_cna)
  base_segments <- CINmetrics::countingBaseSegments(cnvData = cnvData, segmentMean = segmentMean_base_segments, numProbes = numProbes)
  break_points <- CINmetrics::countingBreakPoints(cnvData = cnvData, segmentMean = segmentMean_break_points, numProbes = numProbes)
  fga <- CINmetrics::fga(cnvData = cnvData, segmentMean = segmentMean_fga, numProbes = numProbes, genomeSize = genomeSize_fga)
  cinmetrics <- merge(merge(merge(merge(tai, cna, by="sample_id"), base_segments, by="sample_id"), break_points, by="sample_id"), fga, by="sample_id")
  return(cinmetrics)
}
