
#' Total Aberration Index
#' Total Aberration Index calculation takes the sum of lengths of each segment
#' times its segmentation mean for each sample and divides it by the sum of the
#' lengths of each sample.
#' @param cnvData dataframe containing following columns: Sample, Start, End, Num_Probes, Segment_Mean
#' @param segmentMean numerical value for the minimum segment_mean cutoff/ threshold. Default is 0.2
#' @param numProbes Number of Probes
#' @return Average of lengths weighted by segmentation mean for each unique sample
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
#' Modified Total Aberration Index calculation takes the sum of lengths of each segment
#' times its segmentation mean for each sample and divides it by the sum of the
#' lengths of each sample.
#' @param cnvData dataframe containing following columns: Sample, Start, End, Num_Probes, Segment_Mean
#' @param segmentMean numerical value for the minimum segment_mean cutoff/ threshold. Default is 0.2
#' @param numProbes Number of Probes
#' @return Average of lengths weighted by segmentation mean for each unique sample
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
#' Calculates the number of copy number aberrations
#' Nearly identical to countingBreakPoints, except this one calculates breaks as adjacent segments that have a difference in segment means of >= 0.2.
#' @param cnvData dataframe containing following columns: Sample, Start, End, Num_Probes, Segment_Mean
#' @param segmentMean numerical value for the minimum segment_mean cutoff/ threshold. Default is 0.2
#' @param numProbes Number of Probes
#' @param segmentDistance Segment distance threshold
#' @param minSegSize Minimum segment size
#' @return Number of copy number aberrations between segments
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
#' Function for counting altered base segments
#' The Altered Base Segment calculation takes all the CNV data for a single patient and first filters it for a segmentation mean of > 0.2 and, if specified, the minimum number of probes
#' covering that area. Then, it calculates the sums of the lengths of each segment for a particular patient and outputs that.
#' @param cnvData dataframe containing following columns: Sample, Start, End, Num_Probes, Segment_Mean
#' @param segmentMean numerical value for the minimum segment_mean cutoff/ threshold. Default is 0.2
#' @param numProbes Number of Probes
#' @return Number of Base segments for each unique sample
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
#' The Break Point calculation takes all the CNV data for a single patient and first filters it for segmentation mean of > 0.2 and, if specified, the minimum number of probes
#' covering that area. Then it counts the number of rows of data and multiplies it by 2. This represents the break points at the 5' and 3' ends of each segment.
#' @param cnvData dataframe containing following columns: Sample, Start, End, Num_Probes, Segment_Mean
#' @param segmentMean numerical value for the minimum segment_mean cutoff/ threshold. Default is 0.2
#' @param numProbes Number of Probes
#' @return Number of Break points for each unique sample
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
#' Fraction Genome Altered looks at the fraction of the genome that deviates from a diploid state
#' @param cnvData dataframe containing following columns: Sample, Start, End, Num_Probes, Segment_Mean
#' @param segmentMean numerical value for the minimum segment_mean cutoff/ threshold. Default is 0.2
#' @param numProbes Number of Probes
#' @param genomeSize Size of the genome derived from Affymetrix 6.0 array probe. Default is 2873203431 calculated based on hg38 **excluding sex chromosomes**
#' @return Fraction of the genome altered
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
