% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/intersection_of_intervals.R
\name{intersectionBedList}
\alias{intersectionBedList}
\title{Intersection of bed intervals list}
\usage{
intersectionBedList(bed_list_1, bed_list_2, verbose = FALSE)
}
\arguments{
\item{bed_list_1}{dataframe, bed format dataframe noting one type of the interval mutation}

\item{bed_list_2}{dataframe, bed format dataframe noting one type of the interval mutation}

\item{verbose}{bool, print verbose information?}
}
\value{
bed_df
}
\description{
Intersection of bed intervals list
}
\examples{

bed_df_1 <- data.frame(Chr = c(paste0("chr", c(1,1,1))), Start=c(100, 200, 300), End=c(150, 250, 350))
bed_df_2 <- data.frame(Chr = c(paste0("chr", c(1,1,1))), Start=c(120, 220, 320), End=c(150, 250, 350))
bed_list_1 <- list("sample1" = bed_df_1, "sample2" = bed_df_2)
bed_list_2 <- list("sample1" = bed_df_2, "sample2" = bed_df_1)

intersectionBedList(bed_list_1, bed_list_2, verbose = FALSE)

}
