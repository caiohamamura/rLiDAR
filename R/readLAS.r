#'Reading LiDAR data 
#'
#'@description This function reads and returns values associated with the LAS file format. The LAS file is a public file format for the interchange of LiDAR 3-dimensional point cloud data (American Society of Photogrammetry and Remote Sensing - ASPRS)
#'
#'@usage readLAS(LASfile, short=TRUE)
#'
#'@param LASfile A standard LAS data file (ASPRS) 
#'@param short Logical, if TRUE it will return only a 5-column matrix with information on the returned point x, y, z locations, point intensity and the number of return within an individual discrete-return system laser pulse.
#'@return Returns a matrix listing the information stored in the LAS file.
#'@author Michael Sumner and Carlos Alberto Silva. 
#'@examples
#'
#'#=======================================================================#
#'# Importing LAS file:
#'#=======================================================================#
#'LASfile <- system.file("extdata", "LASexample1.las", package="rLiDAR")
#'
#'# Reading LAS file
#'rLAS<-readLAS(LASfile,short=TRUE)
#'
#'# Summary of the LAS file
#'summary(rLAS)
#'
#'#=======================================================================#
#'# LAS file visualization:
#'#=======================================================================#
#'
#'# 01 Set a single color 
#'
#'col<-"forestgreen"
#'
#'# plot 2D
#'plot(rLAS[,1],rLAS[,2], col=col,xlab="UTM.Easting", ylab="UTM.Northing", main="Single color")
#'
#'# plot 3D
#'library(rgl)
#'plot3d(rLAS[,1:3], col=col, axes=FALSE,xlab="", ylab="", zlab="")
#'axes3d(c("x+", "y-", "z-"))                 # axes
#'grid3d(side=c('x+', 'y-', 'z'), col="gray") # grid
#'title3d(xlab = "UTM.Easting", ylab = "UTM.Northing",zlab = "Height(m)", col="red") # title
#'planes3d(0, 0, -1, 0.001, col="gray", alpha=0.7)   # terrain
#'
#'
#'# 02 Set a color by height 
#'
#'# color ramp
#'myColorRamp <- function(colors, values) {
#'v <- (values - min(values))/diff(range(values))
#'x <- colorRamp(colors)(v)
#'rgb(x[,1], x[,2], x[,3], maxColorValue = 255)
#'}
#'
#'# Color by height
#'col <- myColorRamp(c("blue","green","yellow","red"),rLAS[,3])
#'
#'# plot 2D
#'plot(rLAS[,1], rLAS[,2], col=col, xlab="UTM.Easting", ylab="UTM.Northing", main="Color by height")
#'
#'# plot 3D
#'plot3d(rLAS[,1:3], col=col, axes=FALSE, xlab="", ylab="", zlab="")
#'axes3d(c("x+", "y-", "z-"))                     # axes
#'grid3d(side=c('x+', 'y-', 'z'), col="gray")     # grid
#'title3d(xlab = "UTM.Easting", ylab = "UTM.Northing",zlab = "Height(m)", col="red") # title
#'planes3d(0, 0, -1, 0.001, col="gray",alpha=0.7) # terrain
#'
#'@export
#'@importFrom bitops bitAnd bitShiftR
readLAS = function(LASfile, short = TRUE, fields) 
{
    if (class(short) != "logical") {
        stop("The short parameter is invalid. It must to be a TRUE or FALSE logical statement")
    }
    skip <- 0
    nrows <- NULL
    hd <- publicHeaderDescription()
    pheader <- vector("list", nrow(hd))
    names(pheader) <- hd$Item
    con <- file(LASfile, open = "rb")
    isLASFbytes <- readBin(con, "raw", size = 1, n = 4, endian = "little")
    pheader[[hd$Item[1]]] <- readBin(isLASFbytes, "character", 
        size = 4, endian = "little")
    if (!pheader[[hd$Item[1]]] == "LASF") {
        stop("The LASfile input is not a valid LAS file")
    }
    for (i in 2:nrow(hd)) {
        pheader[[hd$Item[i]]] <- readBin(con, what = hd$what[i], 
            size = hd$Rsize[i], endian = "little", n = hd$n[i])
    }
    close(con)
    `?`(readBin)
    numberPointRecords <- pheader[["Number of point records"]]
    offsetToPointData <- pheader[["Offset to point data"]]
    pointDataRecordLength <- pheader[["Point Data Record Length"]]
    xyzScaleOffset <- cbind(unlist(pheader[c("X scale factor", 
        "Y scale factor", "Z scale factor")]), unlist(pheader[c("X offset", 
        "Y offset", "Z offset")]))
    con <- file(LASfile, open = "rb")
    junk <- readBin(con, "raw", size = 1, n = offsetToPointData)
    if (skip > 0) {
        junk <- readBin(con, "raw", size = 1, n = pointDataRecordLength * 
            skip)
        numberPointRecords <- numberPointRecords - skip
    }
    if (!is.null(nrows)) {
        if (numberPointRecords > nrows) 
            numberPointRecords <- nrows
    }
    if (numberPointRecords < 1) 
        stop("no records left to read")
    allbytes <- matrix(readBin(con, "raw", n = pointDataRecordLength * 
        numberPointRecords, size = 1, endian = "little"), ncol = pointDataRecordLength, 
        nrow = numberPointRecords, byrow = TRUE)
    close(con)
    mm <- matrix(readBin(t(allbytes[, 1:(3 * 4)]), "integer", 
        size = 4, n = 3 * numberPointRecords, endian = "little"), 
        ncol = 3, byrow = TRUE)
    mm[, 1] <- mm[, 1] * xyzScaleOffset[1, 1] + xyzScaleOffset[1, 
        2]
    mm[, 2] <- mm[, 2] * xyzScaleOffset[2, 1] + xyzScaleOffset[2, 
        2]
    mm[, 3] <- mm[, 3] * xyzScaleOffset[3, 1] + xyzScaleOffset[3, 
        2]
    colnames(mm) <- c("X", "Y", "Z")
    Intensity <- readBin(t(allbytes[, 13:14]), "integer", size = 2, 
        n = numberPointRecords, signed = FALSE, endian = "little")
    bytesList <- readBin(t(allbytes[, 15]), "integer", size = 1, 
        n = numberPointRecords, signed = FALSE, endian = "little")
    Classification <- readBin(t(allbytes[, 16]), "integer", 
        size = 1, n = numberPointRecords, signed = FALSE, 
        endian = "little")
       
    return(cbind(mm, Classification))
}