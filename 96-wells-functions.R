### * Functions

### ** checkWellNames(wellNames)

checkWellNames = function(wellNames) {
    #' Check and convert well names to the appropriate format
    #' E.g. A1 -> A01
    #'
    #' @param wellNames String vector (e.g. c("A1", "A2", "B3", "B4"))
    #'
    #' @return String vector
    #'
    o = wellNames
    rows = substr(wellNames, 1, 1)
    stopifnot(all(rows %in% toupper(letters)[1:8]))
    columns = as.integer(substr(wellNames, 2, 10))
    stopifnot(all(columns >= 1 & columns <= 12))
    columns = as.character(columns)
    columns = sapply(columns, function(x) {
        if (nchar(x) == 1) {
            return(paste0("0", x))
        } else {
            return(x)
        }
    })
    return(paste0(rows, columns))
}

### ** drawCircle(x, y, radius, ...)

drawCircle = function(x, y, radius, ...) {
    nSteps = 24
    xc = x
    yc = y
    angle = seq(0, 2*pi, length.out = nSteps + 1)
    angle = angle[1:(length(angle) - 1)]
    x = xc + cos(angle) * radius
    y = yc + sin(angle) * radius
    polygon(x, y, ...)
}

### ** displayWells(x, labels, wells = "well", col = NULL)

displayWells = function(x, labels = NULL, wells = "well", col = NULL) {
    #' Display some labels in the 96-well plate format
    #'
    #' @param x Data frame with wells and labels in columns
    #' @param labels String, name of the column with the labels to display
    #' @param wells String, name of the column with the well names
    #' @param col String, name of the column with the color values
    #'
    par(mar = rep(0, 4))
    plot(0, type = "n", xlim = c(0, 13), ylim = c(0, 9), axes = F,
         xlab = F, ylab = F, asp = 1)
    # Plot row and column names
    rows = toupper(letters)[1:8]
    columns = c(paste0("0", 1:9), as.character(10:12))
    text(0.5, y = 7:0 + 0.5, labels = rows)
    text(1:12 + 0.5, y = 8.5, labels = columns)
    # Draw border
    rect(1, 0, 13, 8)
    # Draw circles
    for (i in 1:8) {
        for (j in 1:12) {
            xc = j + 0.5
            yc = 8 - i + 0.5
            drawCircle(xc, yc, radius = 0.45)
        }
    }
    # Draw labels
    if (is.null(labels)) {
        x$temp.plot.labels = rownames(x)
        labels = "temp.plot.labels"
    }
    if (!is.null(labels)) {
        labels = x[[labels]]
        if (is.numeric(labels[1])) {
            labels = round(labels, 3)
        }
        labels = as.character(labels)
        wells = as.character(x[[wells]])
        rowsX = substr(wells, 1, 1)
        columnsX = substr(wells, 2, 3)
        for (i in 1:length(labels)) {
            rowI = which(rows == rowsX[i])
            columnI = which(columns == columnsX[i])
            if (!is.null(col)) {
                drawCircle(columnI + 0.5, 8 - rowI + 0.5, radius = 0.45,
                           col = x[[col]][i])
            }
            text(columnI + 0.5, 8 - rowI + 0.5, labels = labels[i],
                 col = "black")
        }
    }
}

### *** Test

a = data.frame(well = c("A01", "A02", "C10"),
               sample = c(1, 2, 3))
displayWells(a, "sample")

### ** processPerkinElmerTable(data)

processPerkinElmerTable = function(data) {
    #' Convert the raw data frame from Perkin Elmer to light table
    #'
    absColumns = grepl("Absorbance", names(data))
    a = data[, c("Well", names(data)[which(absColumns)])]
    wl = substr(names(a), 22, 24)[2:ncol(a)]
    names(a) = c("well", paste0("A", wl))
    return(a)
}

### ** loadPerkinElmerFile(data)

loadPerkinElmerFile = function(filename) {
    #' Load the data from a Perkin Elmer text file. The text file will be
    #' loaded until the first line with "End time"
    #'
    a = readLines(filename)
    maxRow = min(grep("End time", a)) - 2
    a = read.table(filename, header = T, sep = "\t", dec = ",",
                   nrows = maxRow)
    return(processPerkinElmerTable(a))
}

### ** mergeRepeats(data)

mergeRepeats = function(data) {
    #' Merge the measurement repeats from a Perkin Elmer table
    #'
    #' @param data Output from loadPerkinElmerFile
    #'
    #' @return List
    #'
    #' @examples
    #' a = loadPerkinElmerFile("2106-05-08-15-51-data.txt")
    #' b = mergeRepeats(a)
    #'
    o = aggregate(. ~ well, data = data, FUN = mean)
    oVar = aggregate(. ~ well, data = data, FUN = var)
    for (i in 2:ncol(data)) {
        data[, i] = as.integer(!is.na(data[, i]))
    }
    oN = aggregate(. ~ well, data = data, FUN = sum)
    return(list(mean = o,
                var = oVar,
                n = oN))
}

### ** rectWells(topleft, bottomright)

rectWells = function(topleft, bottomright, transpose = F) {
    #' Generate a vector with well names within a rectangular region of a plate
    #'
    #' @examples
    #' rectWells("A01", "C06")
    #'
    rows = toupper(letters)[1:8]
    columns = c(paste0("0", 1:9), as.character(10:12))
    minRow = which(rows == substr(topleft, 1, 1))
    maxRow = which(rows == substr(bottomright, 1, 1))
    minCol = which(columns == substr(topleft, 2, 3))
    maxCol = which(columns == substr(bottomright, 2, 3))
    o = vector()
    i = 1
    for (k in minRow:maxRow) {
        for (l in minCol:maxCol) {
            o[i] = paste0(rows[k], columns[l])
            i = i + 1
        }
    }
    if (transpose) {
        o = vector()
        i = 1
        for (k in minCol:maxCol) {
            for (l in minRow:maxRow) {
                o[i] = paste0(rows[l], columns[k])
                i = i + 1
            }
        }
    }
    return(o)
}

### ** plateCoord(data)

plateCoord = function(data) {
    #' Add a colPos and rowCol columns to a data frame with a well column
    #'
    data$colPos = as.integer(substr(data$well, 2, 3))
    data$rowPos = sapply(substr(data$well, 1, 1),
                         function(x) which(toupper(letters)[1:8] == x))
    data$colPos = data$colPos - min(data$colPos)
    data$rowPos = data$rowPos - min(data$rowPos)
    return(data)
}

### ** distribSamples(samples, replicates = 3, wells = NULL, blockLength = NULL)

distribSamples = function(samples, replicates = 3, wells = NULL, blockLength = NULL) {
    #' Distribute replicated samples randomly across a 96-well plate
    #'
    #' @param samples Vector of string containing the sample names
    #' @param replicates Integer, number of replicate to generate for each
    #'   sample name
    #' @param wells Optional, vector containing the wells to be used, for
    #'   example generated using rectWells(). If NULL, uses the whole plate.
    #'   If blockLength is used, this has no effect.
    #' @param blockLength Optional, can be either 2, 4 or 8. If given, tries 
    #'   to optimize distribution for use with multichannel pipette using
    #'   sample blocks of length blockLength
    #'
    #' @return Data frame with the plate plan
    minGrey = 0.9
    maxGrey = 0.6
    stripsInfo = NULL
    if (is.null(blockLength)) {
        if (is.null(wells)) {
            wells = rectWells("A01", "H12")
        }
        samples = rep(samples, each = replicates)
        stopifnot(length(samples) <= length(wells))
        out = data.frame(sample = sample(samples),
                         well = sample(wells, size = length(samples)))
    } else {
        stopifnot(blockLength %in% c(2, 4, 8))
        if (length(samples) %% blockLength != 0) {
            emptySamples = blockLength - length(samples) %% blockLength
        } else {
            emptySamples = 0
        }
        samples = sample(c(samples, rep(NA, emptySamples)))
        naIndices = which(is.na(samples))
        if (length(naIndices) > 0) {
            for (i in 1:length(naIndices)) {
                samples[naIndices[i]] = paste0("NA_sample-", i)
            }
        }
        samplesStrips = split(samples,
                              rep(1:96, each = blockLength)[1:length(samples)])
        samplesStripsRef = samplesStrips
        #greyLevels = as.list(sapply(seq(minGrey, maxGrey, length.out = length(samplesStripsRef)),grey))
        greyLevels = as.list(adjustcolor(RColorBrewer::brewer.pal(length(samplesStripsRef), "Set2"), alpha.f = 0.5)[1:length(samplesStripsRef)])
        stripsInfo = list(strips = samplesStripsRef,
                          colors = greyLevels)
        for (i in 1:length(greyLevels)) {
            greyLevels[[i]] = rep(greyLevels[[i]], blockLength)
        }
        samplesStrips = rep(samplesStrips, replicates)
        greyLevels = rep(greyLevels, replicates)
        for (i in 1:length(samplesStrips)) {
            if (sample(c(1, 2), 1) == 2) {
                samplesStrips[[i]] = rev(samplesStrips[[i]])
            }
        }
        samplesStripsOrdering = sample(1:length(samplesStrips))
        samplesStrips = samplesStrips[samplesStripsOrdering]
        greyLevels = greyLevels[samplesStripsOrdering]
        samples = unlist(samplesStrips)
        stopifnot(length(samples) <= 96)
        out = data.frame(samples,
                         well = rectWells("A01", "H12", transpose = T)[1:length(samples)],
                         col = unlist(greyLevels),
                         stringsAsFactors = F)
    }
    class(stripsInfo) = c("strips", "list")
    outSamples = out$samples
    outSamples[grepl("NA_sample-", outSamples)] = NA
    out$sampleId = as.numeric(factor(outSamples))
    class(out) = c("platePlan", "data.frame")
    out = out[order(out$well), ]
    o = list(plan = out, strips = stripsInfo)
    class(o) = "plateInfo"
    return(o)
}

### ** plot.plateInfo(plateInfo)

plot.plateInfo = function(plateInfo) {
    nRows = ceiling(length(plateInfo[["strips"]][["strips"]])/3)
    layout(matrix(c(1, 1, 1, 2:(nRows*3 + 1)), ncol = 3, byrow = T),
           heights = c(0.5, rep(0.5/nRows, nRows)))
    platePlan = plateInfo[["plan"]]
    if ("col" %in% names(platePlan)) {
        displayWells(platePlan, labels = "sampleId", wells = "well", col = "col")
        plotStrips(plateInfo)
    } else {
        displayWells(platePlan, labels = "sampleId", wells = "well")
    }
}

### ** plotStrips(plateInfo)

plotStrips = function(plateInfo) {
    strips = plateInfo[["strips"]][["strips"]]
    nStrips = length(strips)
    s2id = unique(plateInfo[["plan"]][, c("samples", "sampleId")])
    s2col = unique(plateInfo[["plan"]][, c("samples", "col")])
    for (i in 1:length(strips)) {
        strip = strips[[i]]
        par(mar = c(0, 0, 0, 0))
        plot(0, type = "n", bty = "n", axes = F, xlab = "", ylab = "",
             xlim = c(0, length(strip) + 1),
             ylim = c(0, length(strip) + 1), asp = 1)
        for (j in 1:length(strip)) {
            drawCircle(1, j, radius = 0.40, col = s2col$col[s2col$samples == strip[j]])
            text(1, j, s2id$sampleId[s2id$samples == strip[j]], cex = 1)
            text(2, j, strip[j], pos = 4)
        }
    }
}

### ** savePlan(plateInfo, filename)

savePlan = function(plateInfo, filename) {
    pdf(filename, width = 8.27, height = 11.69, pointsize = 24)
    plot.plateInfo(plateInfo)
    dev.off()
}

### *** Test

samples = c("c0", "c1", "c2", "c3", "c4", "c5", "c6", "w", "DNase+", "DNase-", "EGTA-", "incub90", "incub75", "incub60", "dil2", "dil4", "dil8", "dil16", "dil32")
d = distribSamples(samples, replicates = 3, blockLength = 4)
savePlan(d, "toto.pdf")

### ** analyzePlate(filename, wells, samples, ref)

analyzePlate = function(filename, wells, samples, ref) {
    #' Analyze the results of a Perkin Elmer plate data
    #'
    #' @param filename String, name of the file with the data
    #' @param wells Vector with well names
    #' @param samples Vector with sample names
    #' @param ref Reference sample
    #'
    #' @examples
    #' a = buildPlate("2106-05-08-15-51.txt",
    #'                wells = rectWells("F05", "H12"),
    #'                samples = c(1, 7, 4, NA, 3, 7, 6, 3,
    #'                            2, 5, 1, 2, 7, 5, NA, 6,
    #'                            4, 3, 5, 6, 2, 4, NA, 1),
    #'                ref = 1)
    #'
    data = loadPerkinElmerFile(filename)
    samples = as.factor(samples)
    samples = relevel(samples, ref = ref)
    plate = data.frame(well = wells, sample = samples)
    b = mergeRepeats(data)$mean
    b = na.omit(merge(plate, b))
    b = plateCoord(b)
    m = lm(A570 ~ sample + colPos + rowPos, data = b)
    o = list(plate = b,
             results = m)
    class(o) = "plateData"
    return(o)
}

### ** summary.plateData

summary.plateData = function(x) {
    summary(x$results)
}

### ** plot.plateData

plot.plateData = function(x) {
    par(mfrow = c(2, 1))
    displayWells(x$plate, "sample")
    displayWells(x$plate, "A570")
}
