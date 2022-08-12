# https://rdrr.io/cran/TCGA2STAT/src/R/TCGA2STAT.R
TumorNormalMatch2 <- function(dat){
  temp <- dat
  
  if(is.null(temp)){
    message("Empy data object")
    return(NULL)
  }
  
  if(class(temp) == "data.frame" | class(temp) == "matrix"){
    temp.type <- sapply(colnames(temp), function(s) unlist(strsplit(s, "-"))[4]) # Split the sample name up by dashes. The sample type is the fourth number (e.g. "01A") in the sample name.
    primary.tumor <-  temp[, grep("^01", temp.type)] # The primary tumors are the ones whose type has "01" in it. 
    normal <-  temp[, c(grep("^10", temp.type))] # The blood-derived normals are the ones whose type has "10" in it. 
    # normal <-  temp[, c(grep("^10", temp.type), grep("^11", temp.type), grep("^12", temp.type))] # The blood-derived normals are the ones whose type has "10" in it. 
    
    tum.bcr <- substr(colnames(primary.tumor), 1, 12)
    norm.bcr <- substr(colnames(normal), 1, 12)
    
    matching <- intersect(tum.bcr, norm.bcr)
    
    if(length(matching) == 0){
      message("No matching samples on tumor and normal.")
      return(NULL)
    }
    
    if(length(matching) > 0){
      # Create a table that will match the truncated names to the original names. 3 columns: truncated name, matching tumor, matching normal. This part has to come before the next chunk of code, because the next chunk of code changes the colnames of the primary.tumor and normal data frames to their truncated versions.
      # First, put together all the truncated tumor names and all their original names (not just the samples with matching normals.)
      tumor.dict <- data.frame(BaseName=tum.bcr, FullName=colnames(primary.tumor), stringsAsFactors = F)
      # Next, put together all the truncated normal-sample names and all their original names (not just the sample with matching tumors.)
      normal.dict <- data.frame(BaseName=norm.bcr, FullName=colnames(normal), stringsAsFactors = F)
      # isdup <- function (x) duplicated (x) | duplicated (x, fromLast = TRUE)
      # View(normal.dict[isdup(normal.dict$BaseName),])
      # Some of the patients have matching normal from both blood and solid tissue. (See above line of code and http://www.omnesres.com/tools/tcga/ for the explanation of TCGA sample IDs.) We only want the blood normals. (https://www.biostars.org/p/284806/#285183)
      normal.dict <- normal.dict[grep("10A", normal.dict$FullName),,drop=FALSE]
      # Then subset tumor.dict and normal.dict to include only the samples with matched normal/tumor tissues.
      matching <- intersect(normal.dict$BaseName, tumor.dict$BaseName)
      tumor.dict <- tumor.dict[tumor.dict$BaseName %in% matching,,drop=FALSE]
      normal.dict <- normal.dict[normal.dict$BaseName %in% matching,,drop=FALSE]
      
      # Make sure everything is in the correct order.
      tumor.dict <- tumor.dict[order(tumor.dict$BaseName),,drop=FALSE]
      normal.dict <- normal.dict[order(normal.dict$BaseName),,drop=FALSE]
      # identical(tumor.dict$BaseName, normal.dict$BaseName)
      # identical(tumor.dict$BaseName, matching)
      # Then merge the two tables.
      tumor.normal.dict <- data.frame(BaseName=tumor.dict$BaseName, TumorFullName=tumor.dict$FullName, NormalFullName=normal.dict$FullName, stringsAsFactors = F)
      
      # Subset the primary.tumor and normal tables. 
      primary.tumor.match <- primary.tumor[,colnames(primary.tumor) %in% tumor.normal.dict$TumorFullName, drop=FALSE]
      normal.match <- normal[,colnames(normal) %in% tumor.normal.dict$NormalFullName, drop=FALSE]
      identical(colnames(primary.tumor.match), tumor.normal.dict$TumorFullName)
      identical(colnames(normal.match), tumor.normal.dict$NormalFullName)
      # Rename the columns. 
      colnames(primary.tumor.match) <- substr(colnames(primary.tumor.match), 1, 12)
      colnames(normal.match) <- substr(colnames(normal.match), 1, 12)
      # identical(substr(colnames(primary.tumor.match), 1, 12), substr(colnames(normal.match), 1, 12))
      # identical(colnames(primary.tumor.match), colnames(normal.match))
      
      return(list(primary.tumor=primary.tumor.match, normal=normal.match, dictionary=tumor.normal.dict))
    }
  }
}