#!/usr/bin/Rscript

##############################################################################
# Standalone script for
##############################################################################
# Change to 1 if running as standalone process during
# generation of validation data

args <- commandArgs(TRUE)
if (length(args)==0) {
	write <- FALSE
} else {
	# only write out if argument specified
	write <- as.logical(args[1])
	parentdir <- as.character(args[2])

	# also only reloads packages and sets parentdir if running via Rscript
	suppressMessages(require(dplyr))
	suppressMessages(require(openxlsx))
}

parentdir

# Read DNMs
gonl_dnms <- read.table(paste0(parentdir,
		"/reference_data/DNMs/GoNL_DNMs.txt"),
	header=T, stringsAsFactors=F)
gonl_dnms <- gonl_dnms[,1:5]
names(gonl_dnms) <- c("ID", "CHR", "POS", "REF", "ALT")

itmi_dnms <- read.xlsx(paste0(parentdir,
	"/reference_data/DNMs/goldmann_2016_dnms.xlsx"), sheet=1)
itmi_dnms$ID <- "goldmann"
itmi_dnms <- itmi_dnms[,c(7,1,2,4,5)]
names(itmi_dnms) <- c("ID", "CHR", "POS", "REF", "ALT")
itmi_dnms$CHR <- gsub("chr", "", itmi_dnms$CHR)

dnms_full <- rbind(gonl_dnms, itmi_dnms)
dnms_full$CAT <- paste(dnms_full$REF, dnms_full$ALT, sep="")

dnms_full$Category[dnms_full$CAT=="AC" | dnms_full$CAT=="TG"] <- "AT_CG"
dnms_full$Category[dnms_full$CAT=="AG" | dnms_full$CAT=="TC"] <- "AT_GC"
dnms_full$Category[dnms_full$CAT=="AT" | dnms_full$CAT=="TA"] <- "AT_TA"
dnms_full$Category[dnms_full$CAT=="GA" | dnms_full$CAT=="CT"] <- "GC_AT"
dnms_full$Category[dnms_full$CAT=="GC" | dnms_full$CAT=="CG"] <- "GC_CG"
dnms_full$Category[dnms_full$CAT=="GT" | dnms_full$CAT=="CA"] <- "GC_TA"

dnms_full <- dnms_full %>%
  dplyr::select(ID, CHR, POS, Category) %>%
  arrange(CHR, POS)

if(write){
  categs <- c("AT_CG", "AT_GC", "AT_TA", "GC_AT", "GC_CG", "GC_TA")
  for(i in 1:6){
    cat <- categs[i]
    outdat <- dnms_full %>%
      filter(Category==cat)

    outfile <- paste0(parentdir, "/reference_data/DNMs/GoNL_", cat, ".txt")
    write.table(outdat, outfile, col.names=F, row.names=F, sep="\t", quote=F)
  }
}
