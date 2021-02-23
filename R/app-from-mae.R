#' Load a MultiAssayExperiment object
#'
#' Reads a RDS file containing the MAE object, and returns a list containing a \code{SummarizedExperiment} object with the results from differential expression & geneset enrichment analyses if performed.
#'
#' @param rds.file String - location of RDS file of the MAE object
#'
#' @return A list with a \code{se} component (the \code{SummarizedExperiment} object), and optionally \code{BioTypes}, \code{Contrast.Names} and \code{Statistics} if DEA/GSEA results are present; these contain the biotypes of the genes analysed, the contrasts and the names of the summary statistics of the DEA and GSEA runs, respectively.
#'
#' @examples
#' \dontrun{
#' # MultiAssayExperiment with DEA & GSE data
#' rds.file1 <- "path_to_rds_file1/mae-dea.rds"
#' dea.lst <- Read.MAE_RDS(rds.file1)
#'
#' # Exploratory ShinyNGS app. - MultiAssayExperiment with only raw counts, VSTs etc.
#' rds.file2 <- "path_to_rds_file2/mae-explr.rds"
#' explr.lst <- Read.MAE_RDS(rds.file2)
#' }
#'
#' @export
Read.MAE_RDS <- function(rds.file) {

	if (file.exists(rds.file)) {
		cat("\nLoading _serialised_ 'MultiAssayExperiment' object from \"",
			rds.file,"\" ... ", sep="")
		mae <- readRDS(rds.file)
		app.type <- ifelse(
			"Exploratory.Data" %in% names(experiments(mae)), "Exploratory.Data",
			ifelse("DEA.Results" %in% names(experiments(mae)), "DEA.Results", NA)
		)
		cat("done\nType of ShinyNGS app. = '",app.type,"'\n", sep="")
		print(mae)
		
		if (is.na(app.type)) stop("***NO 'SummarizedExperiment' object found***")
	}
	else {
		cat("\nLocation of file containing 'MultiAssayExperiment' object : \"",rds.file,"\"\n", sep="")
		stop("RDS file ***DOES NOT EXIST***")
	}

	dea.summExpt <- mae[[ app.type ]]
	colData(dea.summExpt) <- colData(mae)

			## -------------------------------------------------------------------- ##
	
	if (! "ensembl_gene_id" %in% colnames(rowData(dea.summExpt))) {
		if ( all(grepl("^ENS[A-Z]*G\\d+$", rowData(dea.summExpt)$geneID, perl=T)) ) {
			rowData(dea.summExpt)[[ "ensembl_gene_id" ]] = as.character( rowData(dea.summExpt)$geneID )
		}
		else { stop("***INCORRECTLY** FORMATTED EnsEMBL identifiers") }

		if (any(duplicated(rowData(dea.summExpt)$ensembl_gene_id))) {
			stop("**DUPLICATED EnsEMBL identifiers***")
		}
	}

	if (! "external_gene_name" %in% colnames(rowData(dea.summExpt))) {
		rowData(dea.summExpt)[[ "external_gene_name" ]] = as.character( rowData(dea.summExpt)$geneName )

		## Genes w/o a symbol, assign the EnsEMBL identifier ##
		Missing_Symbol.idx <- ! grepl("^\\s*\\S.*$", rowData(dea.summExpt)$external_gene_name, perl=T)
		if (sum(Missing_Symbol.idx)>0) {
			rowData(dea.summExpt)$external_gene_name[ Missing_Symbol.idx ] = rowData(dea.summExpt)$ensembl_gene_id[ Missing_Symbol.idx ]
		}
		## Duplicated gene symbols ##
		if (any(duplicated(rowData(dea.summExpt)$external_gene_name))) {
			number.Dups <- length(unique(
				rowData(dea.summExpt)$external_gene_name[ duplicated(rowData(dea.summExpt)$external_gene_name) ]
			))
			cat("\nNo. *DUPLICATED* gene symbols : ",number.Dups,"\n", sep="")
			rowData(dea.summExpt)$external_gene_name = make.unique(
				rowData(dea.summExpt)$external_gene_name, sep="__"
			)
		}
	}

	if ("description" %in% colnames(rowData(dea.summExpt))) {
		## Tidy up the functional descriptions ##
		if ( any(grepl("^\\S.*\\S\\s+\\[\\S.*\\S\\]\\s*$",rowData(dea.summExpt)$description, perl=T)) ) {
			rowData(dea.summExpt)$description = sub(
				"^(\\S.*\\S)\\s+\\[\\S.*\\S\\]\\s*$", "\\1", rowData(dea.summExpt)$description, perl=TRUE
			)
		}
	}

			## -------------------------------------------------------------------- ##

	if (app.type=="Exploratory.Data") {
		return( list(se=dea.summExpt) )
	}
	else {
		bt.Extensions  <- c("all", "pcg", "reg_rna")
		Contr.RegEx <- paste0(
			"^(\\S.*\\S)\\.(",paste(bt.Extensions, collapse="|"),")\\.\\_\\.(\\S+)$"
		)
		dea.res.idx <- grep(Contr.RegEx, colnames(rowData(dea.summExpt)), perl=T)

		BioTypes <- unique(sub(Contr.RegEx, "\\2", colnames(rowData(dea.summExpt))[ dea.res.idx ], perl=T))
		Contrast.Names <- unique(sub(Contr.RegEx, "\\1", colnames(rowData(dea.summExpt))[ dea.res.idx ], perl=T))

		Statistics <- unique(sub(Contr.RegEx, "\\3", colnames(rowData(dea.summExpt))[ dea.res.idx ], perl=T))
		if ("fail.detect.cnt" %in% Statistics) Statistics <- Statistics[ Statistics!="fail.detect.cnt" ]

		return(
			list(se=dea.summExpt, BioTypes=BioTypes, Contrast.Names=Contrast.Names, Statistics=Statistics)
		)
	}
}







