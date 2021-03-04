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

	SummExpt <- mae[[ app.type ]]
	colData(SummExpt) <- colData(mae)

			## -------------------------------------------------------------------- ##
	
	if (! "ensembl_gene_id" %in% colnames(rowData(SummExpt))) {
		if ( all(grepl("^ENS[A-Z]*G\\d+$", rowData(SummExpt)$geneID, perl=T)) ) {
			rowData(SummExpt)[[ "ensembl_gene_id" ]] = as.character( rowData(SummExpt)$geneID )
		}
		else { stop("***INCORRECTLY** FORMATTED EnsEMBL identifiers") }

		if (any(duplicated(rowData(SummExpt)$ensembl_gene_id))) {
			stop("**DUPLICATED EnsEMBL identifiers***")
		}
	}

	if (! "external_gene_name" %in% colnames(rowData(SummExpt))) {
		rowData(SummExpt)[[ "external_gene_name" ]] = as.character( rowData(SummExpt)$geneName )

		## Genes w/o a symbol, assign the EnsEMBL identifier ##
		Missing_Symbol.idx <- ! grepl("^\\s*\\S.*$", rowData(SummExpt)$external_gene_name, perl=T)
		if (sum(Missing_Symbol.idx)>0) {
			rowData(SummExpt)$external_gene_name[ Missing_Symbol.idx ] = rowData(SummExpt)$ensembl_gene_id[ Missing_Symbol.idx ]
		}
		## Duplicated gene symbols ##
		if (any(duplicated(rowData(SummExpt)$external_gene_name))) {
			number.Dups <- length(unique(
				rowData(SummExpt)$external_gene_name[ duplicated(rowData(SummExpt)$external_gene_name) ]
			))
			cat("\nNo. *DUPLICATED* gene symbols : ",number.Dups,"\n\n", sep="")
			rowData(SummExpt)$external_gene_name = make.unique(
				rowData(SummExpt)$external_gene_name, sep="__"
			)
		}
	}

	if ("description" %in% colnames(rowData(SummExpt))) {
		## Tidy up the functional descriptions ##
		if ( any(grepl("^\\S.*\\S\\s+\\[\\S.*\\S\\]\\s*$",rowData(SummExpt)$description, perl=T)) ) {
			rowData(SummExpt)$description = sub(
				"^(\\S.*\\S)\\s+\\[\\S.*\\S\\]\\s*$", "\\1", rowData(SummExpt)$description, perl=TRUE
			)
		}
	}

			## -------------------------------------------------------------------- ##

	species <- if ("species" %in% names(metadata(mae))) {
		metadata(mae)$species
	} else { stop("***NO SPECIES NAME IN MultiAssayExperiment Object***") }

	expt.type <- if ("expt.type" %in% names(metadata(mae))) {
		metadata(mae)$expt.type
	} else { stop("***NO EXPERIMENT TYPE IN MultiAssayExperiment Object***") }

	summ.expt.lst <- list(
		App.Type=app.type, Expt.Type=expt.type, Species=species, se=SummExpt
	)

	if (app.type=="DEA.Results") {

		bt.Extensions  <- c("all", "pcg", "reg_rna")
		Contr.RegEx <- paste0(
			"^(\\S.*\\S)\\.(",paste(bt.Extensions, collapse="|"),")\\.\\_\\.(\\S+)$"
		)
		dea.res.idx <- grep(Contr.RegEx, colnames(rowData(SummExpt)), perl=T)

		BioTypes <- unique(sub(Contr.RegEx, "\\2", colnames(rowData(SummExpt))[ dea.res.idx ], perl=T))
		Contrast.Names <- unique(sub(Contr.RegEx, "\\1", colnames(rowData(SummExpt))[ dea.res.idx ], perl=T))

		Statistics <- unique(sub(Contr.RegEx, "\\3", colnames(rowData(SummExpt))[ dea.res.idx ], perl=T))
		if ("fail.detect.cnt" %in% Statistics) Statistics <- Statistics[ Statistics!="fail.detect.cnt" ]

		summ.expt.lst <- append(summ.expt.lst,
			list(
				BioTypes       = BioTypes,
				Contrast.Names = Contrast.Names,
				Statistics     = Statistics,
				Contrasts      = metadata(mae)$contrasts,
				ESEL.Info      = metadata(mae)$esel.info
			)
		)
		gse.idx <- which(! names(experiments(mae)) %in% app.type)
		if (length(gse.idx)>0) {
			#summ.expt.lst <- append(summ.expt.lst, list(GSEA=experiments(mae)[gse.idx]))
			summ.expt.lst <- append(summ.expt.lst, list( GSEA=metadata(mae)$ShinyNGS.gsea) )
		}
	}
	return(summ.expt.lst)
}


#' Creates an object of the ExploratorySummarizedExperiment class
#' 
#' This is fixes a bug in \code{shinyngs::ExploratorySummarizedExperiment}; the container \code{colData} stores the phenotypic data and if there is only a single factor (i.e. one column), it gets converted to a vector & IS NOT preserved as a dataframe. Hence assignment of 'assays' fails and throws an error because the sample names have been lost! Setting "drop=FALSE" when manipulating this object prevents this behaviour.
#'
#' @param See ?shinyngs::ExploratorySummarizedExperiment
#'
#' @export
ExploratorySummarizedExperiment.SingleFact <- function(
	assays, colData, annotation, idfield,
	labelfield=character(), entrezgenefield=character(), contrast_stats=list(),
	assay_measures=list(), gene_set_analyses=list(), dexseq_results=list(), read_reports=list() )
{
	# https://github.com/pinin4fjords/shinyngs/blob/master/R/ExploratorySummarizedExperiment-class.R

	# Reset NULLs to empty
	if (is.null(entrezgenefield)) entrezgenefield <- character()
	
	# The assays slot of a summarised experiment needs the same dimensions for every matrix
	all_rows <- Reduce(union, lapply(assays, rownames))
	
	add_missing_rows <- function(x) {
	  missing_rows <- all_rows[! all_rows %in% rownames(x)]
	  empty_rows   <- data.frame(matrix(NA, nrow=length(missing_rows), ncol=ncol(x)), row.names=missing_rows)
	  colnames(empty_rows) <- colnames(x)
	  rbind(x, empty_rows)[all_rows, , drop = FALSE]
	}

	# Subset colData to remove any samples not present in the first assay
	
								  ### ***BUG FIX*** ###
	#  If drop isn't set to FALSE, then colData with a single factor i.e. one column    #
	#  gets converted to a vector & IS NOT preserved as a dataframe. Hence assignment   #
	#  of 'assays' fails and throws an error because the sample names have been lost!!! #
	colData <- colData[rownames(colData) %in% colnames(assays[[1]]),,drop=FALSE]
	
	assays <- SimpleList(
	  lapply(assays, function(as) { round(add_missing_rows(as)[,rownames(colData)], 2) })
	)  

	# The same fix for contrast_stats
	if (length(contrast_stats) > 0)
	{
	  contrast_stats <- lapply(contrast_stats,
		function(stats) { lapply(stats, function(test) add_missing_rows(test) ) }
	  )
	}
	
	# Annotations need to be strings
	annotation <- data.frame(
	  lapply(annotation, as.character), stringsAsFactors=FALSE, check.names=FALSE, row.names=rownames(annotation)
	)[all_rows, ]
	
	# Build the object
	sumexp <- SummarizedExperiment(assays = assays, colData = DataFrame(colData))
	mcols(sumexp) <- annotation
	
	new(
	  "ExploratorySummarizedExperiment", sumexp, idfield = idfield, labelfield = labelfield, entrezgenefield = entrezgenefield, assay_measures = assay_measures, 
	  contrast_stats = contrast_stats, gene_set_analyses = gene_set_analyses, dexseq_results = dexseq_results, read_reports = read_reports
	)
} 

#' @export
Excluded_Samples_Check <- function(Samples, Contr.lst) {
	DEA_BioTypes     <- NULL
	excluded.samples <- NULL

	for (l in 1:length(Contr.lst)) {
		dea_Design.lst <- Contr.lst[l]

		## Backwards compatibility of biotypes not previously included in the DEA list ##
		Gene_Bt <- ifelse("DEA.Gene.type" %in% names(dea_Design.lst), dea_Design.lst$DEA.Gene.type, "All Gene Biotypes")
		DEA_BioTypes <- c(DEA_BioTypes, Gene_Bt)

		## Need to account for EXCLUDED samples --> included in DEA object ##
		if ("Excluded.Samples" %in% names(dea_Design.lst)) {
			excluded.samples <- c(excluded.samples, dea_Design.lst$Excluded.Samples)
		}
	}

	if (length(unique(DEA_BioTypes))!=1) {
		cat("\n\n"); stop("There is ***NO UNIQUE DESCRIPTION*** of the biotype of DEA tested genes.")
	}

	if (length(excluded.samples)>0) {
		excluded.samples <- unique(excluded.samples)
		cat("Samples *EXCLUDED* from DEA [",length(excluded.samples),"] : '",
			paste(sort(excluded.samples), collapse="', '"),"'\n\n", sep="")
		Samples <- Samples[ ! Samples %in% excluded.samples ]
	}
	else { cat("NO Samples excluded from DEA.\n\n") }

	return( list(samples = Samples, biotype = unique(DEA_BioTypes)) )
}

## --------------------------------------------------------------------------------------------------------- ##

#' @export
Construct_Exploratory.Ranged.SummExpt <- function(
	Summ.Expt, Contr.Metadata, Biotype, Contrasts, DEA.Stats)
{
	if (! "counts" %in% names(assays(Summ.Expt))) {
		cat("\n"); stop("Assay 'counts' **MUST BE INCLUDED** in assay list")
	} else {
		Sample.IDs <- colnames(Summ.Expt)
		cat("\nSamples in dataset analysed [",length(Sample.IDs),"] :-\n", sep="")
		print(DataFrame(Sample.IDs)); cat("\n")

		Samples.lst <- Excluded_Samples_Check(Sample.IDs, Contr.Metadata[ paste0(Contrasts,".",Biotype) ])
	}


	DEA_rds.lst <- lapply(Contrasts,
		function(c) {
			Results.Names <- paste0(c,".",Biotype,"._.",DEA.Stats)

			dea.res.idx <- which(colnames(rowData(Summ.Expt)) %in% Results.Names)
			dea.res.df  <- rowData(Summ.Expt)[, dea.res.idx]
			colnames(dea.res.df) <- DEA.Stats

			Tested.Idx <- apply(dea.res.df, 1, function(row) all(is.na(row)))
			dea.res.df[! Tested.Idx, ]
		}
	)
	names(DEA_rds.lst) <- Contrasts
	cat("\n\"DifferentialExpression\" list objects per contrast :-\n")
	str(DEA_rds.lst); cat("\n")

	Tested.Genes <- Reduce(union, lapply(DEA_rds.lst, function(x) rownames(x)))
	cat("\nUnion of genes tested in the differential expression analysis :-\n")
	print(DataFrame(Tested.Genes)); cat("\n")

	raw.counts <- round(
		as.matrix(assays(Summ.Expt)$counts[Tested.Genes, Samples.lst$samples])
	)
	mode(raw.counts) <- "integer"

	Assays.lst <- list()
	Assays.lst[[ "counts" ]] = raw.counts

	if ( any(grepl("^VST", names(assays(Summ.Expt)), perl=T, ignore.case=T)) ) {

		vst.Assays.idx <- grep("^VST", names(assays(Summ.Expt)), perl=T, ignore.case=T)
		for (i in vst.Assays.idx) {
			vst.Assay <- names(assays(Summ.Expt))[i]
			Assays.lst[[ vst.Assay ]] = assays(Summ.Expt)[[ vst.Assay ]][ Tested.Genes, Samples.lst$samples ]
		}
	}

	if ("tpm" %in% names(assays(Summ.Expt))) {
		Assays.lst[[ "tpm" ]] = assays(Summ.Expt)[[ "tpm" ]][ Tested.Genes, Samples.lst$samples ]
	}

	AssayNames <- names(Assays.lst)
	vst.idx <- grep("^VST", AssayNames, ignore.case=T)
	Sorted_AssayNames <- c( sort(AssayNames[vst.idx]), sort(AssayNames[-vst.idx]) )

	Feature.Regex <- paste0("\\.\\_\\.(",paste(c(DEA.Stats,"fail.detect.cnt"), collapse="|"),")$")
	Feature.Idx   <- grep(Feature.Regex, colnames(rowData(Summ.Expt)), perl=TRUE)
	Gene.Features <- rowData(Summ.Expt)[ Tested.Genes, -Feature.Idx ] 


	cat("<< Creating 'ExploratorySummarizedExperiment' object >>\n")
	ese <- ExploratorySummarizedExperiment.SingleFact(
		assays     = SimpleList(Assays.lst[ Sorted_AssayNames ]),
		colData    = colData(Summ.Expt)[ Samples.lst$samples, ],
		annotation = Gene.Features,
		idfield = "ensembl_gene_id", labelfield = "external_gene_name"
	)
	cat("DONE.\n")

	return( list(ese=ese, DEA.Genes=Tested.Genes, dea.results=DEA_rds.lst, Biotype=Samples.lst$biotype) )
}


Species_URL <- function(MetaData) {

	if (! "Species" %in% names(MetaData)) {
		cat("\n"); stop("***NO VALID SPECIES NAME IN METADATA*** - e.g. 'Homo sapiens'")
	}

	if ( grepl("^\\s*\\S+\\s+\\S+\\s*$", MetaData$Species, perl=T) ) {
		organism <- sub("^\\s*(\\S{1})\\S+\\s+(\\S+)\\s*$", "\\L\\1\\2", MetaData$Species, perl=T)
		cat("\nEnsEMBL database name : \"",organism,"\" [Species Latin name='",MetaData$Species,"']", sep="")
	}
	else {
		cat("\nSpecies : '",MetaData$Species,"'\n", sep="")
		stop("***SPECIES NAME MUST BE GIVEN AS FULL BINOMINAL LATIN NAME*** - e.g. 'Homo sapiens'")
	}

	SpeciesNames.lst <- list(
		"mmusculus"="Mus_musculus", "hsapiens"="Homo_sapiens", "rnorvegicus"="Rattus_norvegicus"
	)
	if (organism %in% names(SpeciesNames.lst)) {
		cat(" --> URL name : \"",SpeciesNames.lst[[ organism ]],"\"\n", sep="")
		url <- paste0(
			"http://www.ensembl.org/",SpeciesNames.lst[[ organism ]],"/Gene/Summary?db=core;g="
		)
	}
	else {
		cat("\n",paste0("Species supported by EnsEMBL : '",paste(sort(unlist(SpeciesNames.lst)), collapse="', '"),"'\n"),
			sep="")
		stop("***NO VALID SPECIES URL NAME FOUND*** - e.g. 'Homo_sapiens'")
	}
	return(url)
}


Check.Analysis_Variables <- function(Contrast.MetaData, Phenotypes) {
	
	Designs <- lapply(Contrast.MetaData,
		function(x) x[ names(x)=="Design.Formula" ]
	)
	Design.Vars <- sort(unique(
		unlist(lapply(Designs,
			function(design) {
				design <- sub("^\\~\\s*(\\S.*\\S)$", "\\1", design, perl=TRUE)
				Variables <- strsplit(design, "\\s*\\+\\s*", fixed=FALSE, perl=TRUE)
				return(Variables)
			}
		))
	))
	if (any(Design.Vars=="0")) Design.Vars = Design.Vars[ Design.Vars!="0" ]

	if (all(Design.Vars %in% Phenotypes)) {
		cat("\nAll variables present in the experimental designs are also in the manifest:\n'",
			paste(sort(Design.Vars), collapse="', '"),"'\n", sep="")
	} else {
		cat("\nVariables used in the experimental designs absent from the manifest:\n'",
			paste(sort(Design.Vars[! Design.Vars %in% Phenotypes]), collapse="', '"),"'\n", sep="")
		stop("***MISSING DESIGN VARIABLES***")
	}

	design.types <- unlist(
		lapply(Contrast.MetaData, function(x) x[ names(x)=="Design.Type" ])
	)

	atlas <- if (any(design.types=="Atlas")) TRUE else FALSE
	return(atlas)
}


#' Creates a ExploratorySummarizedExperiment list object
#' 
#' Contain one or more ExploratorysummarizedExperiments with the same sets of samples/columns but different feature sets. This facilitates the examination of expression at, for example, both transcript and gene levels or for different gene biotypes in RNA-seq experiments explorted via ‘Shinyngs’
#'
#' @param MAE_Object List containing - minimally - a SummarizedExperiment object, species data derived from and type of transcriptomics analysis
#'
#' @export
Make_ESE_list <- function(MAE_Object) {

	Summ.Expt <- MAE_Object$se

	if (MAE_Object$App.Type=="DEA.Results") {
		ESE_Objs.lst    <- list()
		DEA_Results.lst <- list()

		for (gene.biotype in dea.res.lst$BioTypes) {
			DEA_ESE.lst <- Construct_Exploratory.Ranged.SummExpt(
				Summ.Expt, MAE_Object$Contrasts, gene.biotype, MAE_Object$Contrast.Names, MAE_Object$Statistics
			)

			biotype_desc <- DEA_ESE.lst$Biotype
			DEA_Results.lst[[ gene.biotype ]][[ "desc" ]]        = biotype_desc
			DEA_Results.lst[[ gene.biotype ]][[ "dea.results" ]] = list( DEA_ESE.lst$dea.results )
			DEA_Results.lst[[ gene.biotype ]][[ "DEA.Genes" ]]   = DEA_ESE.lst$DEA.Genes
			ESE_Objs.lst[[ biotype_desc ]] = DEA_ESE.lst$ese
		}
		myesel <- ExploratorySummarizedExperimentList( eses=ESE_Objs.lst )
	}
	else {
		myesel <- ExploratorySummarizedExperimentList(
	        eses = list( "Exploratory Data Analysis" = Read_Exploratory.SummExpt(Summ.Expt) )
	    )
	}
	myesel <- Add_App.Info(MAE_Object$ESEL.Info, myesel)

	myesel@url_roots$ensembl_gene_id = Species_URL(
		MAE_Object[ ! names(MAE_Object) %in% c("se","Contrast.Names","Statistics","Contrasts","GSEA") ]
	)

	if (MAE_Object$App.Type=="DEA.Results") {

		esel.dat.lst <- list(
			esel          = myesel,
			dea.stats     = DEA_Results.lst,
			atlas.designs = Check.Analysis_Variables(
				MAE_Object$Contrasts, colnames(colData(Summ.Expt))
			)
		)

		esel.stats.lst <- Merge_DEA_Statistics(
			esel.dat.lst, MAE_Object$Expt.Type, MAE_Object$Contrasts
		)

		if ("GSEA" %in% names(MAE_Object)) {
			myesel <- Add_GSEA_Results(MAE_Object, esel.stats.lst, DEA_Results.lst)
		}
	}
	else { return(myesel) }
}

## --------------------------------------------------------------------------------------------------------- ##

Add_App.Info <- function(ESEL.Info, esel.obj) {

	Info.Names <- names(ESEL.Info)[ ! names(ESEL.Info) %in% "atlas_levs" ]
	for (Attr in Info.Names) slot(esel.obj, Attr) = ESEL.Info[[ Attr ]]
	
	return(esel.obj)
}


# Coefficients [unlogged fold-changes] calculated from #
# Negative Binominal model (MLE's) or by limma.        #
Model_Coefficients <- function(DEA_Stats.lst, Tested_Gene.IDs) {
	
	cat("<<\"fold_changes\">>\n")

	fcs.mat <- do.call("cbind",
		lapply(names(DEA_Stats.lst),
			function(Contr) {
				fcs <- rep(NA, length(Tested_Gene.IDs))
				names(fcs) <- Tested_Gene.IDs

				LFCs <- DEA_Stats.lst[[ Contr ]]$log2FoldChange
				LFCs[ is.na(LFCs) ] <- 0

				fcs[ rownames(DEA_Stats.lst[[ Contr ]]) ] <- 2^LFCs

				TestedGene.cnt <- sum(!is.na(fcs))
				cat("Contrast '",Contr,"' : No. Genes tested=",
					formatC(TestedGene.cnt, format="d", big.mark=",")," [\"log2FoldChange\"]\n", sep="")
				if (TestedGene.cnt!=length(LFCs)) {
					stop("No. genes tested **DOES NOT EQUAL** number in SummarizedExperiment")
				}
				return(fcs)
			}
		)
	)
	colnames(fcs.mat) <- names(DEA_Stats.lst)
	
	print(DataFrame(fcs.mat, check.names=F)); cat("\n\n")
	return(fcs.mat)
}

## --------------------------------------------------------------------------------------------------------- ##

# Hypothesis tests raw probabilities or adjusted p-values #
Test_Probabilities <- function(DEA_Stats.lst, Tested_Gene.IDs, Prob.Type) {

	pval.desc <- ifelse(Prob.Type=="prob", "P-value", "Q-value/FDR")
	cat("<<\"",pval.desc,"s\">>\n", sep="")

	probs.mat <- do.call("cbind",
		lapply(names(DEA_Stats.lst),
			function(Contr) {
				probs <- rep(NA, length(Tested_Gene.IDs))
				names(probs) <- Tested_Gene.IDs

				P.vals <- DEA_Stats.lst[[ Contr ]][[ Prob.Type ]]
				if ( any(is.na(P.vals)) ) {
					cat("Assigning ",pval.desc," \"NAs\" to 1 [",
						sum(is.na(P.vals))," genes] : Contrast '",Contr,"'\n", sep="")
					P.vals[ is.na(P.vals) ] <- 1
				}
				probs[ rownames(DEA_Stats.lst[[ Contr ]]) ] <- P.vals

				TestedGene.cnt <- sum(!is.na(probs))
				cat("No. Genes tested=",formatC(TestedGene.cnt, format="d", big.mark=","),
					" [\"",Prob.Type,"\"]\n", sep="")
				if (TestedGene.cnt!=length(P.vals)) {
					stop("No. genes tested **DOES NOT EQUAL** number in SummarizedExperiment")
				}
				return(probs)
			}
		)
	)
	colnames(probs.mat) <- names(DEA_Stats.lst)

	print(DataFrame(probs.mat, check.names=F)); cat("\n\n")
	return(probs.mat)
}

## --------------------------------------------------------------------------------------------------------- ##

Check_TestLevs_Present <- function(Contrasts.lst, ComparisonName, SummExpt_Levels, Test_Grp) {

	Levels <- Contrasts.lst[[ ComparisonName ]][[ paste0(Test_Grp,"Levels") ]]
	SummExpt_Levels <- unique(as.character(SummExpt_Levels))

	if (! all(Levels %in% SummExpt_Levels)) {
		lev.idx      <- which( Levels %in% SummExpt_Levels )
		levs_Present <- sort(Levels[ lev.idx ])
		levs_Absent  <- sort(Levels[ -lev.idx ])

		cat("\nAbsent _",toupper(Test_Grp),"_ levels : \"",paste(levs_Absent, collapse="\", "),"\" / [Present : \"",
			paste(levs_Present, collapse="\", "),"\"]\nLevels for \"",Contrasts.lst[[ ComparisonName ]]$Factor,
			"\" in phenotypes/colData() :-\n'",
			paste(sort(SummExpt_Levels), collapse="'\n'"),"'\n\n", sep="")
		stop("***MISSING FACTOR LEVEL(S)***")
	}
	else {
		cat("\nNo missing _",toupper(Test_Grp),"_ levels [Contrast='",ComparisonName,"' / Factor='",
			Contrasts.lst[[ ComparisonName ]]$Factor,"'] : \"",
			paste(sort(Levels), collapse="\", \""),"\".", sep="")
		return(Levels)
	}
}

## --------------------------------------------------------------------------------------------------------- ##

DEA.Contrast.Stats <- function(DEA_rds.lst, Gene.IDs, trxnExpt.type) {
	

	DEA_Results.Vars <- c("log2FoldChange", "prob", "adj.pval")
	DEA_Results.Vars <- c(
		DEA_Results.Vars, ifelse(trxnExpt.type %in% c("Affymetrix","NanoString"), "lfcSE", "shrunken.log2FC")
	)

	cat("\n\"Contrast Stats\" list *str*ucture :-\n")
	dea.Results.lst <- lapply(DEA_rds.lst, function(x) x[, DEA_Results.Vars] )
	str(dea.Results.lst)

	Genes.Identical <- unlist( lapply(dea.Results.lst, function(x) all(rownames(x) %in% Gene.IDs)) )
	cat("\nChecking all the genes tested for each contrast are in the SummarizedExperiment object :-\n")
	print(DataFrame(Genes.Identical)); cat("\n\n")

	if (all(Genes.Identical)) {
		stats.shiny.lst <- list()

		stats.shiny.lst[[ "fold_changes" ]] = Model_Coefficients(dea.Results.lst, Gene.IDs)
		stats.shiny.lst[[ "pvals" ]]        = Test_Probabilities(dea.Results.lst, Gene.IDs, "prob")
		stats.shiny.lst[[ "qvals" ]]        = Test_Probabilities(dea.Results.lst, Gene.IDs, "adj.pval")
	}
	else {
		cat("\n!!CHECK!! Contrasts : '",
			paste(sort(names(Genes.Identical)[ ! Genes.Identical ]),collapse="', '"),"'\n", sep="")
		stop("Tested genes are **ABSENT** from the SummarizedExperiment object")
	}

	return( list(contrast_stats=stats.shiny.lst) )
}

## --------------------------------------------------------------------------------------------------------- ##

DEA.Contrasts <- function(MetaData, BioType, Comparisons, Atlas.Groupings, Phenotypic.Data) {

	cat("\nSample Sheet :-\n")
	print(Phenotypic.Data)

	comp.idx   <- grep(paste0(BioType,"$"), names(MetaData), perl=T)
	comp.names <- names(MetaData)[ comp.idx ]

	cntrs.lst <- lapply(comp.names, function(x) {
		contr <- MetaData[x][[1]]
		contr[ c("Factor","NumeratorLevels","DenominatorLevels") ]
	})
	names(cntrs.lst) <- sub(paste0("\\.",BioType,"$"), "", comp.names, perl=T)
	cat("\nNo. Contrasts : ",length(cntrs.lst),"\n", sep="")
	str(cntrs.lst)

	Contr.lst <- list()
	Contrasts <- lapply(Comparisons,
		function(contr) {
			if (! contr %in% names(cntrs.lst)) stop("***NO ASSOCIATED FACTOR and/or LEVELS for contrast***")

			contr.regex <- paste0("^(\\S+)\\_{2}(\\S+)\\_vs\\_(\\S+)$")
			Factor      <- sub(contr.regex, "\\1", contr, perl=T)
			numerator   <- sub(contr.regex, "\\2", contr, perl=T)
			denominator <- sub(contr.regex, "\\3", contr, perl=T)

			if (denominator=="Continuous") {
				##
				cleaned.fact <- make.names(paste0(DEA_rds.lst[[ contr ]]$Factor,".ContDE"))

				if (! cleaned.fact %in% colnames(Phenotypic.Data)) {
					cat("\nFACTOR ## \"",cleaned.fact,"\" ## __IS NOT PRESENT__ in the phenotypes/colData() :-\n'",
						paste(sort(colnames(Phenotypic.Data)), collapse="'\n'"),"'\n\n", sep="")
					stop("***CONTINUOUS VARIABLE QUARTILE IS MISSING FROM THE SAMPLE SHEET***")
				}
				c(cleaned.fact, "Q1", "Q4")
			}
			else if (Factor==denominator) {
				num.lev     <- DEA_rds.lst[[ contr ]]$NumeratorLevels
				atlasFactor <- Atlas.Groupings[ num.lev ]
				## unname() is _IMPORTANT_ - otherwise it breaks volcano plots etc. ##
				return( c(unname(atlasFactor), "Others", num.lev) )
			} else {
				
				Factor.Name <- make.names(cntrs.lst[[ contr ]]$Factor)

				if (! Factor.Name %in% colnames(Phenotypic.Data)) {
					cat("\nFACTOR ## \"",Factor.Name,"\" ## __IS NOT PRESENT__ in the phenotypes/colData() :-\n'",
						paste(sort(colnames(Phenotypic.Data)), collapse="'\n'"),"'\n\n", sep="")
					stop("***VARIABLE IS MISSING FROM THE SAMPLE SHEET***")
				}

				Num.Levs <- Check_TestLevs_Present(
					cntrs.lst, contr, Phenotypic.Data[[ Factor.Name ]], "Numerator"
				)
				Den.Levs <- Check_TestLevs_Present(
					cntrs.lst, contr, Phenotypic.Data[[ Factor.Name ]], "Denominator"
				)
				cat("\n")
				return( c(Factor.Name, Den.Levs, Num.Levs) )
			}
		}
	)
	names(Contrasts) <- Comparisons

	Contr.lst[[ "contrasts" ]] = Contrasts
	return(Contr.lst)
}

## --------------------------------------------------------------------------------------------------------- ##

Add_Contrast_Levels <- function(FactorName, numerator.levs, denominator.levs, Phenotypic.Data) {

    phenos <- colnames(Phenotypic.Data)
    fact.count  <- sum(grepl(paste0("^",FactorName), phenos, perl=T))
    constr.name <- make.names(paste0(FactorName,"_C",fact.count))

    num.strg <- paste(sort(numerator.levs), collapse="_")
    den.strg <- paste(sort(denominator.levs), collapse="_")

    constr.lst <- list( Contrast=c(constr.name, den.strg, num.strg) )

    constr.lst[[ "Levels" ]] = setNames(
        data.frame(
            unlist(lapply(Phenotypic.Data[[ FactorName ]],
                function(l) {
                    ifelse(l %in% numerator.levs, num.strg,
                        ifelse(l %in% denominator.levs, den.strg, "OMITTED")
                    )
                }
            )), row.names=rownames(Phenotypic.Data)),
        constr.name
    )
    return(constr.lst)
}

## --------------------------------------------------------------------------------------------------------- ##

AtlasDesign.Levels <- function(Atlas.GroupVar, Atlas.Levels, Phenotypes) {

	cat("\nLevels used in \"Atlas\" design for factor \"",Atlas.GroupVar,"\" :-\n",
		paste0("'",paste(sort(Atlas.Levels), collapse="', '"),"'"), "\n\n", sep="")

	Levels <- Phenotypes[[ Atlas.GroupVar ]]
	if (! is.character(Levels)) Levels <- as.character(Levels)
	names(Levels) <- rownames(Phenotypes)

	atlas.levs.df <- data.frame(
		do.call("cbind", lapply(Atlas.Levels,
			function(l) {
				ifelse(Levels %in% l, l,
					ifelse(Levels %in% Atlas.Levels[ Atlas.Levels!=l ], "Others", "OMITTED")
				)
			}
		))
	)

	Atlas_Grp.Names <- make.names(
		paste0(Atlas.GroupVar,"_Grp",1:length(Atlas.Levels))
	)
	names(Atlas_Grp.Names) <- Atlas.Levels
	dimnames(atlas.levs.df) <- list(names(Levels), Atlas_Grp.Names)

	cat("DataFrame of \"Atlas\" contrasts **REQUIRED** for the ShinyNGS app. :-\n", sep="")
	print(DataFrame(atlas.levs.df, check.names=F)); cat("\n")

	return( list(Groups=Atlas_Grp.Names, Levels=atlas.levs.df) )
}


#' Add DEA results ExploratorySummarizedExperiment list 
#'
#' @param ESEL_Object List containing - minimally - a SummarizedExperiment object, species data derived from and type of transcriptomics analysis
#'
#' @export
Merge_DEA_Statistics <- function(ESEL_Object, Expt_Type, Contrasts) {

	myesel <- ESEL_Object$esel

	Atlas_Grp.Names <- if (ESEL_Object$atlas.designs) {
		atlas.pheno.lst <- AtlasDesign.Levels(Default.GroupVar, atlas.levs, colData(mae))
		atlas.pheno.lst$Groups
	} else {
		cat("\nThis is *NOT* an \"Atlas\" design!\n\n", sep="")
		NA
	}

	Comparisons_lst   <- list()

	for (biotype in names(ESEL_Object$dea.stats)) {
		BioType_Desc <- ESEL_Object$dea.stats[[ biotype ]]$desc
		
		if (ESEL_Object$atlas.designs) {
			colData(myesel[[ BioType_Desc ]]) <- cbind(
				colData(myesel[[ BioType_Desc ]]),
				atlas.pheno.lst$Levels[ rownames( colData(myesel[[ BioType_Desc ]]) ), ]
			)
			cat("Modified manifest for '",BioType_Desc,"' :-\n", sep="")
			print( colData(myesel[[ BioType_Desc ]]) ); cat("\n")
		}

		DEA_Object <- DEA.Contrast.Stats(
			ESEL_Object$dea.stats[[ biotype ]]$dea.results[[1]],
			ESEL_Object$dea.stats[[ biotype ]]$DEA.Genes,
			Expt_Type
		)
		cat("\nAdding differential expression analysis results for gene-biotype \"",biotype,"\" :-\n", sep="")
		str(DEA_Object); cat("\n") 

		contrs <- lapply(DEA_Object$contrast_stats, function(x) colnames(x))

		chk.contrs <- unlist(
			lapply(2:length(contrs),
				function(i) all.equal(contrs[i], contrs[1], check.names=F)
			)
		)
		names(chk.contrs) <- names(contrs)[ 2:length(contrs) ]

		Checked.Contrasts <- if (all(chk.contrs)) {
			unlist(contrs[1])
		}
		else {
			mm.contrs <- sort(names(chk.contrs)[ ! chk.contrs ])
			cat("Mismatched contrasts for metric(s) : '",
				paste(mm.contrs, collapse="', '"),"'\n<<Contrasts List>>\n", sep="")
			print(contrs)
			stop("***COMPARISONS ARE NOT IDENTICAL FOR ALL DEA STATISTICS***")
		}

		## Deal with continuous variables in design ##
		Experimental_Variables <- colnames( colData(myesel[[ BioType_Desc ]]) )
		if (all(DEA_Object$group_vars %in% Experimental_Variables))
		{
			for (ASSAY in names(assays(myesel[[ BioType_Desc ]])) ) {
				cat("Adding 'contrast_stats' to \"",ASSAY,"\" ... ", sep="")
				myesel[[ BioType_Desc ]]@contrast_stats[[ ASSAY ]] = DEA_Object$contrast_stats
				cat("done\n")
			}
			cat("Finished adding results of DEA to the ExploratorySummarizedExperimentList\n\n", sep="")
		}
		else {
			Absent.Vars <- sort( DEA_Object$group_vars[ ! DEA_Object$group_vars %in% Experimental_Variables ] )
			cat("\nABSENT from sample sheet but PRESENT in design : \"",paste(Absent.Vars,collapse="\", \""),"\"\n", sep="")
			stop("**MISSING ANALYSIS VARIABLES**")
		}

		comparison.grps.lst <- DEA.Contrasts(
			Contrasts, biotype, Checked.Contrasts,
			Atlas_Grp.Names, colData(myesel[[ BioType_Desc ]])
		)
		Comparisons_lst <- append(Comparisons_lst, comparison.grps.lst$contrasts)
	}


	Contrasts.df <- do.call("rbind", lapply(Comparisons_lst, function(Grp) Grp))
	# Remove duplicated contrasts i.e. rows
	Contrasts.df <- Contrasts.df[! duplicated(Contrasts.df), , drop=FALSE]

	# Turn into a list
	Contrasts.lst <- lapply(rownames(Contrasts.df), function(contr) Contrasts.df[contr,])
	names(Contrasts.lst) <- rownames(Contrasts.df)

	## Setting this enables DEA & GSEA panels in app. ##
	myesel@contrasts <- Contrasts.lst
	cat("\nALL comparisons (contrasts) performed for DEA :-\n")
	print(myesel@contrasts)

	return( list(esel=myesel) )
}

## --------------------------------------------------------------------------------------------------------- ##

Merge.GeneSets <- function(Analysis.GeneSets) {

	if (length(Analysis.GeneSets)==1) {
		return( Analysis.GeneSets[[1]] )
	} else {
		cat("\n\nInternal *str*ucture of 'Gene_Sets' :-\n")
		GeneSet_Classes <- unique(unlist( lapply(Analysis.GeneSets, function(x) names(x)) ))

		Merged.GeneSets.lst <- list()
		# Restructure (reverse) the gene set to biotype mappings
		for (gs.Class in GeneSet_Classes) {
			print(gs.Class)

			for (biotype in names(Analysis.GeneSets)) {
				cat(biotype," : ",sep="")
				print(gs.Class %in% names(Analysis.GeneSets[[ biotype ]]))

				if ( gs.Class %in% names(Analysis.GeneSets[[ biotype ]]) ) {
					Merged.GeneSets.lst[[ gs.Class ]][[ biotype ]] = Analysis.GeneSets[[ biotype ]][[ gs.Class ]]
				}
			}
			cat("\n")
		}
		cat("\n<< Restructured/reversed GENE SET --> BIOTYPE mappings >>\n")
		print( lapply(Merged.GeneSets.lst, function(x) lapply(x, function(y) length(y))) )


		ESEL_GeneSets.lst <- lapply( names(Merged.GeneSets.lst),
			function(gs_Class) {
				gs_Class.terms.lst <- Merged.GeneSets.lst[[ gs_Class ]]

				if (length(gs_Class.terms.lst)>1) {
					terms.lst <- gs_Class.terms.lst[[1]]
					cat("\nNo. terms for gene set class '",gs_Class,"' [Gene Biotype : \"",
						names(gs_Class.terms.lst)[1],"\"] = ",length(terms.lst),"\n", sep="")

					for (l in 2:length(gs_Class.terms.lst))
					{
						Biotype  <- names(gs_Class.terms.lst)[l]
						Appd.IDx <- which(! names(gs_Class.terms.lst[[ l ]]) %in% names(terms.lst))
						cat("No. terms to append for gene biotype \"",Biotype,"\" : ",length(Appd.IDx),"\n", sep="")

						if (length(Appd.IDx)>0) terms.lst <- append(terms.lst, gs_Class.terms.lst[[ l ]][ Appd.IDx ])
					}
					return(terms.lst)
				}
				else {
					cat("\nNO terms to merge for gene set class '",gs_Class,"' & Gene Biotype \"",
						names(gs_Class.terms.lst)[1],"\" (",length(gs_Class.terms.lst[[1]]),")\n", sep="")
					return( gs_Class.terms.lst[[1]] )
				}
			}
		)
		names(ESEL_GeneSets.lst) <- names(Merged.GeneSets.lst)

		cat("\n\nNon-redundant/merged Gene Set sizes :-\n")
		print( lapply(ESEL_GeneSets.lst, function(x) length(x)) )
		return(ESEL_GeneSets.lst)
	}
}

## --------------------------------------------------------------------------------------------------------- ##

Add_GSEA_Results <- function(MAE.Object, ESEL.Data, DEA.Statistics) {

	myesel <- ESEL.Data$esel

	### !!! Important to ASSIGN correctly for platform i.e whether microarray or RNA-seq !!! ###
	Assay.Variable <- ifelse(
		MAE.Object$Expt.Type=="RNA-seq", "counts", "NormVals"
	)

	Gene_Sets <- list()
	for (biotype in names(DEA.Statistics)) {
		BioType_Desc <- DEA.Statistics[[ biotype ]]$desc

		myesel[[ BioType_Desc ]]@gene_set_analyses[[ Assay.Variable ]] = MAE.Object$GSEA[[ biotype ]]$stats
		Gene_Sets[[ BioType_Desc ]] = MAE.Object$GSEA[[ biotype ]]$gene.sets
	}

	myesel@gene_set_id_type = "external_gene_name"
	myesel@gene_sets[[ myesel@gene_set_id_type ]] = Merge.GeneSets(Gene_Sets)

	return(myesel)
}


