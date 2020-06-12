

scatterInput <- function(id, eselist) {
    # Create a namespace function using the provided id
    ns <- NS(id)
    
    Sample.Sheet   <- SummarizedExperiment::colData(eselist[[1]], stringsAsFactors=FALSE, check.names=F)
    Sample.Classes <- unlist(lapply(Sample.Sheet, class))
    Sample.Sheet   <- Sample.Sheet[, Sample.Classes=="factor"]
    
    pheno.idx <- -grep("^(description|qualityFormat|sampleName|bcbio\\.Dirs|\\S+\\.Quantiles)$", colnames(Sample.Sheet), perl=T)
    Phenos <- colnames(Sample.Sheet)[ pheno.idx ]
    
    ColorBrewer_Categories <- c("Diverging","Qualitative")
    
    ColorBrewer.Palettes.lst <- lapply(ColorBrewer_Categories,
      function(x) {
        Data.Nature <- sub("^(\\S{3})\\S+$", "\\L\\1", x , perl=T)
        idx <- grep(paste0("^",Data.Nature), as.character(RColorBrewer::brewer.pal.info$category), perl=T)
        rownames(RColorBrewer::brewer.pal.info)[ idx ]
      }
    )
    names(ColorBrewer.Palettes.lst) <- ColorBrewer_Categories
    
    list(
      selectizeInput(
        ns("gene_symbol"), "", NULL, multiple = TRUE,
        options = list(
          placeholder  = "Please type a value or scroll",
          onInitialize = I('function() { this.setValue(""); }'),
          maxItems     = 5
        )
      ),
      
      selectizeInput(
        inputId=ns("Pheno"), label="Conditions", choices=Phenos,
        options = list(
          placeholder = "Please type a value or scroll", onInitialize = I('function() { this.setValue(""); }')
        )
      ),
      
      shinyWidgets::pickerInput(
        inputId = ns("levels"),
        label   = "Choose condition group(s):",
        choices = "",
        
        options = list(
          `actions-box`= TRUE,
           style = "btn-royal btn-sm", # ?actionBttn & https://github.com/dreamRs/shinyWidgets/issues/74
          `selected-text-format` = "count > 3"
        ),
        multiple = TRUE
      ),
      
      shinyWidgets::awesomeRadio( inputId = ns("AssayType"), choices = "", label = "Quantitation Type"),
      
      selectizeInput(
        inputId=ns("Pt.Colour"), label="Colour Palettes", choices=ColorBrewer.Palettes.lst, selected="Set1"
      )
    )
}


scatterOutput <- function(id) {
    ns <- NS(id)
    plotOutput(ns("bxPlot1"))
}


scatter <- function(input, output, session, eselist) {
    
    ## Select the genes - a maximum of five ##
    observe({ validate( need(input$gene_symbol, FALSE) ) })

    Gene.Annots <- SummarizedExperiment::rowData(eselist[[1]], stringsAsFactors=FALSE, check.names=F)
    gene.Count  <- formatC(nrow(Gene.Annots), format="d", big.mark=",")
    
    gn.symbol.col <- "external_gene_name"
    # See https://shiny.rstudio.com/articles/selectize.html
    updateSelectizeInput(
      session, "gene_symbol",
      label=paste0(prettifyVariablename(gn.symbol.col)," [",gene.Count," genes]"),
      choices=Gene.Annots[[ gn.symbol.col ]], server=TRUE
    )

    # Require the EnsEMBL ids to access the raw counts & TPMs
    Gene.Identifiers <- reactive({
      gn.idx <- which( Gene.Annots[[ gn.symbol.col ]] %in% input$gene_symbol )
      setNames(
        Gene.Annots[[ gn.symbol.col ]][ gn.idx ], Gene.Annots$ensembl_gene_id[ gn.idx ]
      )
    })

    
    ## Select the condition/tratment etc. of interest ##
    observe({ validate(need(input$Pheno, FALSE)) })

    Samples.df  <- SummarizedExperiment::colData(eselist[[1]], stringsAsFactors=FALSE, check.names=F)
    Col.Classes <- unlist(lapply(Samples.df, class))
    Samples.df  <- Samples.df[, Col.Classes=="factor"]
    
    pheno.dat <- do.call("cbind",
      lapply( Samples.df[, -grep("^(description|qualityFormat|sampleName|bcbio\\.Dirs|\\S+\\.Quantiles)$", colnames(Samples.df), perl=T)], as.character)
    )
    rownames(pheno.dat) <- rownames(Samples.df)
    
    observeEvent(input$Pheno, {
      levs <- if (! isTruthy(input$Pheno)) {
        reactiveVal("")
      }
      else { reactive({ names(table(pheno.dat[,input$Pheno])) }) }
      shinyWidgets::updatePickerInput(session, "levels", choices=levs() )
    }, ignoreInit = TRUE)
    # Will _NOT_ run when the ‘observeEvent’/‘eventReactive’ is created #
    # (because ‘ignoreInit = TRUE’), but it will run every other time.  #


    ## Select analysis group(s) - **NEED AT LEAST ONE** ##
    observe({ validate(need(input$levels, FALSE)) })

    assays.lst <- list(`Transcripts per million`="tpm", `Raw counts`="counts")
    shinyWidgets::updateAwesomeRadio(
        session, "AssayType", choices = assays.lst, selected = "tpm"
    )

    Assay_Type <- reactive({ input$AssayType })


    # https://code-examples.net/en/q/211f7c7
    # https://stackoverflow.com/questions/34731975/how-to-listen-for-more-than-one-event-expression-within-a-shiny-eventreactive-ha
    changeData <- reactive({ paste(input$gene_symbol, input$Pheno, input$levels, input$AssayType) })
    
    Plt.Vals.df <- eventReactive(changeData(), {
      req(
        input$gene_symbol, input$Pheno, input$levels,
        all( input$levels %in% unique(pheno.dat[, input$Pheno]) ) # Check ALL are levels of the selected factor
      )

      Samples <- reactive({
        rownames(pheno.dat)[ which(pheno.dat[, input$Pheno] %in% input$levels) ]
      })
      Expr.vals.mat <- t(
        assays(eselist[[1]])[[ Assay_Type() ]][ names(Gene.Identifiers()), Samples(), drop=F ]
      )
      
      setNames( data.frame(
        pheno.dat[ Samples(), input$Pheno, drop=F ],
        Expr.vals.mat[ Samples(), , drop=F ],
        row.names = Samples(), stringsAsFactors=FALSE, check.names=F
      ), c("grps",colnames(Expr.vals.mat)) )
    })

    Bxp.Vals.df <- eventReactive(changeData(), {
      tidyr::gather(Plt.Vals.df(), gene, vals, grep("^ENS",colnames(Plt.Vals.df())) )
    })

    observe({
        cat("The LEVELS of \"",input$Pheno,"\" selected are :-\n", sep="")
        print(input$levels)
        cat("State of \"input$levels\" : '",isTruthy(input$levels),"'\n", sep="")


      output$bxPlot1 <- renderPlot({
        #shinyngs:::Display.Bxp(
        Display.Bxp(
          Plt.Vals.df(), Bxp.Vals.df(), Gene.Identifiers(), Assay_Type(), input$Pheno, input$Pt.Colour,
          names(assays.lst)[ assays.lst==input$AssayType ]
        )
      }, height = 900)

    })
}


Display.Bxp <- function(Bxp_Vals, Vals_to_Plot, Gene_IDs, Assay, Condition, Palette.Name, Y.label) {
    
    Groups <- as.factor( unlist(Bxp_Vals$grps) )
    #cat("Factor levels :-\n")
    #print(table(Groups, useNA="always"))
    #cat("\n")
    
    Leg.Tab <- table(Groups)
    print(Leg.Tab)
    Leg.Strg <- paste0(names(Leg.Tab)," (",Leg.Tab,")")


    # ColorBrewer palettes : 
    # "Diverging" --> 3-11 / "Qualitative" -->  3-[8,9,12]
    
    #Palette <- if (length(Leg.Tab)==1) { "tan3" }
    #else if (length(Leg.Tab)==2)       {  c("gold3","darkgreen") }

    Palette <- if (length(Leg.Tab)<3) {
      RColorBrewer:::brewer.pal(3, Palette.Name)[1:length(Leg.Tab)]
    } else {
      max.colours <- RColorBrewer:::brewer.pal.info[Palette.Name, "maxcolors"]
      
      if (length(Leg.Tab)>max.colours) {
        colorRampPalette( RColorBrewer:::brewer.pal(max.colours, name=Palette.Name), space="Lab")( length(Leg.Tab) )
      }
      else { RColorBrewer:::brewer.pal(length(Leg.Tab), Palette.Name) }
    }
    names(Palette) <- names(Leg.Tab)
    
    #cat("Palette :-\n")
    #print(class(Palette))
    #print(Palette[ as.character(Groups) ])
    

    no.Grps  <- length(Leg.Tab)
    no.Genes <- length(Gene_IDs)
    no.Plots <- no.Grps*no.Genes
    cat("\nNo. of box plots to be displayed : ",no.Plots,"\n", sep="")

    
    # ‘plt’ A vector of the form ‘c(x1, x2, y1, y2)’ giving the coordinates
    #  of the plot region as fractions of the current figure region.
    op.PlotRegion <- par()$plt
    op.PlotRegion[3] = 0.175
    op.PlotRegion[4] = 0.95 

    PlotRegion.coords <- if (no.Plots<4) {
      c(0.375, 0.625, op.PlotRegion[3:4])   # 25% width for three or less plots
    } else if (no.Plots<8) {
      c(0.25, 0.75, op.PlotRegion[3:4])     # 50% for between four & seven
    } else if (no.Plots<11) { 
      c(0.125, 0.875, op.PlotRegion[3:4])   # 75% for eight to ten
    } else { c(0.085, op.PlotRegion[2:4]) }


    op.margins <-par()$mar # ‘c(bottom, left, top, right)’
    op.MarginLine <- par()$mgp # c("axis title", "axis labels", "axis line")

    par( 
      mar = c(12, op.margins[2]+1.85, 2.5, op.margins[4]),
      mgp = c(op.MarginLine[1]+1.75, op.MarginLine[2:3]),
      xpd = TRUE,
      plt = PlotRegion.coords
    )
    
    cat("\nTEST #1 : 'boxplot' formula\n")
    Bxp <- graphics::boxplot(vals ~ gene*grps, data=Vals_to_Plot, ylab=Y.label, las=1, xaxt="n", xlab="")

    print(head(Vals_to_Plot))
    #box("inner", col="red")
    
    Plotting.Symbs <- c(15:18,21:25)
    Pt.Shapes <- setNames( Plotting.Symbs[1:length(Gene_IDs)], names(Gene_IDs))
    #cat("\nTEST #1 : Point Shapes.\n")
    #print(Pt.Shapes)
    
    Plt.Grps <- setNames(Bxp$n, Bxp$names)
    #cat("\nTEST #2 : Groups.\n")
    #print(Plt.Grps)
    
    Point.Coords.df <- do.call("rbind", lapply(names(Plt.Grps),
      function(x) {
        Reg.Ex  <- paste0("^(ENS\\S*G\\d+)\\.(\\S.*\\S)\\s*$")
        gene.id <- sub(Reg.Ex, "\\1", x, perl=T)
        Cond    <- sub(Reg.Ex, "\\2", x, perl=T)
        
        sample.ids <- rownames(Bxp_Vals)[ Bxp_Vals$grps==Cond ]
        
        setNames( data.frame(
          sample.ids,
          rep(which(names(Plt.Grps)==x), length(sample.ids)),
          Bxp_Vals[ sample.ids, gene.id ],
          rep(Pt.Shapes[ gene.id ], length(sample.ids)),
          rep(Palette[Cond], length(sample.ids)),
          row.names=NULL, stringsAsFactors=FALSE, check.names=F
        ), c("Sample","Pt.X_Val", "Pt.Y_Val", "Pt.Shape", "Pt.Colour"))
      }
    ))

    Pt.Size <- ifelse(no.Plots<11, 1.75, ifelse(no.Plots<21, 1.45, 0.975) )

    with(
      Point.Coords.df,
      points(jitter(Pt.X_Val), Pt.Y_Val, pch=Pt.Shape, col=Pt.Colour, cex=Pt.Size)
    )
    #box(col="gold", lwd=6, which="inner")
    #box(col="darkgreen", lwd=4, which="outer")
    
    
    no.key.cols <- ifelse(length(Palette)>15, 3,
      ifelse(length(Palette)>2, 2, length(Palette))
    )

    legend(
      "top", inset=1.01,
      lty=1, lwd=5, title.adj=0, bty="n", 
      legend=Leg.Strg, pch=NA, pt.cex=1.75, title=paste0("Phenotype : '",Condition,"'"), horiz=F,
      col=Palette, ncol=no.key.cols
    )
    
    legend(
      "top", inset=-0.07,
      bty="n", col="snow4", cex=1.5, pt.cex=2.5, text.font=4, text.col="snow4", # "navajowhite4"
      legend=Gene_IDs[ names(Pt.Shapes) ], pch=Pt.Shapes, horiz=T
    )
    
    par(xpd=FALSE)
    if (length(Groups)>1) {
      run.lengths <- rle(Point.Coords.df$Pt.Colour)$lengths
      section.idx <- cumsum(run.lengths)[ -length(run.lengths) ]
      
      grp.divider.pos <- Point.Coords.df$Pt.X_Val[ section.idx ] + 0.5
      abline(v=grp.divider.pos, col="navy", lty=3, lwd=0.50)
    }
    
}
    
    
    
    
    


## --------------------------------------------------------------------------------------------------------- ##

## Shiny modules
## http://zevross.com/blog/2016/04/19/r-powered-web-applications-with-shiny-a-tutorial-and-cheat-sheet-with-40-example-apps/#shiny-modules


## Module UI ##
#  list(
#   plotOutput(ns("plot1"))
#  )

## Module server ##
#  output$plot1 <- renderPlot({
#    ggplot(mtcars, aes(wt, mpg)) + geom_point(color=input$Pt.Colour, size=2)
#  })

## app ui ##
#  h3("This is not part of the module but the plot below was created with a module"),
#  scatterUI("prefix")

## app server ##
#  server <- function(input, output,session) {
#    callModule(scatter, "prefix", "purple")
#  }



