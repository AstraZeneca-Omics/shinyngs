#' The UI input function of the sampleselect module
#' 
#' This module provides controls for selecting matrix columns by sample or group
#' name.
#' 
#' This will generally not be called directly, but by other modules such as the 
#' selectmatrix module.
#' 
#' @param id Submodule namespace
#' @param eselist ExploratorySummarizedExperimentList object containing
#'   ExploratorySummarizedExperiment objects
#' @param getExperiment Reactive expression that returns a
#'   \code{ExploratorySummarizedExperiment} with assays and metadata. Usually a 
#'   result of a user selection
#' @param select_samples Select samples at all? If set to false, a hidden input
#'   indicating the selection of all samples is produced. 
#'   
#' @return output An HTML tag object that can be rendered as HTML using 
#'   as.character()
#'   
#' @keywords shiny
#'   
#' @examples
#' sampleselectInput(ns('heatmap'))

sampleselectInput <- function(id, eselist, getExperiment, select_samples = TRUE) {
    
    ns  <- NS(id)
    ese <- getExperiment()
    
    if (select_samples) {
        
        # If grouping variables have been supplied we can use them to define sample selection
        selectby <- "name"
        if (length(eselist@group_vars) > 0) {
            selectby <- c(selectby, "group")
            default_groupvar <- ese$group_vars[1]
            if (length(eselist@group_vars) > 0) {
                default_groupvar <- eselist@default_groupvar
            }
        }

        inputs <- list(
            h5("Query conditions ('group') or sample identifiers ('name')."),
            selectInput(ns("sampleSelect"), "Select by", selectby, selected = selectby[length(selectby)])
        )
        # Add in group selection if relevant
        if (length(eselist@group_vars) > 0) {
            inputs <- pushToList(inputs,
                selectInput(
                    inputId  = ns("sampleGroupVar"), 
                    label    = "Choose condition of interest:",
                    choices  = structure(eselist@group_vars, names = prettifyVariablename(eselist@group_vars)),
                    selected = eselist@default_groupvar
                )
            )
        }
        
        # We can select by sample in any case
        inputs <- pushToList(inputs,
            list(
                conditionalPanel(
                    condition = paste0("input['", ns("sampleSelect"), "'] == 'name'"),
                    shinyWidgets::pickerInput(
                        inputId = ns("samples"),
                        label   = "",
                        choices = "",
                        options = list(
                            `actions-box`= TRUE,
                            style = "btn-royal btn-sm",
                            `selected-text-format` = "count > 3"
                        ),
                        multiple = TRUE
                    )
                ),
                conditionalPanel(
                    condition = paste0("input['", ns("sampleSelect"), "'] == 'group'"),
                    shinyWidgets::pickerInput(
                        inputId = ns("sampleGroupVal"),
                        label   = "",
                        choices = "",
                        options = list(
                            `actions-box`= TRUE,
                            style = "btn-royal btn-sm",
                            `selected-text-format` = "count > 3"
                        ),
                        multiple = TRUE
                    )
                ),
                uiOutput(ns("groupSamples"))
            )
        )
    } else {
        inputs <- list(hiddenInput(ns("sampleSelect"), "all"))
    }
    
    tagList(inputs)
}

#' The server function of the sampleselect module
#'  
#' This module provides controls for selecting matrix columns by sample or 
#' group name.
#' 
#' This function is not called directly, but rather via callModule() (see 
#' example).
#'
#' @param input Input object
#' @param output Output object
#' @param session Session object
#' @param eselist ExploratorySummarizedExperimentList object containing
#'   ExploratorySummarizedExperiment objects
#' @param getExperiment Reactive expression that returns a
#'   \code{ExploratorySummarizedExperiment} with assays and metadata. Usually a 
#'   result of a user selection
#'
#' @return output A list of reactive functions for interrogating the selected
#' samples/ columns.
#'
#' @keywords shiny
#' 
#' @examples
#' selectSamples <- callModule(sampleselect, 'selectmatrix', getExperiment)

sampleselect <- function(input, output, session, eselist, getExperiment) {
    
    ns <- session$ns
    getSummaryType <- callModule(summarisematrix, "summarise")
    #print(slotNames(eselist))

    output$groupSamples <- renderUI({
        list(
            radioButtons(
                inputId  = ns("autoPlotting"),
                label    = "Automatic Plotting of Data:",
                choices  = list(`'ON'`="Auto", `'OFF' - click the below button to display`="Manual"),
                selected = "Manual"
            ),
            conditionalPanel(
                condition = paste0("input['",ns("autoPlotting"),"'] == 'Manual'"),
                shinyWidgets::actionBttn(
                    inputId = ns("DisplayPlots"),
                    label   = "Analyse Selected Samples",
                    style   = "gradient",
                    color   = "success",
                    block   = TRUE, size = "sm",
                    icon    = icon("stream")
                )
            ),
            hr(),
            summarisematrixInput( ns("summarise") )
        )
    })

    # Render the sampleGroupVal() element based on sampleGroupVar
    observeEvent(input$sampleGroupVar, {
        ese <- getExperiment()

        validate(need(input$sampleGroupVar, FALSE))
        cond.levs.tab <- sort(table(ese[[ input$sampleGroupVar ]]), decreasing=TRUE)

        cond.choices.lst <- as.list( names(cond.levs.tab) )
        names(cond.choices.lst) <- paste0(names(cond.levs.tab)," [",cond.levs.tab,"]")

        no.grps.toShow <- ifelse(length(cond.choices.lst)>2, 2, length(cond.choices.lst))

        shinyWidgets::updatePickerInput(
            session, "sampleGroupVal",
            label    = paste0("Add/remove '",input$sampleGroupVar,"' group(s):"),
            choices  = cond.choices.lst,
            selected = cond.choices.lst[1:no.grps.toShow]
        )

        sampleGrp.lst <-  reactive({
            Grps <- ese[[ input$sampleGroupVar ]]
            names(Grps) <- colnames(ese)

            grp.lst <- lapply(unique(Grps),
                function(x) {
                    # https://github.com/dreamRs/shinyWidgets/issues/58
                    as.list( sort(names(Grps)[ Grps==x ]) )
                }
            )
            names(grp.lst) <- unique(Grps)
            return(grp.lst)
        })

        GrpList.tab  <- unlist( lapply(sampleGrp.lst(), function(x) length(unlist(x))) )
        GrpList.Rank <- order(GrpList.tab, -xtfrm(names(GrpList.tab)), decreasing = T)
        cat("Dropdown list :-\n")
        print(GrpList.tab[GrpList.Rank])

        Selected.ranks   <- if (length(GrpList.Rank)>2) GrpList.Rank[1:2] else GrpList.Rank
        Selected.samples <- Reduce(
            union, unlist( lapply(sampleGrp.lst()[ Selected.ranks ], function(x) unlist(x)) )
        )
        cat("Selected sample IDs :-\n")
        print(Selected.samples)

        shinyWidgets::updatePickerInput(
            session, "samples",
            choices  = sampleGrp.lst()[ GrpList.Rank ],
            selected = Selected.samples,
            label    = paste0("Sample names grouped by '",input$sampleGroupVar,"' :"),
        )
    })


    # Output a reactive so that other modules know whether we've selected by sample or group
    getSampleSelect <- reactive({ input$sampleSelect })
    # Return summary type
    getSampleGroupVar <- reactive({ input$sampleGroupVar })
    
    Samples.lst <- list(
        getSampleGroupVar = getSampleGroupVar, getSummaryType = getSummaryType, getSampleSelect = getSampleSelect
    )

    # Reactive expression for selecting the specified columns
    Samples.lst[[ "selectSamples" ]] = reactive({

        validate(need(!is.null(input$sampleSelect), "Waiting for form to provide sampleSelect"))

        withProgress(message = "Selecting samples", value = 0, {
            ese <- getExperiment()
            sample.sheet <- colData(ese)

            if (! isTruthy(input$autoPlotting) || input$autoPlotting=="Manual") {
                manual.selected <- eventReactive(input$DisplayPlots, {
                    Selected_Sample_IDs(
                        input$sampleSelect, input$sampleGroupVar, input$sampleGroupVal,
                        rownames(sample.sheet), input$samples, length(eselist@group_vars), sample.sheet[[ input$sampleGroupVar ]]
                    )
                }, ignoreNULL=FALSE)
                
                # https://stackoverflow.com/questions/57182912/r-action-button-and-eventreactive-not-working
                observe( manual.selected() )
                return( manual.selected() )
            
            } else {
                Selected_Sample_IDs(
                    input$sampleSelect, input$sampleGroupVar, input$sampleGroupVal,
                    rownames(sample.sheet), input$samples, length(eselist@group_vars), sample.sheet[[ input$sampleGroupVar ]]
                )
            }
        })
    })
    return(Samples.lst)
}

Selected_Sample_IDs <- function(Selection_Type, Group_Variable, Group_Values, Sample_IDs, Selected_IDs, No_Grps, Levels) {

    selected.ids <- if (Selection_Type == "all") {
        Sample_IDs
    } else {
        validate(need(!is.null(Selected_IDs), "Waiting for form to provide samples"))

        if (No_Grps>0) {
            validate(need(!is.null(Group_Variable), FALSE))
        }

        if (Selection_Type == "name") {
            Selected_IDs
        } else {
            samplegroups <- as.character(Levels)
            # Any NA in the colData will become string '' via the inputs, so make sure we consider that when matching
            samplegroups[ is.na(samplegroups) ] <- ""

            Sample_IDs[ samplegroups %in% Group_Values ]
        }
    }
    return(selected.ids)
}




