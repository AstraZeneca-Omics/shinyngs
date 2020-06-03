#' The input function of the colorby module
#' 
#' This module provides a drop-down for picking an RColorBrewer color palette
#' and provides that palette given a reactive which supplied the required 
#' number of colors.
#' 
#' This funcion provides the form elements to control the display
#' 
#' @param id Submodule namespace
#'   
#' @return output An HTML tag object that can be rendered as HTML using 
#'   as.character()
#'   
#' @keywords shiny
#'   
#' @examples
#' colormakerInput('myid')

colormakerInput <- function(id) {
    
    ns <- NS(id)

    ColorBrewer_Categories   <- c("Diverging","Qualitative")
    ColorBrewer.Palettes.lst <- lapply(ColorBrewer_Categories,
        function(x) {
            Data.Nature <- sub("^(\\S{3})\\S+$", "\\L\\1", x , perl=T)
            idx <- grep(paste0("^",Data.Nature), as.character(RColorBrewer::brewer.pal.info$category), perl=T)
            rownames(RColorBrewer::brewer.pal.info)[ idx ]
        }
    )
    names(ColorBrewer.Palettes.lst) <- ColorBrewer_Categories

    selectizeInput(inputId=ns("palette_name"), label="Colour Palettes", choices=ColorBrewer.Palettes.lst, selected="Set1")
}

#' The output function of the colorby module
#' 
#' This module provides a drop-down for picking an RColorBrewer color palette
#' and provides that palette given a reactive which supplied the required 
#' number of colors.
#' 
#' This function is not called directly, but rather via callModule() (see 
#' example).
#' 
#' @param id Submodule namespace
#' @param getNumberCategories A reactive supplying the number of categories 
#' that require a color.
#'   
#' @return output An HTML tag object that can be rendered as HTML using 
#'   as.character()
#'   
#' @keywords shiny
#'   
#' @examples
#' callModule(colormaker, 'myid', getNumberCategories)

colormaker <- function(input, output, session, getNumberCategories) {
    
    getPaletteName <- reactive({
        validate(need(!is.null(input$palette_name), "Waiting for palette"))
        input$palette_name
    })
    
    reactive({
        palette <- getPaletteName()
        n_colors <- getNumberCategories()
        
        makeColorScale(n_colors, palette)
    })
}

#' Make a color palette of a specified length
#' 
#' Given an integer, make a palette with a specified number of colors using
#' palettes from RColorBrewer, and interpolation where necessary.
#'
#' @param ncolors Integer specifying the number of colors
#' @param palette RColorBrewer palette name. (default: 'Set1')
#'
#' @return output Character vector of colors
#' @export
#'
#' @examples
#' makeColorScale(10)
#' [1] '#999999' '#EC83BA' '#B75F49' '#E1C62F' '#FFB716' '#D16948' '#7E6E85' '#48A462' '#4A72A6' '#E41A1C'

makeColorScale <- function(ncolors, palette="Set1") {

    cols <- if (ncolors<3) {
      RColorBrewer::brewer.pal(3, palette)[1:ncolors]
    } else {
      max.colours <- RColorBrewer::brewer.pal.info[palette, "maxcolors"]

      if (ncolors>max.colours) {
        colorRampPalette( RColorBrewer::brewer.pal(max.colours, name=palette), space="Lab")( ncolors )
      }
      else { RColorBrewer::brewer.pal(ncolors, palette) }
    }
} 













