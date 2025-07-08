#####################################################
# Perform personal selection of spatial domain
#####################################################

#' @name Run_Interactive
#' @title Interactive Visualization for Spatial Transcriptomics Data and spot plot with save functionality
#' @details This function generates interactive plots to visualize spatial transcriptomics data. It takes matched spatial coordinates and raw UMI counts to produce customized t-SNE or UMAP plots overlayed on a optional background H&E image.
#' @param matched_data A data frame containing matched spatial coordinates with raw UMI counts.
#' @param clustering_result A data frame containing matched spatial coordinates with raw UMI counts.
#' @param background_img Optional background H&E image.
#' @param reduction_method T-SNE ("tsne") or UMAP ("umap") result data frame obtained from 'Run_Clustering' function. Default set as "tsne".
#' @param point_size Size of point shown in the spatial heatmap. Default set as 1.
#' @examples
#' matched_data <- data.frame(
#' X_coordinate = runif(10, 0, 100),
#' Y_coordinate = runif(10, 0, 100),
#' UMI_count = sample(seq_len(100), 10),
#' spatial_name = paste0("Spot", seq_len(10)),
#' stringsAsFactors = FALSE
#' )
#' clustering_result <- data.frame(
#' TSNE1 = runif(10, -50, 50),
#' TSNE2 = runif(10, -50, 50),
#' spot = paste0("Spot", seq_len(10)),
#' cluster = sample(seq_len(3), 10, replace = TRUE),
#' stringsAsFactors = FALSE
#' )
#' if (interactive()) {
#' Run_Interactive(
#'   matched_data = matched_data, 
#'   clustering_result = clustering_result, 
#'   background_img = NULL, 
#'   reduction_method = "tsne", 
#'   point_size = 1)
#'                                }
#' @return R-shiny interactive webpage.
#' @export
#' @importFrom ggplot2 ggplot scale_fill_brewer aes geom_bar geom_text theme_minimal theme labs ggsave element_text geom_point scale_color_gradient element_blank element_rect xlim ylim scale_color_brewer
#' @importFrom shiny reactive

if (getRversion() >= "2.15.1") {
  utils::globalVariables(
    c("TSNE1", "TSNE2", "UMAP1", "UMAP2", "cluster",
      "X_coordinate", "Y_coordinate", "UMI_count",
      "spatial_name")
  )
}

Run_Interactive <- function(matched_data, clustering_result, background_img = NULL, reduction_method = "tsne", point_size = 1) {
  # Data pre-check
  required_cols_matched <- c("X_coordinate", "Y_coordinate", "UMI_count", "spatial_name")

  # Define required columns based on reduction_method
  if (reduction_method == "tsne") {
    required_cols_reduction <- c("TSNE1", "TSNE2", "spot", "cluster")
  } else if (reduction_method == "umap") {
    required_cols_reduction <- c("UMAP1", "UMAP2", "spot", "cluster")
  } else {
    stop("reduction_method must be either 'tsne' or 'umap'.")
  }

  # Check if required columns are present in matched_data
  if (!all(required_cols_matched %in% colnames(matched_data))) {
    stop(paste("matched_data must contain columns:", paste(required_cols_matched, collapse = ", ")))
  }

  # Check if required columns are present in clustering_result
  if (!all(required_cols_reduction %in% colnames(clustering_result))) {
    stop(paste("clustering_result must contain columns:", paste(required_cols_reduction, collapse = ", ")))
  }

  # Create UMI count spatial data frame
  plot_d <- matched_data %>%
    dplyr::select(X_coordinate, Y_coordinate, UMI_count, spatial_name)

  # Combine clustering_result with plot_d
  plot_d <- plot_d %>%
    dplyr::left_join(clustering_result, by = c("spatial_name" = "spot"))
  plot_d$cluster <- as.factor(plot_d$cluster)


  # Create a highlight key for interactivity
  key <- plotly::highlight_key(plot_d, ~spatial_name)

  # Plot left UMI count spatial plot
  if (!is.null(background_img)) {
    maxX <- dim(background_img)[1]
    maxY <- dim(background_img)[2]
    p1 <- ggplot2::ggplot(mapping = ggplot2::aes(seq_len(maxX), seq_len(maxY))) +
      ggplot2::annotation_raster(background_img, xmin = 1, xmax = maxX, ymin = 1, ymax = maxY)
  } else {
    maxX <- max(plot_d$X_coordinate)
    minX <- min(plot_d$X_coordinate)
    maxY <- max(plot_d$Y_coordinate)
    minY <- min(plot_d$Y_coordinate)
    p1 <- ggplot2::ggplot(mapping = ggplot2::aes(minX:maxX, minX:maxY))
  }

  p1 <- p1 +
    ggplot2::geom_point(
      data = key,
      ggplot2::aes(x = X_coordinate, y = maxY - Y_coordinate, colour = UMI_count, alpha = ifelse(is.null(background_img), 1, 0.85)),
      size = point_size
    ) +
    ggplot2::coord_fixed() +
    ggplot2::scale_color_viridis_c() +
    ggplot2::labs(color = "UMI Count") +
    ggplot2::theme_void()

  # Plot right UMAP or t-SNE plot
  if (reduction_method == "tsne") {
    p2 <- ggplot2::ggplot(key, ggplot2::aes(x = TSNE1, y = TSNE2, colour = as.factor(cluster))) +
      ggplot2::geom_point(size = 1) +
      ggplot2::scale_color_discrete(name = "Cluster") +
      ggplot2::labs(color = "Cluster") +
      ggplot2::theme_minimal()
  } else {
    p2 <- ggplot2::ggplot(key, ggplot2::aes(x = UMAP1, y = UMAP2, colour = as.factor(cluster))) +
      ggplot2::geom_point(size = 1) +
      ggplot2::scale_color_discrete(name = "Cluster") +
      ggplot2::labs(color = "Cluster") +
      ggplot2::theme_minimal()
  }

  # Convert plots to interactive
  p1 <- plotly::ggplotly(p1, source = "A") %>%
    plotly::highlight(on = "plotly_selected", off = "plotly_deselect") %>%
    plotly::config(
      toImageButtonOptions = list(
        format = "svg",
        filename = "UMI_Spatial_Plot",
        width = 800,
        height = 500
      )
    )

  p2 <- plotly::ggplotly(p2, source = "A") %>%
    plotly::highlight(on = "plotly_selected", off = "plotly_deselect") %>%
    plotly::config(
      toImageButtonOptions = list(
        format = "svg",
        filename = "Reduction_Plot",
        width = 800,
        height = 500
      )
    )

  # Create Shiny app for saving functionality
  shiny::shinyApp(
    ui = shiny::fluidPage(
      shiny::fluidRow(
        shiny::column(6, plotly::plotlyOutput("spatial_plot")),
        shiny::column(6, plotly::plotlyOutput("reduction_plot"))
      ),
      shiny::actionButton("add_selection", "Add Selection", icon = shiny::icon("plus")),
      shiny::actionButton("clear_last", "Clear Last Selection", icon = shiny::icon("eraser")),
      shiny::actionButton("reset", "Reset All Selections", icon = shiny::icon("refresh")),  # New reset button
      shiny::actionButton("save", "Save All Selected Regions of Interest (ROI)", icon = shiny::icon("save")),
      shiny::tableOutput("selected_data")
    ),
    server = function(input, output, session) {
      # Reactive expressions for generating the initial plots
      initial_spatial_plot <- reactive({
        plotly::ggplotly(p1, source = "A") %>%
          plotly::highlight(on = "plotly_selected", off = "plotly_deselect") %>%
          plotly::config(
            toImageButtonOptions = list(
              format = "svg",
              filename = "UMI_Spatial_Plot",
              width = 800,
              height = 500
            )
          )
      })

      initial_reduction_plot <- reactive({
        plotly::ggplotly(p2, source = "A") %>%
          plotly::highlight(on = "plotly_selected", off = "plotly_deselect") %>%
          plotly::config(
            toImageButtonOptions = list(
              format = "svg",
              filename = "Reduction_Plot",
              width = 800,
              height = 500
            )
          )
      })

      # Render the initial plots
      output$spatial_plot <- plotly::renderPlotly({ initial_spatial_plot() })
      output$reduction_plot <- plotly::renderPlotly({ initial_reduction_plot() })

      # Save selected subset in reactive list
      selected_regions <- shiny::reactiveValues(data = list())
      selected_spatials <- shiny::reactiveVal(character())

      # Add selection to the list
      shiny::observeEvent(input$add_selection, {
        selected_points <- plotly::event_data("plotly_selected", source = "A") # Get selected data
        if (!is.null(selected_points)) {
          selected_spatials(unique(c(selected_spatials(), selected_points$key))) # Update selected spatials
          selected_subset <- matched_data %>% dplyr::filter(spatial_name %in% selected_spatials()) # Use updated selected_spatials
          selected_regions$data <- list(selected_subset)
        }
      })

      # Clear the last selection
      shiny::observeEvent(input$clear_last, {
        if (length(selected_regions$data) > 0) {
          selected_regions$data <- list()
          selected_spatials(character()) # Clear selected spatials
        }
      })

      # Reset all selections and fully re-render the plots
      shiny::observeEvent(input$reset, {
        selected_regions$data <- list()
        selected_spatials(character())

        # Trigger re-rendering by reassigning the plots
        output$spatial_plot <- plotly::renderPlotly({ initial_spatial_plot() })
        output$reduction_plot <- plotly::renderPlotly({ initial_reduction_plot() })
      })

      # Display combined selected data
      output$selected_data <- shiny::renderTable({
        if (length(selected_regions$data) > 0) {
          combined_selected_data <- do.call(rbind, selected_regions$data)
          return(combined_selected_data)
        }
      })

      # Save all selected data
      shiny::observeEvent(input$save, {
        if (length(selected_regions$data) > 0) {
          combined_selected_data <- do.call(rbind, selected_regions$data)
          assign("selected_ROI", combined_selected_data, envir = .GlobalEnv)
          shiny::showNotification("All selected regions saved to global environment as 'selected_ROI'", type = "message")
        }
      })
    }
  )
}
