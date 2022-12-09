# Shiny app for visualization and exploration of collateral effects networks
# You can run the application by clicking 'Run App' above.
#
#

#Libraries
library(shiny)
library(visNetwork)
library(tidyverse)

#Data
#results from collateral effects analysis
load("data/collateral_effects.Rdata")
all_results <- all_results %>% filter(q < 0.05) %>%
  mutate(effect_type = ifelse(effect_type == "CR", "collateral resistance", "collateral sensitivity"))
#information on antibiotics
meta_antibio_abbr <- read.csv("data/meta_antibio_abbr.csv", header = T, stringsAsFactors = F)
#colors of the antibiotic classes
color_ab_groups <- read.csv("data/color_ab_groups.csv", header = T, stringsAsFactors = F)

##### Define UI ######
ui <- fluidPage(
  # Application title
  titlePanel("Collateral effects network"),

  # Sidebar with inputs for app
  sidebarLayout(
    sidebarPanel(width = 3,
      checkboxGroupInput(
        "CE_type", "Select collateral effect", unique(all_results$effect_type),
        selected = "collateral sensitivity"
      ),

      uiOutput("species"),
      uiOutput("antibiotics"),

      sliderInput(inputId = "effect_size", label = "Minimal absolute log2(fold change)", min = 0, max = 10, step = 0.1, value = 0.5),
      sliderInput(inputId = "balance", label = "Minimal balance between groups", min = 0, max = 1, step = 0.01, value = 0.05),

      actionButton("execute_button", "Create network"),
      htmlOutput("explanation")
    ),

    # Show a plot of the generated network
    mainPanel(width = 9,
      visNetworkOutput("network", height = "700px")
    )
  )
)


##### Define server functionality #####
server <- function(input, output) {
  filterData <- reactiveVal(all_results)

  # Filter data based on input data source
  filtered_df <- reactive({
    filterData()
  })

  # Add button for species selection
  output$species <- renderUI({
    selectInput(
      "species", "Which species do you want to include? (Leave empty to include all)", sort(unique(filtered_df()$species)),
      multiple = TRUE
    )
  })
  # Add button for antibiotics selection
  output$antibiotics <- renderUI({
    selectInput(
      "antibiotics", "Which antibiotics do you want to include? (Leave empty to include all)",
      sort(unique(c(filtered_df()$A,filtered_df()$B))),
      multiple = TRUE
    )

  })

  # Filter based on species and antibiotic selection
  filtered_df_species <- reactive({
    user_chosen <- filtered_df()
    if (!is.null(input$species)) {
      user_chosen <- filtered_df() %>% filter(species %in% input$species)
    }
    if (!is.null(input$antibiotics)) {
      user_chosen <- user_chosen %>% filter(A %in% input$antibiotics & B %in% input$antibiotics)
    }
    user_chosen
  })

  # Create network edges and nodes from chosen filters
  network_data <- eventReactive(input$execute_button, {
    # Create nodes for every antibiotic
    nodes <- data.frame(ab = sort(unique(c(filtered_df_species()$A, filtered_df_species()$B))),
                        stringsAsFactors = FALSE) %>%

      left_join(meta_antibio_abbr %>%
                filter(!is.na(group)) %>%
                rename(ab = abbreviation), by = "ab") %>%
      distinct(across(-antibiotic), .keep_all = TRUE) %>%
      mutate(id = paste0("n", row_number())) %>%
      left_join(color_ab_groups %>% rename(group = ab_group), by = "group")
    # Create edges for every collateral effect relation
    full_edges <- filtered_df_species() %>%
      filter(effect_type %in% input$CE_type) %>%
      filter(abs(effect_size) > input$effect_size) %>%
      filter(balance > input$balance) %>%
      select(A, B, species, effect_size, effect_type) %>%
      left_join(nodes %>% rename(A = ab, to = id) %>% select(A, to), by = "A") %>%
      left_join(nodes %>% rename(B = ab, from = id) %>% select(B, from), by = "B")
    # Combine edges that are of different species
    edges <-  full_edges %>%
      group_by(from, to) %>%
      summarise(weight = mean(effect_size), sources = paste(species, collapse = ","),
                effect_type = ifelse(mean(effect_size) > 0, "#f46e32", "#001158"), .groups = "keep") %>%
      ungroup() %>% distinct()

    # Create elements of network
    vis.nodes <- nodes
    if (is.null(input$antibiotics)) {
      vis.nodes <- nodes %>% filter(id %in% edges$from | id %in% edges$to)
    }
    vis.links <- edges

    # Add options to network e.g. color, labeling etc.
    if (nrow(vis.links) > 0 | nrow(vis.nodes) > 0) {
      vis.nodes$shape  <- "ellipse"
      vis.nodes$label  <- vis.nodes$ab # Node label
      vis.nodes$title  <- vis.nodes$antibiotic

      vis.links$arrows <- "to"
      #if (length(unique(input$species)) != 1) {
        vis.links$title <- vis.links$sources #add labels species
      #}
      vis.links$color <- vis.links$effect_type
      vis.links$width <- abs(vis.links$weight)
      vis.links$arrowStrikethrough <- FALSE
    }

    # Create custom legend data frame
    color_nodes <- color_ab_groups %>% filter(ab_group %in% vis.nodes$group) %>%
      mutate(label = ab_group, value = 13, font.size = 13*3)
    color_edges <- data.frame(color = c("#f46e32", "#001158"),
                              label = c("collateral resistance", "collateral sensitivity"),
                              arrows = c("to", "to"),
                              width = 7, font.size = 13*3, font.align = "top", align = "bottom") %>%
      filter(color %in% edges$effect_type)

    list(nodes = vis.nodes, edges = vis.links, color_nodes = color_nodes,
         color_edges = color_edges)
  })

  # Actual output for the shiny app to show
  output$network <- renderVisNetwork({
    visNetwork(network_data()$nodes, network_data()$edges) %>%
      visOptions(highlightNearest = list(enabled = T, degree = 1),
                 nodesIdSelection = list(enabled = TRUE, main = "Highlight antibiotic")) %>%
      visLegend(position = "right", useGroups = FALSE, addNodes = network_data()$color_nodes,
                addEdges = network_data()$color_edges) %>%
      visIgraphLayout(layout = "layout_with_lgl", smooth = FALSE, randomSeed = 123, type = "full")

  })

  output$explanation <- renderUI({
    HTML("<br><b>Explanation</b><br>
         <ul>
            <li>Choose which antibiotics and species you want to include and press
         'Create network' to create a collateral effect network.</li>
            <li>Hover over the antibiotic to show the full antibiotic name and
                hover over the arrow to show the name(s) of the species.</li>
            <li>Use scrolling to zoom in and out.</li>
            <li>You can drag the nodes around to improve visibility.</li>
         <ul>")
  })


}

# Running the shiny app
shinyApp(ui = ui, server = server)
