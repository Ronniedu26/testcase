# MLDV Latest

if (!require('rsconnect')) {
  install.packages(c('rsconnect', 'shiny', 'pacman', 'ggplot2'))
}

pacman::p_load(
  stringr,
  readr,
  magrittr,
  dplyr,
  plotly,
  data.table,
  plotly,
  shinydashboard,
  shinyWidgets,
  readr,
  ggthemes,
  showtext
)

# Get data
n_protein_data_stats_sig_adj <-
  read_csv("ronnie_ecoli_data_norm.csv", col_names = TRUE)

# Transpose matrix and select rows that include LOG and STAT
invertedProteome <- as.data.frame(t(n_protein_data_stats_sig_adj))
names(invertedProteome) <- invertedProteome[3,]
invertedProteome <- invertedProteome[-1:-4,]
TMTs <-
  invertedProteome[str_which(rownames(invertedProteome), "\\wLOG|\\wSTAT"),]

## above for intensities (relabel)

# Set font to Arial
arial <- list(family = "arial",
              size = 11)

## Housekeeping for E. coli graph
low_cutoff <- -1
p_data_stats_sig_adj <-
  read_csv("ronnie_ecoli_data_norm.csv", col_names = TRUE)
mukha_4cols <- c("gray", "#440154FF", "#55C667FF", "#FDE725FF")

sig_colors_3 <- p_data_stats_sig_adj

## Renaming 1, 2, 3, 0 to respective significance
sig_colors_3$sig_all <-
  str_replace(sig_colors_3$sig_all, "1", "Significant PISA")
sig_colors_3$sig_all <-
  str_replace(sig_colors_3$sig_all, "2", "Significant PROT")
sig_colors_3$sig_all <-
  str_replace(sig_colors_3$sig_all, "3", "Significant PISA and PROT")
sig_colors_3$sig_all <-
  str_replace(sig_colors_3$sig_all, "0", "Not significant")
notSignificant <-
  sig_colors_3 %>% filter(sig_colors_3$sig_all == "Not significant")
sig_colors_3$sig_all[is.na(sig_colors_3$sig_all)] = "Not significant"

PISAsignificant <-
  sig_colors_3 %>% filter(sig_colors_3$sig_all == "Significant PISA")
BOTHsignificant <-
  sig_colors_3 %>% filter(sig_colors_3$sig_all == "Significant PISA and PROT")
PROTsignificant <-
  sig_colors_3 %>% filter(sig_colors_3$sig_all == "Significant PROT")

## Thermal Stability versus Protein Abundance changes graph
FC_PISA_plot_1 <- sig_colors_3 %>%
  arrange(desc(sig_PISA)) %>%
  ggplot(aes(
    x = log_prot,
    y = log_PISA,
    description = Gene,
    color = sig_all
  )) +
  geom_point(alpha = 0.7) +
  xlim(-6, 6) +
  geom_hline(yintercept = 0,
             linetype = "dotted",
             color = "black") +
  geom_vline(xintercept = 0,
             linetype = "dotted",
             color = "black") +
  labs(
    # title = "Thermal Stability Versus Protein Abundance Changes",
    title = "nPISA Versus Protein Abundance Changes",
    x = "Protein Abundance log2 (stationary/log)",
    y = "nPISA Intensity log2 (stationary/log)",
    color = "Significance"
  ) +
  theme_clean() +
  theme(
    # Arial is compatible with most browsers and OS's
    axis.text = element_text(size = 12, family = "arial"),
    axis.title = element_text(family = "arial"),
    panel.background = element_rect(fill = "#f5f5f5"),
    title = element_text(size = 20, family = "arial"),
    legend.text = element_text(size = 10, family = "arial"),
    legend.title = element_text(color = NA, family = "arial"),
    legend.background = element_rect(color = NA, fill = "white", size = 6),
    legend.spacing = unit(0.3, "cm"),
    legend.margin = margin(0.2, 0.2, 0.2, 0.2, "cm"),
    legend.key.size = unit(1, "cm")
  ) +
  scale_colour_manual(values = mukha_4cols)

### GLUCOSE STARVATION SECTION

# Get data
GLU_data_stats_sig_adj <-
  read_csv("ronnie_293T_data_norm.csv", col_names = TRUE)

# Transpose matrix and select rows that include LOG and STAT
GLUinvertedProteome <- as.data.frame(t(GLU_data_stats_sig_adj))
names(GLUinvertedProteome) <- GLUinvertedProteome[3,]
GLUinvertedProteome <- GLUinvertedProteome[-1:-3,]
GLUTMTs <-
  GLUinvertedProteome[str_which(rownames(GLUinvertedProteome), "\\wNoGlu|\\wGlu"),]

# Get data
GLU_p_data_stats_sig_adj <-
  read_csv("ronnie_293T_data_norm.csv", col_names = TRUE)

mukha_4cols <- c("gray", "#440154FF", "#55C667FF", "#FDE725FF")
arial <- list(family = "arial",
              size = 12)

GLU_low_cutoff <- -0.5

sig_colors_2 <- GLU_p_data_stats_sig_adj

sig_colors_2 %>%
  mutate(
    sig_PISA = ifelse(
      log_adj_p_PISA < 0.05 &
        (log_PISA < GLU_low_cutoff |
           log_PISA > 1),
      "1",
      "0"
    ),
    sig_PROT = ifelse(
      log_adj_p_PROT < 0.05 &
        (log_PROT < GLU_low_cutoff |
           log_PROT > 1),
      "2",
      "0"
    )
  )
sig_colors_2$sig_PISA <- as.numeric(sig_colors_2$sig_PISA)
sig_colors_2$sig_PROT <- as.numeric(sig_colors_2$sig_PROT)
sig_colors_2 <- sig_colors_2 %>%
  mutate(sig_all = sig_PISA + sig_PROT)
sig_colors_2$sig_all <- as.character(sig_colors_2$sig_all)

sig_colors_2$sig_all <-
  str_replace(sig_colors_2$sig_all, "1", "Significant PISA")
sig_colors_2$sig_all <-
  str_replace(sig_colors_2$sig_all, "2", "Significant PROT")
sig_colors_2$sig_all <-
  str_replace(sig_colors_2$sig_all, "3", "Significant PISA and PROT")
sig_colors_2$sig_all <-
  str_replace(sig_colors_2$sig_all, "0", "Not significant")
sig_colors_2$sig_all[is.na(sig_colors_2$sig_all)] = "Not significant"
GLUnotSignificant <-
  sig_colors_2 %>% filter(sig_colors_2$sig_all == "Not significant")

GLUPISAsignificant <-
  sig_colors_2 %>% filter(sig_colors_2$sig_all == "Significant PISA")
GLUBOTHsignificant <-
  sig_colors_2 %>% filter(sig_colors_2$sig_all == "Significant PISA and PROT")
GLUPROTsignificant <-
  sig_colors_2 %>% filter(sig_colors_2$sig_all == "Significant PROT")

label_cutoff_1 <- -0.5

FC_PISA_norm_vs_PROT <- sig_colors_2 %>%
  ggplot(aes(
    x = log_PROT,
    y = log_PISA,
    description = Gene,
    color = sig_all
  )) +
  geom_point(alpha = 0.7) +
  geom_hline(yintercept = 0,
             linetype = "dotted",
             color = "black") +
  geom_vline(xintercept = 0,
             linetype = "dotted",
             color = "black") +
  xlim(-4, 4) +
  ylim(-4, 4) +
  labs(
    title = "293T Glucose Starvation Prot vs PISA_norm Fold Changes Plot",
    x = "log2 fold change of PROT",
    y = "log2 fold change of PISA_norm",
    color = "Significance"
  ) +
  theme_clean() +
  theme(
    # axis.text=element_text(size=14),
    axis.text = element_text(size = 12, family = "arial"),
    axis.title = element_text(family = "arial"),
    panel.background = element_rect(fill = "#f5f5f5"),
    title = element_text(size = 20, family = "arial"),
    legend.text = element_text(size = 10, family = "arial"),
    legend.title = element_text(color = NA, family = "arial"),
    legend.background = element_rect(color = NA, fill = "white", size = 6),
    legend.spacing = unit(0.3, "cm"),
    legend.margin = margin(0.2, 0.2, 0.2, 0.2, "cm"),
    legend.key.size = unit(1, "cm")
  ) +
  scale_colour_manual(values = mukha_4cols)

### GLUCOSE STARVATION ENDS

# UI
ui <- dashboardPage(
  # Change to (Protein) Thermal Stability Viewer
  dashboardHeader(title = "Protein Thermostability Viewer"),
  dashboardSidebar(
    tags$head(tags$style(
      HTML(
        '.main-header > .navbar {        margin-left: 300px;      }      .main-header .logo {         width: 300px;   font-family: Source Sans Pro, sans-serif; font-size: 19px;     font-weight: 500;  }      .content-wrapper, .main-footer, .right-side {          margin-left: 300px;      }    '
      )
    )),
    width = 300,
    sidebarMenu(
      menuItem(
        "Discover E. coli",
        tabName = "LOGLOG",
        icon = icon("th", lib = "glyphicon")
      ),
      menuItem(
        "Browse E. coli Proteome",
        tabName = "Search",
        icon = icon("globe", lib = "glyphicon")
      ),
      menuItem(
        "Discover 293T Glucose Starvation",
        tabName = "GlucoseStarvation",
        icon = icon("apple", lib = "glyphicon")
      ),
      menuItem(
        "Browse 293T Glucose Starvation Proteome",
        tabName = "GLUSearch",
        icon = icon("globe", lib = "glyphicon")
      )
    )
  ),
  dashboardBody(tabItems(
    # Boxes need to be put in a fluid row (or column)
    # For E. COLI FIGURE tab
    tabItem(tabName = "LOGLOG",
            fluidRow(
              box(
                title = "Sultonova et al. (2022), Figure 2A)",
                background = "light-blue",
                "Comparison of thermally shifted proteins to their protein abundances in logarithmic and stationary phase E.coli."
              )
            ),
            fluidRow(
              box(title = "Comparison of mean log2 normalized nPISA ratios to mean log2 ratios of protein abundance. Proteins significantly different by nPISA, abundance, or both (Benjamini-Hochberg-adjusted p-value <0.05 and log2 fold change <-1 or > 1) appear as indicated. Click with extreme prejudice.", status = "primary", plotlyOutput("LOGLOG")),
              box(
                title = "Intensities",
                status = "primary",
                plotlyOutput("LOGLOGPlot")
              )
            )),
    #For E. COLI SEARCH
    tabItem(tabName = "Search",
            fluidRow(
              box(
                title = "Sultonova et al. (2022)",
                background = "light-blue",
                "Protein abundances in logarithmic and stationary phase E.coli."
              )
            ),
            fluidRow(box(
              status = "primary",
              selectInput(
                'query',
                'Search proteome by gene symbol',
                choices = names(invertedProteome)[-1]
              )
            )),
            fluidRow(box(
              status = "primary",
              plotlyOutput("foundation")
            ))),
    #FOR GLUCOSE STARVATION
    tabItem(tabName = "GlucoseStarvation",
            fluidRow(
              box(
                title = "Sultonova et al. (2022), Figure 5A)",
                background = "light-blue",
                "Thermally-shifted, but stably-abundant, proteins between starved and unstarved Hek293T cells. "
              )
            ),
            fluidRow(
              box(title = "Scatter plot comparing nPISA versus mean log2 changes in protein abundance. Significant differences in nPISA, abundance, or both (Benjamini-Hochberg-adjusted p-value <0.05 and log2 fold change < -0.5 or > 0.5) are indicated as shown. Click with extreme prejudice.", status = "primary", plotlyOutput("Glucose")),
              box(title = "Intensities",
                  status = "primary",
                  plotlyOutput("GLUPlot"))
            )),
    #FOR GLUCOSE STARVATION SEARCH
    tabItem(tabName = "GLUSearch",
            fluidRow(
              box(
                title = "Sultonova et al. (2022)",
                background = "light-blue",
                "Thermally-shifted, but stably-abundant, proteins between starved and unstarved Hek293T cells."
              )
            ),
            fluidRow(box(
              status = "primary",
              selectInput(
                'GLUquery',
                'Search proteome by gene symbol',
                choices = names(GLUinvertedProteome)[-1]
              )
            )),
            fluidRow(box(
              status = "primary",
              plotlyOutput("GLUfoundation")
            )))
    
  ))
)


server <- function(input, output) {
  ## FOR E. COLI FIGURE tab
  
  LOGLOGtesterZ <- reactive({
    d <- event_data("plotly_click", source = "LOGLOG_gene")
    if (is.null(d))
      return("gnd")
    which <-
      str_remove_all(d["curveNumber"], "curveNumber 1\\s") %>% as.numeric()
    if (which == 0) {
      LOGLOGgeneIndex <-
        str_remove_all(d["pointNumber"], "pointNumber 1\\s") %>% as.numeric() + 1
      LOGLOGgeneSymbol <-
        str_remove_all(notSignificant[LOGLOGgeneIndex , "Gene"], "# A tibble: 1 × 1 Gene <chr> 1\\s")
    } else if (which == 1) {
      LOGLOGgeneIndex <-
        str_remove_all(d["pointNumber"], "pointNumber 1\\s") %>% as.numeric() + 1
      LOGLOGgeneSymbol <-
        str_remove_all(PISAsignificant[LOGLOGgeneIndex , "Gene"], "# A tibble: 1 × 1 Gene <chr> 1\\s")
    } else if (which == 3) {
      LOGLOGgeneIndex <-
        str_remove_all(d["pointNumber"], "pointNumber 1\\s") %>% as.numeric() + 1
      LOGLOGgeneSymbol <-
        str_remove_all(PROTsignificant[LOGLOGgeneIndex , "Gene"], "# A tibble: 1 × 1 Gene <chr> 1\\s")
    } else if (which == 2) {
      LOGLOGgeneIndex <-
        str_remove_all(d["pointNumber"], "pointNumber 1\\s") %>% as.numeric() + 1
      LOGLOGgeneSymbol <-
        str_remove_all(BOTHsignificant[LOGLOGgeneIndex , "Gene"], "# A tibble: 1 × 1 Gene <chr> 1\\s")
    } else {
      return(NULL)
    }
    LOGLOGgeneSymbol
  })
  
  LOGLOGfound <- reactive({
    searching <- dplyr::select(TMTs, as.character(LOGLOGtesterZ()))
    searching <- t(searching)
    searching
  })
  
  output$LOGLOGPlot <- renderPlotly({
    plot_ly(
      x = colnames(LOGLOGfound()),
      y = as.numeric(LOGLOGfound()),
      type = 'bar',
      color = str_extract(colnames(LOGLOGfound()), "[^_]+")
      # marker = list(
      #   color = c(
      #     '#78cdff',
      #     '#78cdff',
      #     '#78cdff',
      #     '#78cdff',
      #     '#3c8dbc',
      #     '#3c8dbc',
      #     '#3c8dbc',
      #     '#3c8dbc',
      #     '#828282',
      #     '#828282',
      #     '#828282',
      #     '#828282',
      #     '#333',
      #     '#333',
      #     '#333',
      #     '#333'
      #   )
      # )
    ) %>%
      layout(
        title = LOGLOGtesterZ(),
        xaxis = list(title = ""),
        yaxis = list(title = "Intensity"),
        # Remove font
        font = arial
      )
  })
  
  output$LOGLOG <- renderPlotly({
    ggplotly(FC_PISA_plot_1, source = "LOGLOG_gene")
  })
  
  ## END OF LOGLOG
  
  
  ### GLUCOSE STARVATION
  
  GLUtesterZ <- reactive({
    d <- event_data("plotly_click", source = "GlucoseGraph")
    if (is.null(d))
      return("GNL2")
    which <-
      str_remove_all(d["curveNumber"], "curveNumber 1\\s") %>% as.numeric()
    if (which == 0) {
      GLUgeneIndex <-
        str_remove_all(d["pointNumber"], "pointNumber 1\\s") %>% as.numeric() + 1
      GLUgeneSymbol <-
        str_remove_all(GLUnotSignificant[GLUgeneIndex , "Gene"], "# A tibble: 1 × 1 Gene <chr> 1\\s")
    } else if (which == 1) {
      GLUgeneIndex <-
        str_remove_all(d["pointNumber"], "pointNumber 1\\s") %>% as.numeric() + 1
      GLUgeneSymbol <-
        str_remove_all(GLUPISAsignificant[GLUgeneIndex , "Gene"], "# A tibble: 1 × 1 Gene <chr> 1\\s")
    } else if (which == 3) {
      GLUgeneIndex <-
        str_remove_all(d["pointNumber"], "pointNumber 1\\s") %>% as.numeric() + 1
      GLUgeneSymbol <-
        str_remove_all(GLUPROTsignificant[GLUgeneIndex , "Gene"], "# A tibble: 1 × 1 Gene <chr> 1\\s")
    } else if (which == 2) {
      GLUgeneIndex <-
        str_remove_all(d["pointNumber"], "pointNumber 1\\s") %>% as.numeric() + 1
      GLUgeneSymbol <-
        str_remove_all(GLUBOTHsignificant[GLUgeneIndex , "Gene"], "# A tibble: 1 × 1 Gene <chr> 1\\s")
    } else {
      return(NULL)
    }
    
    GLUgeneSymbol
  })
  
  GLUfound <- reactive({
    searching <- dplyr::select(GLUTMTs, as.character(GLUtesterZ()))
    searching <- t(searching)
    searching
  })
  
  output$GLUPlot <- renderPlotly({
    plot_ly(
      x = colnames(GLUfound()),
      y = as.numeric(GLUfound()),
      color = str_extract(colnames(GLUfound()), "[^_]+"),
      type = 'bar'
      # marker = list(
      #   color = c(
      #     '#78cdff',
      #     '#78cdff',
      #     '#78cdff',
      #     '#78cdff',
      #     '#78cdff',
      #     '#3c8dbc',
      #     '#3c8dbc',
      #     '#3c8dbc',
      #     '#3c8dbc',
      #     '#3c8dbc',
      #     '#828282',
      #     '#828282',
      #     '#828282',
      #     '#333',
      #     '#333',
      #     '#333'
      #   )
      # )
    ) %>%
      layout(
        title = GLUtesterZ(),
        xaxis = list(title = ""),
        yaxis = list(title = "Intensity"),
        font = arial
      )
  })
  
  output$Glucose <- renderPlotly({
    ggplotly(FC_PISA_norm_vs_PROT, source = "GlucoseGraph")
  })
  
  ### END OF GLUCOSE STARVATION
  
  ## FOR E. COLI SEARCH
  
  findation <- reactive({
    searching <- dplyr::select(TMTs, as.character(input$query))
    searching <- t(searching)
    searching
  })
  
  output$foundation <- renderPlotly({
    plot_ly(
      x = colnames(findation()),
      y = as.numeric(findation()),
      type = 'bar',
      color = str_extract(colnames(findation()), "[^_]+")
      # marker = list(
      #   color = c(
      #     '#78cdff',
      #     '#78cdff',
      #     '#78cdff',
      #     '#78cdff',
      #     '#3c8dbc',
      #     '#3c8dbc',
      #     '#3c8dbc',
      #     '#3c8dbc',
      #     '#828282',
      #     '#828282',
      #     '#828282',
      #     '#828282',
      #     '#333',
      #     '#333',
      #     '#333',
      #     '#333'
      #   )
      # )
    ) %>%
      layout(
        title = input$query,
        xaxis = list(title = ""),
        yaxis = list(title = "Intensity"),
        font = arial
      )
  })
  
  ## END OF E. COLI SEARCH
  
  ## GLUCOSE STARVATION SEARCH
  
  GLUfindation <- reactive({
    searching <- dplyr::select(GLUTMTs, as.character(input$GLUquery))
    searching <- t(searching)
    searching
  })
  
  output$GLUfoundation <- renderPlotly({
    plot_ly(
      x = str_remove(colnames(GLUfindation()), "TMT_|TMT"),
      y = as.numeric(GLUfindation()),
      type = 'bar',
      color = str_extract(colnames(GLUfindation()), "[^_]+")
      # marker = list(
      #   color = c(
      #     '#78cdff',
      #     '#78cdff',
      #     '#78cdff',
      #     '#78cdff',
      #     '#78cdff',
      #     '#3c8dbc',
      #     '#3c8dbc',
      #     '#3c8dbc',
      #     '#3c8dbc',
      #     '#3c8dbc',
      #     '#828282',
      #     '#828282',
      #     '#828282',
      #     '#333',
      #     '#333',
      #     '#333'
      #   )
      # )
    ) %>%
      layout(
        title = input$GLUquery,
        xaxis = list(title = ""),
        yaxis = list(title = "Intensity"),
        font = arial
      )
  })
  
  ## GLUCOSE STARVATION SEARCH ENDS
  
}

shinyApp(ui, server)