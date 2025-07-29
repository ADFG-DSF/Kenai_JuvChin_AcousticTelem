library(shiny)

server <- shinyServer(function(input, output) {
  
  df2 <- reactive({
    df1 <- unfiltered
    
    if(input$f1=="None") df2 <- df1
    if(input$f1=="Event") df2 <- filter_tags_with_3_in_30s(df1, input$n_events, input$t_sec)
    if(input$f1=="Tag Code") df2 <- apply_tagcode_filter(df1, message=FALSE)
    if(input$f1=="Multipath") df2 <- apply_multipath_filter(df1, input$thresh)
    if(input$f1=="Interval") df2 <- filter_tags_with_consistent_3s_intervals(df1, input$period, input$delta)
    
    if(input$f2=="None") df3 <- df2
    if(input$f2=="Event") df3 <- filter_tags_with_3_in_30s(df2, input$n_events, input$t_sec)
    if(input$f2=="Tag Code") df3 <- apply_tagcode_filter(df2, message=FALSE)
    if(input$f2=="Multipath") df3 <- apply_multipath_filter(df2, input$thresh)
    if(input$f2=="Interval") df3 <- filter_tags_with_consistent_3s_intervals(df2, input$period, input$delta)

    if(input$f3=="None") df4 <- df3
    if(input$f3=="Event") df4 <- filter_tags_with_3_in_30s(df3, input$n_events, input$t_sec)
    if(input$f3=="Tag Code") df4 <- apply_tagcode_filter(df3, message=FALSE)
    if(input$f3=="Multipath") df4 <- apply_multipath_filter(df3, input$thresh)
    if(input$f3=="Interval") df4 <- filter_tags_with_consistent_3s_intervals(df3, input$period, input$delta)

    if(input$f4=="None") df5 <- df4
    if(input$f4=="Event") df5 <- filter_tags_with_3_in_30s(df4, input$n_events, input$t_sec)
    if(input$f4=="Tag Code") df5 <- apply_tagcode_filter(df4, message=FALSE)
    if(input$f4=="Multipath") df5 <- apply_multipath_filter(df4, input$thresh)
    if(input$f4=="Interval") df5 <- filter_tags_with_consistent_3s_intervals(df4, input$period, input$delta)

    return(df5)
  })
  
  output$thePlot <- renderPlot(
    # plot(df2()$DateTime[1:100])
    if(length(unique(df2()$TagCode)) < 100) {
      df2() %>%
        ggplot(aes(x=DateTime, y=TagCode, col=TagCode %in% tagcodes_to_keep)) + 
        geom_point()
    }
  )
  
  output$nrow_orig <- renderText(paste("Initial Rows:", nrow(unfiltered)))
  output$ntags_orig <- renderText(paste("Initial n tags:", length(unique(unfiltered$TagCode))))
  output$nrow_filtered <- renderText(paste("Filtered Rows:", nrow(df2())))
  output$ntags_filtered <- renderText(paste("Filtered n tags:", length(unique(df2()$TagCode))))
})
ui <- shinyUI(fluidPage(
  titlePanel("Add Title Here"),
  sidebarLayout(
    sidebarPanel(
      radioButtons("f1", "Filter 1", choices=c("None", "Event", "Multipath", "Tag Code", "Interval"), inline=TRUE),
      radioButtons("f2", "Filter 2", choices=c("None", "Event", "Multipath", "Tag Code", "Interval"), inline=TRUE),
      radioButtons("f3", "Filter 3", choices=c("None", "Event", "Multipath", "Tag Code", "Interval"), inline=TRUE),
      radioButtons("f4", "Filter 4", choices=c("None", "Event", "Multipath", "Tag Code", "Interval"), inline=TRUE),
      
      hr(),
      h5("Event filter:"),
      sliderInput("n_events", "Number of events:", min=1, max=100, value=3),
      sliderInput("t_sec", "in time interval (sec):", min=1, max=100, value=30),
      
      hr(),
      h5("Multipath filter:"),
      sliderInput("thresh", "Threshold:", min=0.01, max=1, value=0.3, step=0.01),
      
      hr(),
      h5("Interval filter:"),
      sliderInput("period", "Period between pings (sec):", min=1, max=100, value=3),
      sliderInput("delta", "tolerance (sec):", min=0, max=3, value=0.5, step=0.05)
      
    ),
    #beginning of main section
    mainPanel(
      textOutput("nrow_orig"),
      textOutput("ntags_orig"),
      textOutput("nrow_filtered"),
      textOutput("ntags_filtered"),
      plotOutput("thePlot", height="600px", width="700px")
    )
  )
))
shinyApp(ui = ui, server = server)
