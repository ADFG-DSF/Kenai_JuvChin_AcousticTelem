library(shiny)

server <- shinyServer(function(input, output) {
  
  options(shiny.maxRequestSize=500*1024^2)
  
  printsize <- function(df) {
    paste0(format(nrow(df), big.mark = ","),
          " rows (", format(n_indiv(df), big.mark = ","), " tags)")
  }
  
  df_in <- reactive({
    req(input$file1)
    clean_but_unfiltered_data(read_clean_vroom(input$file1$datapath))
  })
  
  df1 <- reactive({
    apply_prefilter(df_in(), n=input$n_prefilter)
  })
  
  df2 <- reactive({
    if(input$f1=="None") df2 <- df1()
    if(input$f1=="Event") df2 <- filter_tags_with_3_in_30s_v2(df1(), input$n_events, input$t_sec)
    if(input$f1=="Tag Code") df2 <- apply_tagcode_filter(df1(), message=FALSE)
    if(input$f1=="Multipath") df2 <- apply_multipath_filter(df1(), input$thresh)
    if(input$f1=="Interval") df2 <- filter_tags_with_consistent_3s_intervals(df1(), input$period, input$delta)
    df2
  })
  
  df3 <- reactive({
    if(input$f2=="None") df3 <- df2()
    if(input$f2=="Event") df3 <- filter_tags_with_3_in_30s_v2(df2(), input$n_events, input$t_sec)
    if(input$f2=="Tag Code") df3 <- apply_tagcode_filter(df2(), message=FALSE)
    if(input$f2=="Multipath") df3 <- apply_multipath_filter(df2(), input$thresh)
    if(input$f2=="Interval") df3 <- filter_tags_with_consistent_3s_intervals(df2(), input$period, input$delta)
    df3
  })
  
  df4 <- reactive({
    if(input$f3=="None") df4 <- df3()
    if(input$f3=="Event") df4 <- filter_tags_with_3_in_30s_v2(df3(), input$n_events, input$t_sec)
    if(input$f3=="Tag Code") df4 <- apply_tagcode_filter(df3(), message=FALSE)
    if(input$f3=="Multipath") df4 <- apply_multipath_filter(df3(), input$thresh)
    if(input$f3=="Interval") df4 <- filter_tags_with_consistent_3s_intervals(df3(), input$period, input$delta)
    df4
  })
  
  df_out <- reactive({
    if(input$f4=="None") df5 <- df4()
    if(input$f4=="Event") df5 <- filter_tags_with_3_in_30s_v2(df4(), input$n_events, input$t_sec)
    if(input$f4=="Tag Code") df5 <- apply_tagcode_filter(df4(), message=FALSE)
    if(input$f4=="Multipath") df5 <- apply_multipath_filter(df4(), input$thresh)
    if(input$f4=="Interval") df5 <- filter_tags_with_consistent_3s_intervals(df4(), input$period, input$delta)
    return(df5)
  })
  
  output$thePlot <- renderPlot(
    if(length(unique(df_out()$TagCode)) <= 150) {
      df_out() %>%
        mutate(`Tags in library` = factor(ifelse(TagCode %in% tagcodes_to_keep, "Yes","No"),
               levels = c("Yes", "No"))) %>%
        ggplot(aes(x=DateTime, y=TagCode, col=`Tags in library`)) + 
        geom_point(show.legend = TRUE) +
        scale_color_discrete(drop=FALSE) +
        theme_bw()
    }
  )
  
  output$size_in <- renderText(paste(printsize(df_in()), "   -- -- -- Initial"))
  output$size1 <- renderText(paste(printsize(df1()), "   -- -- -- after Prefilter"))
  output$size2 <- renderText(paste(printsize(df2()), "   -- -- -- after", input$f1))
  output$size3 <- renderText(paste(printsize(df3()), "   -- -- -- after", input$f2))
  output$size4 <- renderText(paste(printsize(df4()), "   -- -- -- after", input$f3))
  output$size_out <- renderText(paste(printsize(df_out()), "   -- -- -- after", input$f4))
})
ui <- shinyUI(fluidPage(
  titlePanel("Add Title Here"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file1", "Choose .csv File",
                accept = c(".csv")),
      sliderInput("n_prefilter", "Prefilter n (tags with at least this many hits)",
                  min=2, max=100, value=2),
      
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
      textOutput("size_in"),
      textOutput("size1"),
      textOutput("size2"),
      textOutput("size3"),
      textOutput("size4"),
      textOutput("size_out"),
      
      plotOutput("thePlot", height="600px", width="700px")
    )
  )
))
shinyApp(ui = ui, server = server)
