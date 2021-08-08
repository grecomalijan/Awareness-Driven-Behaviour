library(shiny)
library(shinydashboard)
library(tidyverse)
library(deSolve)
library(zoo)

#### Data Wrangling ####
datacases <- read_csv("PHcases.csv")%>%
    mutate(date=as.Date(DateRepConf,format="%m/%e/%y")) %>%
    select(date, dailyCases=Count) %>%
    drop_na() %>%
    arrange(date)

datadeaths <- read_csv("PHdeaths.csv") %>%
    mutate(date=as.Date(DateDied,format="%m/%e/%y")) %>%
    select(date, dailyDeaths=Count) %>%
    drop_na() %>%
    arrange(date)

dataSI <- read_csv("PHSI.csv") %>%
    mutate(date=as.Date(Day, format="%m/%e%y")) %>%
    select(date, stringency_index) %>%
    drop_na() %>%
    arrange(date)

dateEpoch <- min(min(datacases$date), min(datadeaths$date))

data <- datacases %>%
    left_join(datadeaths, by='date') %>%
    complete(date=seq.Date(first(date), last(date), by='day'), fill=list(dailyCases=0, dailyDeaths=0)) %>%
    mutate(time=as.numeric(date-dateEpoch)) %>%
    mutate(rollingaverage.dataCases=(rollmean(dailyCases, k=7, fill=NA))) %>%
    mutate(rollingaverage.dataDeaths=(rollmean(dailyDeaths, k=7, fill=NA))) %>%
    left_join(dataSI, by='date')

str_ind <- approxfun(data$stringency_index, rule=2)
str_ind(seq(from=0, to=300, by=1))

#### Model functions ####
ratesFun <- function(t, state, parameters) {
    with(as.list(c(state, parameters)),
         {
             dday <- gammaH*H
             N <- S+E+I+R+H
             SI <- str_ind(t)
             ld <- (1-NPI*(SI/100))
             lambda <- (beta*ld*I/N)/((1+(dday/dcrit)^k)+(D/dtotcrit)^k)
             dS <- -lambda*S
             dE <- lambda*S - mu*E
             dI <- mu*E - gamma*I
             dR <- gamma*I*(1-fd)
             dH <- fd*gamma*I-dday
             dD <- dday
             return(list(c(dS, dE, dI, dR, dH, dD)))
         }
    )
}

runModel <- function(state, ratesFun, parms, times) {
    out <- ode(y = state, times = times, func = ratesFun, parms = parms)
    out
}

#### Shiny UI ####
ui <- fluidPage(
    titlePanel(title=div("Model Exploration - Philippines", img(src="MyImage.png", align = "right", height = "25%", width = "25%"))),
    h6("This is a tool developed in support of a dissertation submitted in partial fulfillment of the requirements for the degree MSc International Health and Tropical Medicine 
            in the University of Oxford, Centre for Tropical Medicine and Global Health. The research placement was conducted in partnership with the Modelling and Simulation Hub, Africa."),
    sidebarLayout(
        sidebarPanel(
            h3("Transmission and Reporting"),
            sliderInput("parm_R0",
                        "R0",
                        min = 0.1,
                        max = 5,
                        value = 2.5,
                        step = 0.1),
            sliderInput("parm_NPI",
                        "NPI",
                        min = 0,
                        max = 1,
                        value = 0.5,
                        step = 0.05),
            sliderInput("parm_report",
                        "report",
                        min = 0,
                        max = 1,
                        value = 0.04,  
                        step = 0.01),
            sliderInput("parm_pa",
                        "pa = Proportion of asymptomatic cases",  
                        min = 0,
                        max = 1,
                        value = 0.3, 
                        step = 0.01),
            h3("Awareness-driven behavior"),
            sliderInput("parm_k",
                        "k",
                        min = 1,
                        max = 10,
                        value = 2.5,
                        step = 0.1),
            sliderInput("parm_dcrit",
                        "dcrit",
                        min = 0,
                        max = 1000,
                        value = 10,
                        step = 10),
            sliderInput("parm_dtotcrit",
                        "dtotcrit",
                        min = 500,
                        max = 20000,
                        value = 10000,
                        step = 100),
            checkboxInput("parm_lt",
                          "Include dtotcrit (Long Term Awareness)?",
                          value = FALSE),
            checkboxInput("parm_st",
                          "Include dcrit (Short Term Awareness)?",
                          value = TRUE),
            h3("Other params"),
            numericInput("initE",
                         "Initial Size of Infection",
                         min = 1,
                         max = 20000,
                         value = 16177),
            sliderInput("parm_inv_mu",
                        "1/mu = Mean latent period (days)",
                        min = 1,
                        max = 20,
                        value = 7,
                        step = 1),
            sliderInput("parm_inv_gamma",
                        "1/gamma = Mean infectious period (days)",
                        min = 1,
                        max = 20,
                        value = 7, 
                        step = 1),
            sliderInput("parm_inv_gammaH",
                        "1/gammaH = Mean time in a hospital before fatality (days)",
                        min = 1,
                        max = 30,
                        value = 13, 
                        step = 1),
            sliderInput("parm_fd",
                        "fd = Infection fatality probability",
                        min = 0,
                        max = 0.5,
                        value = 0.017,
                        step = 0.001)
        ),
        mainPanel(
            tabsetPanel(type='tabs',
                        tabPanel("Fit Deaths",
                                 plotOutput("plotDeathsDaily"),
                                 plotOutput("plotAwareness"),
                                 plotOutput("plotDeathsTotal"),
                                 textOutput("R0"),
                                 plotOutput("plotSI")),
                        tabPanel("Fit Cases",
                                 plotOutput("plotCases")),
                        tabPanel("Plot Compartments",
                                 plotOutput("plotCompartments")),
                        tabPanel("Initial Conditions",
                                 numericInput("initS", "S", 1.108e8),
                                 numericInput("initI", "I", 0),
                                 numericInput("initR", "R", 0),
                                 numericInput("initH", "H", 0),
                                 textOutput("initN")),
                        tabPanel("About This App",
                                 includeMarkdown("README.md")))
        )
    )
)

#### Shiny Server ####
server <- function(input, output) {
    state = reactive({
        c(S = input$initS, E = input$initE, I = input$initI, R = input$initR, H = input$initH, D = 0)
    })
    parms <- reactive({
        c(beta=(input$parm_R0)*(1/input$parm_inv_gamma),
          NPI=input$parm_NPI,
          mu=1/input$parm_inv_mu,
          gamma=1/input$parm_inv_gamma, 
          gammaH=1/input$parm_inv_gammaH,
          fd=input$parm_fd,
          report=input$parm_report,
          # Parameters we can vary in the app:
          dcrit=if(input$parm_st==FALSE){Inf}else{input$parm_dcrit},
          dtotcrit=if(input$parm_lt==FALSE){Inf}else{input$parm_dtotcrit}, # If select not to use totcrit, use Inf instead
          k=input$parm_k,
          pa=input$parm_pa
        )
    })
    times <- reactive({
        data$time[0:300]
    })
    modelOutput <- reactive({
        tbModel <- runModel(
            state(), ratesFun, parms(),
            times = times()
        )
        with(as.list(parms()),{
            tbModel[,] %>%
                as_tibble() %>%
                mutate(modelReportedCases=mu*E*(1-pa)*report) %>%
                left_join(data%>%rename(dataReportedCases=dailyCases), by='time') %>%
                pivot_longer(!c(time,date), names_to = "compartment")
        })
    })
    output$R0 <- renderText({
        beta <- parms()['beta']
        gamma <- parms()['gamma']
        R0 <- beta/gamma
        glue::glue("R0: {R0}")
    })
    output$plotSI <- renderPlot({
        modelOutput() %>%
            filter(compartment%>%str_detect("SI")) %>%
            ggplot(aes(x=date, y=value, color=compartment)) +
            geom_line()
    })
    output$initN <- renderText({
        glue::glue("Population Size: {sum(state())}")
    })
    output$plotCases <- renderPlot({
        modelOutput() %>%
            filter(compartment%>%str_detect("Cases")) %>%
            ggplot(aes(x=date, y=value, color=compartment)) +
            geom_line()
    })
    tbAwareness <- reactive({
        dcrit<-parms()['dcrit']
        dtotcrit<-parms()['dtotcrit']
        k<-parms()['k']
        modelOutput() %>%
            filter(str_detect(compartment, 'D')) %>%
            pivot_wider(c(time, date), names_from = compartment, values_from = value) %>%
            rename(dataDeaths=dailyDeaths) %>%
            mutate(modelDeaths=(D-lag(D)),
                   dataTotalDeaths=cumsum(dataDeaths),
                   awareness=1/(1+(modelDeaths/dcrit)^k+(D/dtotcrit)^k)) %>%
            rename(modelTotalDeaths=D) %>%
            pivot_longer(!c(time,date)) %>%
            drop_na()
    })
    output$plotDeathsDaily <- renderPlot({
        tbAwareness() %>%
            filter(!str_detect(name,'awareness'),
                   !str_detect(name,'Total')) %>%
            ggplot(aes(x=date, y=value, color=name)) +
            geom_line()
    })
    output$plotDeathsTotal <- renderPlot({
        tbAwareness() %>%
            filter(!str_detect(name,'awareness'),
                   str_detect(name,'Total')) %>%
            ggplot(aes(x=date, y=value, color=name)) +
            geom_line()
    })
    output$plotAwareness <- renderPlot({
        tbAwareness() %>%
            filter(str_detect(name,'awareness')) %>%
            ggplot(aes(x=date, y=value, color=name)) +
            geom_line()
    })
    output$plotCompartments <- renderPlot({
        modelOutput() %>%
            filter(str_detect('EIRDHL',compartment)) %>%
            mutate(compartment=factor(compartment, levels=c('S','E','I','R','H','D','L'))) %>%
            ggplot(aes(x=time, y=value, color=compartment)) +
            geom_line()
    })
    output$plotSI <- renderPlot({
        modelOutput() %>%
            filter(compartment%>%str_detect("stringency_index")) %>%
            ggplot(aes(x=date, y=value, color=compartment)) +
            geom_line()
    })
}

# Run the application
shinyApp(ui = ui, server = server)