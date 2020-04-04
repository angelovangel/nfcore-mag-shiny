 #### nfcore mag app ####
 
 # a shiny frontend for the nf-core/mag pipeline
 # https://github.com/nf-core/mag.git
 
 library(shiny)
 library(shinyFiles)
 library(shinyjs)
 library(shinyalert)
 library(dplyr)
 library(processx)
 library(stringr)
 library(digest)
 library(yaml)
 library(shinyFeedback)
 library(pingr) # to check if server has internet
 
 # define reactive to track user counts
 users <- reactiveValues(count = 0)
 
 source("ncct_modal.R", local = FALSE)$value # don't share across sessions, who knows what could happen!
 source("ncct_make_yaml.R")
 
 
 #### ui ####
 ui <- function(x) {
   # I like to mix up R, JS and CSS
   navbarPage(title = tags$button("nf-core/mag",
                               id = "magButton", # events can be listened to as input$cazytableButton in server
                               class = "action-button", #shiny needs this class def to work
                               title = "If you want to start over, just reload the page.",
                               onMouseOver = "this.style.color='#843500'", # old school js
                               onMouseOut = "this.style.color='#D35400'",
                               style = "color: #D35400; font-weight: bold; border: none; background-color: inherit;"),
             
              windowTitle = "nf-core/mag", 
              collapsible = TRUE,
    
    tabPanel("nf-core/mag output",
            # attempts to use external progress bar
            includeCSS("css/custom.css"),
            useShinyFeedback(),
            useShinyjs(),
            useShinyalert(), 
            
            # snackbars begin
            snackbarWarning(id = "tower_snackbar", 
                            message = "Is TOWER_ACCESS_TOKEN available in Sys.getenv() ?"),
            snackbarSuccess("fastp_trimmed", 
                            message = "Default fastp parameters will be used"),
            # snackbars end
            
            shiny::uiOutput("mqc_report_button", inline = TRUE),
            shiny::uiOutput("where_are_my_results", inline = TRUE),
            
            shiny::div(id = "commands_pannel",
              selectizeInput("step1", label = NULL,
                             choices = list(
                               "Short reads only" = list("single_end", "paired_end"),
                               "Short + long reads" = list("hybrid"),
                               "Test run" = list("test", "test_hybrid")
                               ),
                             multiple = FALSE, selected = "paired_end"),
              
              shinyDirButton(id = "fastq_folder", 
                             label = "Select fastq folder", 
                             style = "color: #D35400; font-weight: bold;", 
                             title = "Please select a folder with fastq files", 
                             icon = icon("folder-open")),
              shinyFilesButton(id = "manifest_file", 
                               label = "Choose manifest file",
                               style = "color: #D35400; font-weight: bold;", 
                               title = "Select manifest file for hybrid data (see nf-core/mag docs)",
                               multiple = FALSE, 
                               icon = icon("file")),
              
              actionButton("run", "Run nf-core/mag pipeline", 
                         style = "color: #D35400; font-weight: bold;", 
                         onMouseOver = "this.style.color = '#843500' ", 
                         onMouseOut = "this.style.color = '#D35400' ", 
                         icon = icon("play")),
              actionButton("reset", "Reset", 
                         style = "color: #D35400; font-weight: bold;",
                         onMouseOver = "this.style.color = '#843500' ",
                         onMouseOut = "this.style.color = '#D35400' ", 
                         icon = icon("redo")),
            
            actionButton("more", "More options", 
                         icon = icon("cog"),
                         class = "rightAlign"),
            actionButton("landing_page", "Go to home page", 
                         icon = icon("home"),
                         class = "rightAlign", 
                         onclick ="window.open('http://google.com', '_blank')"),
            
            tags$div(id = "optional_inputs",
              absolutePanel(top = 180, right = 20,
                          textInput(inputId = "reads_pattern", 
                                    label = "Fastq reads pattern", 
                                    value = "*R{1,2}_001.fastq.gz"),
                          tags$hr(),
                          
                          selectizeInput("taxonomy", 
                                         label = "Taxonomic classification", 
                                         choices = c("none", "kraken2", "centrifuge"), 
                                         selected = "none",
                                         multiple = FALSE),
                          tags$hr(),
                          
                          selectizeInput("nxf_profile", 
                                         label = "Select nextflow profile", 
                                         choices = c("docker", "conda"),
                                         selected = "docker", 
                                         multiple = FALSE),
                          tags$hr(),
                          
                          actionButton("ncct", "Enter NCCT project info"),
                          tags$hr(),
                          
                          checkboxInput("tower", "Use Nextflow Tower to monitor run", value = FALSE),
                          tags$hr()
                          
              )
            )
          ),
            
            verbatimTextOutput("stdout")
            
    ),
    tabPanel("Help", 
             includeMarkdown("help.md"))
    
  )
 }
 #### server ####
  server <- function(input, output, session) {
    options(shiny.launch.browser = TRUE, shiny.error=recover)
    
    #----
    # reactive for optional params for nxf, 
    # like tower, optional multiqc config (all nf-core pipes take this), and if hybrid is used
    # set TOWER_ACCESS_TOKEN in ~/.Renviron
    optional_params <- reactiveValues(tower = "", mqc = "", taxonomy = "")
    
    # update user counts at each server call
    isolate({
      users$count <- users$count + 1
    })
    
    # observe changes in users$count and write to log, observers use eager eval
    observe({
      writeLines(as.character(users$count), con = "userlog")
    })
    
    # observer for optional inputs
    hide("optional_inputs")
    observeEvent(input$more, {
      shinyjs::toggle("optional_inputs")
    })
    
    
    # shinyFeeback observers
    
    observe({
      if(input$tower) {
      showSnackbar("tower_snackbar")
      }
    })
    
    
    #----
    # strategy for ncct modal and multiqc config file handling:
    # if input$ncct_ok is clicked, the modal inputs are fed into the ncct_make_yaml() function, which generates
    # a multiqc_config.yml file and saves it using tempfile()
    # initially, the reactive value mqc_config$rv is set to "", if input$ncct_ok then it is set to
    # c("--multiqc_config", mqc_config_temp) and this reactive is given as param to the nxf pipeline
    
    # observer to generate ncct modal
    observeEvent(input$ncct, {
      if(pingr::is_online()) {
        ncct_modal_entries <- yaml::yaml.load_file("https://gist.githubusercontent.com/angelovangel/d079296b184eba5b124c1d434276fa28/raw/ncct_modal_entries")
        showModal( ncct_modal(ncct_modal_entries) )
      } else {
        shinyalert("No internet!", 
                   text = "This feature requires internet connection", 
                   type = "warning")
      }
      
    })
    
    # generate yml file in case OK of modal was pressed
    # the yml file is generated in the app exec env, using temp()
    observeEvent(input$ncct_ok, {
      mqc_config_temp <- tempfile()
      optional_params$mqc <- c("--multiqc_config", mqc_config_temp) 
      ncct_make_yaml(customer = input$customer, 
                     project_id = input$project_id, 
                     ncct_contact = input$ncct_contact, 
                     project_type = input$project_type, 
                     lib_prep = input$lib_prep, 
                     indexing = input$indexing, 
                     seq_setup = input$seq_setup, 
                     ymlfile = mqc_config_temp)
      shinyalert(text = "Project info saved!", type = "info", timer = 1500, showConfirmButton = FALSE)
      removeModal()
    })
    
    
    # generate random hash for multiqc report temp file name
    mqc_hash <- sprintf("%s_%s.html", as.integer(Sys.time()), digest::digest(runif(1)) )
    
    # dir choose management --------------------------------------
    volumes <- c(Home = fs::path_home(), getVolumes()() )
    
    shinyDirChoose(input, "fastq_folder", 
                   roots = volumes, 
                   session = session, 
                   restrictions = system.file(package = "base")) 
    
    # file choose management - manifest
    shinyFileChoose(input, "manifest_file", 
                    roots = volumes, 
                    session = session)
    
    #-----------------------------------
    # The main work of setting args for the nxf call is done here
    # in case the reactive vals are "", then they are not used by nxf
    
    output$stdout <- renderPrint({
    
      # set optional parameters, valid for all CASEs
      optional_params$tower <- if(input$tower) {
        "-with-tower"
      } else {
        ""
      }
      
      optional_params$taxonomy <- case_when(
        input$taxonomy == "kraken2" ~ "--kraken2_db ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/minikraken2_v2_8GB_201904_UPDATE.tgz",
        input$taxonomy == "centrifuge" ~ "--centrifuge_db ftp://ftp.ccb.jhu.edu/pub/infphilo/centrifuge/data/p_compressed+h+v.tar.gz",
        input$taxonomy == "none" ~ "")
    
    # CASE1: -profile test or test_hybrid
      if (input$step1 == "test" | input$step1 == "test_hybrid") {
        
        # hide all that is not needed
        shinyjs::hide(id = "fastq_folder")
        shinyjs::hide(id = "manifest_file")
        shinyjs::hide(id = "reads_pattern")
        
        # wd <<- tempdir() # in case test run is made, just a volatile temp dir is needed to store mqc report etc
        wd <<- file.path(getwd(), "test_wd") # use temp dir in current wd 
        resultsdir <<- file.path(wd, "tests") # set to 'tests' by nf-core/mag configs
        # build nxf command
        nxf_args <<- c("run", "nf-core/mag", "-profile", paste(input$nxf_profile, input$step1, sep = ",") ) 
        
        cat(
          " When running with '-profile test' or 'test_hybrid' there is no need to select a fastq folder, just press run \n",
          "Nextflow command to be executed:\n\n",
          "nextflow", nxf_args 
        )
      
    # CASE2: hybrid data selected
      } else if (input$step1 == "hybrid") {
        
        shinyjs::hide(id = "fastq_folder")
        shinyjs::hide(id = "reads_pattern")
        shinyjs::show(id = "manifest_file")
        
        cat(" Please select a manifest file\n\n")
        
        optional_params$tower <- if(input$tower) {
          "-with-tower"
        } else {
          ""
        }
        # set wd to the dir where the manifest file is
        wd <<- fs::path_dir( parseFilePaths(volumes, input$manifest_file)$datapath )
        resultsdir <<- file.path(wd, 'results')
        
        nxf_args <<- c("run" ,"nf-core/mag",
                       "--manifest", parseFilePaths(volumes, input$manifest_file)$datapath, 
                       "-profile", input$nxf_profile, 
                       optional_params$tower, 
                       optional_params$mqc)
        
        cat(" Nextflow command to be executed:\n\n",
              "nextflow", nxf_args)
          
    # CASE3: SE or PE selected
      } else if (input$step1 == "single_end" | input$step1 == "paired_end") {
        
        shinyjs::hide("manifest_file")
        shinyjs::show("fastq_folder")
        shinyjs::show("reads_pattern")
        
        if(is.numeric(input$fastq_folder)) {
          cat("Please select a fastq folder")
        
        } else {
          nfastq <<- length(list.files(path = parseDirPath(volumes, input$fastq_folder), pattern = "*fast(q|q.gz)$"))
          reads <<- file.path(parseDirPath(volumes, input$fastq_folder),input$reads_pattern)
          wd <<- parseDirPath(volumes, input$fastq_folder) # set wd to where the fastq file are
          resultsdir <<- file.path(wd, 'results')
          
          
          optional_params$tower <- if(input$tower) {
            "-with-tower"
          } else {
            ""
          }
          
          nxf_args <<- c("run", "nf-core/mag", 
                       "--reads", reads, 
                       "-profile", input$nxf_profile,
                       optional_params$taxonomy,
                       optional_params$tower,
                       optional_params$mqc)
            
          cat("Number of fastq files found:\n",
              nfastq, "\n",
              "------------------\n\n",
          
              "Nextflow command to be executed:\n",
              "nextflow", nxf_args, "\n",
              "------------------\n"
            )
        }
      }
    })

    #---
    # real call to nf-core/mag-------
    #----      
    # setup progress bar and callback function to update it
    progress <- shiny::Progress$new(min = 0, max = 1, style = "old")
    
    # callback function, to be called from run() on each chunk of output
    cb_count <- function(chunk, process) {
      counts <- str_count(chunk, pattern = "executor >")
      progress$inc(amount = counts/50)
    }
    
    # using processx to better control stdout
    observeEvent(input$run, {
      # check for consistency
      # important that brackets: T | F & F --> T while (T | F) & F --> F
      if( (input$step1 == "paired_end" | input$step1 == "single_end") & is.integer(input$fastq_folder) ) {
        shinyjs::html(id = "stdout", "\nPlease first select a folder containing the fastq files to be analysed, then press 'Run'...\n", add = FALSE)
      
      } else if ( input$step1 == "hybrid" & is.integer(input$manifest_file) ) {
        shinyjs::html(id = "stdout", "\nPlease select manifest file first", add = FALSE)
      
      } else {
        # set run button color to red?
        shinyjs::disable(id = "commands_pannel")
       
         # change label during run
        shinyjs::html(id = "run", html = "Running... please wait")
        progress$set(message = "Running... please wait ", value = 0)
        on.exit(progress$close() )
        
      # Dean Attali's solution
      # https://stackoverflow.com/a/30490698/8040734
        withCallingHandlers({
          
          shinyjs::html(id = "stdout", "")
          p <- processx::run("nextflow", 
                      args = nxf_args,
                      wd = wd, # wd is hard set in the renderPrint call for stdout above
                      #echo_cmd = TRUE, echo = TRUE,
                      stdout_line_callback = function(line, proc) { message(line) }, 
                      stdout_callback = cb_count,
                      stderr_to_stdout = TRUE, 
                      error_on_status = FALSE
                      )
          }, 
            message = function(m) {
              shinyjs::html(id = "stdout", html = m$message, add = TRUE); 
              runjs("document.getElementById('stdout').scrollTo(0,1e9);") # scroll the page to bottom with each message, 1e9 is just a big number
            }
        )
        if(p$status == 0) {
          # run finished successfully
          shinyjs::enable("commands_pannel")
          shinyjs::html(id = "run", html = "Run nf-core/mag pipeline")
          
          # clean work dir in case run finished ok
          work_dir <- file.path(wd, "work")
          system2("rm", args = c("-rf", work_dir))
          cat("deleted", work_dir, "\n")
          
            
          # copy mqc to www/ to be able to open it, also use hash to enable multiple concurrent users
          # if '-profile test' then outdir is 'tests', otherwise 'results' (set by the nf-core/mag configs)
          mqc_report <- file.path(resultsdir, "MultiQC/multiqc_report.html")
          system2("cp", args = c(mqc_report, paste("www/", mqc_hash, sep = "")) )
          
          
          # render the new action buttons to show report and location of results
          output$mqc_report_button <- renderUI({
            actionButton("mqc", label = "MultiQC report", 
                         icon = icon("th"), 
                         onclick = sprintf("window.open('%s', '_blank')", mqc_hash)
            )
          })
          
          output$where_are_my_results <- renderUI({
            actionButton("results_location_button", 
                         label = "Where are my results?", 
                         icon = icon("question"))
          })
          
          #
          # build js callback string for shinyalert
          js_cb_string <- sprintf("function(x) { if (x == true) {window.open('%s') ;} } ", mqc_hash)
          
          shinyalert("Run finished!", type = "success", 
                   animation = "slide-from-bottom",
                   text = "Pipeline finished, check results folder", 
                   showCancelButton = TRUE, 
                   confirmButtonText = "Open report",
                   callbackJS = js_cb_string, 
                   #callbackR = function(x) { js$openmqc(mqc_url) }
                   )
        } else {
          shinyjs::html(id = "run", html = "Finished with errors")
          shinyjs::enable(id = "commands_pannel")
          shinyjs::disable(id = "run")
          shinyalert("Error!", type = "error", 
                     animation = "slide-from-bottom", 
                     text = "Pipeline finished with errors, press OK to reload the app and try again.", 
                     showCancelButton = TRUE, 
                     callbackJS = "function(x) { if (x == true) {history.go(0);} }"
                     )
        }
      }
      
    })
  #----
  # OBSERVERS
    
  observeEvent(input$results_location_button, {
    shinyalert(title = "The location of the results folder is:", 
               text = resultsdir, 
               type = "info", 
               showCancelButton = FALSE, 
               showConfirmButton = TRUE)
  })
  # ask to start over if title or reset clicked
  #----                     
  observeEvent(input$magButton, {
    shinyalert(title = "",
               type = "warning",
               text = "Start again or stay on page?", 
               html = TRUE, 
               confirmButtonText = "Start again", 
               showCancelButton = TRUE, 
               callbackJS = "function(x) { if (x == true) {history.go(0);} }" # restart app by reloading page
               )
  })
  observeEvent(input$reset, {
    shinyalert(title = "",
               type = "warning",
               text = "Start again or stay on page?", 
               html = TRUE, 
               confirmButtonText = "Start again", 
               showCancelButton = TRUE, 
               # actually, session$reload() as an R callback should also work
               callbackJS = "function(x) { if (x == true) {history.go(0);} }" # restart app by reloading page
      )
    })
  
  #------------------------------------------------------------
  session$onSessionEnded(function() {
    # delete own mqc from www, it is meant to be temp only 
    system2("rm", args = c("-rf", paste("www/", mqc_hash, sep = "")) )
    
    #user management
    isolate({
      users$count <- users$count - 1
      writeLines(as.character(users$count), con = "userlog")
    })
    
  })
   
    
 }
 
 shinyApp(ui, server)
 