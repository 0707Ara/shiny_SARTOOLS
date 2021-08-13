library(shinydashboard)
library(shiny)
library(SARTools)
library(Glimma)
library(biomaRt)
library(gplots)
library(RColorBrewer)
library(dplyr)
library(pheatmap)
library(countToFPKM)
library(biomaRt)

options(shiny.trace = FALSE)


header <- dashboardHeader(title = "RNA Seq")
# sidebar ----
sidebar <- dashboardSidebar(
  # we need to use special function instead of uiOutput
  sidebarMenu(
    menuItem(text = "Count", tabName = "sartools"),
    menuItem(text = "FPKM", tabName = "fpkm")
  )
)

# tab1 
upload_box <- shinydashboard::box(title = "Upload Count File",
                                  status = "info", solidHeader = TRUE, width = 12,
                                  fluidRow(
                                    column(3, h4(icon("upload"), "Upload"))
                                  ),
                                  fluidRow(
                                    column(12, fileInput('file1',
                                                         label = "Choose Count File",
                                                         multiple = TRUE,
                                                         buttonLabel = "Browse or Drop...",
                                                         placeholder = "multiple count files ")
                                           )          
                                  ),
                                  fluidRow(
                                    column(12, fileInput('tele_file', label = "test",
                                              multiple = TRUE,
                                              buttonLabel = "Browse or Drop...",
                                              placeholder = "(multiple) csv or zip"))
                                  )
)
select_box <- shinydashboard::box(title = "Select Control and experiment group",
                                  status = "primary", solidHeader = TRUE, width = 12,
                                  fluidRow(
                                    column(6,
                                           uiOutput("control_checkbox"),
                                           textInput("control_name", "Control group name"),
                                           textOutput("txt")),
                                    column(6, 
                                           uiOutput("experiment_checkbox"),
                                           textInput("experiment_name", "Experiment group name"),
                                           textOutput("txt2"))
                                  )
                                  
)

input_info_box <- shinydashboard::box(title = "Select Control and experiment group",
                                      status = "primary", solidHeader = TRUE, width = 12,
                                      fluidRow(
                                        column(6, 
                                               textInput("workingDirectory", "Set working directory path", "/home/arajo0707/Rat/analysis_rat1")),
                                        column(6, 
                                               textInput("projectName", "Set project name", "example: Rat_rnaSeq "))),
                                      fluidRow(
                                        column(6,
                                               textInput("countDirectory", "Set count directory", "/home/arajo0707/Rat/count")),
                                        column(6, 
                                               textInput("projectAuthor", "Set author name", "example: Kim "))),
                                      fluidRow(
                                        column(2, uiOutput("resultHTML")),
                                        column(2, offset = 0, actionButton("createTarget", "Create Target file")))
                                      
)
heatmap_slider <- shinydashboard::box(title = "Pheatmap",
                                      status = "primary", solidHeader = TRUE, width = 12,
                                      fluidRow(
                                        column(8,
                                        sliderInput("significant_gene_number", "Number of significant genes",
                                                    min = 5, max = 35, value = 20, step = 5,
                                                    width = "80%")),
                                        br(),
                                        br(),
                                        br(),
                                        plotOutput("ttm_heatmap",
                                                   width = "auto", height = "800px"
                                               
                                        ),
                                        dataTableOutput('check_logcpm')
                                        )
                                      
)
#tab2


fpkm_input_info_box <- shinydashboard::box(title = "Type the fpmk directory",
                                      status = "primary", solidHeader = TRUE, width = 12,
                                      fluidRow(
                                        column(6, 
                                               textInput("fpkmDirectory", "Set fpkm directory path", "/home/arajo0707/Rat/analysis_rat1/fpkm")),
                                              )
                                      
)
fpkm_heatmap <- shinydashboard::box(title = "Pheatmap",
                                      status = "primary", solidHeader = TRUE, width = 12,
                                      fluidRow(
                                        column(8,
                                               sliderInput("significant_fpkm_gene_number", "Number of significant genes",
                                                           min = 5, max = 35, value = 20, step = 5,
                                                           width = "80%")),
                                        br(),
                                        br(),
                                        br(),
                                        plotOutput("fpkm_heatmap",
                                                   width = "auto", height = "800px")
                                      )
                                      
)
# body ----
body <- dashboardBody(
  # includeCSS("www/styles.css"),
  
  tabItems(
    tabItem(tabName = "sartools", fluidRow(select_box,
                                           upload_box,
                                           input_info_box,
                                           heatmap_slider)),
    tabItem(tabName = "fpkm",
            fluidRow(
                     fpkm_input_info_box,
                     fpkm_heatmap))
  )  
)

ui <- dashboardPage(header, sidebar, body,skin = "green")

server <- function(input, output, session) {
  # output$file_name_checkbox <- renderPrint({
  #     input$file1$name
  #     checkboxGroupInput("datasetSelector","Specify the datasets to compare:", choices = fileOptions$currentOptions)
  #   })
  targets <-data.frame(matrix(ncol=3,nrow=100, dimnames=list(NULL, c("label","files","condition"))))
  logcpm <- NULL
  
  # filelist <- reactive({input$file1})
  output$contents <- renderTable({
    ##input$file1 structure > name, size, type, datapath.
    # req(input$file1)
    # ext <- tools::file_ext(file$datapath)
    # file <- input$file1
    # print(input)
    print('####################')
    # validate(need(ext == "count", "Please upload a csv file"))
    
    upload = list()
    # read.csv(i)
    # file <- input$file1
    # re_re <- read.csv(
    #   file[1, 'datapath'],
    #   header = input$header,
    #   sep = "\t"
    # )
    # print(re_re)
    # for(i in 1:length(file)){
    #   read.csv(
    #     file[i, 'datapath'],
    #     header = input$header,
    #     sep = "\t"
    #   )
    # }
    # print(upload[[1]])
    return(upload())
  })
  
  # filelist <- reactive ({
  #   filenamelist <- list()
  #   filenamelist <- input$file1$name
  #   print(filenamelist)
  #   return(filenamelist)
  # })
  
  # observeEvent(input$file1, {
  #   fileOptions$currentOptions = list.append(fileOptions$currentOptions, input$file1$datapath)
  # })
  # print(filelist())
  print("checkbox")
  output$control_checkbox<-renderUI({
    checkboxGroupInput("controlSelector","Specify the control datasets:", choices = input$file1$name)
  })
  output$value <- renderText({ input$control_name })
  output$experiment_checkbox<-renderUI({
    checkboxGroupInput("experimentSelector","Specify the experiment datasets:", choices = input$file1$name)
  })
  # output$value <- renderText({ input$experiment_name })
  
  output$txt <- renderText({
    icons <- paste(input$controlSelector, collapse = ", ")
    # paste("You chose",paste(input$controlSelector, collapse = ", "))
    paste("As a control ", paste(input$controlSelector, collapse = ", "))
  })
  output$txt2 <- renderText({
    icons <- paste(input$experimentSelector, collapse = ", ")
    # paste("You chose",paste(input$controlSelector, collapse = ", "))
    paste("As a experiment ", paste(input$experimentSelector, collapse = ", "))
  })
  # print(filelist)
  # print("targetFile")
  
  # targetFile <- reactiveValues(
  #   data.frame(matrix(ncol = length(filelist$name), nrow = 0)),
  #   data.frame(matrix(ncol=length(filelist$name),nrow=0, dimnames=list(NULL, c("label","files","condition"))))
  # )
  # print("dataframe")
  
  ################################################
  ###### Declare new edgeR for MDS function ######
  ################################################
  run.edgeR_New <- function(counts, target, varInt, condRef, batch=NULL, cpmCutoff=1, 
                            normalizationMethod="TMM", pAdjustMethod="BH", ...){
    
    # filtering: select features which contain at least 
    # minReplicates (smallest number of replicates) with
    # at least cpmCutoff counts per million
    minReplicates <- min(table(target[,varInt]))
    fcounts <- counts[rowSums(cpm(counts) >= cpmCutoff) >= minReplicates,]
    cat("Number of features discarded by the filtering:\n")
    cat(nrow(counts)-nrow(fcounts),"\n")
    
    # building dge object
    design <- formula(paste("~", ifelse(!is.null(batch), paste(batch,"+"), ""), varInt))
    dge <- DGEList(counts=fcounts, remove.zeros=TRUE)
    dge$design <- model.matrix(design, data=target)
    cat("\nDesign of the statistical model:\n")
    cat(paste(as.character(design),collapse=" "),"\n")					  
    
    
    # Annotation
    geneid <- rownames(dge)
    ensembl = useMart("ENSEMBL_MART_ENSEMBL")
    ensembl = useDataset("rnorvegicus_gene_ensembl", mart = ensembl)
    geneanno <- getBM(attributes = c('ensembl_gene_id','transcript_length', 'chromosome_name', 'external_gene_name', 'gene_biotype'), filters = 'ensembl_gene_id', values = geneid, mart = ensembl)
    geneanno <- geneanno[order(geneanno$ensembl_gene_id,-geneanno$transcript_length),]
    geneanno <- geneanno[!duplicated(geneanno$ensembl_gene_id),]
    colnames(geneanno) <- c('EnsGid', 'TLength', 'chr', 'GeneName', 'GeneType')
    dge$genes <- geneanno
    
    # normalization
    dge <- calcNormFactors(dge, method=normalizationMethod)
    cat("\nNormalization factors:\n")
    print(dge$samples$norm.factors)
    
    # estimating dispersions
    dge <- estimateGLMCommonDisp(dge, dge$design)
    dge <- estimateGLMTrendedDisp(dge, dge$design)
    dge <- estimateGLMTagwiseDisp(dge, dge$design)
    
    # statistical testing: perform all the comparisons between the levels of varInt
    fit <- glmFit(dge, dge$design, ...)
    cat(paste("Coefficients of the model:",paste(colnames(fit$design),collapse="  ")),"\n")
    colsToTest <- grep(varInt,colnames(fit$design))
    namesToTest <- paste0(gsub(varInt,"",colnames(fit$design)[colsToTest]),"_vs_",condRef)
    results <- list()
    lrtJ <- list()
    # testing coefficients individually (tests againts the reference level)
    for (i in 1:length(colsToTest)){
      cat(paste0("Comparison ",gsub("_"," ",namesToTest[i]),": testing coefficient ",colnames(fit$design)[colsToTest[i]]),"\n")
      lrt <- glmLRT(fit, coef=colsToTest[i])
      results[[namesToTest[i]]] <- topTags(lrt,n=nrow(dge$counts),adjust.method=pAdjustMethod,sort.by="none")$table
      lrtJ[[namesToTest[i]]] <- glmLRT(fit, coef=colsToTest[i])
    }
    # defining contrasts for the other comparisons (if applicable)
    if (length(colsToTest)>=2){
      colnames <- gsub(varInt,"",colnames(fit$design))
      for (comp in combn(length(colsToTest),2,simplify=FALSE)){ 
        contrast <- numeric(ncol(dge$design))
        contrast[colsToTest[comp[1:2]]] <- c(-1,1)
        namecomp <- paste0(colnames[colsToTest[comp[2]]],"_vs_",colnames[colsToTest[comp[1]]])
        cat(paste0("Comparison ",gsub("_"," ",namecomp),": testing contrast (",paste(contrast,collapse=", "),")"),"\n")
        lrt <- glmLRT(fit, contrast=contrast)
        results[[namecomp]] <- topTags(lrt,n=nrow(dge$counts),adjust.method=pAdjustMethod,sort.by="none")$table
        lrtJ[[namecomp]] <- glmLRT(fit, contrast=contrast)
      }
    }
    
    return(list(dge=dge,results=results, lrtJ=lrtJ))
  }
  ################################################
  ###### End of new edgeR for MDS function  ######
  ################################################
  
  
  observeEvent(input$createTarget, {
      ##input$file1 structure > name, size, type, datapath.
      req(input$file1)
      file <- input$file1
      # print(file)
      for(i in 1:nrow(file)){
        targets$label[i] <- strsplit(file$name[i], '[.]')[[1]][1]
        targets$files[i] <- file$name[i]
        
        if(file$name[i] %in% input$controlSelector){
          targets$condition[i] <- input$control_name
        }else{
          targets$condition[i] <- input$experiment_name
        }
      }
      # targetFile <- targetFile[!apply(is.na(data))]
      targets <- na.omit(targets)
      targets <<- na.omit(targets)
      write.csv(targets, file = paste0(input$workingDirectory,"targetFile3.txt"), sep = "\t",
                row.names = FALSE, col.names = TRUE)
      targetFile <- paste0(input$workingDirectory,"targetFiled.txt")
      print("########target first result#########")
      
      print(targets)
      # print(upload[[1]])
      
      ########################end of observation event function#############################
      # rm(list=ls())                          # remove all the objects from the R session
      workDir <- paste0(input$workingDirectory)      # working directory for the R session
      projectName <- paste0(input$projectName)       # name of the project
      author <- paste0(input$projectAuthor)          # author of the statistical analysis/report
      # targetFile <- targets                  # path to the design/target file
      rawDir <- paste0(input$countDirectory)         # path to the directory containing raw counts files(.count file afterHTseq)
      featuresToRemove <- c("alignment_not_unique", "ambiguous", "no_feature", "not_aligned", "too_low_aQual")       
      # names of the features to be removed (specific HTSeq-count information and rRNA for example) NULL if no feature to remove
      
      varInt <- "condition"                                # factor of interest
      condRef <- input$controlSelector                     # reference biological condition
      batch <- NULL                                        # blocking factor: NULL (default) or "batch" for example
      alpha <- 0.05                                        # threshold of statistical significance
      pAdjustMethod <- "BH"                                # p-value adjustment method: "BH" (default) or "BY"
      cpmCutoff <- 1                                       # counts-per-million cut-off to filter low counts
      gene.selection <- "pairwise"                         # selection of the features in MDSPlot
      normalizationMethod <- "TMM"                         # normalization method: "TMM" (default), "RLE" (DESeq) or "upperquartile"
      
      colors <- c("#f38400","#875692", "#be0032","#008856", "#0067a5")
      
      # checking parameters
      print('########  checking parameters  ############')
      checkParameters.edgeR(projectName=projectName,
                            author=author, targetFile=targetFile,
                            rawDir=rawDir, featuresToRemove=featuresToRemove,
                            varInt=varInt, condRef=condRef, batch=batch,
                            alpha=alpha, pAdjustMethod=pAdjustMethod,
                            cpmCutoff=cpmCutoff, gene.selection=gene.selection,
                            normalizationMethod=normalizationMethod, colors=colors)
      target <- targets
      # loading counts
      counts <- loadCountData(target=target, rawDir=rawDir,
                              featuresToRemove=featuresToRemove)
      fcounts <- counts[rowSums(cpm(counts) >= 1) >= 10,]
      indexx <- unique(target$condition)
      
      # description plots
      majSequences <- descriptionPlots(counts=counts, group=target[,varInt], col=colors)
      
      # edgeR analysis
      out.edgeR <- run.edgeR(counts=counts, target=target, varInt=varInt,
                             condRef=condRef, batch=batch, cpmCutoff=cpmCutoff,
                             normalizationMethod=normalizationMethod, pAdjustMethod=pAdjustMethod)  
      
      # MDS + clustering
      exploreCounts(object=out.edgeR$dge, group=target[,varInt],
                    gene.selection=gene.selection, col=colors)
      
      # summary of the analysis (boxplots, dispersions, export table, nDiffTotal, histograms, MA plot)
      summaryResults <- summarizeResults.edgeR(out.edgeR, group=target[,varInt],
                                               counts=counts, alpha=alpha, col=colors)
      
      # generating HTML report
      writeReport.edgeR(target=target, counts=counts, out.edgeR=out.edgeR,
                        summaryResults=summaryResults, majSequences=majSequences,
                        workDir=workDir, projectName=projectName, author=author,
                        targetFile=targetFile, rawDir=rawDir, featuresToRemove=featuresToRemove,
                        varInt=varInt, condRef=condRef, batch=batch, alpha=alpha,
                        pAdjustMethod=pAdjustMethod, cpmCutoff=cpmCutoff,
                        colors=colors, gene.selection=gene.selection,
                        normalizationMethod=normalizationMethod)  
    
    print(input$projectName)
    output$resultHTML<-renderUI({
      workDirectoryMod <- ""
      for (i in 4:length(strsplit(input$workingDirectory, '[/]')[[1]])){
        workDirectoryMod <- paste0(workDirectoryMod,strsplit(input$workingDirectory, '[/]')[[1]][i],"/")
      }
      a("Result HTML",
        target = "_blank",
        href = paste0("/files/",workDirectoryMod,input$projectName,"_report.html"),
        style = "text-decoration: underline;")
    })
    
    aJSH <- run.edgeR_New(counts=counts, target=target, varInt=varInt, condRef=condRef,
                          batch=batch, cpmCutoff=cpmCutoff, normalizationMethod=normalizationMethod,
                          pAdjustMethod=pAdjustMethod)
    
    result <- names(aJSH$results)
    result.filtered <- aJSH$results
    
    
    #filter genews where pvalue >0.05 && logFC>1 
    for (i in result){
      result.filtered[[i]] <- result.filtered[[i]] %>% filter(abs(logFC)>1 & FDR<0.05)
      result.filtered[[i]] <- merge(x = result.filtered[[i]], y = aJSH$dge$genes[,c(1,4)], by.x = "row.names", by.y = "EnsGid" , all.x = TRUE) 
      result.filtered[[i]] <- result.filtered[[i]] %>% select(c(2,5,7,11)) %>% arrange(desc(abs(logFC)))  #selecting specific columns using index; can be replace by column name
    }
    
    
    
  
    #Heatmap 
    logcpm <- cpm(aJSH$dge,prior.count=2,log=TRUE)
    rownames(logcpm) <- aJSH$dge$genes$GeneName
    colnames(logcpm) <- target$label #sample_Name
    
    logcpm <- logcpm[result.filtered[[1]][1:35, "GeneName.x"],]
    logcpm <<- t(scale(t(logcpm)))
    # print out glMDPlot
    ####################### FIGURE OUT HOW TO ACCESS THESE FILES ########################

    
    output$check_logcpm <- renderTable(logcpm)
    
    output$ttm_heatmap <- renderPlot({
      if(!is.null(logcpm)){
        pheatmap(as.matrix(logcpm[1:input$significant_gene_number, ]), color = colorRampPalette(c('#2471A3','white','#C0392B'))(50), fontsize=15, main=paste("AE versus OH (sig n=389)"),  border_color = 'white', show_rownames = T, cutree_cols = 2, cluster_row = T, scale = "row")
      }
      
    })
    
  }) ### observe function of create HTML FIle button ###
  observeEvent(input$significant_gene_number, {
    output$ttm_heatmap <- renderPlot({
      if(!is.null(logcpm)){
        pheatmap(as.matrix(logcpm[1:input$significant_gene_number, ]), color = colorRampPalette(c('#2471A3','white','#C0392B'))(50), fontsize=15, main=paste("AE versus OH (sig n=389)"),  border_color = 'white', show_rownames = T, cutree_cols = 2, cluster_row = T, scale = "row")
      }
      
    })
    
  })
  
}

shinyApp(ui, server)