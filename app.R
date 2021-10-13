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

options(shiny.trace = FALSE)


header <- dashboardHeader(title = "RNA Seq")

# sidebar ----
sidebar <- dashboardSidebar(
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
                                  )
)
select_box <- shinydashboard::box(title = "Select Control and Case group",
                                  status = "primary", solidHeader = TRUE, width = 12,
                                  fluidRow(
                                    column(6,
                                           uiOutput("control_checkbox"),
                                           textInput("control_name", "Control group name"),
                                           textOutput("txt")),
                                    column(6, 
                                           uiOutput("case_checkbox"),
                                           textInput("case_name", "Case group name"),
                                           textOutput("txt2"))
                                  )
                                  
)

input_info_box <- shinydashboard::box(title = "Select Control and Case group",
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
                                        column(2,uiOutput("mds_HTML")),
                                        column(2,uiOutput("glima_HTML")),
                                        column(2, offset = 0, actionButton("createTarget", "Create Target file")))
                                      
)
heatmap_box <- shinydashboard::box(title = "Pheatmap",
                                      status = "primary", solidHeader = TRUE, width = 12,
                                      fluidRow(
                                        column(8,
                                               sliderInput("significant_gene_number", "Number of significant genes",
                                                           min = 5, max = 35, value = 20, step = 5,
                                                           width = "80%")),
                                        br(),
                                        br(),
                                        br(),
                                        column(11,
                                               plotOutput("ttm_heatmap", width = "99%", height = "900px")),
                                        br(),
                                        br(),
                                        br(),
                                        dataTableOutput('check_logcpm')
                                      )
                                      
)

plots_box <- shinydashboard::box(title = "Plots",
                                   status = "primary", solidHeader = TRUE, width = 12,
                                   fluidRow(
                                     br(),
                                     column(6,
                                            br(),
                                            imageOutput("densplot",width = "99%", height = "400px"),
                                            br(),
                                            br(),
                                            ),
                                     br(),
                                     column(9,
                                            br(),
                                            br(),
                                            imageOutput("countsBoxplots",width = "99%", height = "400px"),
                                            br(),
                                            br(),
                                     ),
                                     br(),
                                     column(8,
                                            br(),
                                            br(),
                                            imageOutput("volcanoPlot",width = "99%", height = "900px"),
                                            br(),
                                            br(),
                                     ),
                                     br(),
                                     br(),
                                     br(),
                                     br(),
                                     br(),
                                     br()
                                   )
                                 
)
#tab2
fpkm_upload_box <- shinydashboard::box(title = "Upload Fpkm File",
                                       status = "info", solidHeader = TRUE, width = 12,
                                       fluidRow(
                                         column(12, fileInput('fpkm_file', label = "test",
                                                              multiple = TRUE,
                                                              buttonLabel = "Browse or Drop...",
                                                              placeholder = "(multiple) fpkm count"))
                                       )
)

fpkm_input_info_box <- shinydashboard::box(title = "Type the fpkm directory",
                                           status = "primary", solidHeader = TRUE, width = 12,
                                           fluidRow(
                                             column(6, 
                                                    textInput("fpkmDirectory", "Set fpkm directory path", "/home/arajo0707/Rat/fpkm/")),
                                             column(2, offset = 0, actionButton("create_fpkm_heatmap", "Create fpkm heatmap"))
                                           )
                                           
)
fpkm_heatmap <- shinydashboard::box(title = "FPKM heatmap",
                                    status = "primary", solidHeader = TRUE, width = 12,
                                    fluidRow(
                                      column(8,
                                             sliderInput("significant_fpkm_gene_number", "Number of significant genes",
                                                         min = 5, max = 35, value = 20, step = 5,
                                                         width = "80%")),
                                      br(),
                                      br(),
                                      br(),
                                      column(11,
                                             plotOutput("fpkm_heatmap", width = "auto", height = "900px"))
                                    )
                                    
)

# body ----
body <- dashboardBody(
  tabItems(
    tabItem(tabName = "sartools", fluidRow(
                                           upload_box,
                                           select_box,
                                           input_info_box,
                                           heatmap_box,
                                           plots_box)),
    tabItem(tabName = "fpkm", fluidRow(
      fpkm_input_info_box,
      fpkm_heatmap))
  )  
)

ui <- dashboardPage(header, sidebar, body, skin = "green")

server <- function(input, output, session) {
  
  ####### declare global variables ##############
  targets <-data.frame(matrix(ncol=3,nrow=100, dimnames=list(NULL, c("label","files","condition"))))
  logcpm <- NULL
  counts <- NULL
  mean_insert_size <- NULL
  fpkm_count <- NULL
  result_names <- NULL
  work_dir <- NULL
  ####### global variables ######################
  
  output$control_checkbox<-renderUI({
    checkboxGroupInput("controlSelector","Specify the control datasets:", choices = input$file1$name)
  })
  output$value <- renderText({ input$control_name })
  output$case_checkbox<-renderUI({
    checkboxGroupInput("caseSelector","Specify the case datasets:", choices = input$file1$name)
  })

  
  output$txt <- renderText({
    icons <- paste(input$controlSelector, collapse = ", ")
    paste("As a control ", paste(input$controlSelector, collapse = ", "))
  })
  output$txt2 <- renderText({
    icons <- paste(input$caseSelector, collapse = ", ")
    paste("As a case ", paste(input$caseSelector, collapse = ", "))
  })
  
  
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
    geneanno <- read.table("/home/arajo0707/Rat/analysis_rat1/mart_export.txt",header =TRUE, sep =',',fill = TRUE)
    geneanno <- dplyr::transmute(geneanno, Gene.stable.ID, Gene.name, length = Gene.end..bp. - Gene.start..bp.)
    colnames(geneanno) <- c('geneId', 'geneName', 'length')
    rownames(geneanno) <- geneanno$geneId
    geneanno <- geneanno[rownames(fcounts),]
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
  
  
  #### observe function of create HTML FIle button ###
  observeEvent(input$createTarget, {
    ## input$file1 structure > name, size, type, datapath.
    
    req(input$file1)
    file <- input$file1
    for(i in 1:nrow(file)){
      targets$label[i] <- strsplit(file$name[i], '[.]')[[1]][1]
      targets$files[i] <- file$name[i]
      
      if(file$name[i] %in% input$controlSelector){
        targets$condition[i] <- input$control_name
      }else{
        targets$condition[i] <- input$case_name
      }
    }

    targets <- na.omit(targets)
    targets <<- na.omit(targets)

    ########################end of observation event function#############################
    # rm(list=ls())                          # remove all the objects from the R session
    work_dir <<- paste0(input$workingDirectory)      # working directory for the R session
    print(work_dir)
    setwd(work_dir)
    workDir <- work_dir
    write.csv(targets, file = paste0(input$workingDirectory,"/targetFile.txt"), row.names = FALSE, col.names = FALSE)
    
    targetFile <- paste0(input$workingDirectory,"/targetFile.txt")
    projectName <- paste0(input$projectName)       # name of the project
    author <- paste0(input$projectAuthor)          # author of the statistical analysis/report
    # targetFile <- targets                  # path to the design/target file
    rawDir <- paste0(input$countDirectory)         # path to the directory containing raw counts files(.count file afterHTseq)
    featuresToRemove <- c("alignment_not_unique", "ambiguous", "no_feature", "not_aligned", "too_low_aQual")       
    # names of the features to be removed (specific HTSeq-count information and rRNA for example) NULL if no feature to remove
    
    varInt <- "condition"                                # factor of interest
    condRef <- input$control_name                        # reference biological condition
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
    counts <<- counts
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
    
    # print(input$projectName)
    
    
    output$resultHTML<-renderUI({
      workDirectoryMod <- ""
      for (i in 4:length(strsplit(input$workingDirectory, '[/]')[[1]])){
        workDirectoryMod <- paste0(workDirectoryMod,strsplit(input$workingDirectory, '[/]')[[1]][i],"/")
      }
      work_dir <<- workDirectoryMod
      a("Result HTML",
        target = "_blank",
        href = paste0("/files/",workDirectoryMod,input$projectName,"_report.html"),
        style = "text-decoration: underline;")
    })
    
    
    aJSH <- run.edgeR_New(counts=counts, target=target, varInt=varInt, condRef=condRef,
                          batch=batch, cpmCutoff=cpmCutoff, normalizationMethod=normalizationMethod,
                          pAdjustMethod=pAdjustMethod)
    lcpm <- cpm(out.edgeR$dge, log=TRUE)
    glMDSPlot(lcpm, labels=target$label, groups=target$condition, folder = "EdgeR", html = "MDS_Plot", launch=FALSE)
    result <- names(aJSH$results)
    result.filtered <- aJSH$results
    result_names <<- names(result.filtered)
    
    output$mds_HTML<-renderUI({

      a("MDS HTML",
        target = "_blank",
        href = paste0("/files/",work_dir,"/EdgeR/MDS_Plot.html"),
        style = "text-decoration: underline;")
    })
    output$glima_HTML<-renderUI({

      a("Glima HTML",
        target = "_blank",
        href = paste0("/files/",work_dir,"/EdgeR/", result_names,".html"),
        style = "text-decoration: underline;")
    })
    
    
    result_names <<- names(result.filtered)
    #filter genews where pvalue >0.05 && logFC>1 
    for (i in result){
      glMDPlot(aJSH$lrtJ[[i]], status=decideTests(aJSH$lrtJ[[i]]), main = i, side.main="Ensembl_Gene_ID", counts=aJSH$dge$counts, groups=target$condition, launch=FALSE, folder = "EdgeR", html = i)
    }
    for (i in length(result)){
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
    #
    output$check_logcpm <- renderTable(logcpm)
    
    output$ttm_heatmap <- renderPlot({
      if(!is.null(logcpm)){
        pheatmap(as.matrix(logcpm[1:input$significant_gene_number, ]), color = colorRampPalette(c('#2471A3','white','#C0392B'))(50), fontsize=15, main=paste0(result_names,"( n = ",input$significant_gene_number,")"),  border_color = 'white', show_rownames = T, cutree_cols = 2, cluster_row = T, scale = "row")
      }
      
    })
    output$densplot <- renderImage({
      
      outfile <- paste0(input$workingDirectory,"/figures/densplot.png")
      
      list(src = outfile,
           contentType = 'image/png',
           width = "90%",
           height = "120%")
    }, deleteFile = FALSE)
    
    output$countsBoxplots <- renderImage({
      
      outfile <- paste0(input$workingDirectory,"/figures/countsBoxplots.png")
      
      list(src = outfile,
           contentType = 'image/png',
           width = "90%",
           height = "120%")
    }, deleteFile = FALSE)
    
    output$volcanoPlot <- renderImage({
      
      outfile <- paste0(input$workingDirectory,"/figures/volcanoPlot.png")
      
      list(src = outfile,
           contentType = 'image/png',
           width = "100%",
           height = "100%" )
    }, deleteFile = FALSE)
    
  }) ### observe function of create HTML File button ###

  
  observeEvent(input$significant_gene_number, {
    output$ttm_heatmap <- renderPlot({
      if(!is.null(logcpm)){
        pheatmap(as.matrix(logcpm[1:input$significant_gene_number, ]), color = colorRampPalette(c('#2471A3','white','#C0392B'))(50), fontsize=15, main=paste(result_names,"(n = ",input$significant_gene_number,")"),  border_color = 'white', show_rownames = T, cutree_cols = 2, cluster_row = T, scale = "row")
      }
      
    })
    
  }) ##obeserve function significant gene number
  

  observeEvent(input$create_fpkm_heatmap, {
    target1 <- paste0(targets$label,"_size.txt")
    for (i in 1:length(target1)){
      meanFragment <- read.table(file.path(paste0(input$fpkmDirectory,target1[i])), nrows = 1, header =TRUE, sep ='\t',fill = TRUE)
      mean_insert_size[i] <<- meanFragment$MEAN_INSERT_SIZE
      
    }
    
    rat = useMart("ensembl", dataset = "rnorvegicus_gene_ensembl")
    rat <- useMart('ENSEMBL_MART_ENSEMBL')
    rat <- useDataset('rnorvegicus_gene_ensembl', rat)
    
    gene.annotations <- getBM(mart = rat, attributes=c("ensembl_gene_id", "external_gene_name",  "start_position", "end_position"))
    gene.annotations <- dplyr::transmute(gene.annotations, external_gene_name, ensembl_gene_id, length = end_position - start_position)
    gene.annotations <- gene.annotations %>% dplyr::filter(ensembl_gene_id %in% rownames(counts))
    gene.annotations <- gene.annotations[order(match(gene.annotations$ensembl_gene_id, rownames(counts))),]
    featureLength <- gene.annotations$length
    # meanFragmentLength <- t(meanFragmentLength)
    fpkm_matrix <- fpkm(counts, featureLength, mean_insert_size)
    
    topvar <- 35 #how many genes you want to see
    log10fpkm <- apply((fpkm_matrix + 1), 2, log10) #reduce the size of fpkm by putting in log10
    var_genes <- apply(log10fpkm, 1, var) #select top 30 significant genes
    select_var <- names(sort(var_genes, decreasing = TRUE))[1:topvar]
    fpkm_count <<- fpkm_matrix[select_var,]
    rownames(gene.annotations) = gene.annotations$ensembl_gene_id
    rownames(fpkm_count) <- gene.annotations[select_var,"external_gene_name"]
    fpkmheatmap(fpkm_count, topvar=topvar, showfeaturenames=TRUE, return_log = TRUE)
    
    output$fpkm_heatmap <- renderPlot({
      if(!is.null(fpkm_count)){
        fpkmheatmap(fpkm_count, topvar=input$significant_gene_number, showfeaturenames=TRUE, return_log = TRUE)
      }
      
    })
  })
  
}

shinyApp(ui, server)
