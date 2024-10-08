#' Params
#' @param genelist must be a object of type list() containing an undetermined number of character vectors. Each component of the list is a char vector of gene names. This function will perform gene ontology enrichment using the char vector as test sample and the totality of your organism's genes as population (Universe). 
#' @gene_universe is a vector of characters containing the genes you want to use as universe. This is normally the totality of genes in your species, or the totality of genes with recollected expression in your dataset.
#' @alg is the algorithm of choice, default is 'Classic', but 'elim' is advised
#' @cols is a char vector of the colors used in the barplots
#' 
require(topGO)
require(colorspace)
# getGOs(list(tfs,egs,active_tfs,target_genes,V(graph)$name),toy_universe)
getGOs <- function(
  genelist,
  gene_universe,
  gene2GO,
  alg = "Classic",
  stat = "fisher",
  category = "BP",
  cols = sequential_hcl(7,"PinkYl")[1:6],
  max_terms = NULL
) {
  res <- list(
    GOtable = list(),
    GOplot = list()
  )
  geneID2GO <- gene2GO
  for (j in 1:length(genelist)){
    print(paste0("Starting analysis ",j," of ",length(genelist)))
    x <- genelist[[j]]
    numgenes <- length(x)
    allgenes <- data.frame(
      id = gene_universe)
    allgenes$diffreg <- 0
    allgenes$diffreg[allgenes$id %in% x] <- 1
    genelist_x <- allgenes$diffreg
    names(genelist_x) <- allgenes$id
    signif <- function (allScore) {return(allScore == 1)} # function for TopGO which retrieves the significant genes as those of interest
    # this function should be replaced for something much more efficient. Check other versions of these scripts run by other people
    
    #load data
    #' This is by far the longest and most computationally expensive step.
    #' Can this step be outside the loop and just keep inside the loop 
    #' the bit of which ones to keep?
    GOdata_x <- new( 
      "topGOdata",
      ontology = category,
      allGenes = genelist_x,
      geneSelectionFun = signif,
      annot = annFUN.gene2GO,
      gene2GO = geneID2GO
    )
    
    #test
    res_x <- runTest(GOdata_x, algorithm = alg, statistic = stat)
    
    #table
    allres_x <- GenTable(
      GOdata_x,
      classicFisher = res_x,
      orderBy = "classicFisher", #see if this can be replaced for a static thing that does not change with the test's name
      ranksOf = "classicFisher",
      topNodes = 30,
      numChar = 5000
    )
    allres_x <- allres_x[allres_x$Annotated > 5,]
    allres_x$classicFisher[allres_x$classicFisher=="< 1e-30"] <- 1e-30
    allres_x$classicFisher <-
      parse_number(allres_x$classicFisher)
    allres_x$classicFisher <-
      as.numeric(allres_x$classicFisher)
    
    if(!is.null(max_terms)){
      maxgenesplot <- ifelse(nrow(allres_x) > max_terms, max_terms, nrow(allres_x))
    } else {
      maxgenesplot <- nrow(allres_x)
    }
    
    res$GOtable[[j]] <- allres_x
    
    #plot
    res$GOplot[[j]] <- 
      ggplot(
        allres_x[1:maxgenesplot,], 
        aes(
          x = reorder(Term, -log10(classicFisher)),
          y = Significant/Expected,
          fill = -log10(classicFisher)
        )
      )+
      geom_bar(stat = "identity") +
      scale_fill_gradientn(colors =  sequential_hcl(7,"PinkYl")[1:6], name = "-log10(pvalue)") +
      coord_flip() +
      theme_classic() +
      labs(x = "GO term", y = "Significant/Expected ratio",
           title = names(genelist)[j])+
      geom_text(aes(label=Term, y = max(Significant/Expected)/2))+
      theme(axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank()
      )
  }
  names(res$GOtable) <- names(genelist)
  names(res$GOplot) <- names(genelist)
  return(res)
}
