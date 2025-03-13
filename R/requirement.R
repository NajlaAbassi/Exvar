#' requirement
#'
#' Install required packages based on operating system and selected species
#'
#' @param all Logical. If `FALSE` (default), prompts the user to select a species from a predefined list.
#' If `TRUE`, installs dependencies for all available species (1 to 8) without prompting.
#' @return The function does not return a value but installs necessary packages.
#' @export
#'
#' @examples exvar::requirement()

requirement <- function(all = FALSE){

  #species selection
  if (all == FALSE){
    cat(paste0("These are the species currently supported by exvar, choose the number corresponding to the target specie: \n",
             "[1] Homo sapiens \n",
             "[2] Mus musculus \n",
             "[3] Arabidopsis thaliana \n",
             "[4] Drosophila melanogaster \n",
             "[5] Danio rerio \n",
             "[6] Rattus norvegicus \n",
             "[7] Saccharomyces cerevisiae \n",
             "[8] Caenorhabditis elegans \n"))
  species <- readline("Type the number of the species that you would like to use as a reference: ")
    } else {
    species <- 1:8
    }

  # Detect operating system
  switch(Sys.info()[['sysname']],
         "Linux"={

           ##### TODO use a loop instead of repeating million lines
           message("Detected operating system: Linux")
           #for all species
           message("Checking BiocManager")
           if (!requireNamespace("BiocManager")) install.packages("BiocManager")
           if (requireNamespace("BiocManager")) message("BiocManager installed")
           message("Checking shiny")
           if (!requireNamespace("shiny")) install.packages("shiny")
           if (requireNamespace("shiny")) message("shiny installed")
           message("Checking shinydashboard")
           if (!requireNamespace("shinydashboard")) install.packages("shinydashboard")
           if (requireNamespace("shinydashboard")) message("shinydashboard installed")
           message("Checking DT")
           if (!requireNamespace("DT")) install.packages("DT")
           if (requireNamespace("DT")) message("DT installed")
           message("Checking shinyWidgets")
           if (!requireNamespace("shinyWidgets")) install.packages("shinyWidgets")
           if (requireNamespace("shinyWidgets")) message("shinyWidgets installed")
           message("Checking shinythemes")
           if (!requireNamespace("shinythemes")) install.packages("shinythemes")
           if (requireNamespace("shinythemes")) message("shinythemes installed")
           message("Checking dplyr")
           if (!requireNamespace("dplyr")) install.packages("dplyr")
           if (requireNamespace("dplyr")) message("dplyr installed")
           message("Checking BiocParallel")
           if (!requireNamespace("BiocParallel")) BiocManager::install("BiocParallel")
           if (requireNamespace("BiocParallel")) message("BiocParallel installed")
           message("Checking GenomeInfoDb")
           if (!requireNamespace("GenomeInfoDb")) BiocManager::install("GenomeInfoDb")
           if (requireNamespace("GenomeInfoDb")) message("GenomeInfoDb installed")
           message("Checking GenomicAlignments")
           if (!requireNamespace("GenomicAlignments")) BiocManager::install("GenomicAlignments")
           if (requireNamespace("GenomicAlignments")) message("GenomicAlignments installed")
           message("Checking gmapR")
           if (!requireNamespace("gmapR")) BiocManager::install("gmapR")
           if (requireNamespace("gmapR")) message("gmapR installed")
           message("Checking panelcn.mops")
           if (!requireNamespace("panelcn.mops")) BiocManager::install("panelcn.mops")
           if (requireNamespace("panelcn.mops")) message("panelcn.mops installed")
           message("Checking R.utils")
           if (!requireNamespace("R.utils")) install.packages("R.utils")
           if (requireNamespace("R.utils")) message("R.utils installed")
           message("Checking Rfastp")
           if (!requireNamespace("Rfastp")) BiocManager::install("Rfastp")
           if (requireNamespace("Rfastp")) message("Rfastp installed")
           message("Checking Rsamtools")
           if (!requireNamespace("Rsamtools"))  BiocManager::install("Rsamtools")
           if (requireNamespace("Rsamtools")) message("Rsamtools installed")
           message("Checking VariantTools")
           if (!requireNamespace("VariantTools"))  BiocManager::install("VariantTools")
           if (requireNamespace("VariantTools")) message("VariantTools installed")
           message("Checking tools")
           if (!requireNamespace("tools")) install.packages("tools")
           if (requireNamespace("tools")) message("tools installed")
           message("Checking data.table")
           if (!requireNamespace("data.table")) install.packages("data.table")
           if (requireNamespace("data.table")) message("data.table installed")
           message("Checking readr")
           if (!requireNamespace("readr")) install.packages("readr")
           if (requireNamespace("readr")) message("readr installed")
           message("Checking DESeq2")
           if (!requireNamespace("DESeq2")) BiocManager::install("DESeq2")
           if (requireNamespace("DESeq2")) message("DESeq2 installed")
           message("Checking ggplot2")
           if (!requireNamespace("ggplot2")) install.packages("ggplot2")
           if (requireNamespace("ggplot2")) message("ggplot2 installed")
           message("Checking ComplexHeatmap")
           if (!requireNamespace("ComplexHeatmap")) BiocManager::install("ComplexHeatmap")
           if (requireNamespace("ComplexHeatmap")) message("ComplexHeatmap installed")
           message("Checking clusterProfiler")
           if (!requireNamespace("clusterProfiler")) BiocManager::install("clusterProfiler")
           if (requireNamespace("clusterProfiler")) message("clusterProfiler installed")
           message("Checking AnnotationDbi")
           if (!requireNamespace("AnnotationDbi")) BiocManager::install("AnnotationDbi")
           if (requireNamespace("AnnotationDbi")) message("AnnotationDbi installed")
           message("Checking enrichplot")
           if (!requireNamespace("enrichplot")) BiocManager::install("enrichplot")
           if (requireNamespace("enrichplot")) message("enrichplot installed")
           message("Checking GenomicFeatures")
           if (!requireNamespace("GenomicFeatures"))  BiocManager::install("GenomicFeatures")
           if (requireNamespace("GenomicFeatures")) message("GenomicFeatures installed")
           message("Checking stringr")
           if (!requireNamespace("stringr")) install.packages("stringr")
           if (requireNamespace("stringr")) message("stringr installed")
           message("Checking VariantAnnotation")
           if (!requireNamespace("VariantAnnotation"))  BiocManager::install("VariantAnnotation")
           if (requireNamespace("VariantAnnotation")) message("VariantAnnotation installed")
           message("Checking CNVRanger")
           if (!requireNamespace("CNVRanger")) BiocManager::install("CNVRanger")
           if (requireNamespace("CNVRanger")) message("CNVRanger installed")
           message("Checking Gviz")
           if (!requireNamespace("Gviz")) BiocManager::install("Gviz")
           if (requireNamespace("Gviz")) message("Gviz installed")
           message("Checking AnnotationHub")
           if (!requireNamespace("AnnotationHub")) BiocManager::install("AnnotationHub")
           if (requireNamespace("AnnotationHub")) message("AnnotationHub installed")
           message("Checking regioneR")
           if (!requireNamespace("regioneR")) BiocManager::install("regioneR")
           if (requireNamespace("regioneR")) message("regioneR installed")
           message("Checking ggnewscale")
           if (!requireNamespace("ggnewscale")) install.packages("ggnewscale")
           if (requireNamespace("ggnewscale")) message("ggnewscale installed")
           message("Checking tibble")
           if (!requireNamespace("tibble")) install.packages("tibble")
           if (requireNamespace("tibble")) message("tibble installed")
           ########
         },
         "Windows"={
           message("Detected operating system: Windows")
           message("Checking BiocManager")
           if (!requireNamespace("BiocManager")) install.packages("BiocManager")
           if (requireNamespace("BiocManager")) message("BiocManager installed")
           message("Checking shiny")
           if (!requireNamespace("shiny")) install.packages("shiny")
           if (requireNamespace("shiny")) message("shiny installed")
           message("Checking shinydashboard")
           if (!requireNamespace("shinydashboard")) install.packages("shinydashboard")
           if (requireNamespace("shinydashboard")) message("shinydashboard installed")
           message("Checking DT")
           if (!requireNamespace("DT")) install.packages("DT")
           if (requireNamespace("DT")) message("DT installed")
           message("Checking shinyWidgets")
           if (!requireNamespace("shinyWidgets")) install.packages("shinyWidgets")
           if (requireNamespace("shinyWidgets")) message("shinyWidgets installed")
           message("Checking shinythemes")
           if (!requireNamespace("shinythemes")) install.packages("shinythemes")
           if (requireNamespace("shinythemes")) message("shinythemes installed")
           message("Checking dplyr")
           if (!requireNamespace("dplyr")) install.packages("dplyr")
           if (requireNamespace("dplyr")) message("dplyr installed")
           message("Checking GenomeInfoDb")
           if (!requireNamespace("GenomeInfoDb")) BiocManager::install("GenomeInfoDb")
           if (requireNamespace("GenomeInfoDb")) message("GenomeInfoDb installed")
           message("Checking tools")
           if (!requireNamespace("tools")) install.packages("tools")
           if (requireNamespace("tools")) message("tools installed")
           message("Checking data.table")
           if (!requireNamespace("data.table")) install.packages("data.table")
           if (requireNamespace("data.table")) message("data.table installed")
           message("Checking readr")
           if (!requireNamespace("readr")) install.packages("readr")
           if (requireNamespace("readr")) message("readr installed")
           message("Checking DESeq2")
           if (!requireNamespace("DESeq2")) BiocManager::install("DESeq2")
           if (requireNamespace("DESeq2")) message("DESeq2 installed")
           message("Checking ggplot2")
           if (!requireNamespace("ggplot2")) install.packages("ggplot2")
           if (requireNamespace("ggplot2")) message("ggplot2 installed")
           message("Checking ComplexHeatmap")
           if (!requireNamespace("ComplexHeatmap")) BiocManager::install("ComplexHeatmap")
           if (requireNamespace("ComplexHeatmap")) message("ComplexHeatmap installed")
           message("Checking clusterProfiler")
           if (!requireNamespace("clusterProfiler")) BiocManager::install("clusterProfiler")
           if (requireNamespace("clusterProfiler")) message("clusterProfiler installed")
           message("Checking AnnotationDbi")
           if (!requireNamespace("AnnotationDbi")) BiocManager::install("AnnotationDbi")
           if (requireNamespace("AnnotationDbi")) message("AnnotationDbi installed")
           message("Checking enrichplot")
           if (!requireNamespace("enrichplot")) BiocManager::install("enrichplot")
           if (requireNamespace("enrichplot")) message("enrichplot installed")
           message("Checking stringr")
           if (!requireNamespace("stringr")) install.packages("stringr")
           if (requireNamespace("stringr")) message("stringr installed")
           message("Checking CNVRanger")
           if (!requireNamespace("CNVRanger")) BiocManager::install("CNVRanger")
           if (requireNamespace("CNVRanger")) message("CNVRanger installed")
           message("Checking Gviz")
           if (!requireNamespace("Gviz")) BiocManager::install("Gviz")
           if (requireNamespace("Gviz")) message("Gviz installed")
           message("Checking AnnotationHub")
           if (!requireNamespace("AnnotationHub")) BiocManager::install("AnnotationHub")
           if (requireNamespace("AnnotationHub")) message("AnnotationHub installed")
           message("Checking regioneR")
           if (!requireNamespace("regioneR")) BiocManager::install("regioneR")
           if (requireNamespace("regioneR")) message("regioneR installed")
           message("Checking ggnewscale")
           if (!requireNamespace("ggnewscale")) install.packages("ggnewscale")
           if (requireNamespace("ggnewscale")) message("ggnewscale installed")
           message("Checking tibble")
           if (!requireNamespace("tibble")) install.packages("tibble")
           if (requireNamespace("tibble")) message("tibble installed")
           message("Only packages related to visualisation functions are available for Windows. For full exvar functionality, a Linux operating system is requireNamespaced.")
         },
         "Darwin"={
           message("The package requireNamespace a Linux or Windows system")
         }) # is Darwin correct ?

  #species specific
  for (species in species){
    switch(species,
         "1"={
           message("Checking BSgenome.Hsapiens.UCSC.hg19")
           if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg19")) BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
           if (requireNamespace("BSgenome.Hsapiens.UCSC.hg19")) message("BSgenome.Hsapiens.UCSC.hg19 installed")
           message("Checking BSgenome.Hsapiens.UCSC.hg38")
           if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg38")) BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
           if (requireNamespace("BSgenome.Hsapiens.UCSC.hg38")) message("BSgenome.Hsapiens.UCSC.hg38 installed")
           message("Checking org.Hs.eg.db")
           if (!requireNamespace("org.Hs.eg.db")) BiocManager::install("org.Hs.eg.db")
           if (requireNamespace("org.Hs.eg.db")) message("org.Hs.eg.db installed")
           message("Checking TxDb.Hsapiens.UCSC.hg19.knownGene")
           if (!requireNamespace("TxDb.Hsapiens.UCSC.hg19.knownGene")) BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
           if (requireNamespace("TxDb.Hsapiens.UCSC.hg19.knownGene")) message("TxDb.Hsapiens.UCSC.hg19.knownGene installed")
           message("Checking TxDb.Hsapiens.UCSC.hg38.knownGene")
           if (!requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene")) BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
           if (requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene")) message("TxDb.Hsapiens.UCSC.hg38.knownGene installed")
           message("Checking SNPlocs.Hsapiens.dbSNP144.GRCh37")
           if (!requireNamespace("SNPlocs.Hsapiens.dbSNP144.GRCh37")) BiocManager::install("SNPlocs.Hsapiens.dbSNP144.GRCh37")
           if (requireNamespace("SNPlocs.Hsapiens.dbSNP144.GRCh37")) message("SNPlocs.Hsapiens.dbSNP144.GRCh37 installed")
           message("Checking XtraSNPlocs.Hsapiens.dbSNP144.GRCh37")
           if (!requireNamespace("XtraSNPlocs.Hsapiens.dbSNP144.GRCh37")) BiocManager::install("XtraSNPlocs.Hsapiens.dbSNP144.GRCh37")
           if (requireNamespace("XtraSNPlocs.Hsapiens.dbSNP144.GRCh37")) message("XtraSNPlocs.Hsapiens.dbSNP144.GRCh37 installed")
           message("Checking SNPlocs.Hsapiens.dbSNP151.GRCh38")
           if (!requireNamespace("SNPlocs.Hsapiens.dbSNP151.GRCh38")) BiocManager::install("SNPlocs.Hsapiens.dbSNP151.GRCh38")
           if (requireNamespace("SNPlocs.Hsapiens.dbSNP151.GRCh38")) message("SNPlocs.Hsapiens.dbSNP151.GRCh38 installed")
           message("Checking XtraSNPlocs.Hsapiens.dbSNP144.GRCh38")
           if (!requireNamespace("XtraSNPlocs.Hsapiens.dbSNP144.GRCh38")) BiocManager::install("XtraSNPlocs.Hsapiens.dbSNP144.GRCh38")
           if (requireNamespace("XtraSNPlocs.Hsapiens.dbSNP144.GRCh38")) message("XtraSNPlocs.Hsapiens.dbSNP144.GRCh38 installed")

         },
         "2"={
           message("Checking BSgenome.Mmusculus.UCSC.mm10")
           if (!requireNamespace("BSgenome.Mmusculus.UCSC.mm10")) BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")
           if (requireNamespace("BSgenome.Mmusculus.UCSC.mm10")) message("BSgenome.Mmusculus.UCSC.mm10 installed")
           message("Checking org.Mm.eg.db")
           if (!requireNamespace("org.Mm.eg.db")) BiocManager::install("org.Mm.eg.db")
           if (requireNamespace("org.Mm.eg.db")) message("org.Mm.eg.db installed")
           message("Checking TxDb.Mmusculus.UCSC.mm10.knownGene")
           if (!requireNamespace("TxDb.Mmusculus.UCSC.mm10.knownGene")) BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
           if (requireNamespace("TxDb.Mmusculus.UCSC.mm10.knownGene")) message("TxDb.Mmusculus.UCSC.mm10.knownGene installed")

         },
         "3"={
           message("Checking BSgenome.Athaliana.TAIR.TAIR9")
           if (!requireNamespace("BSgenome.Athaliana.TAIR.TAIR9")) BiocManager::install("BSgenome.Athaliana.TAIR.TAIR9")
           if (requireNamespace("BSgenome.Athaliana.TAIR.TAIR9")) message("BSgenome.Athaliana.TAIR.TAIR9 installed")
           message("Checking org.At.tair.db")
           if (!requireNamespace("org.At.tair.db")) BiocManager::install("org.At.tair.db")
           if (requireNamespace("org.At.tair.db")) message("org.At.tair.db installed")
           message("Checking TxDb.Athaliana.BioMart.plantsmart28")
           if (!requireNamespace("TxDb.Athaliana.BioMart.plantsmart28")) BiocManager::install("TxDb.Athaliana.BioMart.plantsmart28")
           if (requireNamespace("TxDb.Athaliana.BioMart.plantsmart28")) message("TxDb.Athaliana.BioMart.plantsmart28 installed")

         },
         "4"={
           message("Checking BSgenome.Dmelanogaster.UCSC.dm6")
           if (!requireNamespace("BSgenome.Dmelanogaster.UCSC.dm6")) BiocManager::install("BSgenome.Dmelanogaster.UCSC.dm6")
           if (requireNamespace("BSgenome.Dmelanogaster.UCSC.dm6")) message("BSgenome.Dmelanogaster.UCSC.dm6 installed")
           message("Checking org.Dm.eg.db")
           if (!requireNamespace("org.Dm.eg.db")) BiocManager::install("org.Dm.eg.db")
           if (requireNamespace("org.Dm.eg.db")) message("org.Dm.eg.db installed")
           message("Checking TxDb.Dmelanogaster.UCSC.dm6.ensGene")
           if (!requireNamespace("TxDb.Dmelanogaster.UCSC.dm6.ensGene")) BiocManager::install("TxDb.Dmelanogaster.UCSC.dm6.ensGene")
           if (requireNamespace("TxDb.Dmelanogaster.UCSC.dm6.ensGene")) message("TxDb.Dmelanogaster.UCSC.dm6.ensGene installed")

         },
         "5"={
           message("Checking BSgenome.Drerio.UCSC.danRer11")
           if (!requireNamespace("BSgenome.Drerio.UCSC.danRer11")) BiocManager::install("BSgenome.Drerio.UCSC.danRer11")
           if (requireNamespace("BSgenome.Drerio.UCSC.danRer11")) message("BSgenome.Drerio.UCSC.danRer11 installed")
           message("Checking org.Dr.eg.db")
           if (!requireNamespace("org.Dr.eg.db")) BiocManager::install("org.Dr.eg.db")
           if (requireNamespace("org.Dr.eg.db")) message("org.Dr.eg.db installed")
           message("Checking TxDb.Drerio.UCSC.danRer11.refGene")
           if (!requireNamespace("TxDb.Drerio.UCSC.danRer11.refGene")) BiocManager::install("TxDb.Drerio.UCSC.danRer11.refGene")
           if (requireNamespace("TxDb.Drerio.UCSC.danRer11.refGene")) message("TxDb.Drerio.UCSC.danRer11.refGene installed")

         },
         "6"={ # rat
           message("Checking BSgenome.Rnorvegicus.UCSC.rn5")
           if (!requireNamespace("BSgenome.Rnorvegicus.UCSC.rn5")) BiocManager::install("BSgenome.Rnorvegicus.UCSC.rn5")
           if (requireNamespace("BSgenome.Rnorvegicus.UCSC.rn5")) message("BSgenome.Rnorvegicus.UCSC.rn5 installed")
           message("Checking org.Rn.eg.db")
           if (!requireNamespace("org.Rn.eg.db")) BiocManager::install("org.Rn.eg.db")
           if (requireNamespace("org.Rn.eg.db")) message("org.Rn.eg.db installed")
           message("Checking TxDb.Rnorvegicus.UCSC.rn5.refGene")
           if (!requireNamespace("TxDb.Rnorvegicus.UCSC.rn5.refGene"))BiocManager::install("TxDb.Rnorvegicus.UCSC.rn5.refGene")
           if (requireNamespace("TxDb.Rnorvegicus.UCSC.rn5.refGene")) message("TxDb.Rnorvegicus.UCSC.rn5.refGene installed")

         },
         "7"={
           message("Checking BSgenome.Scerevisiae.UCSC.sacCer3")
           if (!requireNamespace("BSgenome.Scerevisiae.UCSC.sacCer3")) BiocManager::install("BSgenome.Scerevisiae.UCSC.sacCer3")
           if (requireNamespace("BSgenome.Scerevisiae.UCSC.sacCer3")) message("BSgenome.Scerevisiae.UCSC.sacCer3 installed")
           message("Checking org.Sc.sgd.db")
           if (!requireNamespace("org.Sc.sgd.db")) BiocManager::install("org.Sc.sgd.db")
           if (requireNamespace("org.Sc.sgd.db")) message("org.Sc.sgd.db installed")
           message("Checking TxDb.Scerevisiae.UCSC.sacCer3.sgdGene")
           if (!requireNamespace("TxDb.Scerevisiae.UCSC.sacCer3.sgdGene")) BiocManager::install("TxDb.Scerevisiae.UCSC.sacCer3.sgdGene")
           if (requireNamespace("TxDb.Scerevisiae.UCSC.sacCer3.sgdGene")) message("TxDb.Scerevisiae.UCSC.sacCer3.sgdGene installed")

         },
         "8"={
           message("Checking BSgenome.Celegans.UCSC.ce11")
           if (!requireNamespace("BSgenome.Celegans.UCSC.ce11")) BiocManager::install("BSgenome.Celegans.UCSC.ce11")
           if (requireNamespace("BSgenome.Celegans.UCSC.ce11")) message("BSgenome.Celegans.UCSC.ce11 installed")
           message("Checking org.Ce.eg.db")
           if (!requireNamespace("org.Ce.eg.db")) BiocManager::install("org.Ce.eg.db")
           if (requireNamespace("org.Ce.eg.db")) message("org.Ce.eg.db installed")
           message("Checking TxDb.Celegans.UCSC.ce11.refGene")
           if (!requireNamespace("TxDb.Celegans.UCSC.ce11.refGene")) BiocManager::install("TxDb.Celegans.UCSC.ce11.refGene")
           if (requireNamespace("TxDb.Celegans.UCSC.ce11.refGene")) message("TxDb.Celegans.UCSC.ce11.refGene installed")

         }
  )
  }
}
