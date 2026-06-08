library("openxlsx")
library("gnumeric")
library("org.Hs.eg.db")
library("jsonlite")
library("processx")
library("pheatmap")
library("seqinr")
library("ape")
library("dplyr")
library("ggplot2")
library("stringr")
library("rentrez")
library("gggenomes")
library("protodeviser")

# INITIAL DIRS ####
# Set top working dir. Trailing "/" *is* required, also for the rest.
top.d <- "/home/petrov/Projects/ecm-evolution/pipeline/"
setwd(top.d)

# Input spreadsheets
table.d <- paste0(top.d, "spreadsheets/")

# BASH scripts
bash.d <- paste0(top.d, "sh/")

# WORK DIRS ####
wrk.d <- paste0(top.d, "work/")

# Downloads, output and temporary files dirs
get.d <- paste0(wrk.d, "downloads/")
out.d <- paste0(wrk.d, "output/")
tmp.d <- paste0(wrk.d, "temp/")

genes.d <- paste0(wrk.d, "genes/")
trees.d <- paste0(wrk.d, "trees/")
blast.d <- paste0(wrk.d, "blast/")
ortho.d <- paste0(wrk.d, "orthol/")
gbloc.d <- paste0(wrk.d, "Gblocks/")
image.d <- paste0(wrk.d, "images/")
jsons.d <- paste0(wrk.d, "json/")

# Create work dirs ####
dir.create(wrk.d)

dir.create(get.d)
dir.create(out.d)
dir.create(tmp.d)

dir.create(genes.d)
dir.create(trees.d)
dir.create(blast.d)
dir.create(ortho.d)
dir.create(gbloc.d)
dir.create(image.d)
dir.create(jsons.d)

# HyPhy work dir. Created by the BASH script upon run.
hyphy.d <- "/var/tmp/hyphy/"

################################################################################

# A function to run an external shell script in xterm, from the bash.d dir.
runSh <- function(sh = NULL){
  system(paste("xterm -e 'bash", paste0(bash.d, sh),  "' &"))
}

# A function to remove accession release version, after the dot, e.g. NM_138455.4 -> NM_138455
# https://www.r-bloggers.com/2024/07/how-to-extract-string-after-a-specific-character-in-r/
rmdot <- function(x = NULL){
  x <- sub("\\..*", "", x)
  return(x)
}

# A function to remove isoform number, after a dash, e.g. Q14766-1 -> Q14766
rmdash <- function(x = NULL){
  x <- sub("-.*", "", x)
  return(x)
}

# A function to get isoform number, after dash, e.g. Q14766-1 -> 1
getiso <- function(x = NULL){
  x <- unlist(strsplit(x, split = "-"))
  x <- x[seq(2, length(x), 2)]
  x <- as.character(paste0(x))
  return(x)
}

# A function to convert sequence from vector to single string, e.g. c("M", "A", "L") - > "MAL"
seqToString <- function(x = NULL){
  x <- gsub(" ", "", paste0(x, collapse = " "))
  return(x)
}

# DOWNLOAD ####
# Get data from GeneCards (https://www.genenames.org/download/) in order to match gene names and (canonical?) UniProt IDs
#download.file(url = "storage.googleapis.com/public-download-files/hgnc/json/json/hgnc_complete_set.json", method = "wget", destfile = paste0(get.d, "hgnc_complete_set.", Sys.Date(), ".json"))

# Accessed: 2026-01-07; name accordingly
hgnc_complete_set.json <- "hgnc_complete_set.2026-01-08.json"

# HSAP SECTION BEGINS ##########################################################

# A function to load matrisome db (Hs_Matrisome_Masterlist_Naba et al_2012.xlsx);
# Get this through browser and save in "downloads", # we are interested just in the core members
# TODO: make the download automatic?
core.matrisome <- function(matrisome.xlsx = NULL, work.d = NULL){
  
  matrisome <- read.xlsx(paste0(work.d, matrisome.xlsx))
  core <- subset(matrisome, matrisome$Matrisome.Division == "Core matrisome")
  core.genes <- data.frame(Matrisome.Division=core$Matrisome.Division,
                           category=core$Matrisome.Category,
                           symbol=core$Gene.Symbol,
                           name=core$Gene.Name,
                           synonyms=core$Synonyms,
                           gene_id=core$HGNC_IDs)
  return(core.genes)
}
core.genes <- core.matrisome(matrisome.xlsx = "Hs_Matrisome_Masterlist_Naba et al_2012.xlsx", work.d = table.d)
core.genes

# A function to import data from GeneCards (hgnc_complete_set.json), process and simplify
simple.genecards <- function(json_input = NULL, work.d = NULL){
  
  # import data from GeneCards and process
  hgnc_complete_set <- fromJSON(txt = paste0(work.d, json_input))
  hgnc_complete_set <- as.data.frame(hgnc_complete_set$response$docs)
  hgnc_complete_set <- flatten(hgnc_complete_set)
  hgnc_complete_set
  
  # create a simple dataframe of the stuff we need
  genecards_subset <- data.frame(symbol = unlist(lapply(hgnc_complete_set$symbol, function(x) ifelse(is.null(x), NA, x))),
                                 name = unlist(lapply(hgnc_complete_set$name, function(x) ifelse(is.null(x), NA, x))),
                                 ensembl_gene_id = unlist(lapply(hgnc_complete_set$ensembl_gene_id, function(x) ifelse(is.null(x), NA, x))),
                                 refseq_accession = unlist(lapply(hgnc_complete_set$refseq_accession, function(x) ifelse(is.null(x), NA, x))),
                                 uniprot_ids = unlist(lapply(hgnc_complete_set$uniprot_ids, function(x) ifelse(is.null(x), NA, x))),
                                 entrez_id = unlist(lapply(hgnc_complete_set$entrez_id, function(x) ifelse(is.null(x), NA, x))),
                                 hgnc_id = gsub("HGNC:", "", unlist(lapply(hgnc_complete_set$hgnc_id, function(x) ifelse(is.null(x), NA, x)))))
  
  return(genecards_subset)
}
genecards_subset <- simple.genecards(json_input = hgnc_complete_set.json, work.d = get.d)
genecards_subset

# A function to match genes from matrisome core against the simplified GeneCards
# table, in order to get the (canonical) UniProt IDs
create.ecm.core <- function(matrisome = NULL, genecards = NULL){
  core.ecm <- data.frame()
  for (r in 1:nrow(matrisome)) {
    
    df.l <- data.frame()
    f <- which(matrisome[r,6] == genecards_subset$hgnc_id)
    
    df.l <- data.frame(symbol = genecards_subset$symbol[f],
                       refseq.cds = genecards_subset$refseq_accession[f], ## This returns mismatching results between the canonical isoforms at NCBI and UniProt sometimes
                       uniprot = genecards_subset$uniprot_ids[f],
                       ensembl = genecards_subset$ensembl_gene_id[f],
                       hgnc_id = genecards_subset$hgnc_id[f],
                       entrez_id = genecards_subset$entrez_id[f],
                       long_name = genecards_subset$name[f],
                       category=matrisome[r,2])
    
    core.ecm <- rbind(core.ecm, df.l)
  }
  return(core.ecm)
}
core.ecm.auto <- create.ecm.core(matrisome = core.genes, genecards = genecards_subset)
core.ecm.auto

# WRITE ####
# Save UniProt IDs to disk (no isoforms version), to be parsed to the "download_uniprot_txt.sh"
# shell script afterwards, which downloads each entry's summary info as a txt.
#write(core.ecm.auto$uniprot, file = paste0(tmp.d, "uniprot.txt"), ncolumns = 1)

# BASH ####
# This downloads UniProt info for each protein in simple TXT format.
#runSh("download_uniprot_txt.sh") # Accessed: 2026-01-08

# A function to import the isoforms listed as "Displayed" (canonical) on UniProt, from the
# downloaded summary file (txt)
import.displayed <- function(displayed.txt = NULL, ecm = core.ecm.auto){
  displayed <- read.table(displayed.txt)
  displayed <- unique(displayed)
  displayed <- data.frame(uniprot = rmdash(displayed[,1]), displayed = displayed[,1], isoform = getiso(displayed[,1]))
  displayed <- merge.data.frame(core.ecm.auto, displayed, all = T)
  displayed$displayed[is.na(displayed$displayed)] <- displayed$uniprot[is.na(displayed$displayed)]
  return(displayed)
}
displayed <- import.displayed(displayed.txt = paste0(tmp.d, "UniProt.displayed.txt"))
displayed

# WRITE ####
# Save the UniProt "displayed" IDs to disk (supposedly canonical isoforms) if found, and the
# "generic" ID if not found. To be parsed to the "download_uniprot.sh" in order to retrieve
# the corresponding FASTA file.
#write(displayed$displayed, file = paste0(tmp.d, "UniProt.canonical.txt"))

# BASH ####
# This downloads UniProt fasta sequences
#runSh("download_uniprot.sh") # Accessed: 2026-01-08

# A function to match the UniProt Displayed isoform to RefSeq IDs, listed in the txt file
# from UniProt.
match.refseq <- function(displayed = NULL){
  refseq <- data.frame()
  for (r in 1:nrow(displayed)) {
    i <- displayed$uniprot[r]
    d <- displayed$displayed[r]
    e <- displayed$entrez_id[r]
    c <- displayed$category[r]

    line.l <- readLines(paste0(get.d, "UniProt/", i, ".txt"))
    
    # get the protein id
    id.l <- head(line.l, n = 1)
    id.l <- strsplit(id.l, " ")[[1]][4]
    
    # get the protein name
    upid <- grep("GN", line.l, value = T)
    upid <- grep("Name=", upid, value = T)
    upid <- strsplit(upid, split = " ")
    upid <- upid[[1]][[4]]
    upid <- gsub("Name=", "", upid)
    upid <- gsub(";", "", upid)

    # get the protein description
    desc.l <- grep("RecName:", line.l, value = T)
    desc.l <- strsplit(desc.l, split = "=")
    desc.l <- gsub(";", "", desc.l[[1]][[2]])
    desc.l <- gsub(" \\{.*", "", desc.l)
    
    # get the isoforms line
    line.l <- grep("RefSeq", line.l, value = T)
    line.l <- data.frame(isoforms = line.l)
    line.l <- strsplit(line.l$isoforms, " ")
    
    # check if displayed isoform is found at all. Thanks UniProt for this
    f <- c()
    for (l in 1:length(line.l)) {
      for (j in 1:length(line.l[[l]])) {
        f.l <- grep(d, line.l[[l]][[j]])
        f <- c(f, f.l)
      }
    }

    df <- data.frame()
    for (l in 1:length(line.l)) {
      pep <- line.l[[l]][[5]]
      cds <- line.l[[l]][[6]]
      
      pep <- sub(".$", "", pep)
      cds <- sub(".$", "", cds)
      
      if (isTRUE(length(line.l[[l]]) == 7)) {
        isoform <- line.l[[l]][[7]]
        isoform <- gsub("\\]", "", isoform)
        isoform <- gsub("\\[", "", isoform)
      } else {
        isoform <- i
      }

      # does the isoform match the one displayed?
      if (isTRUE(isoform == displayed$displayed[r])){
        canonic <- "canonical"
        summary <- "canonical"
        
        # does the isoform NOT match the one displayed, but the displayed is still among the isoforms listed?
      } else if (isTRUE(isoform != displayed$displayed[r]) & isTRUE(length(f) > 0)){
        canonic <- "alternative"
        summary <- "alternative"
        
        # does the isoform NOT match the one displayed, NOT match the unversioned one and the displayed is NOT among the isoforms listed?
      } else if (isTRUE(isoform != displayed$displayed[r]) & isTRUE(isoform != i) & isTRUE(length(f) == 0)){
        canonic <- "missing"
        summary <- "alternative"
        
        # does the isoform NOT match the one displayed, BUT matches the unversioned one and the displayed is NOT among the isoforms listed?
      } else if (isTRUE(isoform != displayed$displayed[r]) & isTRUE(isoform == i) & isTRUE(length(f) == 0)){
        canonic <- "canonical.else"
        summary <- "canonical"
      } else {
        canonic <- "alternative.else"
        summary <- "alternative"
      }
      
      df.l <- data.frame(name = id.l,
                         uniprot = i,
                         entrez_id = e,
                         refseq.pep = pep,
                         refseq.cds = cds,
                         uniprot.isoforms = isoform,
                         isoform = summary,
                         remark = canonic,
                         uniprot.displayed = d,
                         description = desc.l,
                         category = c,
                         common = upid)
      df <- rbind(df, df.l)
    }
    refseq <- rbind(refseq, df)
  }
  return(refseq)
}
refseq <- match.refseq(displayed)
refseq

# Basic check on unique IDs from different DBs; 274
length(unique(displayed$uniprot))
length(unique(refseq$name))
length(unique(refseq$entrez_id))

# Get the RefSeq FASTA files and save to disk. Select whether to use the minor
# release version of the sequence (NP_000079.2) or not (NP_000079). Selecting
# the latter, will download the latest sequence release and name it accordingly
# in the FASTA header with the minor version, e.g. NP_000079 -> NP_000079.3. On
# disk, it will be saved without the minor release version, e.g. NP_000079.fasta.
# The database variable can be either "protein" or "nuccore", depending on what
# ID to use, but will retrieve the CDS in either case, even is "protein" was
# specified. The ret.type can be "fasta_cds_na" or "fasta".
getRefSeq <- function(refseq = NULL, database = "protein", ret.type = "fasta_cds_na", work.d = NULL, min.ver = TRUE){
  dir.create(paste0(work.d, "RefSeq"), showWarnings = F)
  for (r in refseq) {
    if (isTRUE(min.ver)) {
      id.l <- r
    } else {
      id.l <- rmdot(r)
    }
    
    if (isTRUE(file.exists(paste0(work.d, "RefSeq/", id.l, ".fasta")))) {
      print(paste0("SKIP ", id.l))
    } else {
      refseq.na <- entrez_fetch(db = database, id = id.l, rettype = ret.type)
      print(paste0("===> ", id.l))
      write(refseq.na, file = paste0(work.d, "RefSeq/", id.l, ".fasta"))
    }
  }
}
canonical.isoforms <- subset(refseq, refseq$isoform == "canonical")

# Basic check on unique IDs from different DBs; 272
length(unique(canonical.isoforms$uniprot))
length(unique(canonical.isoforms$name))
length(unique(canonical.isoforms$entrez_id))

# DOWNLOAD ####
# Using the RefSeq IDs extracted from UniProt, download the corresponding CDS and
# protein sequence from NCBI. Download point releases, as listed at UniProt and
# the latest (no point). Rerun if download fails, it will skip the files already
# downloaded.
getRefSeq(refseq = canonical.isoforms$refseq.cds, work.d = get.d, min.ver = T, database = "nuccore") # Accessed: 2026-01-07
getRefSeq(refseq = canonical.isoforms$refseq.cds, work.d = get.d, min.ver = F, database = "nuccore") # Accessed: 2026-01-07

getRefSeq(refseq = canonical.isoforms$refseq.pep, work.d = get.d, min.ver = T, ret.type = "fasta", database = "protein") # Accessed: 2026-01-07
getRefSeq(refseq = canonical.isoforms$refseq.pep, work.d = get.d, min.ver = F, ret.type = "fasta", database = "protein") # Accessed: 2026-01-07

# A function to import FASTA sequences into a list, as single string, so we compare
# with what's from UniProt later. Import both point release and latest (no point
# in ID), compare them and save the latest too, if they do not match, because sometimes
# UniProt lists RefSeq ID *.2 but the matching is *.3. This imports CDS, but removes
# the stop codons (no '*' in the FASTA file from UniProt).
importRefSeq <- function(refseq.id = NULL, refseq.d = NULL, compare.point = TRUE){
  uprs.list <- list()
  for (r in 1:nrow(refseq.id)) {
    gene_id <- refseq.id$entrez_id[r]
    
    # Start with non versioned (no point release) sequences. These will be the latest,
    # point release number will be added
    lt.fa <- read.fasta(file = paste0(refseq.d, rmdot(refseq.id$refseq.cds[r]), ".fasta"), seqtype = "DNA", whole.header = T, set.attributes = F, as.string = F)
    aa.lt <- read.fasta(file = paste0(refseq.d, rmdot(refseq.id$refseq.pep[r]), ".fasta"), seqtype = "AA", whole.header = T, set.attributes = F, as.string = F)
    
    # construct the list
    named.lt <- attr(lt.fa,"name")
    named.lt <- strsplit(named.lt, " ")
    named.lt <- named.lt[[1]][[1]]
    named.lt <- gsub("lcl\\|", "", named.lt)
    named.lt <- strsplit(named.lt, "_cds_")
    XM.lt <- named.lt[[1]][[1]]
    
    # TODO: use refseq.id$refseq.pep[r] here?
    XP.lt <- strsplit(named.lt[[1]][[2]], "_")
    XP.lt <- paste0(XP.lt[[1]][[1]], "_", XP.lt[[1]][[2]])
    
    uprs.list[[gene_id]][["metadata"]][["uniprot"]] <- refseq.id$uniprot[r]
    uprs.list[[gene_id]][["metadata"]][["displayed"]] <- refseq.id$uniprot.displayed[r]
    uprs.list[[gene_id]][["metadata"]][["description"]] <- refseq.id$description[r]
    uprs.list[[gene_id]][["metadata"]][["uniprot_name"]] <- refseq.id$name[r]
    uprs.list[[gene_id]][["metadata"]][["category"]] <- refseq.id$category[r]
    uprs.list[[gene_id]][["metadata"]][["common"]] <- refseq.id$common[r]
    uprs.list[[gene_id]][["sequences"]][[XM.lt]][["cds"]] <- lt.fa[[1]]
    uprs.list[[gene_id]][["sequences"]][[XM.lt]][["protein"]] <- aa.lt[[1]]
    
    uprs.list[[gene_id]][["sequences"]][[XM.lt]][["pep.refseq"]] <- XP.lt
    #uprs.list[[gene_id]][["sequences"]][[XM.lt]][["pep.refseq"]] <- refseq.id$refseq.pep[r]
    
    # add point release sequences, as specified. The exact specified version.
    if (isTRUE(compare.point)) {
      rs.fa <- read.fasta(file = paste0(refseq.d, refseq.id$refseq.cds[r], ".fasta"), seqtype = "DNA", whole.header = T, set.attributes = F, as.string = F)
      aa.fa <- read.fasta(file = paste0(refseq.d, refseq.id$refseq.pep[r], ".fasta"), seqtype = "AA", whole.header = T, set.attributes = F, as.string = F)
      
      if (aa.fa[length(aa.fa)] == "*") {
        aa.fa <- aa.fa[-length(aa.fa)]
      }
      
      named.rs <- attr(rs.fa,"name")
      named.rs <- strsplit(named.rs, " ")
      named.rs <- named.rs[[1]][[1]]
      named.rs <- gsub("lcl\\|", "", named.rs)
      named.rs <- strsplit(named.rs, "_cds_")
      XM.rs <- named.rs[[1]][[1]]
      XP.rs <- strsplit(named.rs[[1]][[2]], "_")
      XP.rs <- paste0(XP.rs[[1]][[1]], "_", XP.rs[[1]][[2]])
      
      uprs.list[[gene_id]][["sequences"]][[XM.rs]][["cds"]] <- rs.fa[[1]]
      uprs.list[[gene_id]][["sequences"]][[XM.rs]][["protein"]] <- aa.fa[[1]]
      uprs.list[[gene_id]][["sequences"]][[XM.rs]][["pep.refseq"]] <- XP.rs
      
      if (isTRUE(XM.rs != XM.lt)) {
        uprs.list[[gene_id]][["sequences"]][[XM.lt]][["cds"]] <- lt.fa[[1]]
        uprs.list[[gene_id]][["sequences"]][[XM.lt]][["protein"]] <- aa.lt[[1]]
        uprs.list[[gene_id]][["sequences"]][[XM.lt]][["pep.refseq"]] <- XP.lt
      }
    }
  }
  return(uprs.list)
}
irfs.all <- importRefSeq(refseq.id = canonical.isoforms, refseq.d = paste0(get.d, "RefSeq/"))

# A function to import and store the downloaded UniProt sequences into a list.
import.uniprot.list <- function(UniProt = NULL, work.d = NULL){
  proteins <- list()
  for (u in 1:nrow(UniProt)) {
    rs <- UniProt$entrez_id[u]
    uu <- UniProt$uniprot.displayed[u]
    up <- UniProt$uniprot[u]
    gn.l <- read.fasta(file = paste0(work.d, uu, ".fa"), seqtype = "AA", as.string = F, set.attributes = F, whole.header = T)
    proteins[[rs]] <- gn.l[1]
    proteins[[rs]][["displayed"]] <- uu
    proteins[[rs]][["uniprot"]] <- up
  }
  return(proteins)
}
proteins.uniprot <- import.uniprot.list(UniProt = canonical.isoforms, work.d = paste0(get.d, "UniProt/"))

# A Function to compare UniProt protein sequences against protein sequences from RefSeq
compareSeq <- function(uniprot = proteins.uniprot, refseq = NULL){
  seq_match_df <- data.frame()
  for (gene_id in names(uniprot)) {
    for (rs_id in names(refseq[[gene_id]][["sequences"]])) {
        seq_match <- seqToString(uniprot[[gene_id]][[1]]) == seqToString(refseq[[gene_id]][["sequences"]][[rs_id]][["protein"]])
        seq_match_df.l <- data.frame(gene = gene_id, uniprot = uniprot[[gene_id]][["uniprot"]], refseq = rs_id, matching = seq_match)
        seq_match_df <- rbind(seq_match_df, seq_match_df.l)
    }
  }
  return(seq_match_df)
}
compa <- compareSeq(proteins.uniprot, irfs.all)
subset(compa, compa$matching == F)

# Get only the matching sequences (IDs)
refound <- subset(compa, compa$matching == TRUE)

length(names(irfs.all)) # 272
length(unique(canonical.isoforms$entrez_id)) # 272
length(unique(canonical.isoforms$uniprot)) # 272
length(unique(compa$uniprot)) # 272
length(unique(refound$uniprot)) # 265

# View which sequences did not match and which sequences had no canonical isoform
# listed at UniProt for RefSeq
setdiff(displayed$entrez_id, unique(refound$gene)) # "23145"  "256076" "8785"   "1302"   "80781"  "127731" "375616" "80144"  "54829" 
setdiff(unique(refseq$entrez_id), names(irfs.all)) # "127731" "80144" 

# Compare the UniProt IDs of what was found above, against the original list and
# output the ones still not found.
manual_check <- setdiff(displayed$uniprot, unique(refound$uniprot))
manual_check # "A2VEC9" "A8TX70" "O95460" "P13942" "P39060" "Q5TIE3" "Q6ZWJ8" "Q86XX4" "Q9BXN1"

# A function to make an updated data.frame of core.ecm
make.core.ecm <- function(refound = NULL, refseq = NULL){
  rf.df <- data.frame()
  for (r in 1:nrow(refound)) {
    gn.id <- refound$gene[r]
    rs.id <- refound$refseq[r]
    up.id <- refound$uniprot[r]
    rf.df.l <- data.frame(name = refseq[[gn.id]][["metadata"]][["uniprot_name"]],
                          uniprot = refseq[[gn.id]][["metadata"]][["uniprot"]],
                          entrez_id = gn.id,
                          uniprot.iso = refseq[[gn.id]][["metadata"]][["displayed"]],
                          refseq.pep = refseq[[gn.id]][["sequences"]][[rs.id]][["pep.refseq"]],
                          refseq.cds = names(refseq[[gn.id]][["sequences"]][rs.id]),
                          description = refseq[[gn.id]][["metadata"]][["description"]],
                          category = refseq[[gn.id]][["metadata"]][["category"]],
                          common = refseq[[gn.id]][["metadata"]][["common"]])
    rf.df <- rbind(rf.df, rf.df.l)
  }
  return(rf.df)
}
core.ecm.update <- make.core.ecm(refound = refound, refseq = irfs.all)
core.ecm.update

length(unique(refseq$uniprot)) # 274
length(unique(refound$uniprot)) # 265
length(unique(core.ecm.update$uniprot)) # 265

# A function to generate an updated ecm.core data frame, with columns of single
# CDS/protein identifiers for the matched canonical isoform for each gene. The
# other IDs corresponding to the canonical CDS/protein at NCBI RefSeq are also
# included, in a separate column.
singleID <- function(refound = NULL, core.ecm.update = NULL){
  uniprot.unique <- unique(refound$uniprot)
  single.df <- data.frame()
  for (r in uniprot.unique) {
    up_id <- subset(core.ecm.update, core.ecm.update$uniprot == r)
    
    if (length(up_id$refseq.pep) > 1) {
      pep.alt <- toString(up_id$refseq.pep[-1])
    } else {
      pep.alt <- NA
    }
    
    if (length(up_id$refseq.cds) > 1) {
      cds.alt <- toString(up_id$refseq.cds[-1])
    } else {
      cds.alt <- NA
    }
    
    single.df.l <- data.frame(name = up_id$name[1],
                              uniprot = r,
                              entrez_id = up_id$entrez_id[1],
                              uniprot.iso = up_id$uniprot.iso[1],
                              refseq.pep = up_id$refseq.pep[1],
                              refseq.pep.alt = pep.alt,
                              refseq.cds = up_id$refseq.cds[1],
                              refseq.cds.alt = cds.alt,
                              description = up_id$description[1],
                              category = up_id$category[1],
                              common = up_id$common[1],
                              protodeviser = "UniProt")
    
    single.df <- rbind(single.df, single.df.l)
  }
  return(single.df)
}
single.df <- singleID(refound, core.ecm.update)
nrow(single.df) # 265

# A function to prepare unmatched (UniProt vs RefSeq) genes to be added to the
# final data.frame
prep.unmatched <- function(manual_check = NULL, core.ecm.auto = NULL){
  unmatched <- data.frame()
  for (u in manual_check) {
    un <- subset(core.ecm.auto, core.ecm.auto$uniprot == u)
    unmatched <- rbind(unmatched, un)
  }

  colnames(unmatched) <- c("name", "refseq.cds", "uniprot", "ensembl", "hgnc_id", "entrez_id", "description", "category")
  unmatched$common <- unmatched$name
  unmatched$name <- NA
  return(unmatched)
}
unmatched <- prep.unmatched(manual_check, core.ecm.auto)

# Inspect and check for pseudogenes (RefSeq ID not starting with NP_ or XP_; e.g. NR_).
unmatched

# Exclude the pseudogene (first row)
singleID.pseudogene <- function(pseudogene = NULL){
  df.pseudo <- data.frame("name" = pseudogene$name,
                          uniprot = pseudogene$uniprot,
                          entrez_id = pseudogene$entrez_id,
                          uniprot.iso = NA,
                          refseq.pep = NA,
                          refseq.pep.alt = NA,
                          refseq.cds = pseudogene$refseq.cds,
                          refseq.cds.alt = NA,
                          description = pseudogene$description,
                          category = pseudogene$category,
                          common = pseudogene$common,
                          protodeviser = NA)
  return(df.pseudo)
}
df.pseudo <- singleID.pseudogene(pseudogene = unmatched[1,])
df.pseudo

# A function to extract the protein ID from the corresponding CDS fasta file
getPepID <- function(cds.fa = NULL, refseq.d = NULL){
  fa <- read.fasta(file = paste0(refseq.d, cds.fa, ".fasta"), seqtype = "DNA", whole.header = T, set.attributes = F, as.string = F)
  named <- attr(fa,"name")
  named <- strsplit(named, " ")
  named <- named[[1]][[1]]
  named <- gsub("lcl\\|", "", named)
  named <- strsplit(named, "_cds_")
  XM <- named[[1]][[1]]
  XP <- strsplit(named[[1]][[2]], "_")
  XP <- paste0(XP[[1]][[1]], "_", XP[[1]][[2]])
  return(XP)
}

# Get the unmatched, using the RefSeq ID from GeneCards and import as a list
unmatched[-1,]
getRefSeq(refseq = unmatched[-1,]$refseq.cds, work.d = get.d, min.ver = F, database = "nuccore") # Accessed: 2026-01-08

# Extract the corresponding protein sequences IDs from the CDS fasta header
um.df <- c()
for (r in 1:nrow(unmatched[-1,])) {
  um.df.l <- getPepID(unmatched[-1,]$refseq.cds[r], refseq.d = paste0(get.d, "RefSeq/"))
  um.df <- c(um.df, um.df.l)
}
um.df
unmatched$refseq.pep <- c(NA, um.df)
unmatched

# Get the protein sequences
getRefSeq(refseq = unmatched[-1,]$refseq.pep, work.d = get.d, min.ver = F, ret.type = "fasta", database = "protein") # Accessed: 2026-01-08
irfs.unm <- importRefSeq(refseq.id = unmatched[-1,], refseq.d = paste0(get.d, "RefSeq/"), compare.point = F)

# Turn the "unmatched" list into a data frame, using same columns as the rest
singleID.unmatched <- function(irfs.unm = NULL){
  rf.df <- data.frame()
  for (gn.id in names(irfs.unm)) {
    rf.df.l <- data.frame(name = irfs.unm[[gn.id]][["metadata"]][["uniprot_name"]],
                          uniprot = NA,
                          entrez_id = gn.id,
                          uniprot.iso = NA,
                          refseq.pep = irfs.unm[[gn.id]][["sequences"]][[1]][["pep.refseq"]],
                          refseq.pep.alt = NA,
                          refseq.cds = names(irfs.unm[[gn.id]][["sequences"]][1]),
                          refseq.cds.alt = NA,
                          description = irfs.unm[[gn.id]][["metadata"]][["description"]],
                          category = irfs.unm[[gn.id]][["metadata"]][["category"]],
                          common = irfs.unm[[gn.id]][["metadata"]][["common"]],
                          protodeviser = "NCBI")
    rf.df <- rbind(rf.df, rf.df.l)
  }
  rf.df
}
rf.df <- singleID.unmatched(irfs.unm)

# Merge everything, including the pseudogene, to be saved as a XLSX later.
allgenes <- rbind(single.df, rf.df, df.pseudo)
allgenes

# DOWNLOAD ####
# Download the rest (newest version) with point release
getRefSeq(refseq = allgenes$refseq.cds, work.d = get.d, min.ver = T, database = "nuccore") # Accessed: 2026-01-08
getRefSeq(refseq = allgenes$refseq.pep, work.d = get.d, min.ver = T, ret.type = "fasta", database = "protein") # Accessed: 2026-01-08

# A function to collect all reference sequences from H.sapiens (CDS and protein) into a list
hsap.reference <- function(allgenes = NULL){
  hsap.list <- list()
  allgenes <- subset(allgenes, !is.na(allgenes$protodeviser))
  for (rw in 1:nrow(allgenes)) {
    gene_id <- allgenes$entrez_id[rw]
    cds_id <- allgenes$refseq.cds[rw]
    pep_id <- allgenes$refseq.pep[rw]
    gn.l <- read.fasta(file = paste0(get.d, "RefSeq/", cds_id, ".fasta"), seqtype = "DNA", as.string = F, set.attributes = F, whole.header = T)
    pt.l <- read.fasta(file = paste0(get.d, "RefSeq/", pep_id, ".fasta"), seqtype = "AA", as.string = F, set.attributes = F, whole.header = T)
    hsap.list[[gene_id]][["Homo_sapiens"]][[cds_id]][["cds"]] <- gn.l[[1]]
    hsap.list[[gene_id]][["Homo_sapiens"]][[cds_id]][["protein"]] <- pt.l[[1]]
    hsap.list[[gene_id]][["Homo_sapiens"]][[cds_id]][["pep.refseq"]] <- pep_id
  }
  return(hsap.list)
}
hsap <- hsap.reference(allgenes)
# HSAP SECTION ENDS ############################################################

# QUALITY CHECKS ###############################################################
# A function to remove sequences that have nucleotides other than a, t, g, c (e.g. n)
rmBadNA <- function(seq.list = NULL){
  for (gene_id in names(seq.list)) {
    for (latin_name in names(seq.list[[gene_id]])) {
      for (cds_id in names(seq.list[[gene_id]][[latin_name]])) {
        tes <- seq.list[[gene_id]][[latin_name]][[cds_id]][["cds"]]
        bad <- grep("FALSE", grepl("a|t|g|c", tes))
        if(isTRUE(any(bad))){
          seq.list[[gene_id]][[latin_name]][[cds_id]] <- NULL
          cat("====> Removing: ", cds_id, "\t", gene_id, "\t", latin_name, "\n")
        }
      }
    }
  }
  return(seq.list)
}
hsap.rmBadNA <- rmBadNA(hsap)

# A function to remove sequences that have amino acids other than the standard
# 20; this can check either "protein" or "translated" (needed later).
rmBadAA <- function(seq.list = NULL, aa.seq = NULL){
  for (gene_id in names(seq.list)) {
    for (latin_name in names(seq.list[[gene_id]])) {
      for (cds_id in names(seq.list[[gene_id]][[latin_name]])) {
        tes <- seq.list[[gene_id]][[latin_name]][[cds_id]][[aa.seq]]
        bad <- grep("FALSE", grepl("A|R|N|D|C|Q|E|G|H|I|L|K|M|F|P|S|T|W|Y|V", tes))
        if(isTRUE(any(bad))){
          seq.list[[gene_id]][[latin_name]][[cds_id]] <- NULL
          cat("====> Removing: ", cds_id, "\t", gene_id, "\t", latin_name, "\n")
        }
      }
    }
  }
  return(seq.list)
}
hsap.rmBadNA.rmBadAA <- rmBadAA(hsap.rmBadNA, aa.seq = "protein")

# A function to remove stop codons from CDS
rmStopCodon <- function(cds.list = NULL){
  for (gene_id in names(cds.list)) {
    for (latin_name in names(cds.list[[gene_id]])) {
      for (qseqid in names(cds.list[[gene_id]][[latin_name]])) {
        tes <- cds.list[[gene_id]][[latin_name]][[qseqid]][["cds"]]
        one <- length(tes) - 2
        two <- length(tes) - 1
        end <- length(tes)
        
        if (isTRUE(paste0(tes[c(one, two, end)], collapse = "") == "tag" ) ||
            isTRUE(paste0(tes[c(one, two, end)], collapse = "") == "tga" ) ||
            isTRUE(paste0(tes[c(one, two, end)], collapse = "") == "taa" )) {
          tes <- tes[-c(one, two, end)]
        }
        cds.list[[gene_id]][[latin_name]][[qseqid]][["cds"]] <- tes
      }
    }
  }
  return(cds.list)
}
hsap.rmBadNA.rmBadAA.rmStopCodon <- rmStopCodon(hsap.rmBadNA.rmBadAA)

# Add a translation, to be used for BLAST later
translateHsapCDS <- function(seq.list = NULL){
  for (gene_id in names(seq.list)) {
    for (latin_name in names(seq.list[[gene_id]])) {
      for (cds_id in names(seq.list[[gene_id]][[latin_name]])) {
        seq.list[[gene_id]][[latin_name]][[cds_id]][["translated"]] <- translate(seq.list[[gene_id]][[latin_name]][[cds_id]][["cds"]])
        
        if (seqToString(seq.list[[gene_id]][[latin_name]][[cds_id]][["translated"]]) ==
            seqToString(seq.list[[gene_id]][[latin_name]][[cds_id]][["protein"]])) {
          seq.list[[gene_id]][[latin_name]][[cds_id]][["match"]] <- "YES"
        } else {
          seq.list[[gene_id]][[latin_name]][[cds_id]][["match"]] <- "NO"
        }
        
      }
    }
  }
  return(seq.list)
}
hsap.rmBadNA.rmBadAA.rmStopCodon <- translateHsapCDS(hsap.rmBadNA.rmBadAA.rmStopCodon)
hsap.rmBadNA.rmBadAA.rmStopCodon <- rmBadAA(hsap.rmBadNA.rmBadAA.rmStopCodon, aa.seq = "translated")

# Check if translation matches to protein sequence from database or not (second
# loop should not return anything)
for (gene_id in names(hsap.rmBadNA.rmBadAA.rmStopCodon)) {
  for (latin_name in names(hsap.rmBadNA.rmBadAA.rmStopCodon[[gene_id]])) {
    for (cds_id in names(hsap.rmBadNA.rmBadAA.rmStopCodon[[gene_id]][[latin_name]])) {
      if (hsap.rmBadNA.rmBadAA.rmStopCodon[[gene_id]][[latin_name]][[cds_id]][["match"]] == "YES") {
        print(gene_id)
      }
    }
  }
}

for (gene_id in names(hsap.rmBadNA.rmBadAA.rmStopCodon)) {
  for (latin_name in names(hsap.rmBadNA.rmBadAA.rmStopCodon[[gene_id]])) {
    for (cds_id in names(hsap.rmBadNA.rmBadAA.rmStopCodon[[gene_id]][[latin_name]])) {
      if (hsap.rmBadNA.rmBadAA.rmStopCodon[[gene_id]][[latin_name]][[cds_id]][["match"]] == "NO") {
        print(gene_id)
      }
    }
  }
}


# A function to add amino acid length to ecm core
addLength <- function(all.genes = NULL, pep.hsap = NULL){
  for (r in 1:nrow(all.genes)) {
    gene_id <- all.genes$entrez_id[r]
    
    if ( ! is.null(pep.hsap[[gene_id]])) {
      all.genes$length.aa[r] <- getLength(pep.hsap[[gene_id]][["Homo_sapiens"]][[1]][["protein"]])
    } else {
      all.genes$length.aa[r] <- 0
    }
  }
  return(all.genes)
}
allgenes <- addLength(allgenes, hsap.rmBadNA.rmBadAA.rmStopCodon)

# Add a stamp that protein and translated sequences match (a bit silly, maybe?)
addMatch <- function(all.genes = NULL, pep.hsap = NULL){
  for (r in 1:nrow(all.genes)) {
    gene_id <- all.genes$entrez_id[r]
    all.genes$match[r] <- toString(pep.hsap[[gene_id]][["Homo_sapiens"]][[1]][["match"]])
  }
  return(all.genes)
}
allgenes <- addMatch(allgenes, hsap.rmBadNA.rmBadAA.rmStopCodon)

# Save core genes
core.ecm <- subset(allgenes, !is.na(allgenes$protodeviser))


# OBTAIN HOMOLOGUES (Placentalia) ##############################################

# WRITE ####
# Save the list of genes, to be fed to the BASH script below
#write(core.ecm$entrez_id, file = paste0(out.d, "genes.txt"), ncolumns = 1)

# BASH ####
# Run script to download data, which also checks for empty files and reruns if needed
#runSh("download_orthologues.sh") # Accessed: 2026-01-08; ncbi-datasets v18.14.0

# A function to import all genes and their isoforms for all organisms
import.genes.tsv <- function(entrez_id = NULL, work.d = NULL){
  
  # store everything in a list
  isoforms_list <- list()
  for (gene in entrez_id) {
    df <- read.table(file = paste0(work.d, gene, ".tsv"), header = T, sep = "\t", , quote = "")
    isoforms_list[[gene]] <- df
  }
  return(isoforms_list)
}
isoforms_list <- import.genes.tsv(entrez_id = core.ecm$entrez_id, work.d = genes.d)

# A function to transform the isoforms list into a dataframe
isoforms.df <- function(isoforms_list = NULL){
  genes_species <- data.frame()
  for (gene in 1:length(isoforms_list)) {
    genes_species_local <- unique(isoforms_list[[gene]][, c(2,3,4)])
    genes_species_local$gene <- names(isoforms_list[gene])
    genes_species <- rbind(genes_species, genes_species_local)
  }
  
  rownames(genes_species) <- NULL
  genes_species <- data.frame(symbol = genes_species$Symbol, gene = genes_species$gene, taxid = genes_species$Taxonomic.ID, species = gsub(" ", "_", genes_species$Taxonomic.Name))
  genes_species$taxid <- as.character(genes_species$taxid)
  
  return(genes_species)
}
genes_species <- isoforms.df(isoforms_list = isoforms_list)


# OBTAIN SPECIES TREE ##########################################################

# A function to create a matrix of genes vs species; gs: input dataframe: gene_species;
# rw: rows ("gene", "species", "taxid"), cl: columns ("gene", "species", "taxid").
# NOTE: It can be redone to work directly with numbers, as well.
genes.species.matrix <- function(gs = NULL, rw = NULL, cl = NULL){
  
  if (rw == "gene") {w = 2}
  if (rw == "taxid") {w = 3}
  if (rw == "species") {w = 4}
  if (cl == "gene") {l = 2}
  if (cl == "taxid") {l = 3}
  if (cl == "species") {l = 4}
  
  # empty matrix with defined size
  mtx <- matrix(data = 0, nrow = length(unique(gs[,w])), ncol = length(unique(gs[,l])))
  rownames(mtx) <- unique(gs[,w])
  colnames(mtx) <- unique(gs[,l])
  
  for (r in 1:nrow(gs)) {
    mtx.l <- matrix(data = 0, nrow = length(unique(gs[,w])), ncol = length(unique(gs[,l])))
    rownames(mtx.l) <- unique(gs[,w])
    colnames(mtx.l) <- unique(gs[,l])
    
    row.m <- gs[r,w]
    col.m <- gs[r,l]
    mtx.l[row.m,col.m] <- 1
    mtx <- mtx + mtx.l
  }
  
  return(mtx)
}

# Let's preview genes/species. Set gene names as rows and latin species names as columns.
gen_spe <- genes.species.matrix(genes_species, rw = "gene", cl = "species")
pheatmap(gen_spe, color = c("#eeeeec", "#729fcf"), legend_breaks = c(0, 1), border_color = "white", cluster_cols = TRUE, clustering_distance_rows = "euclidean", fontsize = 6)

# A function to calculate how many genes a taxid "covers", as a 0-1 score, from
# a matrix. In the matrix, genes MUST be as rows, and taxids as columns.
# as.percent.taxid <- function(mtx = NULL){
#   df.taxid <- data.frame()
#   for (i in colnames(mtx)) {
#     pc <- (sum(mtx[,i])/nrow(mtx))
#     df.l <- data.frame(species = i, percent = pc)
#     df.taxid <- rbind(df.taxid, df.l)
#   }
#   return(df.taxid)
# }

# A function to calculate in how many taxids a gene is found, as a 0-1 score, from
# a matrix. In the matrix, genes MUST be as rows, and taxids as columns.
as.percent.gene <- function(mtx = NULL){
  df.gene <- data.frame()
  for (i in rownames(mtx)) {
    pc <- (sum(mtx[i,])/ncol(mtx))
    df.l <- data.frame(gene = i, percent = pc)
    df.gene <- rbind(df.gene, df.l)
  }
  return(df.gene)
}

# A function to calculate in how many taxids a gene is found, as a an absolute number,
# from a matrix. In the matrix, genes MUST be as rows, and taxids as columns.
as.counts.gene <- function(mtx = NULL){
  df.gene <- data.frame()
  for (i in rownames(mtx)) {
    df.l <- data.frame(gene = i, counts = sum(mtx[i,]))
    df.gene <- rbind(df.gene, df.l)
  }
  return(df.gene)
}

# A function to order genes by number of species, from a matrix.
# mtx: input matrix
# as.what = "percent" (0-1) or "counts" (absolute)
orderGenes <- function(mtx = NULL, as.what = NULL){
  
  if (isTRUE(as.what == "percent")) {
    genes.df <- as.percent.gene(mtx)
    genes.df <- genes.df[rev(order(genes.df$percent)),]
  }
  
  if (isTRUE(as.what == "counts")) {
    genes.df <- as.counts.gene(mtx)
    genes.df <- genes.df[rev(order(genes.df$counts)),]
  }
  
  genes.df$gene <- factor(genes.df$gene, levels = genes.df$gene)
  rownames(genes.df) <- NULL
  return(genes.df)
}

# A function to plot genes, calling orderGenes
plotOrderedGenes <- function(mtx = NULL, as.what = NULL){
  genes.df <- orderGenes(mtx, as.what)
  ggplot(data=genes.df, aes(x=gene, y=counts)) + 
    geom_bar(stat="identity", color=NA, fill="blue") + 
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}
plotOrderedGenes(gen_spe, "counts")

# Genes missing in H.sap and M.mus, report here. Nothing missing in H.sap, but these are missing from M.mus:
# "137902" "22798"  "4013"   "3730"   "57829"  "346007" "22932"  "84467"  "64100"  "81578"  "266"    "25878" 
names(which(gen_spe[,"Homo_sapiens"] == 0))
names(which(gen_spe[,"Mus_musculus"] == 0))

# TIMETREE ####
species.pre <- sort(unique(genes_species[,4]))
length(species.pre) # 234

# WRITE ####
# Write species and submit them at TimeTree (https://timetree.org/).
#write(species.pre, file = paste0(tmp.d, "species.pre.txt"), ncolumns = 1)

# GNUMERIC ####
# Inspect the tree/species and create a manually curated spreadsheet
animals <- gnumeric::read.gnumeric.sheet(paste0(table.d, "species.gnumeric"), head = T)

# https://stackoverflow.com/a/67134566
species.tt <- str_subset(animals$Submitted.species.to.TimeTree, pattern = ".+")
length(species.tt) # 228

# WRITE ####
# Write species and submit them at TimeTree (https://timetree.org/).
#write(species.tt, file = paste0(tmp.d, "species.tt.txt"), ncolumns = 1)

# Subset the manually replaced species; TODO: not needed? Use "animals" directly.
##manual <- subset(animals, animals$Source.of.replacement == "Manual")

# A function to process a TimeTree nwk tree (nwk), substituting synonyms back to
# the species names of the sequences. We do not want to mess with the species
# names in the FASTA files and use what was there originally. Keep branches
# length (HyPhy is more stable), and node labels (to be used as a reference when
# comparing the global tree of 221 species to pruned trees), but unroot (pruning
# will likely mess the rooting anyway).
timetree.prep <- function(nwk = NULL, synonyms = NULL){
  timetree <- read.tree(file = nwk)
  is.rooted(timetree)
  timetree <- unroot(timetree)
  
  # Cleaning branches length may lead to a warning in HyPhy, but otherwise it
  # would crash for some genes, e.g. CO1A1 and CO1A2
  timetree$edge.length <- NULL
  timetree$node.label <- NULL
  
  plot(timetree)
  
  for (r in 1:nrow(synonyms)) {
    ori <- synonyms[r,1]
    syn <- synonyms[r,2]
    timetree$tip.label[timetree$tip.label == syn] <- ori
  }
  return(timetree)
}
tree <- timetree.prep(paste0(tmp.d, "species.tt.nwk"), animals) # Accessed: 2026-01-09

# WRITE ####
# Save the processed tree
#write.tree(tree, paste0(out.d, "species.nwk"))

# Use this to remove the sequences of these species, that were not found at TimeTree
tt.remove <- subset(animals, animals$Submitted.species.to.TimeTree == "")[,1]

# A function to extract Latin name of the species from the FASTA header. FASTA
# headers do NOT have the taxid code, that's why we use Latin names.
# .*\\[organism= matches everything up to and including "[organism=".
# (.*?) captures the text between "[organism=" and "]" (non-greedy to stop at the first "]").
# \\].* matches the closing bracket and everything after.
# \\1 keeps only the captured group, i.e., "Canis lupus familiaris".
get.species <- function(x = NULL){
  spl <- names(x)
  spl <- unlist(spl)
  spl <- gsub(".*\\[organism=(.*?)\\].*", "\\1", spl)
  spl <- gsub(" ", "_", spl)
  return(spl)
}

# A function to get RefSeq IDs
get.id <- function(x = NULL){
  spl <- names(x)
  spl <- unlist(strsplit(spl, ":"))
  return(spl[1])
}

# A function to import and store CDS into a list. This calls get.species() and get.id()
import.cds.list <- function(hsap.df = NULL, work.d = NULL){
  genes.cds <- list()
  for (g in 1:length(hsap.df$entrez_id)) {
    gene_id <- as.character(hsap.df$entrez_id[g])
    gn.l <- read.fasta(file = paste0(work.d, gene_id, ".fna"), seqtype = "DNA", as.string = F, set.attributes = F, whole.header = T)
    for (i in 1:length(gn.l)) {
      latin_name <- get.species(gn.l[i])
      seq_id <- get.id(gn.l[i])
      #seq_id <- gn.l[i]
      genes.cds[[gene_id]][[latin_name]][[seq_id]][["cds"]] <- gn.l[[i]]
    }
  }
  return(genes.cds)
}
genes.cds <- import.cds.list(hsap.df = subset(allgenes, !is.na(allgenes$protodeviser)), work.d = genes.d)

# A function to remove species that were not found at TimeTree (tt.remove)
rmMissingSpecies <- function(cds.list = NULL, rm = NULL){
  for (gene_id in names(cds.list)) {
    for (latin_name in rm) {
      cds.list[[gene_id]][[latin_name]] <- NULL
    }
  }
  return(cds.list)
}
genes.cds.rmMissingSpecies <- rmMissingSpecies(genes.cds, tt.remove)

# A function to prepare the resulting genes/species for matrix conversion. Use
# the exact same columns arrangement as in "genes_species" in the beginning.
prep.gen_spe <- function(orthol = NULL){
  mtrx <- data.frame()
  for (gene_id in names(orthol)) {
    for (latin_name in names(orthol[[gene_id]])) {
      df.l <- data.frame(symbol = NA, gene = gene_id, taxid = NA, species = latin_name)
      mtrx <- rbind(mtrx, df.l)
    }
  }
  return(mtrx)
}
gs_filt <- prep.gen_spe(genes.cds.rmMissingSpecies)

# Let's preview again
gen_spe_filt <- genes.species.matrix(gs_filt, rw = "gene", cl = "species")
pheatmap(gen_spe_filt, color = c("#eeeeec", "#729fcf"), legend_breaks = c(0, 1), border_color = "white", cluster_cols = TRUE, clustering_distance_rows = "euclidean", fontsize = 6)
plotOrderedGenes(gen_spe_filt, as.what = "counts")

genes.cds.rmBadNA <- rmBadNA(genes.cds.rmMissingSpecies)

# A function to remove species for which no sequences were left after the above.
rmEmptySpecies <- function(seq.list = NULL){
  for (gene_id in names(seq.list)) {
    for (latin_name in names(seq.list[[gene_id]])) {
      if(isEmpty(seq.list[[gene_id]][[latin_name]])){
        seq.list[[gene_id]][[latin_name]] <- NULL
        cat("====> Removing: ", gene_id, "\t", latin_name, "\n")
      }
    }
  }
  return(seq.list)
}
genes.cds.rmBadNA.rmEmptySpecies <- rmEmptySpecies(genes.cds.rmBadNA)
genes.cds.rmBadNA.rmEmptySpecies.rmStopCodon <- rmStopCodon(genes.cds.rmBadNA.rmEmptySpecies)

# A function to translate to protein from CDS. This is just for quality checks,
# so naming of the pep sequences is the same as for their cds.
translateCDS <- function(cds.list = NULL){
  genes.pep <- list()
  for (gene_id in names(cds.list)) {
    for (latin_name in names(cds.list[[gene_id]])) {
      for (seq_id in names(cds.list[[gene_id]][[latin_name]])) {
        genes.pep[[gene_id]][[latin_name]][[seq_id]][["cds"]] <- cds.list[[gene_id]][[latin_name]][[seq_id]][["cds"]]
        genes.pep[[gene_id]][[latin_name]][[seq_id]][["translated"]] <- translate(cds.list[[gene_id]][[latin_name]][[seq_id]][["cds"]])
      }
    }
  }
  return(genes.pep)
}
genes.translated <- translateCDS(genes.cds.rmBadNA.rmEmptySpecies.rmStopCodon)
genes.translated.rmBadAA <- rmBadAA(genes.translated, aa.seq = "translated")
genes.translated.rmBadAA.rmEmptySpecies <- rmEmptySpecies(genes.translated.rmBadAA)

# A function to remove genes for which only 1 species was left
rmEmptyGenes <- function(seq.list = NULL){
  for (gene_id in names(seq.list)) {
    if(length(names(seq.list[[gene_id]])) == 1){
      seq.list[[gene_id]] <- NULL
    }
  }
  return(seq.list)
}
genes.translated.rmBadAA.rmEmptySpecies.rmEmptyGenes <- rmEmptyGenes(genes.translated.rmBadAA.rmEmptySpecies)

# A function to remove genes for which no sequence for H.sap was found/left
# (there shouldn't be any).
rmEmptyHsap <- function(seq.list = NULL){
  for (gene_id in names(seq.list)) {
    if(! any(grepl("Homo_sapiens", names(seq.list[[gene_id]])))){
      seq.list[[gene_id]] <- NULL
    }
  }
  return(seq.list)
}
genes.translated.rmBadAA.rmEmptySpecies.rmEmptyGenes <- rmEmptyHsap(genes.translated.rmBadAA.rmEmptySpecies.rmEmptyGenes)

# A function to remove H.sap from the original list. This will be used to run
# BLASTP against H.sap.
rmHsap <- function(seq.list = NULL){
  for (gene_id in names(seq.list)) {
    seq.list[[gene_id]][["Homo_sapiens"]] <- NULL
  }
  return(seq.list)
}
genes.nohsap <- rmHsap(genes.translated.rmBadAA.rmEmptySpecies.rmEmptyGenes)

batchFasta <- function(seq.list = NULL, seq.type = NULL){
  seq.local <- list()
  for (gene_id in names(seq.list)) {
    for (latin_name in names(seq.list[[gene_id]])) {
      for (qseqid in names(seq.list[[gene_id]][[latin_name]])) {
        seq.local[[gene_id]][[latin_name]][[qseqid]] <- seq.list[[gene_id]][[latin_name]][[qseqid]][[seq.type]]
      }
    }
  }
  return(seq.local)
}
cds.hsap <- batchFasta(hsap.rmBadNA.rmBadAA.rmStopCodon, "cds")
pep.hsap <- batchFasta(hsap.rmBadNA.rmBadAA.rmStopCodon, "translated") # or use "protein", they are the same for the dataset
cds.nohsap <- batchFasta(genes.nohsap, "cds")
pep.nohsap <- batchFasta(genes.nohsap, "translated")

# A function to save FASTA files from either hsap or nohsap
saveFastaFromList <- function(seq.list, folder.out){
  dir.create(folder.out, showWarnings = F)
  for (gene_id in names(seq.list)) {
    dir.create(paste0(folder.out, "/", gene_id))
    for (latin_name in names(seq.list[[gene_id]])) {
      write.fasta(seq.list[[gene_id]][[latin_name]], names = attr(seq.list[[gene_id]][[latin_name]], "name"), file.out = paste0(folder.out, "/", gene_id, "/", latin_name, ".fasta"), nbchar = 60)
    }
  }
}

# WRITE ####
#saveFastaFromList(pep.hsap, paste0(blast.d, "hsap")) # H.sap canonical sequences, BLAST will run here
#saveFastaFromList(pep.nohsap, paste0(blast.d, "nohsap")) # Orthologues sequences
#write(names(pep.nohsap), file = paste0(out.d, "genes-blast.txt"), ncolumns = 1) # Gene names found in H.sap to a txt first.

# check which gene was excluded
setdiff(names(pep.hsap), names(pep.nohsap)) # 22932

# BASH ####
#runSh("blastPrep.sh") # ncbi-blast-plus, v2.17.0
#runSh("blastRunParallel.sh")

# A function calling find in bash, to check for empty output files
find.empty <- function(blast.dir = NULL){
  run("find", c(blast.dir, "-empty"), stdout = paste0(tmp.d, "empty.txt"))
  empty <- read.table(paste0(tmp.d, "empty.txt"), header = F, sep = "/")
  empty <- empty[, c(10,12)]
  colnames(empty) <- c("gene_id", "latin_name")
  empty$latin_name <- gsub(".csv", "", x = empty$latin_name)
  return(empty)
}
empty <- find.empty(paste0(blast.d, "hsap"))
empty

# A function to remove orthologues that had empty results after BLAST
rmOrtho <- function(ortho.list = NULL, empty = NULL){
  for (r in 1:nrow(empty)) {
    gene_id <- as.character(empty[r,1])
    latin_name <- empty[r,2]
    ortho.list[[gene_id]][[latin_name]] <- NULL
  }
  return(ortho.list)
}
pep.nohsap.rmOrtho <- rmOrtho(pep.nohsap, empty)

# A function to import the BLAST results as dataframes and sort by "bitscore"
import.blast.results <- function(ortho.list = NULL){
  blastResults <- list()
  for (gene_id in names(ortho.list)) {
    for (latin_name in names(ortho.list[[gene_id]])) {
      blast_data <- read.csv(paste0(blast.d, "hsap/", gene_id, "/results/", latin_name, ".csv"), header = F)
      
      # blastp output format:
      # -outfmt="10 qseqid sacc pident length mismatch gapopen qstart qend sstart send evalue bitscore score ppos qlen slen gaps qcovs"
      colnames(blast_data) <- c("qseqid", "sacc", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "score", "ppos", "qlen", "slen", "gaps", "qcovs")
      blastResults[[gene_id]][[latin_name]] <- blast_data
    }
  }
  return(blastResults)
}
blastResults <- import.blast.results(pep.nohsap.rmOrtho)

# A function to filter BLASTP results and create a composite score. Options:
#
# Input results, as a list (blast.list = blastResults)
# E-value (expect value) cutoff (eval = 1e-5)
# Percentage of identical matches cutoff (pide = 60)
# Length ratio (minimum and maximum) between Q/S (ratio_margin = 0.1)
# Maximum allowed absolute start of the Q and S alignments (max_start_abs = 35)
# Maximum allowed relative start of the Q and S alignments (max_start_rel = 0.1)
# Query coverage as a percent (qcover = NULL)
filter.blast.results <- function(blast.list = NULL, eval = NULL, pide = NULL, ratio_margin = NULL, max_start_abs = NULL, max_start_rel = NULL, qcover = NULL){
  filter.list <- list()
  for (gene_id in names(blast.list)) {
    filter.list[[gene_id]] <- list()
    
    for (latin_name in names(blast.list[[gene_id]])) {
      
      blast_data <- blast.list[[gene_id]][[latin_name]]
      
      # filter on E-value and percent identity
      blast_data <- blast_data[blast_data$evalue < eval, ]
      blast_data <- blast_data[blast_data$pident > pide, ]
      
      # filter on query coverage
      blast_data <- blast_data[blast_data$qcovs > qcover, ]
      
      # filter on absolute and relative start position of the alignment
      blast_data <- blast_data[blast_data$qstart <= max_start_abs, ]
      blast_data <- blast_data[blast_data$sstart <= max_start_abs, ]
      blast_data <- blast_data[blast_data$qstart / blast_data$qlen <= max_start_rel, ]
      blast_data <- blast_data[blast_data$sstart / blast_data$slen <= max_start_rel, ]

      # Calculate Query/Subject length ratio
      blast_data$length_ratio <- blast_data$qlen / blast_data$slen
      
      # Calculate composite score, combining Bitscore and percentage identity
      blast_data$composite_score <- blast_data$bitscore * (blast_data$pident / 100)
      
      min_len_ratio <- 1 - ratio_margin
      max_len_ratio <- 1 + ratio_margin

      blast_data <- blast_data[blast_data$length_ratio >= min_len_ratio & blast_data$length_ratio <= max_len_ratio, ]
      
      blast_data_sorted <- blast_data[order(-blast_data$composite_score, 
                                            blast_data$gaps, 
                                            blast_data$gapopen), ]
      
      rownames(blast_data_sorted) <- NULL
      if (nrow(blast_data_sorted) > 0){
        filter.list[[gene_id]][[latin_name]] <- blast_data_sorted
      } else{
        filter.list[[gene_id]][[latin_name]] <- NULL
      }
    }
  }
  return(filter.list)
}
blastResultsFilter <- filter.blast.results(blast.list = blastResults, eval = 1e-5, pide = 75, ratio_margin = 0.1, max_start_abs = 35, max_start_rel = 0.1, qcover = 40)

# A function to collect reference sequence from Hsap first then the best match
# from BLAST per species.
getOrthol <- function(ref.hsap = NULL, cds.list = NULL, blast.hits = NULL){
  
  # start by collecting reference sequences from H.sap
  orthologues <- list()
  for (gene_id in names(ref.hsap)) {
    for (sseq_id in names(ref.hsap[[gene_id]][["Homo_sapiens"]])) {
      orthologues[[gene_id]][["Homo_sapiens"]][[sseq_id]] <- ref.hsap[[gene_id]][["Homo_sapiens"]][[sseq_id]] 
    }
  }
  
  # Now add the rest
  for (gene_id in names(blast.hits)) {
    for (latin_name in names(blast.hits[[gene_id]])) {
      qseqid <- blast.hits[[gene_id]][[latin_name]][1,1]
      orthologues[[gene_id]][[latin_name]][[qseqid]] <- cds.list[[gene_id]][[latin_name]][[qseqid]]
    }
  }
  return(orthologues)
}
orthol.left <- getOrthol(ref.hsap = cds.hsap, cds.list = cds.nohsap, blast.hits = blastResultsFilter)

# A function to remove genes with a number of orthologues below threshold
minOrthol <- function(orthol.list = NULL, orthol.min = NULL){
  for (gene_id in names(orthol.list)) {
    if (length(orthol.list[[gene_id]]) < orthol.min) {
      orthol.list[[gene_id]] <- NULL
    }
  }
  return(orthol.list)
}
orthol <- minOrthol(orthol.left, 1)

# function to translate yet again the CDS left, one per species
translateOrthol <- function(cds.list = NULL){
  pep.list <- list()
  for (gene_id in names(cds.list)) {
    for (latin_name in names(cds.list[[gene_id]])) {
      for (seq_id in names(cds.list[[gene_id]][[latin_name]])) {
        pep.list[[gene_id]][[latin_name]][[seq_id]] <- translate(cds.list[[gene_id]][[latin_name]][[seq_id]])
      }
    }
  }
  return(pep.list)
}
orthol.aa <- translateOrthol(orthol)

# Count the number of orthologues/species for each gene, left after filtering
countSpecies <- function(all.genes = NULL, orthologues = NULL){
  all.genes$species <- 0
  for (r in 1:nrow(all.genes)) {
    gene_id <- all.genes$entrez_id[r]
    num_spe <- length(orthologues[[gene_id]])
    all.genes$species[r] <- num_spe
  }
  return(all.genes)
}
allgenes <- countSpecies(allgenes, orthol)

# WRITE ORTHOLOGUES ENTREZ IDs, more than 1 species ####
# These are the 272 genes that made it
gene_ids <- subset(allgenes, allgenes$species > 1)[,3]
#write(gene_ids, file = paste0(out.d, "genes-orthologues.txt"), ncolumns = 1)

gs_orthol <- prep.gen_spe(orthol.aa)
gen_spe_orthol <- genes.species.matrix(gs_orthol, rw = "gene", cl = "species")
head(gen_spe_orthol)

pheatmap(gen_spe_orthol, color = c("#eeeeec", "#729fcf"), legend_breaks = c(0, 1), border_color = "white", cluster_cols = TRUE, clustering_distance_rows = "euclidean", fontsize = 6)
plotOrderedGenes(gen_spe_orthol, "counts")
names(which(gen_spe_orthol[,"Homo_sapiens"] == 0))
names(which(gen_spe_orthol[,"Mus_musculus"] == 0))

# A function to prune global tree for each gene, leaving just the species where
# the gene sequences are found/left
genesTrees <- function(cds.list = NULL, cds.tree = NULL){
  trees.list <- list()
  for (gene_id in names(cds.list)) {
    species_to_keep <- names(cds.list[[gene_id]])
    pruned_tree <- keep.tip(cds.tree, species_to_keep)
    trees.list[[gene_id]] <- pruned_tree
  }
  return(trees.list)
}

# Prune trees of "final" orthologues list and save into a list
orthol.genesTrees <- genesTrees(orthol, cds.tree = tree)

# WRITE ####
# Save trees to disk
saveTheTrees <- function(trees.list = NULL, tree.path = NULL){
  for (gene_id in names(trees.list)) {
    write.tree(trees.list[[gene_id]], file = paste0(tree.path, gene_id, ".tre"))
  }
}
saveTheTrees(trees.list = orthol.genesTrees, tree.path = trees.d)

# A function to save FASTA files
saveOrtho <- function(genes.list = NULL, folder.out = NULL, seq.type = NULL){
  
  # simplify list
  genes.list.simple <- list()
  for (gene_id in names(genes.list)) {
    for (latin_name in names(genes.list[[gene_id]])) {
      genes.list.simple[[gene_id]][[latin_name]] <- genes.list[[gene_id]][[latin_name]][[1]]
    }
  }
  
  # now save the simplified list
  dir.create(folder.out, showWarnings = F)
  for (gene_id in names(genes.list.simple)) {
    write.fasta(genes.list.simple[[gene_id]], names = attr(genes.list.simple[[gene_id]], "name"), file.out = paste0(folder.out, "/", gene_id, ".", seq.type, ".fasta"), nbchar = 60)
  }
}

# WRITE ####
saveOrtho(genes.list = orthol, folder.out = ortho.d, seq.type = "cds")
saveOrtho(genes.list = orthol.aa, folder.out = ortho.d, seq.type = "pep")

# BASH ####
##runSh("mafftRun.sh")
#runSh("muscle5Run.sh")
#runSh("gblocksRun.sh")
#runSh("pal2nalRun.sh")

# Make a list with species/orthologues refseq df
orthol_xm <- list()
for (gene_id in names(orthol)) {
  df.id <- data.frame()
  for (latin_name in names(orthol[[gene_id]])) {
    df.id.l <- data.frame(Latin_name = latin_name, RefSeq = names(orthol[[gene_id]][[latin_name]]))
    df.id <- rbind(df.id, df.id.l)
  }
  orthol_xm[[gene_id]] <- df.id
}

# Collect codon alignments and save as an object
fasta.codon <- list()
for (gene_id in names(orthol)) {
  fasta.codon[[gene_id]] <- readLines(paste0(ortho.d, gene_id, ".codon.fasta"))
}

# PROTODEVISER ####
# A function to batch generate JSON schemes from core.ecm input. Several color
# schemes can be specified, as a c().
id.JSON.batch <- function(ecm = NULL, color.scheme = NULL){
  proteins.json <- list()
  
  ecm <- subset(ecm, !is.na(ecm$protodeviser))
  
  for (r in 1:nrow(ecm)) {
    gene_id <- ecm$entrez_id[r]
    if (ecm$protodeviser[r] == "UniProt") {
      prot_id <- ecm$uniprot[r]
      data_id <- "UniProt"
    } else if (ecm$protodeviser[r] == "NCBI"){
      prot_id <- ecm$refseq.pep[r]
      data_id <- "NCBI"
    }
    for (c in color.scheme) {
      proteins.json[[c]][[gene_id]] <- id.JSON(input = prot_id, database = data_id, gradient = c)
    }
  }
  return(proteins.json)
}
#proteins.json <- id.JSON.batch(ecm = allgenes, color.scheme = c("rainbow"))

# A function to batch save the JSON schemes from core.ecm input
write.JSON.batch <- function(proteins.json = NULL, work.d = NULL){
  for (c in names(proteins.json)) {
    for (gene_id in names(proteins.json[[c]])) {
      dir.create(paste0(work.d, c), showWarnings = F)
      write(proteins.json[[c]][[gene_id]], file = paste0(work.d, c, "/", gene_id, ".json"))
    }
  }
}
write.JSON.batch(proteins.json, jsons.d)

# Create an object of protein features as a table
proteins.feature.table <- list()
for (gen in gene_ids) {
  proteins.feature.table[[gen]] <- json.TABLE(jsonlite::fromJSON(proteins.json[["rainbow"]][[gen]]))
}

# SEQUENCES LENGTH ####
# Get the amino acid sequence length of the proteins. What will be the maximum
# amino acids length for the plots X-axis?
compareLength <- function(allgenes = NULL, proteins.json = NULL){
  len.df <- subset(allgenes, !is.na(allgenes$protodeviser))
  len.df <- data.frame(entrez_id = len.df$entrez_id, length_translated = len.df$length.aa)
  len.up <- data.frame()
  
  for (gene_id in names(proteins.json[[1]])) {
    fjs <- fromJSON(proteins.json[[1]][[gene_id]])
    ln.l <- fjs$length
    len.up.l <- data.frame(entrez_id = gene_id, length_aa = ln.l)
    len.up <- rbind(len.up, len.up.l)
  }
  
  max(len.up$length_aa)
  mer <- merge(len.df, len.up)
  mer$match <- mer$length_translated == mer$length_aa
  
  return(mer)
}
mer <- compareLength(allgenes, proteins.json)
max(mer$length_aa)

# A function to import H.sap sequence with gaps and Gblocks scores
msaReference <- function(msaFile = NULL, gblocks = NULL, species = NULL){
  msa.prot <- read.fasta(msaFile)
  hsap.prot <- msa.prot[[species]][1:length(msa.prot[[species]])]
  gbl.prot <- read.fasta(gblocks)
  gb <- gbl.prot[["Gblocks"]][1:length(gbl.prot[["Gblocks"]])]
  gb <- sub("[.]", "low", gb)
  gb <- sub("#", "high", gb)
  msa.ref <- data.frame(aa = hsap.prot, msa.quality = gb)
  return(msa.ref)
}

# A function to output the reference sequence without MSA gaps
seqReference <- function(msaFile = NULL, gblocks = NULL, species = NULL){
  seqRefSpecies <- msaReference(msaFile, gblocks, species)
  seqRefSpecies <- data.frame(subset(seqRefSpecies, seqRefSpecies$aa != "-"))
  seqRefSpecies$site <- 1:nrow(seqRefSpecies)
  return(seqRefSpecies)
}

# A function to batch import MSA and unaligned sequence of the reference species
batch.import.ref <- function(gene_ids = NULL, species = "Homo_sapiens"){
  ref.list <- list()
  for (g in gene_ids) {
    ref.list[[g]][["msa"]] <- msaReference(paste0(ortho.d, g, ".muscle5.fasta"),
                                           paste0(gbloc.d, g, ".Gblocks.fa"),
                                           species)
    ref.list[[g]][["seq"]] <- seqReference(paste0(ortho.d, g, ".muscle5.fasta"),
                                           paste0(gbloc.d, g, ".Gblocks.fa"),
                                           species)
  }
  return(ref.list)
}
seq.reference <- batch.import.ref(gene_ids)

# SAVE KEY OBJECTS FROM THE PREPARATION STEPS ####
save(allgenes, file = paste0(out.d, "allgenes.rda"))
save(blastResults, file = paste0(out.d, "blastResults.rda"))
save(cds.hsap, file = paste0(out.d, "cds.hsap.rda"))
save(pep.hsap, file = paste0(out.d, "pep.hsap.rda"))
save(pep.nohsap, file = paste0(out.d, "pep.nohsap.rda"))
save(genes.cds.rmBadNA.rmEmptySpecies.rmStopCodon, file = paste0(out.d, "genes.cds.rmBadNA.rmEmptySpecies.rmStopCodon.rda"))
save(blastResultsFilter, file = paste0(out.d, "blastResultsFilter.rda"))
save(orthol, file = paste0(out.d, "orthol.rda"))
save(orthol.aa, file = paste0(out.d, "orthol.aa.rda"))
save(tree, file = paste0(out.d, "tree.rda"))
save(orthol.genesTrees, file = paste0(out.d, "orthol.genesTrees.rda"))
save(gene_ids, file = paste0(out.d, "gene_ids.rda"))
save(orthol_xm, file = paste0(out.d, "orthol_xm.rda"))
save(fasta.codon, file = paste0(out.d, "fasta.codon.rda"))
save(proteins.json, file = paste0(out.d, "proteins.json.rda"))
save(seq.reference, file = paste0(out.d, "seq.reference.rda"))
save(proteins.feature.table, file = paste0(out.d, "proteins.feature.table.rda"))
save(gen_spe_orthol, file = paste0(out.d, "gen_spe_orthol.rda"))

# LOAD: Import the key objects ####
load(file = paste0(out.d, "allgenes.rda"))
load(file = paste0(out.d, "blastResults.rda"))
load(file = paste0(out.d, "cds.hsap.rda"))
load(file = paste0(out.d, "pep.hsap.rda"))
load(file = paste0(out.d, "pep.nohsap.rda"))
load(file = paste0(out.d, "genes.cds.rmBadNA.rmEmptySpecies.rmStopCodon.rda"))
load(file = paste0(out.d, "blastResultsFilter.rda"))
load(file = paste0(out.d, "orthol.rda"))
load(file = paste0(out.d, "orthol.aa.rda"))
load(file = paste0(out.d, "tree.rda"))
load(file = paste0(out.d, "orthol.genesTrees.rda"))
load(file = paste0(out.d, "gene_ids.rda"))
load(file = paste0(out.d, "orthol_xm.rda"))
load(file = paste0(out.d, "fasta.codon.rda"))
load(file = paste0(out.d, "proteins.json.rda"))
load(file = paste0(out.d, "seq.reference.rda"))
load(file = paste0(out.d, "proteins.feature.table.rda"))
load(file = paste0(out.d, "gen_spe_orthol.rda"))

# BASH: RUN HYPHY ANALYSES #####################################################
#runSh("hyphyRun.sh")
################################################################################

# IMPORT HyPhy RESULTS ####
# For testing, use these, overwriting temporarily gene_ids:
# TODO: comment out, when all analyses are done.
gene_ids <- list.dirs(path = hyphy.d, full.names = F, recursive = F)

# A function to import JSON results
import.json <- function(gene_id = NULL, method = NULL, dir.d = hyphy.d){
  json <- read_json(paste0(dir.d, gene_id, "/", method, "/", gene_id, ".", method, ".json"), simplifyDataFrame = T)
  return(json)
}

# A function to batch import JSON results from both MEME and FUBAR
batch.import.json <- function(gene_ids = NULL, dir.d = hyphy.d){
  json.list <- list()
  for (g in gene_ids) {
    json.list[[g]][["meme"]] <- import.json(g, "meme", dir.d)
    json.list[[g]][["fubar"]] <- import.json(g, "fubar", dir.d)
  }
  return(json.list)
}
json.results <- batch.import.json(gene_ids)


# A function to trim MEME results from JSON file, leaving only p.plus and P.val
trim.meme <- function(gene_id = NULL, json_ls = json.results){
  json <- json_ls[[gene_id]][["meme"]]
  nsites <- json$input$`number of sites`
  
  df <- data.frame()
  for (i in 1:nsites) {
    p_plus <- json$MLE$content[[1]][[i]][[5]]
    P_val <- json$MLE$content[[1]][[i]][[7]]
    if (!is.null(p_plus)) {
      df.i <- data.frame(msa = i, p.plus = p_plus, P.val = P_val)
    }else{
      df.i <- data.frame(msa = i, p.plus = NA, P.val = NA)
    }
    df <- rbind(df, df.i)
  }
  return(df)
}

# A function to trim FUBAR results from JSON file, leaving only beta.alpha and Bayes
trim.fubar <- function(gene_id = NULL, json_ls = json.results){
  json <- json.results[[gene_id]][["fubar"]]
  nsites <- json$input$`number of sites`
  
  df <- data.frame()
  for (i in 1:nsites) {
    beta_alpha <- json$MLE$content[[1]][[i]][[3]]
    P_a_gt_b <- json$MLE$content[[1]][[i]][[4]]
    P_b_gt_a <- json$MLE$content[[1]][[i]][[5]]
    if (!is.null(beta_alpha)) {
      df.i <- data.frame(msa = i, beta.alpha = beta_alpha, P_a.gt.b = P_a_gt_b, P_b.gt.a = P_b_gt_a)
      df <- rbind(df, df.i)}
    else{
      df.i <- data.frame(msa = i, beta.alpha = NA, P_a.gt.b = NA, P_b.gt.a = NA)
      df <- rbind(df, df.i)
    }
  }
  return(df)
}

# A function to batch import trimmed results from both MEME and FUBAR
batch.import.trimmed <- function(gene_ids = NULL, json_ls = json.results){
  trimmed.list <- list()
  for (g in gene_ids) {
    trimmed.list[[g]][["meme"]] <- trim.meme(g)
    trimmed.list[[g]][["fubar"]] <- trim.fubar(g)
  }
  return(trimmed.list)
}
trimmed.hyphy <- batch.import.trimmed(gene_ids)

# HyPhy PROCESSING ####
# A function to process MEME results
processMEME <- function(gene_id = NULL,
                        tr_res = trimmed.hyphy,
                        species = "Homo_sapiens",
                        seq.ref = seq.reference,
                        pvalue = 0.1,
                        pplus = 0.001,
                        quality = "both"){
  
  df <- tr_res[[gene_id]][["meme"]]
  
  msa <- seq.ref[[gene_id]][["msa"]]
  df <- data.frame(msa, df)
  
  df <- subset(df, df$aa != "-")
  df$site <- 1:nrow(df)
  rownames(df) <- df$site
  
  # do we keep sites from bad alignment regions?
  if (quality == "high") {
    bad <- "low"
  }else if (quality == "low"){
    bad = "high"
  }else{
    bad = "both"
  }
  
  selection <- data.frame()
  for (r in 1:nrow(df)) {
    if(is.na(df$p.plus[r])){
      sel <- c(df$msa[r], df$msa.quality[r], NA,  NA, "na")
    }else if(df$P.val[r] <= pvalue & df$p.plus[r] > pplus & df$msa.quality[r] != bad){
      sel <- c(df$msa[r], df$msa.quality[r], df$p.plus[r], df$P.val[r], "pos")
    }else{
      sel <- c(df$msa[r], df$msa.quality[r], NA, NA, "ns")
    }
    selection <- rbind(selection, sel)
  }
  
  colnames(selection) <- c("msa", "quality", "p.plus", "P.value", "selection")
  
  selection$site <- rownames(selection)
  
  seq <- seq.ref[[gene_id]][["seq"]]
  export <- merge.data.frame(selection, seq, by = "site", all.y = T, sort = F)
  export <- data.frame(msa = export$msa, ref = export$aa, quality = export$quality, ref.site = export$site, p.plus = export$p.plus, P.value = export$P.value, episodic = export$sel)
  
  export$p.plus <- as.numeric(export$p.plus)
  export$P.value <- as.numeric(export$P.value)
  export$msa <- as.numeric(export$msa)
  export$ref.site <- as.numeric(export$ref.site)
  
  return(export)
}

# A function to process FUBAR results
processFUBAR <- function(gene_id = NULL,
                         tr_res = trimmed.hyphy,
                         species = "Homo_sapiens",
                         seq.ref = seq.reference,
                         bayes = 0.9,
                         quality = "both"){
  
  df <- tr_res[[gene_id]][["fubar"]]
  
  msa <- seq.reference[[gene_id]][["msa"]]
  df <- data.frame(msa, df)
  
  df <- subset(df, df$aa != "-")
  df$site <- 1:nrow(df)
  rownames(df) <- df$site
  
  if (quality == "high") {
    bad <- "low"
  }else if (quality == "low"){
    bad = "high"
  }else{
    bad = "both"
  }
  
  # loop through the df and filter
  selection <- data.frame()
  for (r in 1:nrow(df)) {
    if(is.na(df$beta.alpha[r])){
      sel <- c(df$msa[r], df$msa.quality[r], NA, NA, NA, "na")
    }else if(df$beta.alpha[r] > 0 & df$P_b.gt.a[r] >= bayes & df$msa.quality[r] != bad){
      sel <- c(df$msa[r], df$msa.quality[r], df$beta.alpha[r], df$P_a.gt.b[r], df$P_b.gt.a[r], "pos")
    }else if(df$beta.alpha[r] < 0 & df$P_a.gt.b[r] >= bayes & df$msa.quality[r] != bad){
      sel <- c(df$msa[r], df$msa.quality[r], df$beta.alpha[r], df$P_a.gt.b[r], df$P_b.gt.a[r], "neg")
    }else{
      sel <- c(df$msa[r], df$msa.quality[r], NA, NA, NA, "ns")
    }
    selection <- rbind(selection, sel)
  }
  
  colnames(selection) <- c("msa", "quality", "beta.alpha", "P_a.gt.b", "P_b.gt.a", "sel")
  selection$site <- rownames(selection)
  
  seq <- seq.reference[[gene_id]][["seq"]]
  export <- merge.data.frame(selection, seq, by = "site", all.y = T, sort = F)
  export <- data.frame(msa = export$msa, ref = export$aa, quality = export$quality, ref.site = export$site, beta.alpha = export$beta.alpha, Bayes.neg = export$P_a.gt.b, Bayes.pos = export$P_b.gt.a, selection = export$sel)
  
  export$beta.alpha <- as.numeric(export$beta.alpha)
  export$Bayes.neg <- as.numeric(export$Bayes.neg)
  export$Bayes.pos <- as.numeric(export$Bayes.pos)
  export$ref.site <- as.numeric(export$ref.site)
  export$msa <- as.numeric(export$msa)
  
  return(export)
}

# A function to batch process HyPhy results, using default filtering parameters.
# This is only needed for global statistics.
processed.hyphy <- function(gene.ids = NULL){
  processed <- list()
  for (g in gene.ids) {
    processed[[g]][["meme"]] <- processMEME(g)
    processed[[g]][["fubar"]] <- processFUBAR(g)
  }
  return(processed)
}
processed <- processed.hyphy(gene_ids)

# A function to calculate global statistics of selection across genes: expressed
# as a percentage of negative, pervasive positive and episodic selected sites,
# compared to protein length from H.sap.
global.selection <- function(processed.hyphy = NULL){
  df <- data.frame()
  for (gene_id in names(processed.hyphy)) {
    m <- processed.hyphy[[gene_id]][["meme"]]
    l <- nrow(m)
    e <- nrow(subset(m, m$episodic == "pos"))
    e <- (e / l) * 100
    
    f <- processed.hyphy[[gene_id]][["fubar"]]
    p <- nrow(subset(f, f$selection == "pos"))
    n <- nrow(subset(f, f$selection == "neg"))
    
    p <- (p / l) * 100
    n <- (n / l) * 100
    
    df.l <- data.frame(entrez_id = gene_id, episodic = e, pervasive = p, negative = n)
    df <- rbind(df, df.l)
  }
  return(df)
}
global.selection.df <- global.selection(processed)

# Update allgenes object, adding selection statistics:
allgenessel <- merge.data.frame(allgenes, global.selection.df, all = T)
allgenessel <- data.frame(common = allgenessel$common,
                          name = allgenessel$name,
                          uniprot = allgenessel$uniprot,
                          uniprot.iso = allgenessel$uniprot.iso,
                          entrez_id = allgenessel$entrez_id,
                          refseq.pep = allgenessel$refseq.pep,
                          refseq.cds = allgenessel$refseq.cds,
                          refseq.pep.alt = allgenessel$refseq.pep.alt,
                          refseq.cds.alt = allgenessel$refseq.cds.alt,
                          length.aa = allgenessel$length.aa,
                          description = allgenessel$description,
                          category = allgenessel$category,
                          species = allgenessel$species,
                          protodeviser = allgenessel$protodeviser,
                          episodic = allgenessel$episodic,
                          pervasive = allgenessel$pervasive,
                          negative = allgenessel$negative)
allgenessel <- allgenessel[order(allgenessel$common), ]
rownames(allgenessel) <- NULL

# A function to save the names of genes, found for more than 1 species. Same as,
# gene_ids, but using gene names.
minAllgenes <- function(allg = allgenes, min.sp = 1){
  ag <- subset(allg, allg$species > min.sp)
  ag <- sort(ag$common)
  return(ag)
}
allnames <- minAllgenes()

# SAVE KEY OBJECTS FROM THE HyPhy POST-RUN STEPS ####
#save(json.results, file = paste0(out.d, "json.results.rda"))
#save(trimmed.hyphy, file = paste0(out.d, "trimmed.hyphy.rda"))
#save(processed, file = paste0(out.d, "processed.rda"))
#save(allgenessel, file = paste0(out.d, "allgenessel.rda"))
#save(allnames, file = paste0(out.d, "allnames.rda"))

# LOAD KEY OBJECTS FROM THE HyPhy POST-RUN STEPS ####
load(paste0(out.d, "json.results.rda"))
load(paste0(out.d, "trimmed.hyphy.rda"))
load(paste0(out.d, "processed.rda"))
load(paste0(out.d, "allgenessel.rda"))
load(paste0(out.d, "allnames.rda"))
