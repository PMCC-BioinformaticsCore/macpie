fetch_normalised_counts<-function(data = NULL, method = NULL, batch = NULL, k = NULL) {

  # Helper function to validate input data
  validate_inputs <- function(data, annotation, batch, k) {
    if (!inherits(data, "Seurat")) {
      stop("Error: argument 'data' must be a Seurat or TidySeurat object.")
    }
    method <- if (is.null(method)) "limma-voom" else method

    batch <- if (is.null(batch)) 1 else batch

    #parameter k for RUVseq-based methods
    k <- if (is.null(k)) 2 else k

    list(data = data, annotation = annotation, batch = batch, k = k)
  }

  # Validate inputs
  validated <- validate_inputs(data, annotation, batch, k)
  annotation <- validated$annotation
  batch <- validated$batch
  k <- validated$k
  data <- validated$data

  #1. raw data
  if(method == "raw"){
    norm_data <- data@assays$RNA$counts
  }

  #2. log normalised (divided by total counts, multiplied by 10000, log1p transformed)
  if(method == "logNorm"){
    lognorm <- NormalizeData(mac, normalization.method = "LogNormalize", scale.factor = 10000)
    norm_data <- lognorm@assays$RNA$data
  }

  #3. CPM - counts per million, not lognormalised
  if(method == "cpm"){
    cpm <- NormalizeData(mac, normalization.method = "RC", scale.factor = 1e6)
    norm_data <- cpm@assays$RNA$data
  }

  #4. CLR - centered log ratio transformation
  if(method == "clr"){
    clr <- NormalizeData(mac, normalization.method = "CLR",margin = 2)
    norm_data <- clr@assays$RNA$data
  }

  #5. SCT from Seurat
  if(method == "SCT"){
    sct <-
      mac %>%
      SCTransform(do.scale = T,
                  return.only.var.genes = F,
                  vars.to.regress = c("percent.mt","Plate_ID"),
                  verbose = FALSE)
    norm_data<-sct@assays$SCT$data
  }

  #5. DESeq2
  if(method == "DEseq2"){
    coldata<-data.frame(batch = batch, condition=data$Treatment_conc)

    dds <- DESeqDataSetFromMatrix(countData = data@assays$RNA$counts,
                                  colData = coldata,
                                  design = ~ condition + batch)
    keep <- rowSums(counts(dds) > 10) >= 2
    dds <- dds[keep, ]
    dds <- estimateSizeFactors(dds)
    norm_data <- counts(dds, normalized=TRUE)
  }

  #6. edgeR
  edgeR<-DGEList(
    counts = data[["raw"]]+1,
    samples = coldata_controls$condition,
    group = coldata_controls$condition
  )
  edgeR <- calcNormFactors(edgeR,methods="TMMwsp")
  edgeR <- estimateDisp(edgeR,
                        design = model.matrix(~coldata_controls$condition+coldata_controls$batch))
  data[["edgeR"]] <- cpm(edgeR, log=FALSE)


  #6. RUVseq normalisation with/without ercc
  filtered <- data[["raw"]]+1
  genes <- rownames(filtered)[!grepl("^ERCC", rownames(filtered))]

  #plotRLE(set1, outline=FALSE)
  #k defines number of sources of variation, two have been chosen for row and column
  set <- newSeqExpressionSet(as.matrix(filtered),
                             phenoData = data.frame(condition=coldata_controls[,2], row.names=colnames(filtered)))
  set <- betweenLaneNormalization(set, which="upper")

  differences <- makeGroups(pData(set)$condition)

  #ruvS
  set3 <- RUVs(set, cIdx=genes, k=2, differences)
  data[["RUVseqS_k2"]]<-normCounts(set3)

  #use genes that haven't been modified
  fit <- glmFit(edgeR, design)
  res <- residuals(fit, type="deviance")
  set4 <- RUVr(set, genes, k=1, res)
  data[["RUVr_TMM_k1"]] <- normCounts(set4)
  set4 <- RUVr(set, genes, k=2, res)
  data[["RUVr_TMM_k2"]] <- normCounts(set4)
  save(data,file=normalized_data_file)

}
