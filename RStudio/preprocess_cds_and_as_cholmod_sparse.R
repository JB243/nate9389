preprocess_cds <- function(cds,
                           method = c('PCA', "LSI"),
                           num_dim = 50,
                           norm_method = c("log", "size_only", "none"),
                           use_genes = NULL,
                           pseudo_count = NULL,
                           scaling = TRUE,
                           verbose = FALSE,
                           build_nn_index = FALSE,
                           nn_control = list()) {

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "method must be one of 'PCA' or 'LSI'")
  method <- match.arg(method)

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(norm_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "norm_method must be one of 'log', 'size_only' or 'none'")
  norm_method <- match.arg(norm_method)

  assertthat::assert_that(assertthat::is.count(num_dim))

  if(!is.null(use_genes)) {
    assertthat::assert_that(is.character(use_genes))
    assertthat::assert_that(all(use_genes %in% row.names(rowData(cds))),
                            msg = paste("use_genes must be NULL, or all must",
                            "be present in the row.names of rowData(cds)"))
  }
  assertthat::assert_that(!is.null(size_factors(cds)),
             msg = paste("You must call estimate_size_factors before calling",
                         "preprocess_cds."))
  assertthat::assert_that(sum(is.na(size_factors(cds))) == 0,
                          msg = paste("One or more cells has a size factor of",
                                      "NA."))

  if(build_nn_index) {
    nn_control <- set_nn_control(mode=1,
                                 nn_control=nn_control,
                                 nn_control_default=get_global_variable('nn_control_annoy_cosine'),
                                 nn_index=NULL,
                                 k=NULL,
                                 verbose=verbose)
  }

  #ensure results from RNG sensitive algorithms are the same on all calls
  set.seed(2016)
  FM <- normalize_expr_data(cds, norm_method, pseudo_count)

  if (nrow(FM) == 0) {
    stop("all rows have standard deviation zero")
  }

  if (!is.null(use_genes)) {
    FM <- FM[use_genes, ]
  }

  fm_rowsums = Matrix::rowSums(FM)
  FM <- FM[is.finite(fm_rowsums) & fm_rowsums != 0, ]

  #
  # Notes:
  #   o  the functions save_transform_models/load_transform_models
  #      expect that the reduce_dim_aux slot consists of a S4Vectors::SimpleList
  #      that stores information about methods with the elements
  #        reduce_dim_aux[[method]][['model']] for the transform elements
  #        reduce_dim_aux[[method]][[nn_method]] for the nn index
  #      and depends on the elements within model and nn_method.
  #
  if(method == 'PCA') {
    cds <- initialize_reduce_dim_metadata(cds, 'PCA')
    cds <- initialize_reduce_dim_model_identity(cds, 'PCA')

    if (verbose) message("Remove noise by PCA ...")

    # Initialize variables
    preproc_res <- NULL
    rotation_matrix <- NULL
    sdev_values <- NULL

    # Determine the dimension of the feature matrix
    dim_FM <- min(dim(FM)) - 1

    # Use irlba if the number of dimensions is less than 20% of the matrix dimension, else use svd
    if (num_dim <= 0.2 * dim_FM) {
        irlba_res <- irlba::irlba(Matrix::t(FM), n = min(num_dim, dim_FM), center = scaling, scale. = scaling)
        preproc_res <- irlba_res$x
        irlba_rotation <- irlba_res$rotation
        rotation_matrix <- irlba_res$rotation
        sdev_values <- irlba_res$sdev
    } else {
        svd_res <- svd(Matrix::t(FM))
        irlba_res = svd_res
        preproc_res <- svd_res$u[, 1:num_dim] %*% diag(svd_res$d[1:num_dim])
        irlba_rotation <- svd_res$v[, 1:num_dim]
        rotation_matrix <- svd_res$v[, 1:num_dim]
        sdev_values <- svd_res$d[1:num_dim]
    }

    row.names(preproc_res) <- colnames(cds)
    SingleCellExperiment::reducedDims(cds)[[method]] <- as.matrix(preproc_res)

    row.names(rotation_matrix) <- rownames(FM)


    # we need svd_v downstream so
    # calculate gene_loadings in cluster_cells.R
    cds@reduce_dim_aux[['PCA']][['model']][['num_dim']] <- num_dim
    cds@reduce_dim_aux[['PCA']][['model']][['norm_method']] <- norm_method
    cds@reduce_dim_aux[['PCA']][['model']][['use_genes']] <- use_genes
    cds@reduce_dim_aux[['PCA']][['model']][['pseudo_count']] <- pseudo_count
    cds@reduce_dim_aux[['PCA']][['model']][['svd_v']] <- irlba_rotation
    cds@reduce_dim_aux[['PCA']][['model']][['svd_sdev']] <- irlba_res$sdev
    cds@reduce_dim_aux[['PCA']][['model']][['svd_center']] <- irlba_res$center
    cds@reduce_dim_aux[['PCA']][['model']][['svd_scale']] <- irlba_res$svd_scale
    # Note that prop_var_expl is the fraction of variance explained by the retained
    # PCs, not the fraction of total variance.
    cds@reduce_dim_aux[['PCA']][['model']][['prop_var_expl']] <- irlba_res$sdev^2 / sum(irlba_res$sdev^2)

    matrix_id <- get_unique_id(SingleCellExperiment::reducedDims(cds)[['PCA']])
    counts_identity <- get_counts_identity(cds)

    cds <- set_reduce_dim_matrix_identity(cds, 'PCA',
                                          'matrix:PCA',
                                          matrix_id,
                                          counts_identity[['matrix_type']],
                                          counts_identity[['matrix_id']],
                                          'matrix:PCA',
                                          matrix_id)
    cds <- set_reduce_dim_model_identity(cds, 'PCA',
                                         'matrix:PCA',
                                         matrix_id,
                                         'none',
                                         'none')

    if( build_nn_index ) {
      nn_index <- make_nn_index(subject_matrix=SingleCellExperiment::reducedDims(cds)[[method]], nn_control=nn_control, verbose=verbose)
      cds <- set_cds_nn_index(cds=cds, reduction_method=method, nn_index=nn_index, verbose=verbose)
    }
    else
      cds <- clear_cds_nn_index(cds=cds, reduction_method=method, nn_method='all')

  } else if(method == "LSI") {
    cds <- initialize_reduce_dim_metadata(cds, 'LSI')
    cds <- initialize_reduce_dim_model_identity(cds, 'LSI')

#    preproc_res <- tfidf(FM)
    tfidf_res <- tfidf(FM)
    preproc_res <- tfidf_res[['tf_idf_counts']]

    num_col <- ncol(preproc_res)
    irlba_res <- irlba::irlba(Matrix::t(preproc_res),
                              nv = min(num_dim,min(dim(FM)) - 1))

    preproc_res <- irlba_res$u %*% diag(irlba_res$d)
    row.names(preproc_res) <- colnames(cds)
    SingleCellExperiment::reducedDims(cds)[[method]] <- as.matrix(preproc_res)

    irlba_rotation = irlba_res$v
    row.names(irlba_rotation) = rownames(FM)
    cds@reduce_dim_aux[['LSI']][['model']][['num_dim']] <- num_dim
    cds@reduce_dim_aux[['LSI']][['model']][['norm_method']] <- norm_method
    cds@reduce_dim_aux[['LSI']][['model']][['use_genes']] <- use_genes
    cds@reduce_dim_aux[['LSI']][['model']][['pseudo_count']] <- pseudo_count
    cds@reduce_dim_aux[['LSI']][['model']][['log_scale_tf']] <- tfidf_res[['log_scale_tf']]
    cds@reduce_dim_aux[['LSI']][['model']][['frequencies']] <- tfidf_res[['frequencies']]
    cds@reduce_dim_aux[['LSI']][['model']][['scale_factor']] <- tfidf_res[['scale_factor']]
    cds@reduce_dim_aux[['LSI']][['model']][['col_sums']] <- tfidf_res[['col_sums']]
    cds@reduce_dim_aux[['LSI']][['model']][['row_sums']] <- tfidf_res[['row_sums']]
    cds@reduce_dim_aux[['LSI']][['model']][['num_cols']] <- tfidf_res[['num_cols']]
    cds@reduce_dim_aux[['LSI']][['model']][['svd_v']] <- irlba_rotation
    cds@reduce_dim_aux[['LSI']][['model']][['svd_sdev']] <- irlba_res$d/sqrt(max(1, num_col - 1))

    # we need svd_v downstream so
    # calculate gene_loadings in cluster_cells.R

    matrix_id <- get_unique_id(SingleCellExperiment::reducedDims(cds)[['LSI']])
    counts_identity <- get_counts_identity(cds)

    cds <- set_reduce_dim_matrix_identity(cds, 'LSI',
                                          'matrix:LSI',
                                          matrix_id,
                                          counts_identity[['matrix_type']],
                                          counts_identity[['matrix_id']],
                                          'matrix:LSI',
                                          matrix_id)
    cds <- set_reduce_dim_model_identity(cds, 'LSI',
                                         'matrix:LSI',
                                         matrix_id,
                                         'none',
                                         'none')

    if( build_nn_index ) {
      nn_index <- make_nn_index(subject_matrix=SingleCellExperiment::reducedDims(cds)[[method]], nn_control=nn_control, verbose=verbose)
      cds <- set_cds_nn_index(cds=cds, reduction_method=method, nn_index=nn_index, verbose=verbose)
    }
    else
      cds <- clear_cds_nn_index(cds=cds, reduction_method=method, nn_method='all')
  }

  if(!is.null(cds@reduce_dim_aux[['Aligned']]) && !is.null(cds@reduce_dim_aux[['Aligned']][['model']][['beta']])) {
    cds@reduce_dim_aux[['Aligned']][['model']][['beta']] <- NULL
  }

  cds
}


# Helper function to normalize the expression data prior to dimensionality
# reduction
normalize_expr_data <- function(cds,
                                norm_method = c("log", "size_only", "none"),
                                pseudo_count = NULL) {
  norm_method <- match.arg(norm_method)

  FM <- SingleCellExperiment::counts(cds)

  # If we're going to be using log, and the user hasn't given us a
  # pseudocount set it to 1 by default.
  if (is.null(pseudo_count)){
    if(norm_method == "log")
      pseudo_count <- 1
    else
      pseudo_count <- 0
  }

  if (norm_method == "log") {
    # If we are using log, normalize by size factor before log-transforming

    FM <- Matrix::t(Matrix::t(FM)/size_factors(cds))

    if (pseudo_count != 1 || is_sparse_matrix(SingleCellExperiment::counts(cds)) == FALSE){
      FM <- FM + pseudo_count
      FM <- log2(FM)
    } else {
      FM@x = log2(FM@x + 1)
    }

  } else if (norm_method == "size_only") {
    FM <- Matrix::t(Matrix::t(FM)/size_factors(cds))
    FM <- FM + pseudo_count
  }
  return (FM)
}

# Andrew's tfidf
tfidf <- function(count_matrix, frequencies=TRUE, log_scale_tf=TRUE,
                  scale_factor=100000, block_size=2000e6) {
  # Use either raw counts or divide by total counts in each cell
  if (frequencies) {
    # "term frequency" method
    col_sums <- Matrix::colSums(count_matrix)
    tf <- Matrix::t(Matrix::t(count_matrix) / col_sums)
  } else {
    # "raw count" method
    col_sums <- NA
    tf <- count_matrix
  }

  # Either TF method can optionally be log scaled
  if (log_scale_tf) {
    if (frequencies) {
      tf@x = log1p(tf@x * scale_factor)
    } else {
      tf@x = log1p(tf@x * 1)
    }
  }

  # IDF w/ "inverse document frequency smooth" method
  num_cols <- ncol(count_matrix)
  row_sums <- Matrix::rowSums(count_matrix > 0)
  idf = log(1 + num_cols / row_sums)

  # Try to just to the multiplication and fall back on delayed array
  # TODO hopefully this actually falls back and not get jobs killed in SGE
  tf_idf_counts = tryCatch({
    tf_idf_counts = tf * idf
    tf_idf_counts
  }, error = function(e) {
    print(paste("TF*IDF multiplication too large for in-memory, falling back",
                "on DelayedArray."))
    options(DelayedArray.block.size=block_size)
    DelayedArray:::set_verbose_block_processing(TRUE)

    tf = DelayedArray::DelayedArray(tf)
    idf = as.matrix(idf)

    tf_idf_counts = tf * idf
    tf_idf_counts
  })

  rownames(tf_idf_counts) = rownames(count_matrix)
  colnames(tf_idf_counts) = colnames(count_matrix)
  tf_idf_counts = methods::as(tf_idf_counts, "sparseMatrix")
  return(list(tf_idf_counts=tf_idf_counts, frequencies=frequencies, log_scale_tf=log_scale_tf, scale_factor=scale_factor, col_sums=col_sums, row_sums=row_sums, num_cols=num_cols))
}


# Functions for matrix and model identities.


initialize_counts_metadata <- function(cds) {
  assertthat::assert_that(methods::is(cds, 'cell_data_set'),
                          msg=paste('cds parameter is not a cell_data_set'))

  if(is.null(SingleCellExperiment::int_metadata(cds))) {
    SingleCellExperiment::int_metadata(cds) <- list()
  }
  SingleCellExperiment::int_metadata(cds)[['counts_metadata']] <- list()

  return(cds)
}


set_counts_identity <- function(cds, matrix_type, matrix_id) {
  assertthat::assert_that(methods::is(cds, 'cell_data_set'),
                          msg=paste('cds parameter is not a cell_data_set'))

  assertthat::assert_that(!is.null(SingleCellExperiment::int_metadata(cds)[['counts_metadata']]),
    msg = paste0('call initialize_counts_metadata() before set_counts_identity()'))                  

  SingleCellExperiment::int_metadata(cds)[['counts_metadata']][['identity']][['matrix_type']] <- matrix_type
  SingleCellExperiment::int_metadata(cds)[['counts_metadata']][['identity']][['matrix_id']] <- matrix_id

  return(cds)
}


get_counts_identity <- function(cds) {
  assertthat::assert_that(methods::is(cds, 'cell_data_set'),
                          msg=paste('cds parameter is not a cell_data_set'))

  if(is.null(SingleCellExperiment::int_metadata(cds)[['counts_metadata']])) {
    initialize_counts_metadata(cds)
  }

  if(!is.null(SingleCellExperiment::int_metadata(cds)[['counts_metadata']][['identity']][['matrix_id']])) {
    matrix_id <- SingleCellExperiment::int_metadata(cds)[['counts_metadata']][['identity']][['matrix_id']]
  } else {
    matrix_id <- 'none'
  }
  if(!is.null(SingleCellExperiment::int_metadata(cds)[['counts_metadata']][['identity']][['matrix_type']])) {
    matrix_type <- SingleCellExperiment::int_metadata(cds)[['counts_metadata']][['identity']][['matrix_type']]
  } else {
    matrix_type <- 'matrix:counts'
  }

  return(list(matrix_id=matrix_id, matrix_type=matrix_type))
}


# Note: SingleCellExperiment::int_metadata(cds) requires a list
initialize_reduce_dim_metadata <- function(cds, reduction_method=c('PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP')) {
  assertthat::assert_that(methods::is(cds, 'cell_data_set'),
                          msg=paste('cds parameter is not a cell_data_set'))

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(reduction_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "reduction_method must be one of 'PCA', 'LSI', 'Aligned', 'tSNE', or 'UMAP'")
  reduction_method <- match.arg(reduction_method)

  if(is.null(SingleCellExperiment::int_metadata(cds))) {
    SingleCellExperiment::int_metadata(cds) <- list()
  }
  if(is.null(SingleCellExperiment::int_metadata(cds)[['reduce_dim_metadata']])) {
    SingleCellExperiment::int_metadata(cds)[['reduce_dim_metadata']] <- list()
  }
  SingleCellExperiment::int_metadata(cds)[['reduce_dim_metadata']][[reduction_method]] <- list()
  SingleCellExperiment::int_metadata(cds)[['reduce_dim_metadata']][[reduction_method]][['identity']] <- list()

  return(cds)
}


set_reduce_dim_matrix_identity <- function(cds, reduction_method=c('PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP'),
                                           matrix_type, matrix_id,
                                           prev_matrix_type, prev_matrix_id,
                                           model_type, model_id) {
  assertthat::assert_that(methods::is(cds, 'cell_data_set'),
                          msg=paste('cds parameter is not a cell_data_set'))

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(reduction_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "reduction_method must be one of 'PCA', 'LSI', 'Aligned', 'tSNE', or 'UMAP'")
  reduction_method <- match.arg(reduction_method)

  if(is.null(SingleCellExperiment::int_metadata(cds)[['reduce_dim_metadata']][[reduction_method]])) {
    cds <- initialize_reduce_dim_metadata(cds=cds, reduction_method=reduction_method)
  }

  SingleCellExperiment::int_metadata(cds)[['reduce_dim_metadata']][[reduction_method]][['identity']][['matrix_type']] <- matrix_type
  SingleCellExperiment::int_metadata(cds)[['reduce_dim_metadata']][[reduction_method]][['identity']][['matrix_id']] <- matrix_id
  SingleCellExperiment::int_metadata(cds)[['reduce_dim_metadata']][[reduction_method]][['identity']][['prev_matrix_type']] <- prev_matrix_type
  SingleCellExperiment::int_metadata(cds)[['reduce_dim_metadata']][[reduction_method]][['identity']][['prev_matrix_id']] <- prev_matrix_id
  SingleCellExperiment::int_metadata(cds)[['reduce_dim_metadata']][[reduction_method]][['identity']][['model_type']] <- model_type
  SingleCellExperiment::int_metadata(cds)[['reduce_dim_metadata']][[reduction_method]][['identity']][['model_id']] <- model_id

  return(cds)
}


get_reduce_dim_matrix_identity <- function(cds, reduction_method=c('PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP')) {
  assertthat::assert_that(methods::is(cds, 'cell_data_set'),
                          msg=paste('cds parameter is not a cell_data_set'))

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(reduction_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "reduction_method must be one of 'PCA', 'LSI', 'Aligned', 'tSNE', or 'UMAP'")
  reduction_method <- match.arg(reduction_method)

  if(is.null(SingleCellExperiment::int_metadata(cds)[['reduce_dim_metadata']][[reduction_method]])) {
    cds <- initialize_reduce_dim_metadata(cds=cds, reduction_method=reduction_method)
    return(list(identity_exists=FALSE))
  }

  return(list(identity_exists=TRUE,
              matrix_type=SingleCellExperiment::int_metadata(cds)[['reduce_dim_metadata']][[reduction_method]][['identity']][['matrix_type']],
              matrix_id=SingleCellExperiment::int_metadata(cds)[['reduce_dim_metadata']][[reduction_method]][['identity']][['matrix_id']],
              prev_matrix_type=SingleCellExperiment::int_metadata(cds)[['reduce_dim_metadata']][[reduction_method]][['identity']][['prev_matrix_type']],
              prev_matrix_id=SingleCellExperiment::int_metadata(cds)[['reduce_dim_metadata']][[reduction_method]][['identity']][['prev_matrix_id']],
              model_type=SingleCellExperiment::int_metadata(cds)[['reduce_dim_metadata']][[reduction_method]][['identity']][['model_type']],
              model_id=SingleCellExperiment::int_metadata(cds)[['reduce_dim_metadata']][[reduction_method]][['identity']][['model_id']]))
}


initialize_reduce_dim_model_identity <- function(cds, reduction_method=c('PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP')) {
  assertthat::assert_that(methods::is(cds, 'cell_data_set'),
                          msg=paste('cds parameter is not a cell_data_set'))

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(reduction_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "reduction_method must be one of 'PCA', 'LSI', 'Aligned', 'tSNE', or 'UMAP'")
  reduction_method <- match.arg(reduction_method)

  if(is.null(cds@reduce_dim_aux[[reduction_method]])) {
    cds@reduce_dim_aux[[reduction_method]] <- S4Vectors::SimpleList()
  }
  cds@reduce_dim_aux[[reduction_method]][['model']] <- S4Vectors::SimpleList()
  cds@reduce_dim_aux[[reduction_method]][['model']][['identity']] <- S4Vectors::SimpleList()

  return(cds)
}


set_reduce_dim_model_identity <- function(cds, reduction_method=c('PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP'),
                                          model_type, model_id,
                                          prev_model_type, prev_model_id,
                                          model_path='none') {
  assertthat::assert_that(methods::is(cds, 'cell_data_set'),
                          msg=paste('cds parameter is not a cell_data_set'))

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(reduction_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "reduction_method must be one of 'PCA', 'LSI', 'Aligned', 'tSNE', or 'UMAP'")
  reduction_method <- match.arg(reduction_method)

  if(is.null(cds@reduce_dim_aux[[reduction_method]][['model']][['identity']])) {
    cds <- initialize_reduce_dim_model_identity(cds=cds, reduction_method=reduction_method)
  }

  cds@reduce_dim_aux[[reduction_method]][['model']][['identity']][['model_type']] <- model_type
  cds@reduce_dim_aux[[reduction_method]][['model']][['identity']][['model_id']] <- model_id
  cds@reduce_dim_aux[[reduction_method]][['model']][['identity']][['prev_model_type']] <- prev_model_type
  cds@reduce_dim_aux[[reduction_method]][['model']][['identity']][['prev_model_id']] <- prev_model_id
  cds@reduce_dim_aux[[reduction_method]][['model']][['identity']][['model_path']] <- model_path

  global_variable_name <- list(PCA='reduce_dim_pca_model_version',
                               LSI='reduce_dim_lsi_model_version',
                               Aligned='reduce_dim_aligned_model_version',
                               tSNE='reduce_dim_tsne_model_version',
                               UMAP='reduce_dim_umap_model_version')

  cds@reduce_dim_aux[[reduction_method]][['model']][['identity']][['model_version']] <- get_global_variable(global_variable_name[[reduction_method]])

  return(cds)
}


get_reduce_dim_model_identity <- function(cds, reduction_method=c('PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP')) {
  assertthat::assert_that(methods::is(cds, 'cell_data_set'),
                          msg=paste('cds parameter is not a cell_data_set'))

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(reduction_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "reduction_method must be one of 'PCA', 'LSI', 'Aligned', 'tSNE', or 'UMAP'")
  reduction_method <- match.arg(reduction_method)

  if(is.null(cds@reduce_dim_aux[[reduction_method]][['model']][['identity']])) {
    cds <- initialize_reduce_dim_model_identity(cds=cds, reduction_method=reduction_method)
    return(list(identity_exists=FALSE))
  }

  return(list(identity_exists=TRUE,
              model_type=cds@reduce_dim_aux[[reduction_method]][['model']][['identity']][['model_type']],
              model_id=cds@reduce_dim_aux[[reduction_method]][['model']][['identity']][['model_id']],
              prev_model_type=cds@reduce_dim_aux[[reduction_method]][['model']][['identity']][['prev_model_type']],
              prev_model_id=cds@reduce_dim_aux[[reduction_method]][['model']][['identity']][['prev_model_id']],
              model_path=cds@reduce_dim_aux[[reduction_method]][['model']][['identity']][['model_path']],
              version=cds@reduce_dim_aux[[reduction_method]][['model']][['identity']][['model_version']]))
}


set_model_identity_path <- function(cds, reduction_method=c('PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP'), model_path='none') {
  assertthat::assert_that(methods::is(cds, 'cell_data_set'),
                          msg=paste('cds parameter is not a cell_data_set'))

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(reduction_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "reduction_method must be one of 'PCA', 'LSI', 'Aligned', 'tSNE', or 'UMAP'")
  reduction_method <- match.arg(reduction_method)

  if(is.null(cds@reduce_dim_aux[[reduction_method]][['model']][['identity']])) {
    cds <- initialize_reduce_dim_model_identity(cds=cds, reduction_method=reduction_method)
  }

  if(!is.null(cds@reduce_dim_aux[[reduction_method]][['model']][['identity']][['model_path']])) {
    cds@reduce_dim_aux[[reduction_method]][['model']][['identity']][['model_path']] <- model_path
  } else {
    cds@reduce_dim_aux[[reduction_method]][['model']][['identity']][['model_path']] <- 'none'
  }

  return(cds)
}


identity_as_string <- function( object_id ) {
  if(is.character(object_id)) {
    object_id_string <- object_id
  }
  else
  if(!is.null(object_id[['checksum']])) {
    object_id_string <- sprintf('%s  dim: %s', object_id[['checksum']], paste(object_id[['dim']], collapse=' '))
  }
  else {
    object_id_string <- object_id
  }

  return( object_id_string)
}


#' @title Report matrix and model identity information.
#'
#' @description Write the cell_data_set matrix and model 
#'   identity information to stdout.
#'
#' @details A matrix identity is a checksum that is stored 
#'   in the cell_data_set when a reduced dimension matrix is 
#'   created and when certain functions read count matrices
#'   into the cell_data_set, such as load_mm_data(). At the 
#'   same time, the same checksum is stored as the model 
#'   identity in order to link the model to its matrix.
#'   
#'   Additionally, Monocle3 stores the identity of the matrix
#'   from which the matrix was made. For example, in the case
#'   of a UMAP reduced dimension matrix made from a PCA
#'   reduced dimension matrix, the cell_data_set has the
#'   identities of both the UMAP and the PCA matrices. The
#'   UMAP identity is stored as 'matrix_id' and the PCA as
#'   'prev_matrix_id'. Similarly, the model and the previous 
#'   model identities are stored as 'model_id' and 
#'   'prev_model_id'. This allows one to trace a matrix
#'   to its origin, which may be helpful when a cell_data_set
#'   is partially reprocessed; for example, if preprocess_cds()
#'   is re-run but reduce_dimension() is not. Also, it may be
#'   helpful when transform models are loaded with the
#'   load_transform_models() function, in which case
#'   the matrix and model identities will differ.
#'   
#'   The identity of the model used to transform a matrix
#'   is stored with the matrix identity information as 
#'   'model_id'. Ordinarily, the matrix 'matrix_id' and
#'   'model_id' and the corresponding model 'model_id' will
#'   have the same string value. However, they differ when
#'   the preprocess_transform() and reduce_dim_transform()
#'   functions are used to transform a matrix.
#' 
#'   Notes:
#'   * Certain file and directory paths may be stored in the
#'   cell_data_set as identifiers.
#'
#'   * Checksums are calculated using the digest function in
#'   the digest package. The matrix dimensions are stored
#'   with the checksum.
#'
#'   * Matrix transformations such as subsetting and row and
#'   or column reordering do not affect the matrix identity.
#'
#'   * The matrix identity string is stored in the internal
#'   metadata slot of the cell_data_set and the model
#'   identity string is stored in the model object in the
#'   cds@reduce_dim_aux slot of the cell_data_set.
#'
#' @param cds the cell_data_set to use.
#'
#' @return Write identity information to stdout.
#'
#' @examples
#'   \donttest{
#'     cds <- load_a549()
#'     cds <- preprocess_cds(cds)
#'     cds <- reduce_dimension(cds)
#'     identity_table(cds)
#'   }
#'
#' @export
identity_table <- function(cds) {
  assertthat::assert_that(methods::is(cds, 'cell_data_set'),
                          msg=paste('cds parameter is not a cell_data_set'))

  reduction_methods <- c('PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP')
  counts_identity <- get_counts_identity(cds)
  write(sprintf('Count matrix identity'), stdout())
  write(sprintf('  matrix_id: %s', identity_as_string(counts_identity[['matrix_id']])), stdout())
  write(sprintf('  matrix_type: %s', counts_identity[['matrix_type']]), stdout())
  write('', stdout())

  write(sprintf('Reduced dimension matrix identity'), stdout())
  for(reduction_method in reduction_methods) {
    matrix_identity <- get_reduce_dim_matrix_identity(cds, reduction_method)
    if(!is.null(SingleCellExperiment::reducedDims(cds)[[reduction_method]]) && matrix_identity[['identity_exists']]) {
      write(sprintf('  %s', reduction_method), stdout())
      write(sprintf('    %s\t%s', 'matrix_type', matrix_identity[['matrix_type']]), stdout())
      write(sprintf('    %s\t%s', 'matrix_id', identity_as_string(matrix_identity[['matrix_id']])), stdout())
      write(sprintf('    %s\t%s', 'prev_matrix_type', matrix_identity[['prev_matrix_type']]), stdout())
      write(sprintf('    %s\t%s', 'prev_matrix_id', identity_as_string(matrix_identity[['prev_matrix_id']])), stdout())
      write(sprintf('    %s\t%s', 'model_type', matrix_identity[['model_type']]), stdout())
      write(sprintf('    %s\t%s', 'model_id', identity_as_string(matrix_identity[['model_id']])), stdout())
    } else {
       write(sprintf('  %s\n    none', reduction_method), stdout())
    }
    write('', stdout())
  }

  write(sprintf('Reduced dimension model identity'), stdout())
  for(reduction_method in reduction_methods) {
    model_identity <- get_reduce_dim_model_identity(cds, reduction_method)
    if(model_identity[['identity_exists']]) {
      write(sprintf('  %s', reduction_method), stdout())
      write(sprintf('    %s\t%s', 'model_type', model_identity[['model_type']]), stdout())
      write(sprintf('    %s\t%s', 'model_id', identity_as_string(model_identity[['model_id']])), stdout())
      write(sprintf('    %s\t%s', 'prev_model_type', model_identity[['prev_model_type']]), stdout())
      write(sprintf('    %s\t%s', 'prev_model_id', identity_as_string(model_identity[['prev_model_id']])), stdout())
      write(sprintf('    %s\t%s', 'model_path', model_identity[['model_path']]), stdout())
      write(sprintf('    %s\t%s', 'model_version', model_identity[['model_version']]), stdout())
    } else {
       write(sprintf('  %s\n    none', reduction_method), stdout())
    }
    write('', stdout())
  }
}

# Test whether a matrix is one of our supported sparse matrices
is_sparse_matrix <- function(x){
  class(x) %in% c("dgCMatrix", "dgTMatrix", "lgCMatrix")
}

#' Function to calculate size factors for single-cell RNA-seq data
#'
#' @param cds The cell_data_set
#' @param round_exprs A logic flag to determine whether or not the expression
#'   value should be rounded
#' @param method A string to specify the size factor calculation approach.
#'   Options are "mean-geometric-mean-total" (default),
#'   "mean-geometric-mean-log-total".
#'
#' @return Updated cell_data_set object with a new colData column called
#'   'Size_Factor'.
#'
#' @examples
#'   \donttest{
#'     cds <- load_a549()
#'     colData(cds)[['Size_Factor']] <- NULL
#'     cds <- estimate_size_factors(cds)
#'   }
#'
#' @export
estimate_size_factors <- function(cds,
                                  round_exprs=TRUE,
                                  method=c("mean-geometric-mean-total",
                                           'mean-geometric-mean-log-total'))
{
  method <- match.arg(method)
  if(any(Matrix::colSums(SingleCellExperiment::counts(cds)) == 0)) {
    warning("Your CDS object contains cells with zero reads. ",
                  "This causes size factor calculation to fail. Please remove ",
                  "the zero read cells using ",
                  "cds <- cds[,Matrix::colSums(exprs(cds)) != 0] and then ",
                  "run cds <- estimate_size_factors(cds)")
    return(cds)
  }
  if (is_sparse_matrix(SingleCellExperiment::counts(cds))){
    size_factors(cds) <- estimate_sf_sparse(SingleCellExperiment::counts(cds),
                                            round_exprs=round_exprs,
                                            method=method)
  }else{
    size_factors(cds) <- estimate_sf_dense(SingleCellExperiment::counts(cds),
                                           round_exprs=round_exprs,
                                           method=method)
  }
  return(cds)
}

# Estimate size factors for each column, given a sparseMatrix from the Matrix
# package
estimate_sf_sparse <- function(counts,
                               round_exprs=TRUE,
                               method="mean-geometric-mean-total"){
  if (round_exprs)
    counts <- round(counts)

  if(method == 'mean-geometric-mean-total') {
    cell_total <- Matrix::colSums(counts)
    sfs <- cell_total / exp(mean(log(cell_total)))
  }else if(method == 'mean-geometric-mean-log-total') {
    cell_total <- Matrix::colSums(counts)
    sfs <- log(cell_total) / exp(mean(log(log(cell_total))))
  }

  sfs[is.na(sfs)] <- 1
  sfs
}

# Estimate size factors for each column, given a matrix
estimate_sf_dense <- function(counts,
                              round_exprs=TRUE,
                              method="mean-geometric-mean-total"){

  CM <- counts
  if (round_exprs)
    CM <- round(CM)
  if(method == "mean-geometric-mean-log-total") {
    cell_total <- apply(CM, 2, sum)
    sfs <- log(cell_total) / exp(mean(log(log(cell_total))))
  } else if(method == 'mean-geometric-mean-total') {
    cell_total <- apply(CM, 2, sum)
    sfs <- cell_total / exp(mean(log(cell_total)))
  }

  sfs[is.na(sfs)] <- 1
  sfs
}

sparse_apply <- function(Sp_X, MARGIN, FUN, convert_to_dense, ...){
  if (convert_to_dense){
    if (MARGIN == 1){
      Sp_X <- Matrix::t(Sp_X)
      res <- lapply(colnames(Sp_X), function(i, FUN, ...) {
        FUN(as.matrix(Sp_X[,i]), ...)
      }, FUN, ...)
    }else{
      res <- lapply(colnames(Sp_X), function(i, FUN, ...) {
        FUN(as.matrix(Sp_X[,i]), ...)
      }, FUN, ...)
    }
  }else{
    if (MARGIN == 1){
      Sp_X <- Matrix::t(Sp_X)
      res <- lapply(colnames(Sp_X), function(i, FUN, ...) {
        FUN(Sp_X[,i], ...)
      }, FUN, ...)
    }else{
      res <- lapply(colnames(Sp_X), function(i, FUN, ...) {
        FUN(Sp_X[,i], ...)
      }, FUN, ...)
    }
  }

  return(res)

}

#' @noRd
split_rows <- function (x, ncl) {
  lapply(parallel::splitIndices(nrow(x), ncl),
         function(i) x[i, , drop = FALSE])
}

#' @noRd
split_cols <- function (x, ncl) {
  lapply(parallel::splitIndices(ncol(x), ncl),
         function(i) x[, i, drop = FALSE])
}

#' @noRd
sparse_par_r_apply <- function (cl, x, FUN, convert_to_dense, ...) {
  par_res <- do.call(c, parallel::clusterApply(cl = cl,
                                                   x = split_rows(x,
                                                                  length(cl)),
                                     fun = sparse_apply, MARGIN = 1L,
                                     FUN = FUN,
                                     convert_to_dense=convert_to_dense, ...),
                     quote = TRUE)
  names(par_res) <- row.names(x)
  par_res
}

#' @noRd
sparse_par_c_apply <- function (cl = NULL, x, FUN, convert_to_dense, ...) {
  par_res <- do.call(c, parallel::clusterApply(cl = cl,
                                                   x = split_cols(x,
                                                                  length(cl)),
                                     fun = sparse_apply, MARGIN = 2L,
                                     FUN = FUN,
                                     convert_to_dense=convert_to_dense, ...),
                     quote = TRUE)
  names(par_res) <- colnames(x)
  par_res
}


#' Multicore apply-like function for cell_data_set
#'
#' mc_es_apply computes the row-wise or column-wise results of FUN, just like
#' esApply. Variables in colData from cds are available in FUN.
#'
#' @param cds A cell_data_set object.
#' @param MARGIN The margin to apply to, either 1 for rows (samples) or 2 for
#'   columns (features).
#' @param FUN Any function.
#' @param required_packages A list of packages FUN will need. Failing to
#'   provide packages needed by FUN will generate errors in worker threads.
#' @param convert_to_dense Whether to force conversion of a sparse matrix to a
#'   dense one before calling FUN.
#' @param reduction_method character, the method used to reduce dimension.
#'   Default "UMAP".
#' @param ... Additional parameters for FUN.
#' @param cores The number of cores to use for evaluation.
#'
#' @importFrom Biobase multiassign
#' @return The result of with(colData(cds) apply(counts(cds)), MARGIN, FUN, ...))
mc_es_apply <- function(cds, MARGIN, FUN, required_packages, cores=1,
                        convert_to_dense=TRUE,
                        reduction_method="UMAP", ...) {
  parent <- environment(FUN)
  if (is.null(parent))
    parent <- emptyenv()
  e1 <- new.env(parent=parent)
  coldata_df = as.data.frame(colData(cds))
  tryCatch({
    coldata_df$cluster = clusters(cds, reduction_method)[colnames(cds)]
    coldata_df$partition = partitions(cds, reduction_method)[colnames(cds)]
  }, error = function(e) {} )

  tryCatch({
    coldata_df$pseudotime = pseudotime(cds)
  }, error = function(e) {} )

  Biobase::multiassign(names(as.data.frame(coldata_df)),
                       as.data.frame(coldata_df), envir=e1)
  environment(FUN) <- e1


  platform <- Sys.info()[['sysname']]

  # Temporarily disable OpenMP threading in functions to be run in parallel
  old_omp_num_threads = as.numeric(Sys.getenv("OMP_NUM_THREADS"))
  if (is.na(old_omp_num_threads)){
    old_omp_num_threads = 1
  }
  RhpcBLASctl::omp_set_num_threads(1)

  # Temporarily set the number of threads the BLAS library can use to be 1
  old_blas_num_threads = as.numeric(Sys.getenv("OPENBLAS_NUM_THREADS"))
  if (is.na(old_omp_num_threads)){
    old_blas_num_threads = 1
  }
  RhpcBLASctl::blas_set_num_threads(1)

  # Note: use outfile argument to makeCluster for debugging
  if (platform == "Windows")
    cl <- parallel::makeCluster(cores)
  if (platform %in% c("Linux", "Darwin"))
    cl <- parallel::makeCluster(cores, type="FORK")

  cleanup <- function(){
    parallel::stopCluster(cl)
    RhpcBLASctl::omp_set_num_threads(old_omp_num_threads)
    RhpcBLASctl::blas_set_num_threads(old_blas_num_threads)
  }
  on.exit(cleanup)

  if (is.null(required_packages) == FALSE){
    parallel::clusterCall(cl, function(pkgs) {
      options(conflicts.policy =
                list(error = FALSE,
                     warn = FALSE,
                     generics.ok = TRUE,
                     can.mask = c("base", "methods", "utils",
                                  "grDevices", "graphics",
                                  "stats"),
                     depends.ok = TRUE))
      for (req in pkgs) {
        suppressMessages(library(req, character.only=TRUE, warn.conflicts=FALSE, quietly=TRUE,
                verbose=FALSE))
      }
    }, required_packages)
  }

  if (MARGIN == 1){
    suppressWarnings(res <- sparse_par_r_apply(cl, SingleCellExperiment::counts(cds), FUN,
                                               convert_to_dense, ...))
  }else{
    suppressWarnings(res <- sparse_par_c_apply(cl, SingleCellExperiment::counts(cds), FUN,
                                               convert_to_dense, ...))
  }

  res
}

#' @importFrom Biobase multiassign
smart_es_apply <- function(cds, MARGIN, FUN, convert_to_dense,
                           reduction_method="UMAP", ...) {
  parent <- environment(FUN)
  if (is.null(parent))
    parent <- emptyenv()
  e1 <- new.env(parent=parent)
  coldata_df = as.data.frame(colData(cds))
  tryCatch({
    coldata_df$cluster = clusters(cds, reduction_method)[colnames(cds)]
    coldata_df$partition = partitions(cds, reduction_method)[colnames(cds)]
    coldata_df$pseudotime = pseudotime(cds)
  }, error = function(e) {} )
  Biobase::multiassign(names(as.data.frame(coldata_df)),
                       as.data.frame(coldata_df), envir=e1)
  environment(FUN) <- e1

  if (is_sparse_matrix(SingleCellExperiment::counts(cds))){
    res <- sparse_apply(SingleCellExperiment::counts(cds), MARGIN, FUN, convert_to_dense, ...)
  } else {
    res <- pbapply::pbapply(SingleCellExperiment::counts(cds), MARGIN, FUN, ...)
  }

  if (MARGIN == 1)
  {
    names(res) <- row.names(cds)
  }else{
    names(res) <- colnames(cds)
  }

  res
}


#' Principal Components Analysis
#'
#' Efficient computation of a truncated principal components analysis of a
#' given data matrix using an implicitly restarted Lanczos method from the
#' \code{\link{irlba}} package.
#'
#' @param x a numeric or complex matrix (or data frame) which provides the data
#'   for the principal components analysis.
#' @param retx a logical value indicating whether the rotated variables should
#'   be returned.
#' @param center a logical value indicating whether the variables should be
#'   shifted to be zero centered. Alternately, a centering vector of length
#'   equal the number of columns of \code{x} can be supplied.
#' @param scale. a logical value indicating whether the variables should be
#'   scaled to have unit variance before the analysis takes place. The default
#'   is \code{FALSE} for consistency with S, but scaling is often advisable.
#'   Alternatively, a vector of length equal the number of columns of \code{x}
#'   can be supplied.
#'
#'   The value of \code{scale} determines how column scaling is performed
#'   (after centering). If \code{scale} is a numeric vector with length equal
#'   to the number of columns of \code{x}, then each column of \code{x} is
#'   divided by the corresponding value from \code{scale}. If \code{scale} is
#'   \code{TRUE} then scaling is done by dividing the (centered) columns of
#'   \code{x} by their standard deviations if \code{center=TRUE}, and the root
#'   mean square otherwise.  If \code{scale} is \code{FALSE}, no scaling is done.
#'   See \code{\link{scale}} for more details.
#' @param n integer number of principal component vectors to return, must be
#'   less than \code{min(dim(x))}.
#' @param ... additional arguments passed to \code{\link{irlba}}.
#'
#' @return
#' A list with class "prcomp" containing the following components:
#' \itemize{
#'    \item{sdev} {the standard deviations of the principal components (i.e.,
#'      the square roots of the eigenvalues of the covariance/correlation
#'      matrix, though the calculation is actually done with the singular
#'      values of the data matrix).}
#'   \item{rotation} {the matrix of variable loadings (i.e., a matrix whose
#'     columns contain the eigenvectors).}
#'   \item {x} {if \code{retx} is \code{TRUE} the value of the rotated data
#'     (the centered (and scaled if requested) data multiplied by the
#'     \code{rotation} matrix) is returned. Hence, \code{cov(x)} is the
#'     diagonal matrix \code{diag(sdev^2)}.}
#'   \item{center, scale} {the centering and scaling used, or \code{FALSE}.}
#' }
#'
#' @note
#' The signs of the columns of the rotation matrix are arbitrary, and so may
#' differ between different programs for PCA, and even between different builds
#' of R.
#'
#' NOTE DIFFERENCES WITH THE DEFAULT \code{\link{prcomp}} FUNCTION! The
#' \code{tol} truncation argument found in \code{prcomp} is not supported. In
#' place of the truncation tolerance in the original function, the
#' \code{prcomp_irlba}  function has the argument \code{n} explicitly giving
#' the number of principal components to return. A warning is generated if the
#' argument \code{tol} is used, which is interpreted differently between the
#' two functions.
#'
#' @examples
#' \dontrun{
#'   set.seed(1)
#'   x  <- matrix(rnorm(200), nrow=20)
#'   p1 <- irlba::prcomp_irlba(x, n=3)
#'   summary(p1)
#'
#'   # Compare with
#'   p2 <- prcomp(x, tol=0.7)
#'   summary(p2)}
#'
#' @seealso \code{\link{prcomp}}
sparse_prcomp_irlba <- function(x, n = 3, retx = TRUE, center = TRUE,
                                scale. = FALSE, ...)
{
  a <- names(as.list(match.call()))
  ans <- list(scale=scale.)
  if ("tol" %in% a)
    warning("The `tol` truncation argument from `prcomp` is not supported by
            `prcomp_irlba`. If specified, `tol` is passed to the `irlba`
            function to control that algorithm's convergence tolerance. See
            `?prcomp_irlba` for help.")
  orig_x <- x
  if (!methods::is(x, "DelayedMatrix")) {
    x = DelayedArray::DelayedArray(x)
  }

  args <- list(A=orig_x, nv=n)
  if (is.logical(center))
  {
    if (center) args$center <- DelayedMatrixStats::colMeans2(x)
  } else args$center <- center
  if (is.logical(scale.))
  {
    if (is.numeric(args$center))
    {
      scale. <- sqrt(DelayedMatrixStats::colVars(x))
      if (ans$scale) ans$totalvar <- ncol(x)
      else ans$totalvar <- sum(scale. ^ 2)
    } else
    {
      if (ans$scale)
      {
        scale. <-
          sqrt(DelayedMatrixStats::colSums2(x ^ 2) / (max(1, nrow(x) - 1L)))
        ans$totalvar <-
          sum(sqrt(DelayedMatrixStats::colSums2(t(t(x)/scale.) ^ 2) /
                     (nrow(x) - 1L)))
      } else
      {
        ans$totalvar <-
          sum(DelayedMatrixStats::colSums2(x ^ 2) / (nrow(x) - 1L))
      }
    }
    if (ans$scale) args$scale <- scale.
  } else
  {
    args$scale <- scale.
    ans$totalvar <-
      sum(sqrt(DelayedMatrixStats::colSums2(t(t(x)/scale.) ^ 2) /
                 (nrow(x) - 1L)))
  }
  if (!missing(...)) args <- c(args, list(...))

  s <- do.call(irlba::irlba, args=args)
  ans$sdev <- s$d / sqrt(max(1, nrow(x) - 1))
  ans$rotation <- s$v
  colnames(ans$rotation) <- paste("PC", seq(1, ncol(ans$rotation)), sep="")
  ans$center <- args$center
  ans$svd_scale <- args$scale
  if (retx)
  {
    ans <- c(ans, list(x = sweep(s$u, 2, s$d, FUN=`*`)))
    colnames(ans$x) <- paste("PC", seq(1, ncol(ans$rotation)), sep="")
  }
  class(ans) <- c("irlba_prcomp", "prcomp")
  ans
}


#' Detects genes above minimum threshold.
#'
#' @description For each gene in a cell_data_set object, detect_genes counts
#' how many cells are expressed above a minimum threshold. In addition, for
#' each cell, detect_genes counts the number of genes above this threshold that
#' are detectable. Results are added as columns num_cells_expressed and
#' num_genes_expressed in the rowData and colData tables respectively.
#'
#' @param cds Input cell_data_set object.
#' @param min_expr Numeric indicating expression threshold
#' @return Updated cell_data_set object
#' @export
#' @examples
#' \dontrun{
#'    cds <- detect_genes(cds, min_expr=0.1)
#' }
detect_genes <- function(cds, min_expr=0){
  assertthat::assert_that(methods::is(cds, "cell_data_set"))
    assertthat::assert_that(is.numeric(min_expr))

  rowData(cds)$num_cells_expressed <- Matrix::rowSums(SingleCellExperiment::counts(cds) > min_expr)
  colData(cds)$num_genes_expressed <- Matrix::colSums(SingleCellExperiment::counts(cds) > min_expr)

  cds
}

#' Return a size-factor normalized and (optionally) log-transformed expression
#' matrix
#'
#' @param cds A CDS object to calculate normalized expression matrix from.
#' @param norm_method String indicating the normalization method. Options are
#'   "log" (Default), "binary" and "size_only".
#' @param pseudocount A pseudocount to add before log transformation. Ignored
#'   if norm_method is not "log". Default is 1.
#' @return Size-factor normalized, and optionally log-transformed, expression
#'   matrix.
#'
#' @examples
#'   \donttest{
#'     cds <- load_a549()
#'     normalized_matrix <- normalized_counts(cds)
#'   }
#'
#' @export
normalized_counts <- function(cds,
                              norm_method=c("log", "binary", "size_only"),
                              pseudocount=1){
  norm_method = match.arg(norm_method)
  norm_mat = SingleCellExperiment::counts(cds)
  if (norm_method == "binary"){
    # The '+ 0' coerces the matrix to type numeric. It's possible
    # to use 'as.numeric(norm_mat > 0)' but the matrix
    # attributes disappear...            
    norm_mat = (norm_mat > 0) + 0
    if (is_sparse_matrix(norm_mat)){
      norm_mat = methods::as(norm_mat, "dgCMatrix")
    }
  }
  else {
    if (is_sparse_matrix(norm_mat)){
      norm_mat@x = norm_mat@x / rep.int(size_factors(cds), diff(norm_mat@p))
      if (norm_method == "log"){
        if (pseudocount == 1){
          norm_mat@x = log10(norm_mat@x + pseudocount)
        }else{
          stop("Pseudocount must equal 1 with sparse expression matrices")
        }
      }
    }else{
      norm_mat = Matrix::t(Matrix::t(norm_mat) / size_factors(cds))
      if (norm_method == "log"){
          norm_mat@x <- log10(norm_mat + pseudocount)
      }
    }
  }
  return(norm_mat)
}



#' Combine a list of cell_data_set objects
#'
#' This function will combine a list of cell_data_set objects into a new
#' cell_data_set object.
#'
#' @param cds_list List of cds objects to be combined.
#' @param keep_all_genes Logical indicating what to do if there is a mismatch
#'   in the gene sets of the CDSs. If TRUE, all genes are kept and cells from
#'   CDSs missing a given gene will be filled in with zeroes. If FALSE, only
#'   the genes in common among all of the CDSs will be kept. Default is TRUE.
#' @param cell_names_unique Logical indicating whether all of the cell IDs
#'   across all of the CDSs are unique. If FALSE, the CDS name is appended to
#'   each cell ID to prevent collisions. These cell IDs are used as count matrix
#'   column names and colData(cds) row names. Cell names stored in other
#'   cds locations are not modified so you will need to modify them manually
#'   for consistency. Default is FALSE.
#' @param sample_col_name A string to be the column name for the colData column
#'   that indicates which original cds the cell derives from. Default is
#'   "sample".
#' @param keep_reduced_dims Logical indicating whether to keep the reduced
#'   dimension matrices. Default is FALSE.
#'
#' @return A combined cell_data_set object.
#' @export
#'
combine_cds <- function(cds_list,
                        keep_all_genes = TRUE,
                        cell_names_unique = FALSE,
                        sample_col_name = "sample",
                        keep_reduced_dims = FALSE) {

  assertthat::assert_that(is.list(cds_list),
                          msg=paste("cds_list must be a list."))

  assertthat::assert_that(all(sapply(cds_list, class) == "cell_data_set"),
                          msg=paste("All members of cds_list must be",
                                    "cell_data_set class."))
  assertthat::assert_that(is.character(sample_col_name))

  if (sample_col_name == "sample" &
      any(sapply(cds_list, function(cds) "sample" %in% names(colData(cds))))) {
    warning("By default, the combine_cds function adds a column called ",
                   "'sample' which indicates which initial cds a cell came ",
                   "from. One or more of your input cds objects contains a ",
                   "'sample' column, which will be overwritten. We recommend ",
                   "you rename this column or provide an alternative column ",
                   "name using the 'sample_col_name' parameter.")
  }
  assertthat::assert_that(!any(sapply(cds_list, function(cds)
    sum(is.na(names(colData(cds)))) != 0)),
                          msg = paste0("One of the input CDS' has a colData ",
                                       "column name that is NA, please ",
                                       "remove or rename that column before ",
                                       "proceeding."))
  assertthat::assert_that(!any(sapply(cds_list, function(cds)
    sum(is.na(names(rowData(cds)))) != 0)),
    msg = paste0("One of the input CDS' has a rowData ",
                 "column name that is NA, please ",
                 "remove or rename that column before ",
                 "proceeding."))
  num_cells <- sapply(cds_list, ncol)
  if(sum(num_cells == 0) != 0) {
    message("Some CDS' have no cells, these will be skipped.")
    cds_list <- cds_list[num_cells != 0]
  }
  if(length(cds_list) == 1) return(cds_list[[1]])

  assertthat::assert_that(is.logical(keep_all_genes))
  assertthat::assert_that(is.logical(cell_names_unique))

  list_named <- TRUE
  if(is.null(names(cds_list))) {
    list_named <- FALSE
  }

  exprs_list <- list()
  fd_list <- list()
  pd_list <- list()
  gene_list <- c()
  overlap_list <- c(row.names(fData(cds_list[[1]])))
  pdata_cols <- c()
  fdata_cols <- c()
  all_cells <- c()

  for(cds in cds_list) {
    gene_list <-  c(gene_list, row.names(fData(cds)))
    overlap_list <- intersect(overlap_list, row.names(fData(cds)))
    if (!keep_all_genes) {
      gene_list <- overlap_list
    }

    pdata_cols <- c(pdata_cols, names(pData(cds)))
    fdata_cols <- c(fdata_cols, names(fData(cds)))
    all_cells <- c(all_cells, row.names(pData(cds)))
  }

  gene_list <- unique(gene_list)
  if(length(overlap_list) == 0) {
    if (keep_all_genes) {
      warning("No genes are shared amongst all the CDS objects.")
    } else {
      stop("No genes are shared amongst all the CDS objects. To generate ",
                 "a combined CDS with all genes, use keep_all_genes = TRUE")
    }
  }
  pdata_cols <- unique(pdata_cols)
  fdata_cols <- unique(fdata_cols)
  if (sum(duplicated(all_cells)) != 0 & cell_names_unique) {
    stop("Cell names are not unique across CDSs - cell_names_unique ",
               "must be FALSE.")
  }
  all_cells <- unique(all_cells)
  for(i in seq(1, length(cds_list), 1)) {
    pd <- as.data.frame(pData(cds_list[[i]]))
    exp <- exprs(cds_list[[i]])
    exp <- exp[intersect(row.names(exp), gene_list),, drop=FALSE]
    if (!cell_names_unique) {
      if(list_named) {
        row.names(pd) <- paste(row.names(pd), names(cds_list)[[i]], sep="_")
        pd[,sample_col_name] <- names(cds_list)[[i]]
      } else {
        row.names(pd) <- paste(row.names(pd), i, sep="_")
        pd[,sample_col_name] <- i
      }
      colnames(exp) <- row.names(pd)
    } else {
      if(list_named) {
        pd[,sample_col_name] <- names(cds_list)[[i]]
      } else {
        pd[,sample_col_name] <- i
      }
    }
    not_in <- pdata_cols[!pdata_cols %in% names(pd)]
    for (n in not_in) {
      pd[,n] <- NA
    }

    fd <- as.data.frame(fData(cds_list[[i]]))
    fd <- fd[intersect(row.names(fd), gene_list),, drop=FALSE]
    not_in <- fdata_cols[!fdata_cols %in% names(fd)]
    for(col in names(fd)) {
      if(methods::is(fd[,col], "factor")) {
        fd[,col] <- as.character(fd[,col])
      }
    }
    for (n in not_in) {
      fd[,n] <- NA
    }
    not_in_g <- gene_list[!gene_list %in% row.names(fd)]

    if (length(not_in_g) > 0) {
      not_in_g_df <- as.data.frame(matrix(NA, nrow = length(not_in_g), ncol=ncol(fd)))
      row.names(not_in_g_df) <- not_in_g
      names(not_in_g_df) <- names(fd)
      fd <- rbind(fd, not_in_g_df)

      extra_rows <- Matrix::Matrix(0, ncol=ncol(exp),
                                   sparse=TRUE,
                                   nrow=length(not_in_g))
      row.names(extra_rows) <- not_in_g
      colnames(extra_rows) <- colnames(exp)
      exp <- rbind(exp, extra_rows)
      exp <- exp
    }


    exprs_list[[i]] <- exp[gene_list, , drop=FALSE]
    fd_list[[i]] <- fd[gene_list, , drop=FALSE]
    pd_list[[i]] <- pd

  }
  all_fd <- array(NA,dim(fd_list[[1]]),dimnames(fd_list[[1]]))

  for (fd in fd_list) {
    for (j in colnames(fd)) {
      col_info <- all_fd[,j]
      col_info[is.na(col_info)] <- fd[is.na(col_info),j]
      col_info[col_info != fd[,j]] <- "conf"
      all_fd[,j] <- col_info
    }
  }

  confs <- sum(all_fd == "conf", na.rm=TRUE)

  if (confs > 0) {
   warning("When combining rowData, conflicting values were found - ",
                  "conflicts will be labelled 'conf' in the combined cds ",
                  "to prevent conflicts, either change conflicting values to ",
                  "match, or rename columns from different cds' to be unique.")
  }
  #all_fd <- do.call(cbind, fd_list)
  all_fd <- all_fd[,fdata_cols, drop=FALSE]

  all_pd <- do.call(rbind, pd_list)
  all_exp <- do.call(cbind, exprs_list)

  all_exp <- all_exp[row.names(all_fd), row.names(all_pd), drop=FALSE]

  new_cds <- new_cell_data_set(all_exp, cell_metadata = all_pd, gene_metadata = all_fd)

  if(keep_reduced_dims) {
    for(red_dim in names(SingleCellExperiment::reducedDims(cds_list[[1]]))) {
      reduced_dims_list <- list()
      for(j in seq(1, length(cds_list), 1)) {
        reduced_dims_list[[j]] <- SingleCellExperiment::reducedDims(cds_list[[j]])[[red_dim]]
      }
      SingleCellExperiment::reducedDims(new_cds)[[red_dim]] <- do.call(rbind, reduced_dims_list, quote=FALSE)
      # The following should not happen; the accessor appears to ensure the
      # correct row order.
      if(!identical(rownames(SingleCellExperiment::reducedDims(new_cds)[[red_dim]]), rownames(all_pd))) {
        stop('Mis-ordered reduced matrix rows.')
      }
    }
  }
  new_cds
}

#' Clear CDS slots
#'
#' Function to clear all CDS slots besides colData, rowData and expression data.
#'
#' @param cds cell_data_set to be cleared
#'
#' @return A cell_data_set with only expression, rowData and colData present.
#' @export
#'
clear_cds_slots <- function(cds) {
  cds@reduce_dim_aux <- S4Vectors::SimpleList()
  cds@principal_graph_aux <- S4Vectors::SimpleList()
  cds@principal_graph <- S4Vectors::SimpleList()
  cds@clusters <- S4Vectors::SimpleList()
  SingleCellExperiment::reducedDims(cds) <- S4Vectors::SimpleList()
  cds
}


add_citation <- function(cds, citation_key) {
  citation_map <- list(
    UMAP = c("UMAP", "McInnes, L., Healy, J. & Melville, J. UMAP: Uniform Manifold Approximation and Projection for dimension reduction. Preprint at https://arxiv.org/abs/1802.03426 (2018)."),
    MNN_correct = c("MNN Correct", "Haghverdi, L. et. al. Batch effects in single-cell RNA-sequencing data are corrected by matching mutual nearest neighbors. Nat. Biotechnol. 36, 421-427 (2018). https://doi.org/10.1038/nbt.4091"),
    partitions = c("paritioning", c("Levine, J. H., et. al. Data-driven phenotypic dissection of AML reveals progenitor-like cells that correlate with prognosis. Cell 162, 184-197 (2015). https://doi.org/10.1016/j.cell.2015.05.047",
                                  "Wolf, F. A. et. al. PAGA: graph abstraction reconciles clustering with trajectory inference through a topology preserving map of single cells. Genome Biol. 20, 59 (2019). https://doi.org/10.1186/s13059-019-1663-x")),
    clusters = c("clustering", "Levine, J. H. et. al. Data-driven phenotypic dissection of AML reveals progenitor-like cells that correlate with prognosis. Cell 162, 184-197 (2015). https://doi.org/10.1016/j.cell.2015.05.047"),
    leiden = c("leiden", "Traag, V.A., Waltman, L. & van Eck, N.J. From Louvain to Leiden: guaranteeing well-connected communities. Scientific Reportsvolume 9, Article number: 5233 (2019). https://doi.org/10.1038/s41598-019-41695-z" )
  )
  if (is.null(S4Vectors::metadata(cds)$citations) | citation_key == "Monocle") {
    S4Vectors::metadata(cds)$citations <- data.frame(method = c("Monocle", "Monocle", "Monocle"),
                                          citations = c("Trapnell C. et. al. The dynamics and regulators of cell fate decisions are revealed by pseudotemporal ordering of single cells. Nat. Biotechnol. 32, 381-386 (2014). https://doi.org/10.1038/nbt.2859",
                                                        "Qiu, X. et. al. Reversed graph embedding resolves complex single-cell trajectories. Nat. Methods 14, 979-982 (2017). https://doi.org/10.1038/nmeth.4402",
                                                        "Cao, J. et. al. The single-cell transcriptional landscape of mammalian organogenesis. Nature 566, 496-502 (2019). https://doi.org/10.1038/s41586-019-0969-x"))
  }
  S4Vectors::metadata(cds)$citations <- rbind(S4Vectors::metadata(cds)$citations,
                                   data.frame(method = citation_map[[citation_key]][1],
                                              citations = citation_map[[citation_key]][2]))
  cds
}

#' Access citations for methods used during analysis.
#'
#' @param cds The cds object to access citations from.
#'
#' @return A data frame with the methods used and the papers to be cited.
#' @export
#'
#' @examples {
#'   \dontrun{
#'      get_citations(cds)
#'   }
#' }
get_citations <- function(cds) {
  message("Your analysis used methods from the following recent work. ",
                "Please cite them wherever you are presenting your analyses.")
  if(is.null(S4Vectors::metadata(cds)$citations)) {
    cds <- add_citation(cds, "Monocle")
  }
  S4Vectors::metadata(cds)$citations
}


# Make a unique identifier string.
get_unique_id <- function(object=NULL) {

  if(!is.null(object)) {
    object_dim <- dim(object)
    object_checksum <- digest::digest(object)
    if(!is.null(object_dim))
      object_id <- list(checksum=object_checksum, dim=object_dim)
    else
      object_id <- list(checksum=object_checksum, dim=length(object))
  }
  else {
    id_count <- get_global_variable('id_count')
    rtime <- as.numeric(Sys.time()) * 100000 + id_count
    object_id <- openssl::md5(as.character(rtime))
    id_count <- id_count + 1
    set_global_variable('id_count', id_count)
  }

  return(object_id)
}


# Get current time stamp.
get_time_stamp <- function() {
  time_stamp <- format(Sys.time(), "%Y%m%d:%H%M%S" )
  return(time_stamp)
}


# Manage parallel processing for matrix multiplication.
matrix_multiply_multicore <- function(mat_a, mat_b, cores=1L) {
  assertthat::assert_that(is.matrix(mat_a) || is_sparse_matrix(mat_a) || methods::is(mat_a, 'DelayedMatrix'),
    msg=paste0('mat_a must be either a matrix or a sparse matrix'))
  assertthat::assert_that(is.matrix(mat_b) || is_sparse_matrix(mat_b) || methods::is(mat_b, 'DelayedMatrix'),
    msg=paste0('mat_b must be either a matrix or a sparse matrix'))

  if(cores > 1) {
    omp_num_threads <- get_global_variable('omp_num_threads')
    blas_num_threads <- get_global_variable('blas_num_threads')
  
    RhpcBLASctl::omp_set_num_threads(1L)
    RhpcBLASctl::blas_set_num_threads(1L)
  
    DelayedArray::setAutoBPPARAM(BPPARAM=BiocParallel::MulticoreParam(workers=as.integer(cores)))
  
    mat_c <- mat_a %*% mat_b
  
    DelayedArray::setAutoBPPARAM(BPPARAM=BiocParallel::SerialParam())
  
    RhpcBLASctl::omp_set_num_threads(as.integer(omp_num_threads))
    RhpcBLASctl::blas_set_num_threads(as.integer(blas_num_threads))
  } else {
    mat_c <- mat_a %*% mat_b
  }

  return(mat_c)
}


# Return the call stack as a character vector.
get_call_stack <- function () {
  cv<-as.vector(sys.calls())
  lcv <- length(cv)
  n <- lcv - 1

  ocv <- vector()
  for(i in seq(1,n,1)) {
    elem <- stringr::str_split(as.character(cv[i]), '[(]', n=2)[[1]][[1]]
    ocv <- c(ocv, elem)
  }

  return(ocv)
}


# Return the call stack as a character string.
get_call_stack_as_string <- function() {
  cs <- get_call_stack()
  scs <- ''
  for(i in seq(length(cs)-1)) {
    csep <- ifelse(i == 1, '', ' => ')
    scs <- sprintf("%s%s%s()", scs, csep, cs[[i]])
  }

  return(scs)
}


# Return the name of object as a string
object_name_to_string <- function( object ) {
  str <- deparse(substitute(object))

  return( str )
}


stop_no_noise <- function() {
  opt <- options(show.error.messages = FALSE)
  on.exit(options(opt))
  stop()
}


# number of tasks per block for multicore processing
# task vector limits
#   vbeg <- c(0,cumsum(tasks_per_block(11,3)))[1:3]+1
#   vend <- cumsum(tasks_per_block(11,3))
tasks_per_block <- function(ntask=NULL, nblock=NULL) {
  tasks_block <- rep(trunc(ntask/nblock), nblock)
  remain <- ntask %% nblock
  if(remain)
    for(i in seq(remain)) 
      tasks_block[i] <- tasks_block[i] + 1
  return(tasks_block)
}


#
# Initialize Monocle3 timer.
#
tick <- function(msg="") {
  set_global_variable('monocle3_timer_t0', Sys.time())
  set_global_variable('monocle3_timer_msg', msg)
}

#
# Return time elapsed since call to tick.
#
tock <- function() {
  t1 <- Sys.time()
  t0 <- get_global_variable('monocle3_timer_t0')
  msg <- get_global_variable('monocle3_timer_msg')
  if(length(msg) > 0) {
    message(sprintf('%s %.2f seconds.',msg, t1 - t0))
  }
  else {
    return(t1 - t0)
  }
}

#
# Set up a global-variable-like environment.
#

set_global_variable <- function(variable_name, value) {
  assign(variable_name, value, envir=._._global_variable_env_._.)
}


# Return value of variable_name. If variable_name is NULL, return a list
# of all global variables.
get_global_variable <- function(variable_name=NULL) {
  value <- tryCatch({
                      v <- get('guard_element', envir=._._global_variable_env_._.) 
                      if(v != 'sanity_check') stop()
                      v
                    }, error=function(msg) {
                      message('Global variable storage is compromised.')
                      return(NA)
                    })

  if(is.na(value)) {
    return(NA)
  }

  if(!is.null(variable_name)) {
    value <- tryCatch({
                        get(variable_name, envir=._._global_variable_env_._.)
                      }, error=function(msg) {
                        message("\'", variable_name, "\'", ' is not a global variable.')
                        return(NA)
                      })
  }
  else {
    value=list()
    variable_names <- ls(envir=._._global_variable_env_._.)
    for(variable_name in variable_names) {
      value[[variable_name]] <- get(variable_name, envir=._._global_variable_env_._.)
    }
  }

  return(value)
}


# The ._._global_variable_env_._. environment stores global
# objects. Use the set_global_variable and get_global_variable
# functions to access them.
._._global_variable_env_._. <- new.env(parent=emptyenv())

# Define some global variables.
.onLoad <- function(libname, pkgname) {
  # A value used to ensure that this is the Monocle3
  # global variables.
  set_global_variable('guard_element', 'sanity_check')

  # Counter used to ensure that matrix time stamps generated
  # by the get_unique_id() function are distinct. This value
  # is incremented each time that get_unique_id() is called
  # during an R session.
  set_global_variable('id_count', 1)

  # Object version numbers. These are used by the Monocle3
  # to recognize objects so we suggest that you not change
  # them.
  set_global_variable('reduce_dim_pca_model_version', 1)
  set_global_variable('reduce_dim_lsi_model_version', 1)
  set_global_variable('reduce_dim_aligned_model_version', 1)
  set_global_variable('reduce_dim_tsne_model_version', 1)
  set_global_variable('reduce_dim_umap_model_version', 1)
  set_global_variable('monocle_objects_version', 1)
  set_global_variable('transform_models_version', 1)
  set_global_variable('monocle3_annoy_index_version', 2)
  set_global_variable('monocle3_hnsw_index_version', 1)
  set_global_variable('monocle3_timer_t0', 0)
  set_global_variable('monocle3_timer_msg', "")

  # Default nn_control list for functions that do not need
  # an index, which is all but the label transfer functions.
  set_global_variable('nn_control_annoy_euclidean', list(method='annoy', metric='euclidean', n_trees=50, M=48, ef_construction=200, ef=150, grain_size=1, cores=1))

  # Default nn_control list for functions that need an index,
  # which are the label transfer functions.
  set_global_variable('nn_control_annoy_cosine', list(method='annoy', metric='cosine', n_trees=50, M=48, ef_construction=200, ef=150, grain_size=1, cores=1))

  # Watching preprocess_cds() it appears that R uses OMP_NUM_THREADS
  # threads if OMP_NUM_THREADS > 1 and OPENBLAS_NUM_THREADS is NA.
  # After setting, RhpcBLASctl::blas_set_num_threads(1L) and then
  # RhpcBLASctl::blas_set_num_threads(NA), preprocess_cds() uses
  # one thread. So I assume that the number of desired threads for
  # both OMP_NUM_THREADS and OPENBLAS_NUM_THREADS is the greater
  # of the two.
  omp_num_threads <- as.numeric(Sys.getenv('OMP_NUM_THREADS'))
  blas_num_threads <- as.numeric(Sys.getenv('OPENBLAS_NUM_THREADS'))

  if(is.na(omp_num_threads))
    omp_num_threads <- 1L
  if(is.na(blas_num_threads))
    blas_num_threads <- 1L

  if(omp_num_threads < blas_num_threads)
    omp_num_threads <- blas_num_threads
  else
    blas_num_threads <- omp_num_threads
  
  set_global_variable('omp_num_threads', omp_num_threads)
  set_global_variable('blas_num_threads', blas_num_threads)

  # for travis
  Sys.setenv('TESTTHAT_MAX_FAILS' = Inf)
}

# Functions that support nearest neighbors use.


# Check whether nn index exists and is consistent with matrix and parameters.
# This function is not in use currently and may fall into disrepair.
check_cds_nn_index_is_current <- function(cds, reduction_method=c('PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP'), nn_control=list(), verbose=FALSE) {

  assertthat::assert_that(methods::is(cds, 'cell_data_set'),
                          msg=paste('cds parameter is not a cell_data_set'))

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(reduction_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "reduction_method must be one of 'PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP'")
  reduction_method <- match.arg(reduction_method)

  assertthat::assert_that(
    !is.null(SingleCellExperiment::reducedDims(cds)[[reduction_method]]),
    msg = paste0('Data has not been processed with',
                 ' chosen reduction_method: ',
                 reduction_method))

  # We check the index build parameters.
  nn_control <- set_nn_control(mode=1, nn_control=nn_control, nn_control_default=list(), nn_index=NULL, k=NULL, verbose=verbose)

  if(verbose) {
    message('check_cds_nn_index_is_current:')
    report_nn_control('  nn_control: ', nn_control=nn_control)
  }

  nn_method <- nn_control[['method']]

  if(is.null(cds@reduce_dim_aux[[reduction_method]]) ||
     is.null(cds@reduce_dim_aux[[reduction_method]][['nn_index']]) ||
     is.null(cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]])) {
    if(verbose) {
      message('check_cds_nn_index_is_current: FALSE')
    }
    return(FALSE)
  }

  if(nrow(SingleCellExperiment::reducedDims(cds)[[reduction_method]]) != cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['nn_index']][['nrow']]) {
    if(verbose) {
      message('check_cds_nn_index_is_current: FALSE')
    }
    return(FALSE)
  }

  if(ncol(SingleCellExperiment::reducedDims(cds)[[reduction_method]]) != cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['nn_index']][['ncol']]) {
    if(verbose) {
      message('check_cds_nn_index_is_current: FALSE')
    }
    return(FALSE)
  }

  if(cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['matrix_id']] != get_reduce_dim_matrix_identity(cds, reduction_method)[['matrix_id']]) {
    if(verbose) {
      message('check_cds_nn_index_is_current: FALSE')
    }
    return(FALSE)
  }

  if(nn_method == 'annoy') {
    if(cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['nn_index']][['metric']] != nn_control[['metric']]) {
      if(verbose) {
        message('check_cds_nn_index_is_current: FALSE')
      }
      return(FALSE)
    }

    if(cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['nn_index']][['n_trees']] != nn_control[['n_trees']]) {
        if(verbose) {
          message('check_cds_nn_index_is_current: FALSE')
        }
       return(FALSE)
    }
  }
  else
  if(nn_method == 'hnsw') {
    if(cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['nn_index']][['metric']] != nn_control[['metric']]) {
      if(verbose) {
        message('check_cds_nn_index_is_current: FALSE')
      }
      return(FALSE)
    }

    if(cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['nn_index']][['M']] != nn_control[['M']]) { 
      if(verbose) {
        message('check_cds_nn_index_is_current: FALSE')
      }
      return(FALSE)
    }

    if(cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['nn_index']][['ef_construction']] != nn_control[['ef_construction']]) { 
      if(verbose) {
        message('check_cds_nn_index_is_current: FALSE')
      }
      return(FALSE)
    }
  }
  else
    stop('check_cds_nn_index_is_current: unsupported nearest neighbor index type \'', nn_method, '\'')

  if(verbose) {
    message('check_cds_nn_index_is_current: TRUE')
  }

  return(TRUE)
}


# Clear nearest neighbor index.
# Note: remove cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]]
#       because the objects stored in it are related; that is, if the object
#       cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['nn_index']]
#       is removed or changed, the other objects in
#       cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]] need to be
#       updated.
clear_cds_nn_index <- function(cds, reduction_method=c('PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP'), nn_method=c('annoy', 'hnsw', 'all')) {
  assertthat::assert_that(methods::is(cds, 'cell_data_set'),
                          msg=paste('cds parameter is not a cell_data_set'))

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(reduction_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "reduction_method must be one of 'PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP'")
  reduction_method <- match.arg(reduction_method)

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(nn_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "nn_method must be one of 'annoy', 'hnsw', or 'all'")

  nn_method <- match.arg(nn_method)

  if(nn_method == 'annoy') {
    if(!is.null(cds@reduce_dim_aux[[reduction_method]][['nn_index']]) &&
       !is.null(cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]])) {
      cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]] <- NULL
    }
  }
  else
  if(nn_method == 'hnsw') {
    if(!is.null(cds@reduce_dim_aux[[reduction_method]][['nn_index']]) &&
       !is.null(cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]])) {
      cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]] <- NULL
    }
  }
  else
  if(nn_method == 'all') {
    if(!is.null(cds@reduce_dim_aux[[reduction_method]][['nn_index']])) {
      cds@reduce_dim_aux[[reduction_method]][['nn_index']] <- S4Vectors::SimpleList()
    }
  }
  else
    stop('clear_cds_nn_index: unsupported nearest neighbor index type \'', nn_method, '\'')

  return(cds)
}


# Check whether nearest neighbor search information exists
# for reduction_method. This is not useful at this time. The
# intention was to store information about the search
# parameters used in cluster_cells() and use them in for
# nearest neighbor searches in later function calls. However,
# this does not seem useful at this time.
# This function is not in use currently and may fall into disrepair.
check_cds_nn_search_exists <- function(cds, reduction_method = c("UMAP", "tSNE", "PCA", "LSI", "Aligned"), search_id, verbose=FALSE)
{
  assertthat::assert_that(methods::is(cds, "cell_data_set"))
  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(reduction_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "reduction_method must be one of 'UMAP', 'tSNE', 'PCA', 'LSI', 'Aligned'")
  reduction_method <- match.arg(reduction_method)

  if(is.null(SingleCellExperiment::reducedDims(cds)[[reduction_method]])) {
    return(FALSE)
  }

  if(is.null(cds@clusters[[reduction_method]])) {
    return(FALSE)
  }

  if(is.null(cds@reduce_dim_aux[[reduction_method]][['nn_search']][[search_id]])) {
    return(FALSE)
  }

  return(TRUE)
}


select_nn_parameter_value <- function(parameter, nn_control, nn_control_default, default_value) {
  if(!is.null(nn_control[[parameter]])) {
    return(nn_control[[parameter]])
  }
  else
  if(!is.null(nn_control_default[[parameter]])) {
    return(nn_control_default[[parameter]])
  }
  return(default_value)
}


#' set annoy search_k parameter
#' precedence for setting search_k
#'   o  nn_control[['search_k']]
#'   o  mode == 2 and the following values exist/set
#'        o  nn_index parameter
#'        o  n_trees value in nn_index parameter object
#'        o  k parameter
#'   o  nn_control[['n_trees']] and k parameter or k default
#'   o  nn_control_default[['search_k']]
#'   o  nn_control_default[['n_trees']] and k parameter or k default
#'   o  n_trees default and k parameter or k default
#' @noRd
select_annoy_search_k <- function(mode, nn_control, nn_control_default, nn_index, k, default_n_trees, default_k) {
  if(!is.null(k)) {
    use_k <- k
    src_k <- 'parameter'
  } 
  else {
    use_k <- default_k
    src_k <- 'default'
  }

  # Sanity test.
  if(mode == 2 &&
     !is.null(nn_index) &&
     is.null(nn_index[['n_trees']])) {
    stop('set_nn_control: unexpected condition: found old version of reduce_dim_aux')
  }

  if(!is.null(nn_control[['search_k']])) {
    return(nn_control[['search_k']])
  }
  else
  if(mode == 2 &&
     !is.null(nn_index) &&
     !is.null(nn_index[['n_trees']]) &&
     !is.null(k)) {
    n_trees <- nn_index[['n_trees']]
    return(2 * n_trees * k)
  }
  else
  if(!is.null(nn_control[['n_trees']])) {
    return(2 * nn_control[['n_trees']] * use_k)
  }
  else
  if(!is.null(nn_control_default[['search_k']])) {
    return(nn_control_default[['search_k']])
  }
  else
  if(!is.null(nn_control_default[['n_trees']])) {
    return(2 * nn_control_default[['n_trees']] * use_k)
  }

  return(2 * default_n_trees * use_k)
}



#' @title Verify and set nearest neighbor parameter list.
#'
#' @description Verifies the listed parameter values
#'   that will be passed to the nearest neighbor function
#'   given by nn_control\[\['method'\]\]. Unspecified
#'   values are set to default values. To see the default
#'   values, call the function with
#'   nn_control=list(show_values=TRUE).
#'
#' @section Optional nn_control parameters:
#' \describe{
#'   \item{method}{The method used to find nearest neighbor points.
#'      The available methods are 'nn2', 'annoy', and 'hnsw'.
#'      Detailed information about each method can be found on the
#'      WWW sites:
#'      https://cran.r-project.org/web/packages/RANN/,
#'      https://cran.r-project.org/web/packages/RcppAnnoy/index.html,
#'      and https://cran.rstudio.com/web/packages/RcppHNSW/index.html.}
#'   \item{metric}{The distance metric used by the nearest neighbor
#'      functions.  Annoy accepts 'euclidean', 'cosine', 'manhattan',
#'      and 'hamming'. HNSW accepts 'euclidean', 'l2', 'cosine', and
#'      'ip'. RANN uses 'euclidean'.}
#'   \item{n_trees}{The annoy index build parameter that affects the build
#'      time and index size. Larger values give more accurate results,
#'      longer build times, and larger indexes.}
#'   \item{search_k}{The annoy index search parameter that affects the
#'      search accuracy and time. Larger values give more accurate
#'      results and longer search times. Default is 2 * n_trees * k.
#'      In order to set search_k, the following conditions are tested
#'      and the first TRUE condition is used: nn_control\[\['search_k'\]\]
#'      exists; mode=2 and nn_index and k are not NULL;
#'      nn_control\[\['n_trees'\]\] exists;
#'      nn_control_default\[\['search_k'\]\] exists;
#'      nn_control_default\[\['n_trees'\]\] exists. If none of those is
#'      TRUE, the fallback default n_trees value is used. If the
#'      set_nn_control k parameter value is not NULL, it is used;
#'      otherwise, the default is used.}
#'   \item{M}{The HNSW index build parameter that affects the
#'      search accuracy and memory requirements. Larger values give
#'      more accurate search results and increase the index memory
#'      use.}
#'   \item{ef_construction}{The HNSW index build parameter that affects
#'      the search accuracy and index build time. Larger values give more 
#'      accurate search results and longer build times. Default is 200.}
#'   \item{ef}{The HNSW index search parameter that affects the search
#'      accuracy and search time. Larger values give more accurate results
#'      and longer search times. ef must be greater than or equal to k.}
#'   \item{grain_size}{The HNSW parameter that gives the
#'      minimum amount of work to do per thread.}
#'   \item{cores}{The annoy and HNSW parameter that gives the number of
#'      threads to use for the annoy index search and for the HNSW index
#'      build and search.}
#'   \item{show_values}{A logical value used to show the
#'      nearest neighbor parameters to use, and then exit the function.
#'      When show_values=TRUE is the only nn_control value, the
#'      parameters are the defaults for the function. Each function
#'      that calls set_nn_control may have its own nn_control_default
#'      list.}
#' }
#'
#' @param mode the nearest neighbor operation for which the nn_control
#'   list will be used. 1=make index, 2=search index, and 3=both make
#'   and search index. Required parameter.
#' @param nn_control an optional list of parameters passed to
#'   the nearest neighbor function specified by nn_control\[\['method'\]\].
#'   If a value is not given in nn_control, the value in nn_control_default
#'   is used. If neither is given, a fallback default is assigned.
#' @param nn_control_default an optional nn_control list to use when
#'   a parameter is not given in nn_control.
#' @param nn_index an nn_index. This may be used to look up parameters
#'   that were used to make an index. For example, the default search_k
#'   parameter depends on the n_trees values used to make the index.
#'   The default is NULL.
#' @param k integer the number of desired nearest neighbor points to
#'   return from a search. k is used to
#'   * set the annoy search_k parameter when search_k is not given in
#'     nn_control or nn_control_default where search_k is 2 * n_trees * k.
#'   * test the hnsw ef parameter, which must be at least as large
#'     as k.
#'
#'   k is ignored for index builds and does not give the number
#'   of nearest neighbors to return for a search.
#' @param verbose a boolean indicating whether to emit verbose output.
#'
#' @return an updated nn_control list.
set_nn_control <- function(mode, nn_control=list(), nn_control_default=list(), nn_index=NULL, k=NULL, verbose=FALSE) {

  default_method <- 'annoy'
  default_metric <- 'euclidean'
  default_k <- 25
  default_n_trees <- 50
  default_M <- 48
  default_ef_construction <- 200
  default_ef <- 150
  default_grain_size <- 1
  default_cores <- 1
  default_annoy_random_seed <- 42

  assertthat::assert_that(methods::is(nn_control, "list"))
  assertthat::assert_that(methods::is(nn_control_default, "list"))

  allowed_control_parameters <- c('method',
                                  'metric',
                                  'n_trees',
                                  'search_k',
                                  'M',
                                  'ef_construction',
                                  'ef',
                                  'grain_size',
                                  'cores',
                                  'show_values',
                                  'annoy_random_seed')

  assertthat::assert_that(all(names(nn_control) %in% allowed_control_parameters),
                          msg = "set_nn_control: unknown variable in nn_control")
  assertthat::assert_that(all(names(nn_control_default) %in% allowed_control_parameters),
                          msg = "set_nn_control: unknown variable in nn_control_default")

  assertthat::assert_that(assertthat::is.count(mode) && mode >= 1 && mode <= 3,
                          msg = paste0("set_nn_control: invalid mode value. Mode must be an integer with the value 1, 2, or 3."))

  nn_control_out <- list()

  nn_control_out[['method']] <- select_nn_parameter_value('method', nn_control, nn_control_default, default_method)

  if(nn_control_out[['method']] == 'nn2') {
    # Do nothing for nn2. Call report_nn_control when verbose <- TRUE
  } else
  if(nn_control_out[['method']] == 'annoy') {
    nn_control_out[['metric']] <- select_nn_parameter_value('metric', nn_control, nn_control_default, default_metric)
    assertthat::assert_that(nn_control_out[['metric']] %in% c('euclidean', 'cosine', 'manhattan', 'hamming'),
                            msg=paste0("set_nn_control: nearest neighbor metric for annoy must be one of 'euclidean', 'cosine', 'manhattan', or 'hamming'"))
    if(bitwAnd(mode, 1)) {
      nn_control_out[['n_trees']] <- select_nn_parameter_value('n_trees', nn_control, nn_control_default, default_n_trees)
      nn_control_out[['annoy_random_seed']] <- select_nn_parameter_value('annoy_random_seed', nn_control, nn_control_default, default_annoy_random_seed)
      assertthat::assert_that(assertthat::is.count(nn_control_out[['n_trees']]))
    }

    if(bitwAnd(mode, 2)) {
      nn_control_out[['search_k']] <- select_annoy_search_k(mode, nn_control, nn_control_default, nn_index, k, default_n_trees, default_k)
      nn_control_out[['grain_size']] <- select_nn_parameter_value('grain_size', nn_control, nn_control_default, default_grain_size)
      nn_control_out[['cores']] <- select_nn_parameter_value('cores', nn_control, nn_control_default, default_cores)
      assertthat::assert_that(assertthat::is.count(nn_control_out[['search_k']]))
      assertthat::assert_that(assertthat::is.count(nn_control_out[['grain_size']]))
      assertthat::assert_that(assertthat::is.count(nn_control_out[['cores']]))
    }

  } else
  if(nn_control_out[['method']] == 'hnsw') {
    nn_control_out[['metric']] <- select_nn_parameter_value('metric', nn_control, nn_control_default, default_metric)
    assertthat::assert_that(nn_control_out[['metric']] %in% c('euclidean', 'l2', 'cosine', 'ip'),
                            msg=paste0("set_nn_control: nearest neighbor metric for HNSW must be one of 'euclidean', 'l2', 'cosine', or 'ip'"))
    nn_control_out[['grain_size']] <- select_nn_parameter_value('grain_size', nn_control, nn_control_default, default_grain_size)
    nn_control_out[['cores']] <- select_nn_parameter_value('cores', nn_control, nn_control_default, default_cores)
    assertthat::assert_that(assertthat::is.count(nn_control_out[['grain_size']]))
    assertthat::assert_that(assertthat::is.count(nn_control_out[['cores']]))

    if(bitwAnd(mode, 1)) {
      nn_control_out[['M']] <- select_nn_parameter_value('M', nn_control, nn_control_default, default_M)
      nn_control_out[['ef_construction']] <- select_nn_parameter_value('ef_construction', nn_control, nn_control_default, default_ef_construction)
      assertthat::assert_that(assertthat::is.count(nn_control_out[['M']]))
      assertthat::assert_that(assertthat::is.count(nn_control_out[['ef_construction']]))
      assertthat::assert_that(nn_control_out[['M']] >= 2,
                              msg=paste0('set_nn_control: HNSW nearest neighbor M index build parameter must be >= 2'))
    }

    if(bitwAnd(mode, 2)) {
      assertthat::assert_that(assertthat::is.count(k),
                              msg=paste0('set_nn_control: parameter k must be set for method=\'hnsw\' and modes 2 and 3.'))
      nn_control_out[['ef']] <- select_nn_parameter_value('ef', nn_control, nn_control_default, default_ef)
      assertthat::assert_that(assertthat::is.count(nn_control_out[['ef']]))
      assertthat::assert_that(nn_control_out[['ef']] >= k,
                              msg=paste0('set_nn_control: HNSW nearest neighbor ef index search parameter must be >= k (',k,')'))
    }
  }
  else
    stop('set_nn_control: unsupported nearest neighbor method \'', nn_control_out[['method']], '\'')

  if(verbose) {
    cs <- get_call_stack_as_string()
    message('set_nn_control: call stack: ', cs)
    report_nn_control('  nn_control: ', nn_control_out)
  }

  if(!is.null(nn_control[['show_values']]) && nn_control[['show_values']] == TRUE)
  {
    report_nn_control('  nn_control: ', nn_control=nn_control_out)
    stop_no_noise()
  }

  return(nn_control_out)
}


# Report nn_control list values.
report_nn_control <- function(label=NULL, nn_control) {
  indent <- ''
  if(!is.null(label)) {
    indent <- '  '
  }

  message(ifelse(!is.null(label), label, ''))

  message(indent, '  method: ', ifelse(!is.null(nn_control[['method']]), nn_control[['method']], as.character(NA)))
  message(indent, '  metric: ', ifelse(!is.null(nn_control[['metric']]), nn_control[['metric']], as.character(NA)))

  if(is.null(nn_control[['method']])) {
    return()
  }

  if(nn_control[['method']] == 'nn2') {
    message(indent, '  nn2 has no parameters')
  }
  else
  if(nn_control[['method']] == 'annoy') {
    message(indent, '  n_trees: ', ifelse(!is.null(nn_control[['n_trees']]), nn_control[['n_trees']], as.character(NA)))
    message(indent, '  search_k: ', ifelse(!is.null(nn_control[['search_k']]), nn_control[['search_k']], as.character(NA)))
    message(indent, '  cores: ', ifelse(!is.null(nn_control[['cores']]), nn_control[['cores']], as.character(NA)))
    message(indent, '  grain_size: ', ifelse(!is.null(nn_control[['grain_size']]), nn_control[['grain_size']], as.character(NA)))
  }
  else
  if(nn_control[['method']] == 'hnsw') {
    message(indent, '  M: ', ifelse(!is.null(nn_control[['M']]), nn_control[['M']], as.character(NA)))
    message(indent, '  ef_construction: ', ifelse(!is.null(nn_control[['ef_construction']]), nn_control[['ef_construction']], as.character(NA)))
    message(indent, '  ef: ', ifelse(!is.null(nn_control[['ef']]), nn_control[['ef']], as.character(NA)))
    message(indent, '  cores: ', ifelse(!is.null(nn_control[['cores']]), nn_control[['cores']], as.character(NA)))
    message(indent, '  grain_size: ', ifelse(!is.null(nn_control[['grain_size']]), nn_control[['grain_size']], as.character(NA)))
  }
  else
    stop('report_nn_control: unsupported nearest neighbor method \'', nn_control[['method']], '\'')
}


new_annoy_index <- function(metric, ndim) {
  nn_class <- switch( metric,
                      cosine = RcppAnnoy::AnnoyAngular,
                      euclidean = RcppAnnoy::AnnoyEuclidean,
                      hamming = RcppAnnoy::AnnoyHamming,
                      manhattan = RcppAnnoy::AnnoyManhattan,
                      stop('unsupported annoy metric ', metric)
                    )
  nn_index <- methods::new(nn_class, ndim)
  return(nn_index)
}


#' @title Make a nearest neighbor index.
#'
#' @description Make a nearest neighbor index from the
#'   subject_matrix using either the default nearest
#'   neighbor method or the method specified in the
#'   nn_control list parameter. The function returns
#'   the index.
#'
#' @param subject_matrix the matrix used to build the
#'   index.
#' @param nn_control a list of parameters used to make the
#'   nearest neighbor index. See the set_nn_control help
#'   for details.
#' @param verbose a boolean indicating whether to emit verbose output.
#'
#' @return a nearest neighbor index.
#'
#' @examples
#'   \donttest{
#'     cds <- load_a549()
#'     cds <- preprocess_cds(cds)
#'     nn_index <- make_nn_index(SingleCellExperiment::reducedDims(cds)[['PCA']])
#'   }
#'
#' @importFrom utils packageVersion
#' @export
make_nn_index <- function(subject_matrix, nn_control=list(), verbose=FALSE) {
  assertthat::assert_that(methods::is(subject_matrix, 'matrix') ||
                          is_sparse_matrix(subject_matrix),
    msg=paste0('make_nn_matrix: the subject_matrix object must be of type matrix'))

  # We check the index build parameters.
  nn_control_default <- get_global_variable('nn_control_annoy_euclidean')
  nn_control <- set_nn_control(mode=1, nn_control=nn_control, nn_control_default=nn_control_default, nn_index=NULL, k=NULL, verbose=verbose)

  if(verbose) {
    message('make_nn_index:')
    report_nn_control('  nn_control: ', nn_control)
    tick('make_nn_index: build time')
  }

  nn_method <- nn_control[['method']]
  metric = nn_control[['metric']]

  num_row <- nrow(subject_matrix)
  num_col <- ncol(subject_matrix)
  if(!is.null(rownames(subject_matrix)))
    checksum_rownames <- digest::digest(sort(rownames(subject_matrix)))
  else
    checksum_rownames <- NA_character_

  if(nn_method == 'nn2') {
    stop('make_nn_index is not valid for method nn2')
  } else
  if(nn_method == 'annoy') {
    monocle3_annoy_index_version <- get_global_variable('monocle3_annoy_index_version')
    annoy_index <- new_annoy_index(metric, num_col)
    annoy_random_seed <- nn_control[['annoy_random_seed']]
    annoy_index$setSeed(annoy_random_seed)
    n_trees <- nn_control[['n_trees']]
    if(num_row > 0 ) {
      for(i in 1:num_row)
        annoy_index$addItem(i-1, subject_matrix[i,])
      annoy_index$build(n_trees)
    }
    annoy_index_version <- packageVersion('RcppAnnoy')
    nn_index <- list(method='annoy', annoy_index=annoy_index, version=monocle3_annoy_index_version, annoy_index_version=annoy_index_version, metric=metric, n_trees=n_trees, nrow=num_row, ncol=num_col, checksum_rownames=checksum_rownames, annoy_random_seed=annoy_random_seed)
  }
  else
  if(nn_method == 'hnsw') {
    monocle3_hnsw_index_version <- get_global_variable('monocle3_hnsw_index_version')
    M <- nn_control[['M']]
    ef_construction <- nn_control[['ef_construction']]
    hnsw_index <- RcppHNSW::hnsw_build(X=subject_matrix,
                                       distance=metric,
                                       M=M,
                                       ef=ef_construction,
                                       verbose=verbose,
                                       n_threads=nn_control[['cores']],
                                       grain_size=nn_control[['grain_size']])
    hnsw_index_version <- packageVersion('RcppHNSW')
    nn_index <- list(method='hnsw', hnsw_index=hnsw_index, version=monocle3_hnsw_index_version, hnsw_index_version=hnsw_index_version, metric=metric, M=M, ef_construction=ef_construction, nrow=num_row, ncol=num_col, checksum_rownames=checksum_rownames)
  }
  else
    stop('make_nn_index: unsupported nearest neighbor index type \'', nn_method, '\'')

  if(verbose) {
    tock()
  }

  return(nn_index)
}


#' @title Set a nearest neighor index in the cell_data_set.
#'
#' @description Store the given nearest neighbor index
#'   in the cell_data_set. The reduction_method parameter
#'   tells set_cds_nn_index where in the cell_data_set to
#'   store the index.
#'
#' @param cds a cell_data_set in which to store the nearest
#'   neighbor index.
#' @param reduction_method a string giving the reduced
#'   dimension matrix used to make the nn_index nearest
#'   neighbor index, and determines where the index is
#'   stored in the cell_data_set.
#' @param nn_index a nearest neighbor index to store in cds.
#' @param verbose a boolean indicating whether to emit verbose output.
#'
#' @return a cell_data_set with the stored index.
#'
#' @export
set_cds_nn_index <- function(cds, reduction_method=c('UMAP', 'PCA', 'LSI', 'Aligned', 'tSNE'), nn_index, verbose=FALSE) {
  assertthat::assert_that(methods::is(cds, 'cell_data_set'),
                          msg=paste('cds parameter is not a cell_data_set'))

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(reduction_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "reduction_method must be one of 'PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP'")
  reduction_method <- match.arg(reduction_method)

  assertthat::assert_that(!is.null(SingleCellExperiment::reducedDims(cds)[[reduction_method]]),
                          msg = paste0("When reduction_method = '", reduction_method,
                                      "' the cds must have been processed for it.",
                                      " Please run the required processing function",
                                      " before this one."))

  nn_method <- nn_index[['method']]

  if(nn_method == 'nn2') {
    stop('set_cds_nn_index is not valid for method nn2')
  } else
  if(nn_method == 'annoy') {
    cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]] <- S4Vectors::SimpleList()
    cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['nn_index']] <- nn_index
    cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['matrix_id']] <- get_reduce_dim_matrix_identity(cds, reduction_method)[['matrix_id']]
  }
  else
  if(nn_method == 'hnsw') {
    cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]] <- S4Vectors::SimpleList()
    cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['nn_index']] <- nn_index
    cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['matrix_id']] <- get_reduce_dim_matrix_identity(cds, reduction_method)[['matrix_id']]
  }
  else
    stop('set_cds_nn_index: unsupported nearest neighbor index type \'', nn_method, '\'')

  return(cds)
}


#' @title Make and store a nearest neighbor index in the cell_data_set.
#'
#' @description Make a nearest neighbor index from the specified
#'   reduction_method matrix in the cell_data_set using either the
#'   default nearest neighbor method or the method specified in the
#'   nn_control list parameter, and store the index in the
#'   cell_data_set. This function returns a cell_data_set.
#'
#' @param cds a cell_data_set with the reduced dimension matrix from
#'   which to make the nearest neighbor index and with which the index
#'   is stored.
#' @param reduction_method a string giving the reduced dimension matrix
#'   to use for making the nn_index nearest neighbor index.
#'   Note: distances in tSNE space reflect spatial differences poorly
#'   so using nearest neighbors with it may be meaningless.
#' @param nn_control a list of parameters to use for making the nearest
#'   neighbor index. See the set_nn_control help for details.
#' @param verbose a boolean indicating whether to emit verbose output.
#'
#' @return a cell_data_set with the stored index.
#'
#' @examples
#'   \donttest{
#'     cds <- load_a549()
#'     cds <- preprocess_cds(cds)
#'     cds <- make_cds_nn_index(cds, 'PCA')
#'   }
#'
#' @export
make_cds_nn_index <- function(cds, reduction_method=c('UMAP', 'PCA', 'LSI', 'Aligned', 'tSNE'), nn_control=list(), verbose=FALSE) {
  assertthat::assert_that(methods::is(cds, 'cell_data_set'),
                          msg=paste('cds parameter is not a cell_data_set'))

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(reduction_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "reduction_method must be one of 'PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP'")
  reduction_method <- match.arg(reduction_method)

  assertthat::assert_that(!is.null(SingleCellExperiment::reducedDims(cds)[[reduction_method]]),
                          msg = paste0("When reduction_method = '", reduction_method,
                                      "' the cds must have been processed for it.",
                                      " Please run the required processing function",
                                      " before this one."))

  nn_control_default <- get_global_variable('nn_control_annoy_euclidean')
  nn_control <- set_nn_control(mode=1, nn_control=nn_control, nn_control_default=nn_control_default, nn_index=NULL, k=NULL, verbose=verbose)

  nn_method <- nn_control[['method']]

  if(nn_method == 'nn2') {
    stop('make_cds_nn_index is not valid for method nn2')
  }

  reduced_matrix <- SingleCellExperiment::reducedDims(cds)[[reduction_method]]
  nn_index <- make_nn_index(subject_matrix=reduced_matrix, nn_control=nn_control, verbose=verbose)
  cds <- set_cds_nn_index(cds=cds, reduction_method=reduction_method, nn_index=nn_index, verbose=verbose)

  return(cds)
}


# Return the nn_index that was made from the reduction_method
# reduced dimension matrix.
get_cds_nn_index <- function(cds, reduction_method=c('UMAP', 'PCA', 'LSI', 'Aligned', 'tSNE'), nn_method, verbose=FALSE) {

  assertthat::assert_that(methods::is(cds, 'cell_data_set'),
                          msg=paste('cds parameter is not a cell_data_set'))

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(reduction_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "reduction_method must be one of 'PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP'")
  reduction_method <- match.arg(reduction_method)

  assertthat::assert_that(
    !is.null(cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]]),
    msg = paste("There is no nearest neighbor index",
                "for reduction_method",
                reduction_method,
                "and nearest neighbor method",
                nn_method))

  assertthat::assert_that(nn_method %in% c('annoy', 'hnsw'),
    msg = paste0('get_cds_nn_index: unsupported nearest neighbor index type \'',
                nn_method, '\'.'))

  if(nn_method == 'annoy') {
    nn_index <- cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['nn_index']]
  }
  else
  if(nn_method == 'hnsw') {
    nn_index <- cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['nn_index']]
  }
  else
    stop('get_cds_nn_index: unsupported nearest neighbor index type \'', nn_method, '\'')

  return(nn_index)
}


# Check that the annoy index exists.
# Returns logical TRUE if the index exists.
test_annoy_index <- function(nn_index, verbose=FALSE) {
  res <- TRUE
  if(is.null(nn_index[['annoy_index']])) {
    if(!verbose) {
      cs <- get_call_stack_as_string()
      message('test_annoy_index: the annoy nearest neighbor does not exist\ncall stack: ', cs)
    }
    else {
      message('test_annoy_index: the annoy nearest neighbor does not exist.')
    }
    return(FALSE)
  }

  tryCatch( {
    dist_res <- nn_index[['annoy_index']]$getDistance(0,1)
  },
  error=function(emsg) {
    if(!verbose) {
      cs <- get_call_stack_as_string()
      message('test_annoy_index: the annoy nearest neighbor does not exist\ncall stack: ', cs)
    }
    else {
      message('test_annoy_index: the annoy nearest neighbor does not exist.')  
    }
    res <<- FALSE
  } )

  return(res)
}


# Check that the hnsw index exists.
# Returns logical TRUE if the index exists.
test_hnsw_index <- function(nn_index, verbose=FALSE) {
  res <- TRUE
  if(is.null(nn_index[['hnsw_index']])) {
    if(!verbose) {
      cs <- get_call_stack_as_string()
      message('test_hnsw_index: the hnsw nearest neighbor does not exist\ncall stack: ', cs)
    }
    else {
      message('test_hnsw_index: the hnsw nearest neighbor does not exist.')
    }
    return(FALSE)
  }

  tryCatch( {
    size_res <- nn_index[['hnsw_index']]$size()
  },
  error=function(emsg) {
    if(!verbose) {
      cs <- get_call_stack_as_string()
      message('test_hnsw_index: the hnsw nearest neighbor does not exist\ncall stack: ', cs)
    }
    else {
      message('test_hnsw_index: the hnsw nearest neighbor does not exist.')
    }
    res <<- FALSE
  } )

  return(res)
}


search_nn_annoy_index <- function(query_matrix, nn_index, metric, k, search_k, beg_row_index, end_row_index) {
  assertthat::assert_that(beg_row_index <= end_row_index,
                          msg=paste0('search_nn_annoy_index: beg_row_index must be <= end_row_index'))

  if(!is.null(nn_index[['version']]))
    annoy_index <- nn_index[['annoy_index']]
  else
  if(!is.null(nn_index[['type']]))
    annoy_index <- nn_index[['ann']]
  else
    annoy_index <- nn_index
  nrow <- end_row_index - beg_row_index + 1
  idx <- matrix(nrow=nrow, ncol=k)
  dists <- matrix(nrow=nrow, ncol=k)
  num_bad <- 0
  offset <- beg_row_index - 1
  for(i in seq(1, nrow)) {
    nn_list <- annoy_index$getNNsByVectorList(query_matrix[i+offset,], k, search_k, TRUE)
    if(length(nn_list$item) != k)
      num_bad <- num_bad + 1
    idx[i,] <- nn_list$item
    dists[i,] <- nn_list$distance
  }
  if(metric == 'cosine') {
    dists <- 0.5 * dists * dists
  }
  return(list(idx=idx+1, dists=dists, num_bad=num_bad))
}


#' @title Search a nearest neighbor index.
#'
#' @description Search a nearest neighbor index for cells near
#' those in the query_matrix.
#'
#' @param query_matrix a reduced dimension matrix used to find the
#'   nearest neighbors in the index nn_index.
#' @param nn_index a nearest_neighbor index.
#' @param k an integer for the number of nearest neighbors to return for
#'   each cell. Default is 25.
#' @param nn_control a list of parameters used to search the nearest
#'  neighbor index. See the set_nn_control help for details. Note: the
#'  default annoy search_k parameter value is set to the default value
#'  of 2 * n_trees * k. It does not know the value of n_trees that was
#'  used to build the annoy index so if a non-default n_trees value
#'  was used to build the index, you may need to set search_k in
#'  nn_control list when you run search_nn_index.
#' @param verbose a boolean indicating whether to emit verbose output.
#'
#' @return a list list(nn.idx, nn.dists) where nn.idx is
#'  a matrix of nearest neighbor indices and nn.dists is a matrix of
#'  the distance between the index given by the row number and the index
#'  given in nn.idx. If the same reduced dim matrix is used to make the
#'  index and search the index, the index given by the row number should
#'  be in the row, usually in the first column.
#'
#' @examples
#'   \donttest{
#'     cds <- load_a549()
#'     cds <- preprocess_cds(cds)
#'     nn_index <- make_nn_index(SingleCellExperiment::reducedDims(cds)[['PCA']])
#'     nn_res <- search_nn_index(SingleCellExperiment::reducedDims(cds)[['PCA']], nn_index, 10)
#'   }
#'
#' @export
search_nn_index <- function(query_matrix, nn_index, k=25, nn_control=list(), verbose=FALSE) {
  assertthat::assert_that(methods::is(query_matrix, 'matrix') ||
                          is_sparse_matrix(query_matrix),
    msg=paste0('make_nn_matrix: the query_matrix object must be of type matrix'))

  assertthat::assert_that(assertthat::is.count(k))

  nn_control_default <- get_global_variable('nn_control_annoy_euclidean')
  nn_control <- set_nn_control(mode=2, nn_control=nn_control, nn_control_default=nn_control_default, nn_index=nn_index, k=k, verbose=verbose)

  nn_method <- nn_control[['method']]

  assertthat::assert_that(nn_method %in% c('annoy', 'hnsw'),
    msg = paste0('search_nn_index: unsupported nearest neighbor index type \'',
                nn_method, '\'.'))

  if(verbose) {
    message('search_nn_index:')
    message('  k: ', k)
    report_nn_control('  nn_control: ', nn_control)
    tick('search_nn_index: search time:')
  }

  k <- min(k, nrow(query_matrix))

  cores <- nn_control[['cores']]

  if(nn_method == 'nn2') {
    stop('search_nn_index is not valid for method nn2')
  } else
  if(nn_method == 'annoy') {
    if(!test_annoy_index(nn_index=nn_index, verbose=verbose)) {
      stop_no_noise()
    }

    # notes:
    #   o  set list names to nn.idx and nn.dists for compatibility with
    #      RANN::nn2()
    #
    num_row <- nrow(query_matrix)
    idx <- matrix(nrow=num_row, ncol=k)
    dists <- matrix(nrow=num_row, ncol=k)
    metric <- nn_control[['metric']]
    search_k <- nn_control[['search_k']]

    if(cores <= 1) {
      nn_res <- search_nn_annoy_index(query_matrix=query_matrix, nn_index=nn_index, metric=metric, k=k, search_k=search_k, beg_row_index=1, end_row_index=num_row)
      if(nn_res[['num_bad']])
        stop('annoy was unable to find ', k, ' nearest neighbors for ', nn_res[["num_bad"]], ' rows. You may need to increase the n_trees and/or search_k parameter values.')
      nn_res <- list(nn.idx=nn_res[['idx']], nn.dists=nn_res[['dists']])
    }
    else {
      omp_num_threads <- get_global_variable('omp_num_threads')
      blas_num_threads <- get_global_variable('blas_num_threads')

      RhpcBLASctl::omp_set_num_threads(1L)
      RhpcBLASctl::blas_set_num_threads(1L)

      if(cores > num_row)
        cores <- num_row
      tasks <- tasks_per_block(num_row, cores)
      beg_block <- c(0,cumsum(tasks))[1:cores]+1
      end_block <- cumsum(tasks)
      nn_blocks <- list()
      inplan <- future::plan()
      future::plan(future::multicore, workers=cores)
      on.exit(future::plan(inplan), add=TRUE)
      for(iblock in seq(cores)) {
        nn_blocks[[iblock]] <- future::future( { search_nn_annoy_index(query_matrix=query_matrix, nn_index=nn_index, metric=metric, k=k, search_k=search_k, beg_row_index=beg_block[[iblock]], end_row_index=end_block[[iblock]]) })
      }
      tot_bad <- 0
      for(iblock in seq(cores)) {
        nn_res <- future::value(nn_blocks[[iblock]])
        if(nrow(nn_res[['idx']]) != end_block[[iblock]] - beg_block[[iblock]] + 1) {
          stop('bad row count in nn_res')
        }
        idx[beg_block[[iblock]]:end_block[[iblock]],] <- nn_res[['idx']]
        dists[beg_block[[iblock]]:end_block[[iblock]],] <- nn_res[['dists']]
        tot_bad <- tot_bad + nn_res[['num_bad']]
      }
      if(tot_bad)
        stop('annoy was unable to find ', k, ' nearest neighbors for ', tot_bad, ' rows. You may need to increase the n_trees and/or search_k parameter values.')
      nn_res <- list(nn.idx=idx, nn.dists=dists)

      RhpcBLASctl::omp_set_num_threads(as.integer(omp_num_threads))
      RhpcBLASctl::blas_set_num_threads(as.integer(blas_num_threads))
    }
  }
  else
  if(nn_method == 'hnsw') {
    if(!test_hnsw_index(nn_index=nn_index, verbose=verbose)) {
      stop_no_noise()
    }

    assertthat::assert_that(nn_control[['ef']] >= k,
      msg=paste0('search_nn_index: ef must be >= k'))

    tmp <- RcppHNSW::hnsw_search(X=query_matrix,
                                  ann=nn_index[['hnsw_index']],
                                  k=k,
                                  ef=nn_control[['ef']],
                                  verbose=verbose,
                                  n_threads=cores,
                                  grain_size=nn_control[['grain_size']])
    # The RcppHNSW documentation says that the L2 metric is the square
    # of the Euclidean distance but my tests indicate that the L2 metric
    # returns sqrt(sum((SingleCellExperiment::reducedDims(cds)[['UMAP']][i,] - SingleCellExperiment::reducedDims(cds)[['UMAP']][j,])^2)).
    # Additionally, nn2, annoy, and hnsw return the same distance values for the euclidean
    # metric.
    nn_res <- list(nn.idx=tmp[['idx']], nn.dists=tmp[['dist']])
  }
  else
    stop('search_nn_index: unsupported nearest neighbor index type \'', nn_method, '\'')

  if(verbose) {
    tock()
  }

  return(nn_res)
}


#' @title Search a nearest neighbor index that is stored in
#'   the cds.
#'
#' @description Search a nearest neighbor index for cells near
#' those in the query_matrix.
#'
#' @param query_matrix a reduced dimension matrix used to find the
#'   nearest neighbors in the index nn_index.
#' @param cds a cell_data_set in which the nearest neighbor index
#'   is stored.
#' @param reduction_method a string giving the reduced
#'   dimension matrix used to make the nearest neighbor index, and
#'   determines where the index is stored in the cell_data_set.
#'   Note: distances in tSNE space reflect spatial differences poorly
#'   so using nearest neighbors with it may be meaningless.
#' @param k an integer for the number of nearest neighbors to return for
#'   each cell. Default is 25.
#' @param nn_control a list of parameters used to make and search
#'   the nearest neighbors indexes. See the set_nn_control help
#'   for additional details. Note that if nn_control\[\['search_k'\]\]
#'   is not defined, transfer_cell_labels will try to use
#'   search_k <- 2 * n_trees * k where n_trees is the value used
#'   to build the index. The default metric is cosine for
#'   reduction_methods PCA, LSI, and Aligned, and is euclidean for
#'   reduction_methods tSNE and UMAP.
#' @param verbose a boolean indicating whether to emit verbose output.
#'
#' @return a list list(nn.idx, nn.dists) where nn.idx is
#'  a matrix of nearest neighbor indices and nn.dists is a matrix of
#'  the distance between the index given by the row number and the index
#'  given in nn.idx. If the same reduced dim matrix is used to make the
#'  index and search the index, the index given by the row number should
#'  be in the row, usually in the first column.
#'
#' @examples
#'   \donttest{
#'     cds <- load_a549()
#'     cds <- preprocess_cds(cds)
#'     cds <- make_cds_nn_index(cds, 'PCA')
#'     nn_res <- search_cds_nn_index(SingleCellExperiment::reducedDims(cds)[['PCA']], cds, 'PCA', 10)
#'   }
#'
#' @export
search_cds_nn_index <- function(query_matrix, cds, reduction_method=c('UMAP', 'PCA', 'LSI', 'Aligned', 'tSNE'), k=25, nn_control=list(), verbose=FALSE) {
  assertthat::assert_that(methods::is(query_matrix, 'matrix') ||
                          is_sparse_matrix(query_matrix),
    msg=paste0('make_nn_matrix: the query_matrix object must be of type matrix'))
  
  assertthat::assert_that(methods::is(cds, 'cell_data_set'),
                          msg=paste('cds parameter is not a cell_data_set'))

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(reduction_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "reduction_method must be one of 'PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP'")
  reduction_method <- match.arg(reduction_method)

  assertthat::assert_that(!is.null(SingleCellExperiment::reducedDims(cds)[[reduction_method]]),
                          msg = paste0("When reduction_method = '", reduction_method,
                                      "' the cds must have been processed for it.",
                                      " Please run the required processing function",
                                      " before this one."))

  assertthat::assert_that(assertthat::is.count(k))

  if(reduction_method == 'tSNE' || reduction_method == 'UMAP')
    nn_control_default <- get_global_variable('nn_control_annoy_euclidean')
  else
    nn_control_default <- get_global_variable('nn_control_annoy_cosine')

  # Use set_nn_control to find nn method, which we need in order to select the correct index,
  # and we may need the index to set the nn_control[['search_k']].
  nn_control_tmp <- set_nn_control(mode=2,
                                   nn_control=nn_control,
                                   nn_control_default=nn_control_default,
                                   nn_index=NULL,
                                   k=k,
                                   verbose=verbose)
  nn_index <- get_cds_nn_index(cds=cds, reduction_method=reduction_method, nn_control_tmp[['method']], verbose=FALSE)

  nn_control <- set_nn_control(mode=2,
                               nn_control=nn_control,
                               nn_control_default=nn_control_default,
                               nn_index=nn_index,
                               k=k,
                               verbose=verbose)
  nn_res <- search_nn_index(query_matrix=query_matrix,
                            nn_index=nn_index,
                            k=k,
                            nn_control=nn_control,
                            verbose=verbose)

  return(nn_res)
}

# Notice that this function uses the available SingleCellExperiment::reducedDims(cds)[[reduction_method]] matrix
# for verifying the query matrix -- be certain that the search (and index build) were run
# on this matrix so call this function immediately after running the search.
# Notes:
#   o  this function was written to store nn search information in the cds when
#      we intended to re-use an index for later processing stages but
#      it appears that there is little or no opportunity to re-use
#      indices so we are not storing search information at this time.
# This function is not in use currently and may fall into disrepair.
set_cds_nn_search <- function(cds, reduction_method=c('UMAP', 'PCA', 'LSI', 'Aligned', 'tSNE'), search_id, ef, k, nn_control=list(), verbose=TRUE) {
  assertthat::assert_that(methods::is(cds, 'cell_data_set'),
                          msg=paste('cds parameter is not a cell_data_set'))

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(reduction_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "reduction_method must be one of 'PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP'")
  reduction_method <- match.arg(reduction_method)

  assertthat::assert_that(!is.null(SingleCellExperiment::reducedDims(cds)[[reduction_method]]),
                          msg = paste0("When reduction_method = '", reduction_method,
                                      "' the cds must have been processed for it.",
                                      " Please run the required processing function",
                                      " before this one."))
  assertthat::assert_that(!is.null(search_id),
                          msg = paste0('You must give a search_id string.'))

  assertthat::assert_that(!is.null(k),
                          msg = paste0('You must give a k value.'))

  reduced_matrix <- SingleCellExperiment::reducedDims(cds)[[reduction_method]]

  nn_control_default <- get_global_variable('nn_control_annoy_euclidean')
  nn_control <- set_nn_control(mode=2, nn_control=nn_control, nn_control_default=nn_control_default, nn_index=NULL, k=k, verbose=verbose)

  nn_method <- nn_control[['method']]

  if(is.null(cds@reduce_dim_aux[[reduction_method]][['nn_search']])) {
    cds@reduce_dim_aux[[reduction_method]][['nn_search']] <- S4Vectors::SimpleList()
  }

  if(nn_method == 'nn2') {
    cds@reduce_dim_aux[[reduction_method]][['nn_search']][[search_id]] <- S4Vectors::SimpleList()
    cds@reduce_dim_aux[[reduction_method]][['nn_search']][[search_id]][['method']] <- nn_method
    cds@reduce_dim_aux[[reduction_method]][['nn_search']][[search_id]][['k']] <- k
    cds@reduce_dim_aux[[reduction_method]][['nn_search']][[search_id]][['nrow']] <- nrow(reduced_matrix)
    cds@reduce_dim_aux[[reduction_method]][['nn_search']][[search_id]][['ncol']] <- ncol(reduced_matrix)
  }
  else
  if(nn_method == 'annoy') {
    cds@reduce_dim_aux[[reduction_method]][['nn_search']][[search_id]] <- S4Vectors::SimpleList()
    cds@reduce_dim_aux[[reduction_method]][['nn_search']][[search_id]][['method']] <- nn_method
    cds@reduce_dim_aux[[reduction_method]][['nn_search']][[search_id]][['search_k']] <- nn_contol[['search_k']]
    cds@reduce_dim_aux[[reduction_method]][['nn_search']][[search_id]][['k']] <- k
    cds@reduce_dim_aux[[reduction_method]][['nn_search']][[search_id]][['nrow']] <- nrow(reduced_matrix)
    cds@reduce_dim_aux[[reduction_method]][['nn_search']][[search_id]][['ncol']] <- ncol(reduced_matrix)
    cds@reduce_dim_aux[[reduction_method]][['nn_search']][[search_id]][['grain_size']] <- nn_control[['grain_size']]
    cds@reduce_dim_aux[[reduction_method]][['nn_search']][[search_id]][['cores']] <- nn_control[['cores']]
  }
  else
  if(nn_method == 'hnsw') {
    cds@reduce_dim_aux[[reduction_method]][['nn_search']][[search_id]] <- S4Vectors::SimpleList()
    cds@reduce_dim_aux[[reduction_method]][['nn_search']][[search_id]][['method']] <- nn_method
    cds@reduce_dim_aux[[reduction_method]][['nn_search']][[search_id]][['ef']] <- ef
    cds@reduce_dim_aux[[reduction_method]][['nn_search']][[search_id]][['k']] <- k
    cds@reduce_dim_aux[[reduction_method]][['nn_search']][[search_id]][['nrow']] <- nrow(reduced_matrix)
    cds@reduce_dim_aux[[reduction_method]][['nn_search']][[search_id]][['ncol']] <- ncol(reduced_matrix)
    cds@reduce_dim_aux[[reduction_method]][['nn_search']][[search_id]][['grain_size']] <- nn_control[['grain_size']]
    cds@reduce_dim_aux[[reduction_method]][['nn_search']][[search_id]][['cores']] <- nn_control[['cores']]
  }
  else {
    stop('set_cds_nn_search: unsupported nearest neighbor index type \'', nn_method, '\'')
  }

  return(cds)
}


# Get nearest neighbor search information stored in the cds.
# This function is not in use currently and may fall into disrepair.
get_cds_nn_search <- function(cds, reduction_method=c('UMAP', 'PCA', 'LSI', 'Aligned', 'tSNE'), search_id, verbose=TRUE) {

  assertthat::assert_that(methods::is(cds, 'cell_data_set'),
                          msg=paste('cds parameter is not a cell_data_set'))

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(reduction_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "reduction_method must be one of 'PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP'")
  reduction_method <- match.arg(reduction_method)
  
  assertthat::assert_that(!is.null(SingleCellExperiment::reducedDims(cds)[[reduction_method]]),
                          msg = paste0("When reduction_method = '", reduction_method,
                                      "' the cds must have been processed for it.",
                                      " Please run the required processing function",
                                      " before this one."))
  assertthat::assert_that(!is.null(search_id),
                          msg = paste0('You must give a search_id string.'))

  assertthat::assert_that(!is.null(cds@reduce_dim_aux[[reduction_method]][['nn_search']][[search_id]]),
                          msg = paste0("Stopping because there is no information in the cds for search_id = '",
                                       search_id, "'."))

  nn_search <- cds@reduce_dim_aux[[reduction_method]][['nn_search']][[search_id]]

  nn_method <- nn_search[['method']]

  if(nn_method == 'nn2') {
    nn_control <- list(method=nn_method,
                       grain_size=nn_search[['grain_size']],
                       cores=nn_search[['cores']])
  }
  else
  if(nn_method == 'annoy') {
    nn_control <- list(method=nn_method,
                       search_k=nn_search[['search_k']],
                       grain_size=nn_search[['grain_size']],
                       cores=nn_search[['cores']])
  }
  else
  if(nn_method == 'hnsw') {
    nn_control <- list(method=nn_method,
                       ef=nn_search[['ef']],
                       grain_size=nn_search[['grain_size']],
                       cores=nn_search[['cores']])
  }
  else {
    stop('get_cds_nn_search: unsupported nearest neighbor search type \'', nn_method, '\'')
  }
                     
  return(nn_control)
}


# Get index build and search information that's stored in the cds.
# This function is not used at this time. It may fall into disrepair.
get_cds_nn_control <- function(cds, reduction_method=c('UMAP', 'PCA', 'LSI', 'Aligned', 'tSNE'), search_id, verbose=TRUE) {

  assertthat::assert_that(methods::is(cds, 'cell_data_set'),
                          msg=paste('cds parameter is not a cell_data_set'))

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(reduction_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "reduction_method must be one of 'PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP'")
  reduction_method <- match.arg(reduction_method)

  assertthat::assert_that(!is.null(SingleCellExperiment::reducedDims(cds)[[reduction_method]]),
                          msg = paste0("When reduction_method = '", reduction_method,
                                      "' the cds must have been processed for it.",
                                      " Please run the required processing function",
                                      " before this one."))
  assertthat::assert_that(!is.null(search_id),
                          msg = paste0('You must give a search_id string.'))

  assertthat::assert_that(!is.null(cds@reduce_dim_aux[[reduction_method]][['nn_search']][[search_id]]),
                          msg = paste0("Stopping because there is no information in the cds for search_id = '",
                                       search_id, "'."))

  nn_search <- cds@reduce_dim_aux[[reduction_method]][['nn_search']][[search_id]]

  nn_method <- nn_search[['method']]

  if(nn_method == 'nn2') {
    nn_control <- list(method=nn_method,
                       grain_size=nn_search[['grain_size']],
                       cores=nn_search[['cores']])
  }
  else
  if(nn_method == 'annoy') {
    nn_control <- list(method=nn_method,
                       metric=cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['nn_index']][['metric']],
                       n_trees=cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['nn_index']][['n_trees']],
                       search_k=nn_search[['search_k']],
                       grain_size=nn_search[['grain_size']],
                       cores=nn_search[['cores']])
  }
  else
  if(nn_method == 'hnsw') {
    nn_control <- list(method=nn_method,
                       metric=cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['nn_index']][['metric']],
                       M=cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['nn_index']][['M']],
                       ef_construction=cds@reduce_dim_aux[[reduction_method]][['nn_index']][[nn_method]][['nn_index']][['ef_construction']],
                       ef=nn_search[['ef']],
                       grain_size=nn_search[['grain_size']],
                       cores=nn_search[['cores']])
  }
  else {
    stop('get_cds_nn_control: unsupported nearest neighbor search type \'', nn_method, '\'')
  }

  if(verbose) {
    cs <- get_call_stack_as_string()
    message('get_cds_nn_control: (',cs,')')
    report_nn_control('  nn_control: ', nn_control=nn_control)
  }

  return(nn_control)
}


#' @title Search a subject matrix for nearest neighbors to a query_matrix.
#'
#' @description Make a nearest neighbors index using the subject matrix and
#' search it for nearest neighbors to the query_matrix.
#'
#' @param subject_matrix a matrix used to build a nearest neighbor index.
#' @param query_matrix a matrix used to search the subject_matrix
#'   nearest neighbor index.
#' @param k an integer for the number of nearest neighbors to return for
#'   each cell. Default is 25.
#' @param nn_control a list of parameters used to make and search the nearest
#'  neighbor index. See the set_nn_control help for details.
#' @param verbose a boolean indicating whether to emit verbose output.
#'
#' @return a list list(nn.idx, nn.dists) where nn.idx is
#'  a matrix of nearest neighbor indices and nn.dists is a matrix of
#'  the distance between the index given by the row number and the index
#'  given in nn.idx. If the query_matrix is the same as the subject
#'  matrix, the index given by the row number should be in the row, usually
#' in the first column.
#'
#' @export
search_nn_matrix <- function(subject_matrix, query_matrix, k=25, nn_control=list(), verbose=FALSE) {
  assertthat::assert_that(methods::is(subject_matrix, 'matrix') ||
                          is_sparse_matrix(subject_matrix),
    msg=paste0('search_nn_matrix: the subject_matrix object must be of type matrix'))

  assertthat::assert_that(methods::is(query_matrix, 'matrix') ||
                          is_sparse_matrix(query_matrix),
    msg=paste0('search_nn_matrix: the query_matrix object must be of type matrix'))

  nn_control_default <- get_global_variable('nn_control_annoy_euclidean')
  nn_control <- set_nn_control(mode=3, nn_control=nn_control, nn_control_default=nn_control_default, nn_index=NULL, k=k, verbose=verbose)

  method <- nn_control[['method']]
  k <- min(k, nrow(subject_matrix))

  if(verbose) {
    message('search_nn_matrix:')
    message('  k: ', k)
    report_nn_control('  nn_control: ', nn_control)
    tick('search_nn_matrix: search_time')
  }


  if(method == 'nn2') {
    nn_res <- RANN::nn2(subject_matrix, query_matrix, k, searchtype = "standard")
  }
  else {
    nn_index <- make_nn_index(subject_matrix, nn_control=nn_control, verbose=verbose)
    nn_res <- search_nn_index(query_matrix=query_matrix, nn_index=nn_index, k=k, nn_control=nn_control, verbose=verbose)
  }

  if(verbose)
    tock()

  return(nn_res)
}


# Check whether the neighbor in the first column == row index.
# The hnsw_search help page states
#   Every item in the dataset is considered to be a neighbor of itself,
#   so the first neighbor of item i should always be i itself.
check_nn_col1 <- function(mat) {
  irow <- seq(nrow(mat))
  if(any(mat[,1] != irow)) {
    return(FALSE)
  }
  return(TRUE)
}


# When the search query matrix is the same as the
# index matrix, we expect that the row index
# exists in the row. Count the cases for which
# this is not true and the distances are not
# all zero.
count_nn_missing_self_index <- function(nn_res, verbose=FALSE) {
  idx <- nn_res[['nn.idx']]
  dst <- nn_res[['nn.dists']]
  len <- length(idx[1,])
  num_missing <- 0

  if(nrow(idx) == 0)
    return(0)

  for (irow in 1:nrow(idx)) {
    vidx <- idx[irow,]
    vdst <- dst[irow,]
    if(vidx[[1]] != irow) {
      dself <- FALSE
      dzero <- TRUE
      if(length(vidx) == 1) {
        num_missing <- num_missing + 1
        next
      }
      for(i in seq(1, len, 1)) {
        if(vidx[[i]] == irow) {
          dself = TRUE
          break
        }
        if(vdst[[i]] > .Machine$double.xmin) {
          dzero <- FALSE
        }
      }

      if(dself == FALSE && dzero == FALSE) {
        num_missing <- num_missing + 1
      }
    }
  }

  if(verbose) {
    message('count_nn_missing_self_index:')
    message('  ', num_missing, ' out of ', nrow(nn_res[['nn.idx']]), ' rows are missing the row index')
    message('  \'recall\': ', formatC((as.double(nrow(nn_res[['nn.idx']]) - num_missing) / as.double(nrow(nn_res[['nn.idx']]))) * 100.0), "%")
  }

  return(num_missing)
}


# For searches in which the index and search sets
# are the same, sort nearest neighbor results so
# that row index point is in the first column. If the
# row index point was missed, shift the results to the
# right one place and store the row index point in the
# first element.
#
# This is required because approximate searches may
# not store the row index point in the first column if
# there is more than one point with zero distance,
# or, the function may miss the row index point.
#
# Notes:
#   o  we assume that the self index distance is zero
#   o  do not use this if the search point set differs
#      from the index build set.
#
swap_nn_row_index_point <- function(nn_res, verbose=FALSE) {

  if(verbose) {
    count_nn_missing_self_index(nn_res, verbose)
  }

  # Skip if there is no need to swap indices.
  if(check_nn_col1(nn_res$nn.idx))
    return(nn_res)

  idx <- nn_res[['nn.idx']]
  dst <- nn_res[['nn.dists']]

  diagnostics <- FALSE

  num_no_recall <- 0

  if(nrow(idx) == 0 )
    return(nn_res) 

  for (irow in 1:nrow(idx)) {
    vidx <- idx[irow,]
    vdst <- dst[irow,]
    if(vidx[[1]] != irow) {
      if(diagnostics) {
        message('swap_nn_row_index_point: adjust nn matrix row: ', irow)
        message('swap_nn_row_index_point: idx row pre fix: ', paste(vidx, collapse=' '))
      }

      if(length(vidx) == 1) {
        num_no_recall <- num_no_recall + 1
        next
      }

      if(vidx[[2]] == irow) {
        vidx[[2]] <- vidx[[1]]
        vidx[[1]] <- irow
      } else {
        match <- FALSE
        for(i in seq(2, length(vidx), 1)) {
          if(vidx[[i]] == irow) {
            vidx[[i]] <- vidx[[1]]
            vidx[[1]] <- irow
            match <- TRUE
            break
          }
        }
        if(!match) {
          if(diagnostics) {
            message('swap_nn_row_index_point: dst row pre fix: ', paste(vdst, collapse=' '))
          }
          for(i in seq(1, length(vidx), 1)) {
            if(vdst[[i]] > .Machine$double.xmin) {
              num_no_recall <- num_no_recall + 1
            }
          }
          vidx <- c(irow, vidx[1:(length(vidx)-1)])
          vdst <- c(0, vdst[1:(length(vdst)-1)])
          dst[irow,] <- vdst
          if(diagnostics) {
            message('swap_nn_row_index_point: dst row post fix: ', paste(vdst, collapse=' '))
          }
        }
      }
      idx[irow,] <- vidx
      if(diagnostics) {
        message('swap_nn_row_index_point: idx row post fix: ', paste(vidx, collapse=' '))
      }
    }
  }

  if(num_no_recall > 0) {
    frac_recall <- (nrow(idx)-num_no_recall) / nrow(idx)
    format_recall <- sprintf('%3.1f', frac_recall * 100.0)
    message('The search result is expected to include the query row value (self)\n',
            'because the NN index includes the query objects; however, this search result\n',
            'is missing ', num_no_recall, ' self values (recall: ', format_recall, '%). Monocle3 has added the self\n',
            'values to the first column of the search result in order to allow further\n',
            'analysis -- but it is missing important nearest neighbors so you need to\n',
            'increase the sensitivity for making and/or searching the index.')
  }
  nn_res_out <- list(nn.idx=idx, nn.dists=dst)

  return(nn_res_out)
}