# #
#
# ann_acc <- function(results, ground_truth) {
#   results <- results$idx - 1
#   ground_truth <- t(ground_truth$neighbors)
#
#   n_row_intersect(results, ground_truth) / prod(dim(ground_truth))
# }
#
# n_row_intersect <- function(A, B) {
#   nr <- nrow(A)
#   agreement <- 0
#   for (i in 1:nr) {
#     agreement <- agreement + length(intersect(A[i, ], B[i, ]))
#   }
#   agreement
# }
#
# ann_bench <- function(nn_data, M = 16, ef_construction = 200,
#                       ef = 10, verbose = TRUE) {
#   tstart <- Sys.time()
#   hnsw_res <- knn_xy(X = t(nn_data$train), Y = t(nn_data$test),
#                 distance = nn_data$metric,
#                 k = 100, M = M, ef_construction = ef_construction, ef = ef,
#                 verbose = verbose)
#   tend <- Sys.time()
#
#   acc <- ann_acc(results = hnsw_res, ground_truth = nn_data)
#   list(
#     nn = hnsw_res,
#     time = tend - tstart,
#     acc = acc
#   )
# }
#
# knn_xy <- function(X, Y, k = 10, distance = "euclidean",
#                    M = 16, ef_construction = 200, ef = ef_construction,
#                    verbose = FALSE) {
#   if (!is.matrix(X)) {
#     stop("X must be matrix")
#   }
#   if (!is.matrix(Y)) {
#     stop("Y must be matrix")
#   }
#   stopifnot(ncol(X) == ncol(Y))
#
#   nr <- nrow(X)
#   max_k <- nr
#   if (k > max_k) {
#     stop("k cannot be larger than ", max_k)
#   }
#   distance <- match.arg(distance, c("l2", "euclidean", "cosine", "ip"))
#
#   ann <- hnsw_build(X = X, distance = distance, M = M, ef = ef_construction,
#                     verbose = verbose)
#   hnsw_search(X = Y, ann = ann, k = k,
#               verbose = verbose)
# }
#
#
# rand_nn_idx <- function(nrow, ncol, include_self = FALSE) {
#   stopifnot(ncol <= nrow)
#   t(sapply(1:nrow, function(i) { sort(sample((1:nrow)[-i], ncol, replace = FALSE))}))
# }
#
#
#
# download_ann_bench <- function(url, destfile = dtemp(ext = ".hdf5"),
#                                base_url = "http://ann-benchmarks.com/",
#                                verbose = FALSE) {
#   url <- paste0(base_url, url)
#   tsmessage("Downloading ", url, " to ", destfile)
#   utils::download.file(url = url, destfile = destfile, quiet = !verbose,
#                        mode = "wb")
#
#   tsmessage("Processing ", destfile)
#   res <- hdf5file_to_ann(destfile)
#
#   tsmessage("Cleaning up")
#   unlink(destfile)
#
#   res
# }
#
# hdf5_to_ann <- function(fileh5) {
#   res <- list()
#   res$distances <- fileh5[["distances"]][,]
#   res$neighbors <- fileh5[["neighbors"]][,]
#   res$train <- fileh5[["train"]][,]
#   res$test <- fileh5[["test"]][,]
#
#   res
# }
#
# hdf5file_to_ann <- function(filename) {
#   fileh5 <- hdf5r::H5File$new(filename, mode = "r")
#   res <- hdf5_to_ann(fileh5)
#   fileh5$close_all()
#
#   res
# }
#
# dtemp <- function(tdir = "D:\\temp", ext = "") {
#   tempfile(tmpdir = tdir, fileext = ext)
# }
#
#
# load_dataset <- function(file, dir = "E:/dev/R/datasets", verbose = FALSE) {
#
#   for (i in seq_along(file)) {
#     f <- file[i]
#     if (!endsWith(f, ".Rda")) {
#       f <- paste0(f, ".Rda")
#     }
#     f <- file.path(dir, f)
#     if (verbose) {
#       message("Loading ", f)
#     }
#     load(file = f, envir = .GlobalEnv)
#   }
# }
#
# loadnndata <- function(dir = "E:/dev/R/datasets", verbose = TRUE) {
#   load_dataset(c("fashion_nn", "glove50_nn", "mnist_nn", "nytimes_nn"),
#                verbose = verbose)
# }
#
# annoy_build <- function(X, metric = "euclidean", n_trees = 50,
#                         verbose = FALSE) {
#   nr <- nrow(X)
#   nc <- ncol(X)
#
#   if (metric == "cosine") {
#     ann <- methods::new(RcppAnnoy::AnnoyAngular, nc)
#   }
#   else if (metric == "manhattan") {
#     ann <- methods::new(RcppAnnoy::AnnoyManhattan, nc)
#   }
#   else {
#     ann <- methods::new(RcppAnnoy::AnnoyEuclidean, nc)
#   }
#
#   tsmessage("Building Annoy index with metric = ", metric)
#   progress <- Progress$new(max = nr, display = verbose)
#
#   # Add items
#   for (i in 1:nr) {
#     ann$addItem(i - 1, X[i, ])
#     progress$increment()
#   }
#
#   # Build index
#   ann$build(n_trees)
#
#   ann
# }
#
# annoy_search <- function(X, k, ann,
#                          search_k = 100 * k,
#                          verbose = FALSE) {
#   if (!is.matrix(X)) {
#     stop("X must be matrix")
#   }
#   nr <- nrow(X)
#
#   max_k <- nr
#   if (k > max_k) {
#     stop("k cannot be larger than ", max_k)
#   }
#
#   idx <- matrix(nrow = nr, ncol = k)
#   dist <- matrix(nrow = nr, ncol = k)
#   tsmessage("Searching Annoy index")
#   nr <- nrow(X)
#   search_progress <- Progress$new(max = nr, display = verbose)
#   idx <- matrix(nrow = nr, ncol = k)
#   dist <- matrix(nrow = nr, ncol = k)
#   for (i in 1:nr) {
#     res <- ann$getNNsByVectorList(X[i, ], k, search_k, TRUE)
#     if (length(res$item) != k) {
#       stop(
#         "search_k/n_trees settings were unable to find ", k,
#         " neighbors for item ", i
#       )
#     }
#     idx[i, ] <- res$item
#     dist[i, ] <- res$distance
#     search_progress$increment()
#   }
#
#   list(idx = idx + 1, dist = dist)
# }
#
#
# annoy_bench <- function(nn_data,
#                         n_trees = 50, search_k = 100 * 100,
#                         verbose = TRUE) {
#   tstart <- Sys.time()
#   annoy_res <- annoy_xy(X = t(nn_data$train), Y = t(nn_data$test),
#                      distance = nn_data$metric,
#                      k = 100, n_trees = n_trees, search_k = search_k,
#                      verbose = verbose)
#   tend <- Sys.time()
#
#   acc <- ann_acc(results = annoy_res, ground_truth = nn_data)
#   list(
#     nn = annoy_res,
#     time = tend - tstart,
#     acc = acc
#   )
# }
#
# annoy_xy <- function(X, Y, k = 10, distance = "euclidean",
#                      n_trees = 50, search_k = 100 * k,
#                      verbose = FALSE) {
#   if (!is.matrix(X)) {
#     stop("X must be matrix")
#   }
#   if (!is.matrix(Y)) {
#     stop("Y must be matrix")
#   }
#   stopifnot(ncol(X) == ncol(Y))
#
#   nr <- nrow(X)
#   max_k <- nr
#   if (k > max_k) {
#     stop("k cannot be larger than ", max_k)
#   }
#   distance <- match.arg(distance, c("euclidean", "cosine", "manhattan", "hamming"))
#
#   ann <- annoy_build(X = X, metric = distance, n_trees = n_trees,
#                     verbose = verbose)
#   annoy_search(X = Y, ann = ann, k = k, verbose = verbose)
# }
#
# readme1 <- function() {
#   library(RcppHNSW)
#   data <- as.matrix(iris[, -5])
#
#   # Create a new index using the L2 (squared Euclidean) distance
#   # nr and nc are the number of rows and columns of the data to be added, respectively
#   # ef and M determines speed vs accuracy trade off
#   # You must specify the maximum number of items to add to the index when it
#   # is created. But you can increase this number: see the next example
#   M <- 16
#   ef <- 200
#   ann <- new(HnswL2, ncol(data), nrow(data), M, ef)
#
#   # Add items to index
#   for (i in 1:nrow(data)) {
#     ann$addItem(data[i, ])
#   }
#
#   # Find 4 nearest neighbors of row 1
#   # indexes are in res$item, distances in res$distance
#   # set include_distances = TRUE to get distances as well as index
#   res <- ann$getNNsList(data[1, ], k = 4, include_distances = TRUE)
#
#   # function interface returns results for all rows in nr x k matrices
#   all_knn <- RcppHNSW::hnsw_knn(data, k = 4, distance = "l2")
#   list(all_knn, res)
# }
#
# readme2 <- function() {
#   library("RcppHNSW")
#   set.seed(12345)
#
#   dim <- 16
#   num_elements <- 100000
#
#   # Generate sample data
#   data <- matrix(stats::runif(num_elements * dim), nrow = num_elements)
#
#   # Split data into two batches
#   data1 <- data[1:(num_elements / 2), ]
#   data2 <- data[(num_elements / 2 + 1):num_elements, ]
#
#   # Create index
#   M <- 16
#   ef <- 10
#   # Set the initial index size to the size of the first batch
#   p <- new(HnswL2, dim, num_elements / 2, M, ef)
#
#   message("Adding first batch of ", nrow(data1), " elements")
#   p$addItems(data1)
#
#   # Query the elements for themselves and measure recall:
#   idx <- p$getAllNNs(data1, k = 1)
#   message("Recall for the first batch: ", formatC(mean(idx == 1:nrow(data1))))
#
#   filename <- "first_half.bin"
#   # Serialize index
#   p$save(filename)
#
#   # Reinitialize and load the index
#   rm(p)
#   message("Loading index from ", filename)
#   # Increase the total capacity, so that it will handle the new data
#   p <- new(HnswL2, dim, filename, num_elements)
#
#   message("Adding the second batch of ", nrow(data2), " elements")
#   p$addItems(data2)
#
#   # Query the elements for themselves and measure recall:
#   idx <- p$getAllNNs(data, k = 1)
#   # You can get distances with:
#   # res <- p$getAllNNsList(data, k = 1, include_distances = TRUE)
#   # res$dist contains the distance matrix, res$item stores the indexes
#
#   message("Recall for two batches: ", formatC(mean(idx == 1:num_elements)))
#
#   res <- hnsw_search(data, p, k = 1, ef = 10, verbose = TRUE)
#   message("Recall for two batches: ", formatC(mean(res$idx == 1:num_elements)))
#
#   unlink(filename)
# }
