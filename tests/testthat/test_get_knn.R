library(RcppHNSW)
context("hnsw_knn")

expect_error(hnsw_knn(uiris), "(?i)matrix")
expect_error(hnsw_knn(ui10, M = 1), "M cannot")
res <- hnsw_knn(ui10, k = 4)
expect_equal(res$idx, self_nn_index4, check.attributes = FALSE)
expect_equal(res$dist, self_nn_dist4, check.attributes = FALSE, tol = 1e-6)

res <- hnsw_knn(ui10, k = 4, distance = "l2")
expect_equal(res$idx, self_nn_index4, check.attributes = FALSE)
expect_equal(res$dist, self_nn_dist4^2, check.attributes = FALSE, tol = 1e-6)

res <- hnsw_knn(ui10, k = 1)
expect_is(res$idx, "matrix")
expect_is(res$dist, "matrix")
expect_equal(res$idx, self_nn_index4[, 1], check.attributes = FALSE)
expect_equal(res$dist, self_nn_dist4[, 1], check.attributes = FALSE, tol = 1e-6)

# test byrow = TRUE
res <- hnsw_knn(t(ui10), k = 4, byrow = FALSE)
expect_equal(t(res$idx), self_nn_index4, check.attributes = FALSE)
expect_equal(t(res$dist), self_nn_dist4, check.attributes = FALSE, tol = 1e-6)
