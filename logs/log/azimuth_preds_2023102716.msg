Error in adata[[md]][[i]][] : object of type 'environment' is not subsettable
Calls: LoadFileInput ... LoadH5AD -> LoadMetadata -> sapply -> lapply -> FUN
Traceback:
6: FUN(X[[i]], ...)
5: lapply(X = X, FUN = FUN, ...)
4: sapply(X = col.names, FUN = function(i) {
       x <- adata[[md]][[i]][]
       if (i %in% factor.cols) {
           x <- factor(x = x, levels = adata[[md]][["__categories"]][[i]][])
       }
       return(x)
   }, simplify = FALSE, USE.NAMES = TRUE)
3: LoadMetadata(md = "obs")
2: LoadH5AD(path = path)
1: LoadFileInput(path = ABS_DATASET_PATH)

