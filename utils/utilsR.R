# Generate hash table
lookup.table <- function(x,y){
  stopifnot(length(x) == length(y))
  e <- new.env(hash=TRUE, size=length(x))
  for (i in 1:length(x)){
    assign(x[[i]],y[[i]], envir=e)
  }
  e
}
