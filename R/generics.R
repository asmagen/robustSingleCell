setGeneric("cache", function(object) standardGeneric("cache"))
setMethod("cache", signature(object = "robustSingleCell"),
          function(object) object@cache)