# Make the route ID for the model

makeID <- function(x){
    paste(x[-length(x)], x[-1])
}

# Need to define each reach as A1 - A2, not by receiver at A2
tars = ls(pattern = "routeSac")
new_vects <- sapply(tars, function(x)
    apply(get(x), 2, makeID), simplify = FALSE)

levs <- unique(unlist(new_vects))
tars <- ls(pattern = "yolo201[2|3]_vector")
new_vects <- sapply(tars, function(x)
   makeID(get(x)), simplify = FALSE)

levs <- unique(c(levs, unlist(new_vects)))

levs

# Reaches
reachSac12 <- apply(routeSac12, 2, function(x){
    match(makeID(x), levs)
})

reachSac13 <- apply(routeSac13, 2, function(x){
    match(makeID(x), levs)
})

reachYolo12 <- match(makeID(yolo2012_vector), levs)
reachYolo13 <- match(makeID(yolo2013_vector), levs)

# Reach length
sac2012_RouteA_reaches <- c(sac2012_RouteA_reaches[1:5],
                            "39" = 0,
                            sac2012_RouteA_reaches[6:8])
sac2012_RouteC_reaches <- c(sac2012_RouteC_reaches[1:5],
                            "39" = 0,
                            sac2012_RouteC_reaches[6:8])


sac2013_RouteB_reaches <- c(sac2013_RouteB_reaches[1:11],
                            "39" = 0,
                            sac2013_RouteB_reaches[12:14])
sac2013_RouteC_reaches <- c(sac2013_RouteC_reaches[1:11],
                            "39" = 0,
                            sac2013_RouteC_reaches[12:14])


tars <- ls(pattern = "_reaches")

reach_length <- lapply(tars, function(x){
    #    browser()
    tmp <- get(x)
    idx <- match(makeID(names(tmp)), levs)
    structure(tmp[-1], names = makeID(names(tmp)))
})
reach_length <- do.call(c, reach_length)
reach_length <- reach_length[levs]
