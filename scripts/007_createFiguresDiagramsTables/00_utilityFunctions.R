
# - Quick modification to ggsn::scalebar ----
# define distance in "m"

scalebar2 = function (data = NULL, location = "bottomright", dist = NULL, 
          dist_unit = NULL, transform = NULL, dd2km = NULL, model = NULL, 
          height = 0.02, st.dist = 0.02, st.bottom = TRUE, st.size = 5, 
          st.color = "black", box.fill = c("black", "white"), box.color = "black", 
          border.size = 1, x.min = NULL, x.max = NULL, y.min = NULL, 
          y.max = NULL, anchor = NULL, facet.var = NULL, facet.lev = NULL, 
          ...) 
{
  if (is.null(data)) {
    if (is.null(x.min) | is.null(x.max) | is.null(y.min) | 
        is.null(y.max)) {
      stop("If data is not defined, x.min, x.max, y.min and y.max must be.")
    }
    data <- data.frame(long = c(x.min, x.max), lat = c(y.min, 
                                                       y.max))
  }
  if (is.null(transform)) {
    stop("transform should be logical.")
  }
  if (any(class(data) %in% "sf")) {
    xmin <- sf::st_bbox(data)["xmin"]
    xmax <- sf::st_bbox(data)["xmax"]
    ymin <- sf::st_bbox(data)["ymin"]
    ymax <- sf::st_bbox(data)["ymax"]
  }
  else {
    xmin <- min(data$long)
    xmax <- max(data$long)
    ymin <- min(data$lat)
    ymax <- max(data$lat)
  }
  if (location == "bottomleft") {
    if (is.null(anchor)) {
      x <- xmin
      y <- ymin
    }
    else {
      x <- as.numeric(anchor["x"])
      y <- as.numeric(anchor["y"])
    }
    direction <- 1
  }
  if (location == "bottomright") {
    if (is.null(anchor)) {
      x <- xmax
      y <- ymin
    }
    else {
      x <- as.numeric(anchor["x"])
      y <- as.numeric(anchor["y"])
    }
    direction <- -1
  }
  if (location == "topleft") {
    if (is.null(anchor)) {
      x <- xmin
      y <- ymax
    }
    else {
      x <- as.numeric(anchor["x"])
      y <- as.numeric(anchor["y"])
    }
    direction <- 1
  }
  if (location == "topright") {
    if (is.null(anchor)) {
      x <- xmax
      y <- ymax
    }
    else {
      x <- as.numeric(anchor["x"])
      y <- as.numeric(anchor["y"])
    }
    direction <- -1
  }
  if (!st.bottom) {
    st.dist <- y + (ymax - ymin) * (height + st.dist)
  }
  else {
    st.dist <- y - (ymax - ymin) * st.dist
  }
  height <- y + (ymax - ymin) * height
  if (dist_unit == "m") {
    dist <- dist/1000
    dist_unit0 <- "m"
    dist_unit <- "km"
  }
  if (!is.null(dd2km)) {
    if (dd2km) {
      transform <- TRUE
    }
    cat("dd2km is deprecated. Use ggsn::transform instead.")
  }
  if (transform) {
    break1 <- maptools::gcDestination(lon = x, lat = y, bearing = 90 * 
                                        direction, dist = dist, dist.units = dist_unit, model = model)[1, 
                                                                                                       1]
    break2 <- maptools::gcDestination(lon = x, lat = y, bearing = 90 * 
                                        direction, dist = dist * 2, dist.units = dist_unit, 
                                      model = model)[1, 1]
  }
  else {
    if (location == "bottomleft" | location == "topleft") {
      if (exists("dist_unit0") | (!exists("dist_unit0") & 
                                  dist_unit == "km")) {
        break1 <- x + dist * 1000
        break2 <- x + dist * 2000
      }
      else if (dist_unit == "m") {
        break1 <- x + dist * 1
        break2 <- x + dist * 2
      }
      else if (dist_unit == "nm") {
        break1 <- x + dist * 1852
        break2 <- x + dist * 1852 * 2
      }
      else if (dist_unit == "mi") {
        break1 <- x + dist * 1609.34
        break2 <- x + dist * 1609.34 * 2
      }
      else {
        break1 <- x + dist
        break2 <- x + dist
      }
    }
    else {
      if (exists("dist_unit0") | (!exists("dist_unit0") & 
                                  dist_unit == "km")) {
        break1 <- x - dist * 1000
        break2 <- x - dist * 2000
      }
      else if (dist_unit == "m") {
        break1 <- x + dist * 1
        break2 <- x + dist * 2
      }
      else if (dist_unit == "nm") {
        break1 <- x - dist * 1852
        break2 <- x - dist * 1852 * 2
      }
      else if (dist_unit == "mi") {
        break1 <- x - dist * 1609.34
        break2 <- x - dist * 1609.34 * 2
      }
      else {
        break1 <- x - dist
        break2 <- x - dist
      }
    }
  }
  box1 <- data.frame(x = c(x, x, rep(break1, 2), x), y = c(y, 
                                                           height, height, y, y), group = 1)
  box2 <- data.frame(x = c(rep(break1, 2), rep(break2, 2), 
                           break1), y = c(y, rep(height, 2), y, y), group = 1)
  if (!is.null(facet.var) & !is.null(facet.lev)) {
    for (i in 1:length(facet.var)) {
      if (any(class(data) == "sf")) {
        if (!is.factor(data[, facet.var[i]][[1]])) {
          data[, facet.var[i]] <- factor(data[, facet.var[i]][[1]])
        }
        box1[, facet.var[i]] <- factor(facet.lev[i], 
                                       levels(data[, facet.var[i]][[1]]))
        box2[, facet.var[i]] <- factor(facet.lev[i], 
                                       levels(data[, facet.var[i]][[1]]))
      }
      else {
        if (!is.factor(data[, facet.var[i]])) {
          data[, facet.var[i]] <- factor(data[, facet.var[i]])
        }
        box1[, facet.var[i]] <- factor(facet.lev[i], 
                                       levels(data[, facet.var[i]]))
        box2[, facet.var[i]] <- factor(facet.lev[i], 
                                       levels(data[, facet.var[i]]))
      }
    }
  }
  if (exists("dist_unit0")) {
    legend <- cbind(text = c(0, dist * 1000, dist * 2000), 
                    row.names = NULL)
  }
  else {
    legend <- cbind(text = c(0, dist, dist * 2), row.names = NULL)
  }
  gg.box1 <- geom_polygon(data = box1, aes(x, y), fill = utils::tail(box.fill, 
                                                                     1), color = utils::tail(box.color, 1), size = border.size)
  gg.box2 <- geom_polygon(data = box2, aes(x, y), fill = box.fill[1], 
                          color = box.color[1], size = border.size)
  x.st.pos <- c(box1[c(1, 3), 1], box2[3, 1])
  if (location == "bottomright" | location == "topright") {
    x.st.pos <- rev(x.st.pos)
  }
  label <- NULL
  if (exists("dist_unit0")) {
    legend2 <- cbind(data[1:3, ], x = unname(x.st.pos), y = unname(st.dist), 
                     label = paste0(legend[, "text"], c("", "", "km")))
  }
  else {
    legend2 <- cbind(data[1:3, ], x = unname(x.st.pos), y = unname(st.dist), 
                     label = paste0(legend[, "text"], c("", "", dist_unit)))
  }
  if (!is.null(facet.var) & !is.null(facet.lev)) {
    for (i in 1:length(facet.var)) {
      if (any(class(data) == "sf")) {
        legend2[, facet.var[i]] <- factor(facet.lev[i], 
                                          levels(data[, facet.var[i]][[1]]))
      }
      else {
        legend2[, facet.var[i]] <- factor(facet.lev[i], 
                                          levels(data[, facet.var[i]]))
      }
    }
  }
  else if (!is.null(facet.var) & is.null(facet.lev)) {
    facet.levels0 <- unique(as.data.frame(data)[, facet.var])
    facet.levels <- unlist(unique(as.data.frame(data)[, facet.var]))
    legend2 <- do.call("rbind", replicate(length(facet.levels), 
                                          legend2, simplify = FALSE))
    if (length(facet.var) > 1) {
      facet.levels0 <- expand.grid(facet.levels0)
      legend2[, facet.var] <- facet.levels0[rep(row.names(facet.levels0), 
                                                each = 3), ]
    }
    else {
      legend2[, facet.var] <- rep(facet.levels0, each = 3)
    }
  }
  gg.legend <- geom_text(data = legend2, aes(x, y, label = label), 
                         size = st.size, color = st.color, ...)
  return(list(gg.box1, gg.box2, gg.legend))
}
