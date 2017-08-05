# =====================================================================
# Auxiliary funtions for lattice graphics
#                                                        Eduardo Junior
#                                                    edujrrib@gmail.com
#                                                            2017-08-01
# =====================================================================

#-----------------------------------------------------------------------
# Lattice configuration
library(lattice)
library(latticeExtra)
library(grid)

## http://www.magesblog.com/2013/04/how-to-change-alpha-value-of-colours-in.html
add.alpha <- function(col, alpha = 1){
    apply(sapply(col, col2rgb)/255, 2,
          function(x) rgb(x[1], x[2], x[3], alpha = alpha))
}

## Define Colors
mycol <- colorRampPalette(c("gray80", "gray10"))(4)
myreg <- colorRampPalette(c("gray90", "gray20"))(100)

## Trellis graphical style.
ps1 <- list(
    superpose.symbol = list(
        col = mycol, pch = 1,
        fill = add.alpha(mycol, alpha = 0.4)),
    box.rectangle = list(col = 1, fill = c("gray70")),
    box.umbrella = list(col = 1, lty = 1),
    box.dot = list(pch = "|"),
    dot.symbol = list(col = 1, pch = 19),
    dot.line = list(col = "gray50", lty = 3),
    plot.symbol = list(col = 1),
    plot.line = list(col = 1),
    plot.polygon = list(col = "gray95"),
    superpose.line = list(col = mycol, lty = 1),
    superpose.polygon = list(col = mycol),
    strip.background = list(col = c("gray90", "gray70")),
    regions = list(col = myreg),
    par.sub.text = list(
        font = 1, just = "left", cex = 0.9,
        x = grid::unit(10, "mm"))
    )

# Space economic settings
ac <- list(pad1 = 0.5, pad2 = 0.5, tck = 0.5)
ps2 <- list(
    layout.widths = list(
        left.padding = 0.25,
        right.padding = -1,
        ylab.axis.padding = 0),
    layout.heights = list(
        bottom.padding = 0.25,
        top.padding = 0,
        axis.xlab.padding = 0,
        xlab.top = 0),
    axis.components = list(
        bottom = ac, top = ac,
        left = ac, right = ac)
)

trellis.par.set(lattice:::updateList(ps1, ps2), warn = FALSE)

#=======================================================================
# Lattice function for graphics

#-------------------------------------------
# Plot profile in lattice
xyprofile <- function(prof, conf = c(0.9, 0.95, 0.99),
                      namestrip = NULL, subset = 4,
                      scales.x = "free", ...) {
    #-------------------------------------------
    conf <- conf[order(conf, decreasing = TRUE)]
    da <- subset(prof, abs(z) <= subset)
    #-------------------------------------------
    fl <- levels(da$param)
    if (!is.null(namestrip)) {
        fl <- namestrip
    }
    xyplot(abs(z) ~ focal | param,
           data = da,
           layout = c(NA, 1),
           # ylab = expression(abs(z)~~(sqrt(~Delta~"deviance"))),
           ylab = expression(sqrt(2~(L(hat(phi))-L(phi)))),
           scales = list(x = scales.x),
           type = c("l", "g"),
           strip = strip.custom(
               factor.levels = fl),
           panel = function(x, y, subscripts, ...) {
               conf <- c(0.9, 0.95, 0.99)
               hl <- sqrt(qchisq(conf, 1))
               #-------------------------------------------
               mle <- x[y == 0]
               xl <- x[x < mle]; yl <- y[x < mle]
               xr <- x[x > mle]; yr <- y[x > mle]
               #-------------------------------------------
               funleft <- approxfun(x = yl, y = xl)
               funright <- approxfun(x = yr, y = xr)
               vxl <- funleft(hl)
               vxr <- funright(hl)
               vz <- c(hl, hl)
               ##-------------------------------------------
               panel.xyplot(x, y, ...)
               panel.arrows(c(vxl, vxr), 0, c(vxl, vxr), vz,
                            code = 1, length = 0.1, lty = 2,
                            col = "gray40")
               panel.segments(vxl, vz, vxr, vz, lty = 2,
                              col = "gray40")
               panel.abline(h = 0, v = mle, lty = 3)
               panel.text(x = rep(mle, 2), y = vz + 0.1,
                          labels = paste(conf*100, "%"),
                          col = "gray20")
           }, ...)
}

#-------------------------------------------
# Panel for envelope bands
panel.cbH <- function(x, y, ly, uy,
                      subscripts, cty,
                      col.line = plot.line$col,
                      lwd = plot.line$lwd,
                      desloc = NULL,
                      fill = 1, alpha = 0.1, length = 0.05, ...) {
    plot.line <- trellis.par.get("plot.line")
    if (is.null(desloc)) {
        desloc <- rep(0, length(uy))
    }
    y <- as.numeric(y)
    x <- as.numeric(x)
    or <- order(x)
    ly <- as.numeric(ly[subscripts])
    uy <- as.numeric(uy[subscripts])
    xo <- x[or]
    yo <- y[or]
    lyo <- ly[or]
    uyo <- uy[or]
    desl <- desloc[subscripts]
    if (cty == "bands") {
        panel.polygon(c(xo, rev(xo)), c(lyo, rev(uyo)), col = fill,
                      alpha = alpha, border = NA)
        panel.lines(xo, lyo, lty = list(...)$lty, lwd = 0.5,
                    col = col.line)
        panel.lines(xo, uyo, lty = list(...)$lty, lwd = 0.5,
                    col = col.line)
    }
    if (cty == "bars") {
        panel.arrows(xo + desl, lyo, xo + desl, uyo, length = length,
                     code = 3, angle = 90, col = col.line, lwd = lwd)
    }
    panel.xyplot(x + desl, y, subscripts = subscripts,
                 col.line = col.line, lwd = lwd, ...)
}

prepanel.cbH <- function(y, ly, uy, subscripts) {
    ly <- as.numeric(ly[subscripts])
    uy <- as.numeric(uy[subscripts])
    y <- as.numeric(y[subscripts])
    list(ylim = range(y, uy, ly, finite = TRUE))
}

#-------------------------------------------
# Panel for boxplot with `whisker.width` argument
my.panel.bwplot <- function(x, y, box.ratio = 1, box.width = box.ratio/(1 +
    box.ratio), horizontal = TRUE, pch = box.dot$pch, col = box.dot$col,
    alpha = box.dot$alpha, cex = box.dot$cex, font = box.dot$font,
    fontfamily = box.dot$fontfamily, fontface = box.dot$fontface,
    fill = box.rectangle$fill, varwidth = FALSE, notch = FALSE,
    notch.frac = 0.5, ..., levels.fos = if (horizontal) sort(unique(y)) else sort(unique(x)),
    stats = boxplot.stats, coef = 1.5, do.out = TRUE,
    identifier = "bwplot", whisker.width) {
    if (all(is.na(x) | is.na(y)))
        return()
    x <- as.numeric(x)
    y <- as.numeric(y)
    box.dot <- trellis.par.get("box.dot")
    box.rectangle <- trellis.par.get("box.rectangle")
    box.umbrella <- trellis.par.get("box.umbrella")
    plot.symbol <- trellis.par.get("plot.symbol")
    fontsize.points <- trellis.par.get("fontsize")$points
    if (!notch)
        notch.frac <- 0
    if (horizontal) {
        blist <- tapply(x, factor(y, levels = levels.fos), stats,
            coef = coef, do.out = do.out)
        blist.stats <- t(sapply(blist, "[[", "stats"))
        blist.out <- lapply(blist, "[[", "out")
        blist.height <- box.width
        if (varwidth) {
            maxn <- max(table(y))
            blist.n <- sapply(blist, "[[", "n")
            blist.height <- sqrt(blist.n/maxn) * blist.height
        }
        blist.conf <- if (notch)
            t(sapply(blist, "[[", "conf")) else blist.stats[, c(2, 4), drop = FALSE]
        xbnd <- cbind(blist.stats[, 3], blist.conf[, 2], blist.stats[,
            4], blist.stats[, 4], blist.conf[, 2], blist.stats[,
            3], blist.conf[, 1], blist.stats[, 2], blist.stats[,
            2], blist.conf[, 1], blist.stats[, 3])
        ytop <- levels.fos + blist.height/2
        ybot <- levels.fos - blist.height/2
        ybnd <- cbind(ytop - notch.frac * blist.height/2, ytop,
            ytop, ybot, ybot, ybot + notch.frac * blist.height/2,
            ybot, ybot, ytop, ytop, ytop - notch.frac * blist.height/2)
        xs <- cbind(xbnd, NA_real_)
        ys <- cbind(ybnd, NA_real_)
        panel.polygon(t(xs), t(ys), lwd = box.rectangle$lwd,
            lty = box.rectangle$lty, col = fill, alpha = box.rectangle$alpha,
            border = box.rectangle$col, identifier = paste(identifier,
                "box", sep = "."))
        panel.segments(
            c(blist.stats[, 2], blist.stats[, 4]),
            rep(levels.fos, 2),
            c(blist.stats[, 1], blist.stats[, 5]),
            rep(levels.fos, 2),
            col = box.umbrella$col,
            alpha = box.umbrella$alpha, lwd = box.umbrella$lwd,
            lty = box.umbrella$lty, identifier = paste(identifier,
                "whisker", sep = "."))
        panel.segments(
            c(blist.stats[, 1], blist.stats[, 5]),
            levels.fos - blist.height * whisker.width/2,
            c(blist.stats[, 1], blist.stats[, 5]),
            levels.fos + blist.height * whisker.width/2,
            col = box.umbrella$col, alpha = box.umbrella$alpha,
            lwd = box.umbrella$lwd, lty = box.umbrella$lty,
            identifier = paste(identifier,
                "cap", sep = "."))
        if (all(pch == "|")) {
            mult <- if (notch)
                1 - notch.frac else 1
            panel.segments(blist.stats[, 3], levels.fos - mult *
                blist.height/2, blist.stats[, 3], levels.fos +
                mult * blist.height/2, lwd = box.rectangle$lwd,
                lty = box.rectangle$lty, col = box.rectangle$col,
                alpha = alpha, identifier = paste(identifier,
                  "dot", sep = "."))
        } else {
            panel.points(x = blist.stats[, 3], y = levels.fos,
                pch = pch, col = col, alpha = alpha, cex = cex,
                fontfamily = fontfamily, fontface = lattice:::chooseFace(fontface,
                  font), fontsize = fontsize.points, identifier = paste(identifier,
                  "dot", sep = "."))
        }
        panel.points(x = unlist(blist.out), y = rep(levels.fos,
            sapply(blist.out, length)), pch = plot.symbol$pch,
            col = plot.symbol$col, alpha = plot.symbol$alpha,
            cex = plot.symbol$cex, fontfamily = plot.symbol$fontfamily,
            fontface = lattice:::chooseFace(plot.symbol$fontface, plot.symbol$font),
            fontsize = fontsize.points, identifier = paste(identifier,
                "outlier", sep = "."))
    } else {
        blist <- tapply(y, factor(x, levels = levels.fos), stats,
            coef = coef, do.out = do.out)
        blist.stats <- t(sapply(blist, "[[", "stats"))
        blist.out <- lapply(blist, "[[", "out")
        blist.height <- box.width
        if (varwidth) {
            maxn <- max(table(x))
            blist.n <- sapply(blist, "[[", "n")
            blist.height <- sqrt(blist.n/maxn) * blist.height
        }
        blist.conf <- if (notch)
            sapply(blist, "[[", "conf") else t(blist.stats[, c(2, 4), drop = FALSE])
        ybnd <- cbind(blist.stats[, 3], blist.conf[2, ], blist.stats[,
            4], blist.stats[, 4], blist.conf[2, ], blist.stats[,
            3], blist.conf[1, ], blist.stats[, 2], blist.stats[,
            2], blist.conf[1, ], blist.stats[, 3])
        xleft <- levels.fos - blist.height/2
        xright <- levels.fos + blist.height/2
        xbnd <- cbind(xleft + notch.frac * blist.height/2, xleft,
            xleft, xright, xright, xright - notch.frac * blist.height/2,
            xright, xright, xleft, xleft, xleft + notch.frac *
                blist.height/2)
        xs <- cbind(xbnd, NA_real_)
        ys <- cbind(ybnd, NA_real_)
        panel.polygon(t(xs), t(ys), lwd = box.rectangle$lwd,
            lty = box.rectangle$lty, col = fill, alpha = box.rectangle$alpha,
            border = box.rectangle$col, identifier = paste(identifier,
                "box", sep = "."))
        panel.segments(rep(levels.fos, 2), c(blist.stats[, 2],
            blist.stats[, 4]), rep(levels.fos, 2), c(blist.stats[,
            1], blist.stats[, 5]), col = box.umbrella$col, alpha = box.umbrella$alpha,
            lwd = box.umbrella$lwd, lty = box.umbrella$lty, identifier = paste(identifier,
                "whisker", sep = "."))
        panel.segments(levels.fos - blist.height/2, c(blist.stats[,
            1], blist.stats[, 5]), levels.fos + blist.height/2,
            c(blist.stats[, 1], blist.stats[, 5]), col = box.umbrella$col,
            alpha = box.umbrella$alpha, lwd = box.umbrella$lwd,
            lty = box.umbrella$lty, identifier = paste(identifier,
                "cap", sep = "."))
        if (all(pch == "|")) {
            mult <- if (notch)
                1 - notch.frac else 1
            panel.segments(levels.fos - mult * blist.height/2,
                blist.stats[, 3], levels.fos + mult * blist.height/2,
                blist.stats[, 3], lwd = box.rectangle$lwd, lty = box.rectangle$lty,
                col = box.rectangle$col, alpha = alpha, identifier = paste(identifier,
                  "dot", sep = "."))
        } else {
            panel.points(x = levels.fos, y = blist.stats[, 3],
                pch = pch, col = col, alpha = alpha, cex = cex,
                fontfamily = fontfamily, fontface = lattice:::chooseFace(fontface,
                  font), fontsize = fontsize.points, identifier = paste(identifier,
                  "dot", sep = "."))
        }
        panel.points(x = rep(levels.fos, sapply(blist.out, length)),
            y = unlist(blist.out), pch = plot.symbol$pch, col = plot.symbol$col,
            alpha = plot.symbol$alpha, cex = plot.symbol$cex,
            fontfamily = plot.symbol$fontfamily, fontface = lattice:::chooseFace(plot.symbol$fontface,
                plot.symbol$font), fontsize = fontsize.points,
            identifier = paste(identifier, "outlier", sep = "."))
    }
}

#-------------------------------------------
# Panel for segplots divide by groups
centfac <- function (group, space = NULL) {
    stopifnot(is.factor(group))
    if (is.null(space)) {
        space <- 0.5/nlevels(group)
    }
    d <- 2 * ((as.integer(group) - 1)/(nlevels(group) - 1)) -
        1
    return(space * d)
}

panel.groups.segplot <- function (x, y, z, centers, groups, gap = NULL,
                                  data, subscripts,  ...) {
    if (!missing(data)) {
        data <- eval(data, envir = parent.frame())
        groups <- data[, deparse(substitute(groups))]
    }
    stopifnot(is.factor(groups))
    stopifnot(length(groups) == length(z))
    if (is.null(gap)) {
        gap <- 0.5/nlevels(groups)
    }
    d <- 2 * ((as.numeric(groups) - 1)/(nlevels(groups) - 1)) -
        1
    z <- as.numeric(z) + gap * d
    panel.segplot(x, y, z, centers = centers,
                  subscripts = subscripts, ...)
}
