# R code from paper: Interactive Data Visualization and Exploration Using the
# Loon R Package
#
# Author: Adrian Waddell
# Date: October 10, 2016
#
# Note:
# The R markdown and PDF paper that put this code into conctext can be found
# in the git repository here:
# https://github.com/waddella/phuse2016_adverse_events


# Generating Adverse Events Data

set.seed(1)

aeterms <- c(
  'HEADACHE', 'CHRONIC BACK PAIN', 'NOSE BLEEDING RIGHT NOSTRIL',
  'PROBLEMS OF HYPOTENSION', 'LOOSE STOOL', 'ABDOMINAL DISCOMFORT',
  'DIARRHEA', 'ABDOMINAL FULLNESS DUE TO GAS', 'NAUSEA (INTERMITTENT)',
  'WEAKNESS', 'HYPOTENSIVE'
)
aesevns <- c(1, 2, 1, 1, 3, 2, 2, 1, 1, 1, 3)

normalize <- function(x)x/sum(x)
weightsA <- normalize(dlnorm(seq(0, 5, length.out = 25), meanlog = 3))
weightsB <- normalize(dlnorm(seq(0, 5, length.out = 25)))

l.aae <- Map(function(id, ARM) {

  TRTSDT <- as.Date("2016-01-01", "%Y-%m-%d") + sample(0:365, 1, replace = TRUE)
  TRTEDT <- TRTSDT + 200 + rbinom(1, 60, 0.5)
  DISCDEAT <- sample(c(TRUE, FALSE), 1, prob=if(ARM=="ARM A") c(.3,.7) else c(.15,.85))

  n_ae <- sample(1:25, 1, prob=if(ARM == "ARM A") weightsA else weightsB)
  i <- sample(1:length(aeterms), size=n_ae, replace=TRUE, prob=c(6,rep(1,10))/16)
  ASTDT <- sample(seq(TRTSDT, TRTEDT-1, by=1), n_ae, replace = TRUE)
  ADURN <- sample(1:18, size=n_ae, replace=TRUE)
  AENDT <- ASTDT + ADURN
  ii <- order(ASTDT)
  if(DISCDEAT & any(AENDT>TRTEDT)) {
    AENDT[AENDT>TRTEDT] <- TRTEDT
    ADURN <- as.numeric(AENDT - ASTDT)
  }

  list(
    USUBJID = id,
    SEX = sample(c('F', 'M'), 1),
    AGE = 20 + rbinom(1, size=40, prob=0.7),
    ARM = ARM, DISCDEAT = DISCDEAT,
    TRTSDT = TRTSDT, TRTEDT = TRTEDT,
    aes = list(
      AESEQ = 1:n_ae, AETERM = aeterms[i],
      AESEVN = aesevns[i], ASTDT = ASTDT[ii],
      AENDT = AENDT[ii], ADURN = ADURN[ii]
    )
  )
}, seq(1, 300, by=1), sample(rep(c('ARM A', 'ARM B'), 150), replace = FALSE))

aae <- Reduce(rbind, Map(as.data.frame, l.aae))
names(aae) <- gsub("aes.", "", names(aae), fixed = TRUE)

head(aae, 3)

View(aae)


## Visualize this adverse events data with loon

library(loon)

age <- sapply(l.aae, function(x)x$AGE)
naes <- sapply(l.aae, function(x)length(x$aes$AESEQ))

p <- l_plot(x=age, y=naes, ylabel="Number of Adverse Events")

# this can be done on the inspector directly
l_move_jitter(p, "all")
l_scaleto_world(p)
p['glyph'] <- 'ccircle'

## Encode information
p['glyph'] <- ifelse(sapply(l.aae, function(x)x$SEX)=='F', 'ccircle', 'csquare')
p['color'] <- ifelse(sapply(l.aae, function(x)x$ARM)=='ARM A', 'tan', 'steelblue')


## ------------------------------------------------------------------------
t.label <- unlist(Map(function(x) {
  t.x <- table(x$aes$AETERM)
  paste(c(
    paste0('Patient ', x$USUBJID,':'),
    apply(cbind(t.x, names(t.x)), 1, function(x)paste(x, collapse = ' '))
  ), collapse = '\n')
}, l.aae))

l_configure(p, itemLabel=t.label, showItemLabels=TRUE)


## ------------------------------------------------------------------------
h <- l_hist(
  x = sapply(l.aae, function(x)sum(x$aes$AETERM == 'HEADACHE')),
  xlabel = 'Number of Headaches per Patient',
  yshows = 'density',
  showScales = TRUE,
  binwidth = 1
)


## this can be done on the inspectory directly
p['linkingGroup'] <- "aes"
l_configure(h, linkingGroup="aes", sync="pull", showStackedColors=TRUE)
p['selected'] <- naes > 15 & age > 46



## ------------------------------------------------------------------------
createAEplot <- function() {
  pae <- l_plot(showItemLabels=TRUE, xlabel="Treatment Relative Day", showScales=TRUE)

  rectHeight <- 4
  y <- 0

  scale01 <- function(x) {
    dx <- diff(range(x))
    if (dx == 0) rep(0, length(x)) else (x-min(x))/dx
  }

  draw_patient <- function(x) {
    patient_label <- paste("Patient", x$USUBJID)
    g <- l_layer_group(pae, label=patient_label)
    l_layer_rectangle(
      pae, parent=g,
      x = c(1, x$TRTEDT - x$TRTSDT + 1), y = c(y, y+rectHeight),
      color = if(x$DISCDEAT) "lemonchiffon1" else "gray80",
      linecolor = "",
      itemLabel = paste("Treatment Period for Patient", x$USUBJID)
    )
    l_layer_text(pae, parent=g, text=patient_label, x=1, y=y+rectHeight,
                 justify='left', anchor='nw', color="black")

    if (length(x$aes$AESEQ)>0) {
      xcoords <- Map(function(t0, t1) as.numeric(c(t0, t1)-x$TRTSDT+1),
                     x$aes$ASTDT, x$aes$AENDT)
      tE <- list()
      ycoords <- Map(function(j)c(j, j), scale01(unlist(Map(function(xc) {
        k <- vapply(tE, function(tEi) xc[2]>tEi, logical(1))
        i <- if (any(k)) which(k)[1] else length(tE) + 1
        tE[[i]] <<- xc[2]
        i
      }, xcoords))) * rectHeight + y)
      col <- ifelse(x$aes$AESEVN == 3, "orangered", "dodgerblue4")
      if (length(xcoords) == 1)
        l_layer_line(pae, parent=g, x=xcoords[[1]], y=ycoords[[1]],
                     itemLabel=x$aes$AETERM, linewidth=2, color=col)
      else
        l_layer_lines(pae, parent=g, x=xcoords, y=ycoords,
                      itemLabel=x$aes$AETERM, linewidth=2, color=col)
    }
    y <<- y + rectHeight + 3
    l_scaleto_world(pae)
  }
  list(
    updateWith = function(selected) {
      y <<- 0
      for (layer in l_layer_getChildren(pae, "root"))
        if (layer != "model") l_layer_expunge(pae, layer)
      if (sum(selected)>0)
        Map(function(x)draw_patient(x), l.aae[selected])
    },
    widget = pae
  )
}



showAEs <- createAEplot()

l_bind_state(p, "selected", function() {
  showAEs$updateWith(p['selected'])
})

## to create the plot at the end of the paper
p['selected'] <- naes > 23 & sapply(l.aae, function(x)x$ARM) == 'ARM A' &
  sapply(l.aae, function(x)x$SEX) == 'M'
