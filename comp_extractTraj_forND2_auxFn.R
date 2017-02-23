require(data.table)
require(ggplot2)

myCheckDigits <- function(x){ grepl('^[-]?[0-9]+[.]?[0-9]*$' , x) }

myCheckLogical <- function(x){ grepl('^TRUE$|^FALSE$' , x) }

myConvertStringListToTypes <- function(in.l){ 
  # convert strings with digits to numeric
  # uses logical indexing: http://stackoverflow.com/questions/42207235/replace-list-elements-by-name-with-another-list
  loc.l = myCheckDigits(in.l)
  in.l[loc.l] = lapply(in.l[loc.l], as.numeric)
  
  # convert strings with TRUE/FALSE to logical
  loc.l = myCheckLogical(in.l)
  in.l[loc.l] = lapply(in.l[loc.l], as.logical)
  
  return(in.l) 
}

myTrajExtr = function(in.dt,
                      in.max.break = 1, 
                      in.met.series = 'Metadata_Series', 
                      in.met.t = 'Metadata_T', 
                      in.met.tracklabel = 'TrackObjects_Label') {
  
  loc.dt = copy(in.dt)
  
  # TrackObjects assigns the same label for different cells from IdentifyPrimaryObjects
  # The following aggregation makes sure there's a unique TrackObjects_Label
  # for every Site at every time point: it takes the mean intensity of duplicated labels.
  # Make sure it makes sense!!!
  # Roughly 10% of objects affected
  
  loc.dt.nuc.aggr = loc.dt[, .(Intensity_MeanIntensity_Ratio = mean(Intensity_MeanIntensity_Ratio)),  by = c(in.met.series, in.met.t, in.met.tracklabel)]
  
  # cells from different sites have the same TrackObjects_Label
  # make it unique accross the experiment and pad numbers with zeros
  loc.dt.nuc.aggr[, TrackObjects_Label_uni := paste(sprintf("%02d", get(in.met.series)), sprintf("%04d", get(in.met.tracklabel)), sep = "_")]
  
  ####
  ## Allow for single-breaks in the middle of the trajectory
  ## Single-breaks usually result from a single frame with lost focus
  ## no interpolation
  
  loc.n.tp = max(loc.dt.nuc.aggr[, in.met.t, with = FALSE])
  
  # Select only cells with at least 1st and last timepoint measured
  # Metadata_T starts at 0, thus cells with T=0 & T=last measured 
  # have max(Metadata_T) - min(Metadata_T) = n.tp
  
  loc.dt.tmp1 = loc.dt.nuc.aggr[, .(Metadata_T.span = max(get(in.met.t)) - min(get(in.met.t))), by = TrackObjects_Label_uni][Metadata_T.span == loc.n.tp]
  
  
  # With cells selected above,
  # select tcourses with at most 1 break.
  # Calculation based on differences between consecutive Metadata_T.
  # Metadata_T increases by 1, thus a single-frame break would result in a difference of 2.
  
  loc.dt.tmp2 = loc.dt.nuc.aggr[TrackObjects_Label_uni %in% loc.dt.tmp1$TrackObjects_Label_uni]
  loc.dt.tmp2[, Metadata_T.diff := c(NA, diff(get(in.met.t))), by = TrackObjects_Label_uni]
  loc.dt.tmp2 = loc.dt.tmp2[Metadata_T.diff <= in.max.break + 1]
  
  # Selected trajectories with at most 1-frame break
  loc.out = loc.dt.nuc.aggr[TrackObjects_Label_uni %in% loc.dt.tmp2$TrackObjects_Label_uni]
  return(loc.out)
}


# Returns original dt with an additional column with normalized quantity
# the column to be normalised is given by 'in.meas.col'
# Normalisation is based on part of the trajectory;
# this is defined by in.rt.in and max, and the column with time in.rt.col.
# Additional parameters:
# in.by.cols - character vector with 'by' columns to calculate normalisation per group
#              if NULL, no grouping is done
# in.robust - whether robust measures should be used (median instead of mean, mad instead of sd)
# in.type - type of normalization: z.score or mean (fold change w.r.t. mean)

myNorm = function(in.dt,
                  in.meas.col,
                  in.rt.col = 'RealTime',
                  in.rt.min = 10,
                  in.rt.max = 20,
                  in.by.cols = NULL,
                  in.robust = TRUE,
                  in.type = 'z.score') {
  loc.dt <- copy(in.dt) # copy so as not to alter original dt object w intermediate assignments
  
  if (is.null(in.by.cols)) {
    if (in.robust)
      loc.dt.pre.aggr = loc.dt[get(in.rt.col) > in.rt.min & get(in.rt.col) < in.rt.max, .(meas.md = median(get(in.meas.col)),
                                                                                          meas.mad = mad(get(in.meas.col)))]
    else
      loc.dt.pre.aggr = loc.dt[get(in.rt.col) > in.rt.min & get(in.rt.col) < in.rt.max, .(meas.md = mean(get(in.meas.col)),
                                                                                          meas.mad = sd(get(in.meas.col)))]
    
    loc.dt = cbind(loc.dt, loc.dt.pre.aggr)
  }  else {
    if (in.robust)
      loc.dt.pre.aggr = loc.dt[get(in.rt.col) > in.rt.min & get(in.rt.col) < in.rt.max, .(meas.md = median(get(in.meas.col)),
                                                                                          meas.mad = mad(get(in.meas.col))), by = in.by.cols]
    else
      loc.dt.pre.aggr = loc.dt[get(in.rt.col) > in.rt.min & get(in.rt.col) < in.rt.max, .(meas.md = mean(get(in.meas.col)),
                                                                                          meas.mad = sd(get(in.meas.col))), by = in.by.cols]
    
    loc.dt = merge(loc.dt, loc.dt.pre.aggr, by = in.by.cols)
  }
  
  
  if (in.type == 'z.score') {
    loc.dt[, meas.norm := (get(in.meas.col) - meas.md) / meas.mad]
  } else {
    loc.dt[, meas.norm := (get(in.meas.col) / meas.md)]
  }
  
  loc.dt[, c('meas.md', 'meas.mad') := NULL]
  return(loc.dt)
}


# add a rescaled column by shifting a selected column to a baseline (in.y.min - in.y.h) and setting height to in in.y.h
# the column to shift is in.meas.col
# the new column has the same name as the original but with a suffix in.meas.resc.col
myResc = function(in.dt,
                  in.meas.col,
                  in.meas.resc.col = '.resc',
                  in.y.min,
                  in.y.h) {
  
  loc.dt = copy(in.dt)
  loc.dt[, paste0(in.meas.col, in.meas.resc.col) := (in.y.min - in.y.h) +
           (in.y.min - (in.y.min - in.y.h)) * (get(in.meas.col) - min(get(in.meas.col))) /
           (max(get(in.meas.col)) - min(get(in.meas.col)))]
  
  return(loc.dt)
}


myGgplotTraj = function(dt.arg,
                        x.arg,
                        y.arg,
                        group.arg,
                        facet.arg = NULL,
                        summary.arg = FALSE,
                        facet.ncol.arg = 2,
                        line.col.arg = NULL,
                        xlab.arg = "Time",
                        ylab.arg = "Fl. int.",
                        plotlab.arg = "",
                        dt.stim.arg = NULL,
                        stim.x.arg,
                        stim.y.arg,
                        maxrt.arg = 60,
                        xaxisbreaks.arg = 10,
                        xlim.arg = NULL,
                        ylim.arg = NULL) {
  
  p.tmp = ggplot(dt.arg,
                 aes_string(x = x.arg,
                            y = y.arg,
                            group = group.arg))
  
  if (is.null(line.col.arg))
    p.tmp = p.tmp + geom_line(alpha = 0.25, size = 0.25)
  else
    p.tmp = p.tmp + geom_line(aes_string(colour = line.col.arg), alpha = 0.5, size = 0.5)
  
  
  if (summary.arg)
    p.tmp = p.tmp + 
      stat_summary(aes_string(y = y.arg, group = 1), 
                   geom="ribbon", 
                   fun.data=mean_cl_boot, 
                   colour='red', alpha=0.5,
                   group = 1) +
      stat_summary(
        aes_string(y = y.arg, group = 1),
        fun.y = mean,
        colour = 'red',
        linetype = 'solid',
        size = 1,
        geom = "line",
        group = 1
      )
  
  
  
  if (!is.null(facet.arg))
    p.tmp = p.tmp + 
      facet_wrap(as.formula(paste("~", facet.arg)),
                 ncol = facet.ncol.arg,
                 scales = "free_x", 
                 drop = FALSE)
  
  if(!is.null(dt.stim.arg)) {
    p.tmp = p.tmp + geom_line(data = dt.stim.arg,
                              aes_string(x = stim.x.arg, y = stim.y.arg),
                              colour = 'blue',
                              size = 1,
                              group = 1) 
  }
  
  if (!is.null(ylim.arg))
    p.tmp = p.tmp +
      coord_cartesian(ylim = ylim.arg)
  
  if (!is.null(xlim.arg))
    p.tmp = p.tmp +
      coord_cartesian(xlim = xlim.arg)
  
  p.tmp = p.tmp + 
    scale_x_continuous(breaks = seq(0, maxrt.arg, xaxisbreaks.arg)) +
    xlab(paste0(xlab.arg, "\n")) +
    ylab(paste0("\n", ylab.arg)) +
    ggtitle(plotlab.arg) +
    theme_bw(base_size = 18, base_family = "Helvetica") +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      axis.line.x = element_line(color = "black", size = 0.25),
      axis.line.y = element_line(color = "black", size = 0.25),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      strip.text.x = element_text(size = 14, face = "bold"),
      strip.text.y = element_text(size = 14, face = "bold"),
      strip.background = element_blank(),
      legend.key = element_blank(),
      legend.key.height = unit(1, "lines"),
      legend.key.width = unit(2, "lines"),
      legend.position = "top"
    )
  
  return(p.tmp)
}

myGgplotTrajRibbon = function(dt.arg,
                              x.arg,
                              y.arg,
                              group.arg, 
                              xlab.arg = "Time",
                              ylab.arg = "Fl. int.",
                              plotlab.arg = "") {
  p.tmp = ggplot(dt.arg, aes_string(x = x.arg, group = group.arg)) +
    geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = 'grey70', alpha = 0.5) +
    geom_line(aes_string(y = y.arg, colour = group.arg)) +
    scale_color_manual(values = rhg_cols, name = "") +
    xlab(paste0(xlab.arg, "\n")) +
    ylab(paste0("\n", ylab.arg)) +
    ggtitle(plotlab.arg) +
    theme_bw(base_size = 18, base_family = "Helvetica") +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      axis.line.x = element_line(color = "black", size = 0.25),
      axis.line.y = element_line(color = "black", size = 0.25),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      strip.text.x = element_text(size = 14, face = "bold"),
      strip.text.y = element_text(size = 14, face = "bold"),
      strip.background = element_blank(),
      legend.key = element_blank(),
      legend.key.height = unit(1, "lines"),
      legend.key.width = unit(2, "lines"),
      legend.position = "top"
    ) 
  return(p.tmp)
}
