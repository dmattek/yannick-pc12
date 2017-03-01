####
## Analyzes Image.csv & Nuclei.csv
## Extracts single-cell trajectories

# Execute in a directory where output.mer is located

# loads required libraries
require(xlsx)
require(gridExtra)
require(Hmisc) # for smean.cl.boot
require(RCurl)


auxScript = getURL("https://raw.githubusercontent.com/dmattek/yannick-pc12/master/comp_extractTraj_forND2_auxFn.R", ssl.verifypeer = FALSE)
eval(parse(text = auxScript))


####
## Main parameters

# names of two files with paramaters for analysis and experiment
# these files should be located one folder up of cp.out from which the script is executed
s.par.plot = 'plotFormat.xlsx'
s.par.exp  = 'experimentDescription.xlsx'

####
## Read parameters from the file

df.par = read.xlsx(
  paste0('../', s.par.plot),
  sheetIndex = 1,
  header = FALSE,
  as.data.frame = TRUE,
  colClasses = rep("character", 2),
  stringsAsFactors = FALSE
)

# convert data frame with parameters to a list 
l.par = split(df.par[, 2], df.par[, 1])

# convert strings with digits to numeric and strings with TRUE/FALSE to logical
l.par = myConvertStringListToTypes(l.par)


#####
## Read experiment description
dt.exp = as.data.table(read.xlsx(paste0('../', s.par.exp), 
                                 sheetIndex = 1, 
                                 header = TRUE, 
                                 as.data.frame = TRUE, 
                                 stringsAsFactors=FALSE))

# sometimes an NA column appears at the end; remove
dt.exp = dt.exp[, names(dt.exp)[!(names(dt.exp) %like% 'NA')], with = FALSE]

# sometimes an NA row appears at the end; remove
dt.exp = dt.exp[dt.exp[, !Reduce(`&`, lapply(.SD, is.na))]]

# create new column with experimental condition for each site and channel
dt.exp[, Stim_All_S  := sprintf('S%02d: %s %s %s',  Position, Stim_Treat, Stim_Conc, Stim_Duration)]
dt.exp[, Stim_All_Ch := sprintf('Ch%02d: %s %s %s',  Channel,  Stim_Treat, Stim_Conc, Stim_Duration)]

# Paramater: number of channels
c.n.devch = length(unique((dt.exp$Channel)))

# Parameter: number of sites per channel
c.n.SperCh = length(unique((dt.exp$Position))) / c.n.devch

# Acquisition frequency (minutes)
# Used to assign real time
c.acq.freq =  median(dt.exp$Acquisition_frequency_min)

# Time before which a 5-time point average is calculated for normalization of time courses
c.t.norm = median(dt.exp$Equilibration_min)



####
## Read data

## Load Nuclei files using custom file reading function
dt.nuc = fread(paste(l.par$dir.out, l.par$files.nuc, sep = '/'))

# load Image files
dt.im = fread(paste(l.par$dir.out, l.par$files.img, sep = '/'))



## assignment of columns:
## s.met.site
## s.met.time
## s.tracklabel
## s.meas


# deciding whether FOV are in Image_Metadata_Site or Image_Metadata_Series
# assign to Image_Metadata_Site or Image_Metadata_Series
b.error = FALSE
if (length(names(dt.nuc)[names(dt.nuc) %like% 'Site']) == 1)
  s.met.site = names(dt.nuc)[names(dt.nuc) %like% 'Site'] else 
    if (length(names(dt.nuc)[names(dt.nuc) %like% 'Series']) == 1)
      s.met.site = names(dt.nuc)[names(dt.nuc) %like% 'Series'] else {
        print("Error extracting Site/Series metadata!")
        b.error = TRUE
      }

stopifnot(!b.error)

# assign to Image_Metadata_T or Image_Metadata_Time
s.met.time = names(dt.nuc)[names(dt.nuc) %like% 'Metadata_T']

# assign to TrackObjects_Label
s.met.trackabel = names(dt.nuc)[names(dt.nuc) %like% 'TrackObjects_Label']

# assign measurement of nuclei
s.meas.nuc = names(dt.nuc)[names(dt.nuc) %like% 'Mean.*Ratio']

# assign measurement of the image
s.meas.img = names(dt.im)[names(dt.im) %like% 'MeanIntensity']

# experiment name column (used in plot titles)
s.exp.name.col = names(dt.im)[names(dt.im) %like% 'Exp']

# if more than 1 experiments in a file, use the 1st name only
s.exp.name = unique(dt.im[, get(s.exp.name.col)])
if(length(s.exp.name) > 1)
  s.exp.name = s.exp.name[[1]]




# Create a table with real time assigned to to Metadata_T
# Assumption: frequency of data acquisition = c.acq.freq
# Metadata_T starts at 0

dt.tass = data.table(Metadata_T = unique(dt.nuc[, get(s.met.time)]), 
                     RealTime = seq(0, max(dt.nuc[, get(s.met.time)])*c.acq.freq, c.acq.freq))


## Nuclei
# obtain data table with cells with trajectories spanning entire experiment
dt.nuc.sel = myTrajExtr(in.dt = dt.nuc, 
                        in.max.break = 1, 
                        in.met.series = s.met.site, 
                        in.met.t = s.met.time, 
                        in.met.tracklabel = s.met.trackabel, 
                        in.aggr.cols = 'Intensity_MeanIntensity_Ratio')

# add real time column
dt.nuc.sel = merge(dt.nuc.sel, dt.tass)

# add experimental description to Nuclei
dt.nuc.sel = merge(dt.nuc.sel, dt.exp[, .(Position, Stim_All_Ch, Stim_All_S)], by.x = s.met.site, by.y = 'Position')


## Image
# add Real time to Image data
dt.im.mer = merge(dt.im, dt.tass, by = s.met.time)

# add experimental description to Image
dt.im.mer = merge(dt.im.mer, dt.exp[, .(Position, Stim_All_Ch, Stim_All_S)], by.x = s.met.site, by.y = 'Position')

## Rescale pulse channel between c.plot.y.min & c.plot.y.max
dt.im.mer = myResc(
  in.dt = dt.im.mer,
  in.meas.col = s.meas.img,
  in.y.min = as.numeric(l.par$plot.y.min),
  in.y.h = as.numeric(l.par$plot.stim.y.height)
)


# create pulse averages per channel
dt.im.mer.aggr = dt.im.mer[, .(mean.pulse.resc = mean(get(paste0(s.meas.img, '.resc')))), by = .(Stim_All_Ch, RealTime)]


# number of cells per site
dt.ncells.persite = setorder(dt.nuc.sel[get(s.met.time) == 1, .N, by = Stim_All_S], Stim_All_S)
dt.ncells.perch   = setorder(dt.nuc.sel[get(s.met.time) == 1, .N, by = Stim_All_Ch], Stim_All_Ch)


# creating output in wide format for Jan
l.nuc.sel.wide = lapply(split(dt.nuc.sel, by = 'Stim_All_Ch'), function(x)  {
  dcast(data = x, formula = RealTime ~ TrackObjects_Label_uni, value.var = s.meas.nuc)
})

l.nuc.sel.wide.pulse = lapply(names(l.nuc.sel.wide), function(x) {loc.dt = merge(l.nuc.sel.wide[[x]], dt.im.mer.aggr[Stim_All_Ch %in% x, c('RealTime', 'mean.pulse.resc'), with = FALSE], by = 'RealTime')})
names(l.nuc.sel.wide.pulse) = names(l.nuc.sel.wide)


#####
## Plotting 

p.out = list()

# Plot: raw trajectories per site
# p.out$imRatio_corrIllum_notNorm_perSite = myGgplotTraj(dt.nuc.sel, x.arg = 'RealTime', y.arg = s.meas.nuc, group.arg = paste0(s.met.trackabel, '_uni'),
#              ylim.arg = c(c.plot.y.min, c.plot.y.max),
#              xaxisbreaks.arg = c.plot.x.interv,
#              maxrt.arg =  max(dt.tass$RealTime),
#              summary.arg = TRUE,
#              facet.arg = 'Stim_All_S', facet.ncol.arg = c.plot.facets.ncol.s,
#              xlab.arg = "\nTime [min]", ylab.arg = "Fluorescence ratio (x1000; from Illum.Corr)\n", plotlab.arg = paste0('Experiment: ', s.exp.name))



# Plot: raw dextrin per site
# p.out$imMnFlInt_dextrin_perSite = myGgplotTraj(dt.im.mer, x.arg = 'RealTime', y.arg = s.meas.img, group.arg = s.met.site, 
#                                                        xaxisbreaks.arg = c.plot.x.interv, 
#                                                        maxrt.arg =  max(dt.tass$RealTime), 
#                                                        summary.arg = TRUE,
#                                                        facet.arg = 'Stim_All_S', facet.ncol.arg = c.plot.facets.ncol.s, 
#                                                        xlab.arg = "\nTime [min]", ylab.arg = "Mean fl.int.; whole image dextrin\n", plotlab.arg = paste0('Experiment: ', s.exp.name))



# Plot: raw trajectories per site with dextrin
p.out$imRatio_wDextrin_notNorm_perSite = myGgplotTraj(dt.arg = dt.nuc.sel, x.arg = 'RealTime', y.arg = s.meas.nuc, group.arg = paste0(s.met.trackabel, '_uni'), 
                                                      dt.stim.arg = dt.im.mer, stim.x.arg = 'RealTime', stim.y.arg = paste0(s.meas.img, '.resc'),
                                                      ylim.arg = c(l.par$plot.y.min - l.par$plot.stim.y.height, l.par$plot.y.max), 
                                                      xaxisbreaks.arg = l.par$plot.x.interv,
                                                      maxrt.arg =  max(dt.tass$RealTime), 
                                                      summary.arg = TRUE,
                                                      facet.arg = 'Stim_All_S', facet.ncol.arg = l.par$plot.facets.ncol.site, 
                                                      xlab.arg = "\nTime [min]", ylab.arg = "Fluorescence ratio (x1000; from Illum.Corr)\n", plotlab.arg = paste0('Experiment: ', s.exp.name))




# Plot: raw trajectories with dextrin per channel
p.out$imRatio_wDextrin_notNorm_perCh = myGgplotTraj(dt.arg = dt.nuc.sel, x.arg = 'RealTime', y.arg = s.meas.nuc, group.arg = paste0(s.met.trackabel, '_uni'), 
                                                    dt.stim.arg = dt.im.mer.aggr, stim.x.arg = 'RealTime', stim.y.arg = 'mean.pulse.resc',
                                                    ylim.arg = c(l.par$plot.y.min - l.par$plot.stim.y.height, l.par$plot.y.max), 
                                                    xaxisbreaks.arg = l.par$plot.x.interv, 
                                                    maxrt.arg =  max(dt.tass$RealTime), 
                                                    summary.arg = TRUE,
                                                    facet.arg = 'Stim_All_Ch', facet.ncol.arg = l.par$plot.facets.ncol.channel, 
                                                    xlab.arg = "\nTime [min]", ylab.arg = "Fluorescence ratio (x1000; from Illum.Corr)\n", plotlab.arg = paste0('Experiment: ', s.exp.name))








## Population averages per channel
dt.nuc.sel.mean = dt.nuc.sel[, as.list(smean.cl.boot(get(s.meas.nuc), B = 1000)), by = .(Stim_All_Ch, RealTime)]

p.out$imRatio_popMeanCI_notNorm_perCh = myGgplotTrajRibbon(dt.nuc.sel.mean, 
                   x.arg = "RealTime", 
                   y.arg = "Mean", 
                   group.arg = "Stim_All_Ch", 
                   xlab.arg = "\nTime [min]",
                   ylab.arg = "Mean\n",
                   plotlab.arg = paste0('Experiment: ', s.exp.name))


## Normalization
dt.nuc.sel.norm = myNorm(in.dt = dt.nuc.sel, 
                         in.meas.col = 'Intensity_MeanIntensity_Ratio', 
                         in.rt.min = l.par$plot.norm.skip,
                         in.rt.max = c.t.norm, 
                         in.type = 'z.score')

# rescale dextrin again according to new normalization
dt.im.mer = myResc(
  in.dt = dt.im.mer,
  in.meas.col = s.meas.img,
  in.y.min = -5,
  in.y.h = 4
)

# create pulse averages per channel
dt.im.mer.aggr = dt.im.mer[, .(mean.pulse.resc = mean(get(paste0(s.meas.img, '.resc')))), by = .(Stim_All_Ch, RealTime)]


p.out$imRatio_wDextrin_rZscore_perCh = myGgplotTraj(dt.arg = dt.nuc.sel.norm, x.arg = 'RealTime', y.arg = 'meas.norm', group.arg = paste0(s.met.trackabel, '_uni'), 
                                                    dt.stim.arg = dt.im.mer.aggr, stim.x.arg = 'RealTime', stim.y.arg = 'mean.pulse.resc',
                                                    xaxisbreaks.arg = l.par$plot.x.interv, 
                                                    maxrt.arg =  max(dt.tass$RealTime), 
                                                    summary.arg = TRUE,
                                                    facet.arg = 'Stim_All_Ch', facet.ncol.arg = l.par$plot.facets.ncol.channel, 
                                                    xlab.arg = "\nTime [min]", ylab.arg = "rZ-score (from Illum.Corr)\n", plotlab.arg = paste0('Experiment: ', s.exp.name))



#####
## Save o files
if (l.par$plot.save) {
  
  # Create directory for plots in the currenty working directory
  ifelse(!dir.exists(file.path(".", l.par$dir.plot)), dir.create(file.path(".", l.par$dir.plot)), FALSE)
  
    # all plots
  lapply(names(p.out),
         function(x)
           ggsave(
             filename = paste0(l.par$dir.plot, '/', x, ".pdf"),
             plot = p.out[[x]],
             width = l.par$plot.width, 
             height = l.par$plot.height
           ))
  
  # tables
  pdf(paste0(l.par$dir.plot, '/', "tab_ncells_perSite.pdf"),
      height = 10,
      width = 10)
  grid.table(dt.ncells.persite)
  dev.off()
  
  pdf(paste0(l.par$dir.plot, '/', "tab_ncells_perChannel.pdf"),
      height = 10,
      width = 10)
  grid.table(dt.ncells.perch)
  dev.off()
  
  
  # Create directory for outputs in wide format
  ifelse(!dir.exists(file.path(".", 'data.wide')), dir.create(file.path(".", 'data.wide')), FALSE)
  lapply(names(l.nuc.sel.wide.pulse), function(x) write.csv(l.nuc.sel.wide.pulse[[x]], file = paste0('data.wide/', gsub(' ', '_', gsub(':|\'|/', '', x)), ".csv"), row.names = FALSE) )

  # Create directory for outputs in long format
  ifelse(!dir.exists(file.path(".", 'data.long')), dir.create(file.path(".", 'data.long')), FALSE)
  write.csv(dt.nuc.sel, 'data.long/tCoursesSelected.csv', row.names = FALSE)
  write.csv(dt.im.mer, 'data.long/stimPulses.csv', row.names = FALSE)
  
}

