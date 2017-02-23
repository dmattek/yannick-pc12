require(data.table)

myFreadNuc = function(fileIn) {
  
  # read the rest of the output (except first two rows)
  dt.loc = fread(fileIn)
  
  s.head = names(dt.loc)
  # remove duplicated columns
  dt.loc = dt.loc[, s.head[!duplicated(s.head)], with = FALSE]
  
  # remove unnecesary columns
  s.cols = c(
    'Metadata_C',
    'Metadata_ChannelName',
    'Metadata_ColorFormat',
    'Metadata_FileLocation',
    'Metadata_Frame',
    'Metadata_Plate',
    'Metadata_Site',
    'Metadata_SizeC',
    'Metadata_SizeT',
    'Metadata_SizeX',
    'Metadata_SizeY',
    'Metadata_SizeZ',
    'Metadata_Z',
    'Metadata_Well'
  )
  dt.loc[, (s.cols) := NULL]
  return(dt.loc)
}

myFreadImg = function(fileIn) {
  
  dt.loc = fread(fileIn)

  # make a new column with experiment name based on Metadata_FileLocation
  dt.loc[, Metadata_Experiment := gsub(".*\\/(.*).nd2", "\\1", Metadata_FileLocation)]

  # remove unnecessary columns
  dt.loc[, (c('Metadata_FileLocation', 
              'Group_Index', 
              'Group_Number')) := NULL]
  return(dt.loc)
}


# define filenames for Images and Nuclei
s.files.nuc.core = 'Nuclei.csv'
s.files.img.core = 'Image.csv'

# define directory with merged output
s.dir.out = "output"
s.dir.out.mer = "output.mer"

# search subdirectories for csv files with Nuclei data
s.files.nuc = list.files(
  path = paste0(s.dir.out, '/.'),
  pattern = s.files.nuc.core,
  recursive = TRUE,
  full.names = TRUE
)

s.files.img = list.files(
  path = paste0(s.dir.out, '/.'),
  pattern = s.files.img.core,
  recursive = TRUE,
  full.names = TRUE
)


## Load Nuclei files using custom file reading function
#dt.nuc = do.call(rbind, lapply(s.files.nuc, myFreadNuc))
dt.img = do.call(rbind, lapply(s.files.img, myFreadImg))


# Create directory for merged output in the currenty working directory
ifelse(!dir.exists(file.path(".", s.dir.out.mer)), dir.create(file.path(".", s.dir.out.mer)), FALSE)

#write.csv(file = paste(s.dir.out.mer, s.files.nuc.core, sep = '/'), x = dt.nuc, row.names = FALSE)
#write.csv(file = paste(s.dir.out.mer, s.files.img.core, sep = '/'), x = dt.img, row.names = FALSE)

