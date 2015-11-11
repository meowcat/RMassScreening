
dir <- "C:/Daten/Algae/Bioconcentration/20151105 screener package test/"

inputDir <- "C:/Input_data/Algae/20151023 100 rescreening"
files <- list.files(inputDir, ".mzXML", full.names=TRUE)

setwd(dir)
outputDir <- paste0(dir, "results/")

# enviPick all files, pos/neg separately, into outputDir. This saves both RData and csv output
batchPick(files, outputDir)


# generate suspect mass list
compounds <- read.csv("input/compounds.csv", stringsAsFactors = FALSE)
reactions <- read.csv("input/reactions.csv", stringsAsFactors = FALSE)
suspects <- combineReactions(compounds, reactions)
save(suspects, file="results/suspects.RData")

# build profile container for positive mode
profiles.pos <- fillProfiles(outputDir, files, "+")
# perform profiling
profiles.pos <- computeProfiles(profiles.pos, dmass=3, dret=60)
save(profiles.pos, file = "results/profiles.pos.RData")

# assign sampleIDs to files and merge the info into the profile table
sampleList <- read.csv("input/sampleList.csv")
sampleList <- assignSamples(files, sampleList)
profiles.pos <- assignProfiles(profiles.pos, sampleList)

# screen profiles
hits.pos <- screenProfiles(profiles.pos, suspects, "+", 2)
save(hits.pos, file = "results/hits.pos.RData")

# build time trends
sampleAssignment <- read.csv("input/sampleAssignment.csv", stringsAsFactors = FALSE)
sampleGroups <- read.csv("input/sampleGroups.csv", stringsAsFactors = FALSE)
summaryTables <- groupSummaries(
  profiles.pos, sampleList, sampleAssignment, sampleGroups, groupBy="time", groups=c("t0", "t1", "t2", "t3")
)

# compute additional values for each sample group, at will
summaryTables <- lapply(summaryTables, function(df)
{
  within(df, {
    mean = rowMeans(cbind(t1,t2,t3))
  })
})

# merge to the "tt.total" table
totalTable <- mergeGroups(profiles.pos, groupList, summaryTables)

timepoints <- data.frame(
  name = c("t0", "t1", "t2", "t3"),
  t=c(0,1,3,5)
  )

runViewer(totalTable, hits.pos,timepoints ,sampleGroups)

