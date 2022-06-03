## MMJ Sleep Study

# load libraries
library(edf)
library(tidyverse)
library(RColorBrewer)
library(eegkit)
library(reshape2)
library(signal)
library(psd)
library(bspec)
library(flux)
library(e1071)
library(ggplot2)
library(factoextra)
library(cluster)

# load overall master excel sheet with participant info
anondatacsv <- read.csv("~/MJsleepstudy/anondatacsv.csv")
# filter for usable thc use and clock time
filtdat = anondatacsv %>% 
  dplyr::filter(THCMinutesSleep != '.') %>% 
  dplyr::filter(Lightsoffclocktime != '.')
# filter for no THC use
nothc = filtdat  %>% 
  dplyr::filter(THCMinutesSleep == '0')
# filter for THC use
yesthc = filtdat %>% 
  dplyr::filter(THCMinutesSleep != '0') 
# change to numeric and choose a middle value of THC use
yesthc$perTHCSleepMODG = as.numeric(yesthc$perTHCSleepMODG)
yesthcfil = yesthc %>% 
  dplyr::filter(perTHCSleepMODG < 0.7 & perTHCSleepMODG > 0.2) %>% 
  dplyr::sample_n(23)

# EEG sampling rate
samprate = 200
n2int = samprate * 60 * 5

testint = 200 * 10
# generate values to pull
teststart = sample(c(1:60000), 6)
testend = teststart + testint
print(teststart); print(testend)

# six 10 second interval samples
teststart = c(34027, 17043, 38171, 47211, 56521, 26842)
testend = c(36027, 19043, 40171, 49211, 58521, 28842)

# sensor names
sensorv = c('A1', 'O2', 'C4', 'C3', 'O1')
# initialize output frame
ff = data.frame(id = character(),
                smoke = character(),
                time = character(),
                timestart = character(),
                delta = numeric(),
                theta = numeric(),
                alpha1 = numeric(),
                alpha2 = numeric(),
                beta = numeric(),
                gamma = numeric(),
                tensecinv = numeric(),
                sensor = character())
# subject ID
ins = c('NJT0884')
# smoke statues
smoke = 'no'
# load EDF file
eegmerge = read.edf(paste('~/MJsleepstudy/usable/', ins,'/', ins, '.edf.edf', 
                          sep = ''))

# set times for early/late N2 sleep
tstart = eegmerge[[1]]$timestamp.start; print(tstart)
tend = eegmerge[[1]]$timestamp.stop; print(tend)
earlystart = c('2017-03-13 23:35:00 MDT')
latestart = c('2017-03-14 05:35:00 MDT')

# begin analysis loop
{
  # starting sample point
  tearst = abs(as.numeric(difftime(tstart, earlystart, units = 'secs'))) * 200
  tearen = tearst + n2int - 1
  # ending sample point
  tlatst = abs(as.numeric(difftime(tstart, latestart, units = 'secs'))) * 200
  tlaten = tlatst + n2int - 1
  # create list
  tf_eeg = list()
  # select EEG leads only, for early N2
  tf_eeg$early = data.frame(
    Time = eegmerge$signal$EEG_A1_A2$t[tearst:tearen],
    EEG_A1_A2 = eegmerge$signal$EEG_A1_A2$data[tearst:tearen],
    EEG_O2_A1 = eegmerge$signal$EEG_O2_A1$data[tearst:tearen],
    EEG_C4_A1 = eegmerge$signal$EEG_C4_A1$data[tearst:tearen],
    EEG_C3_A2 = eegmerge$signal$EEG_C3_A2$data[tearst:tearen],
    EEG_O1_A2 = eegmerge$signal$EEG_O1_A2$data[tearst:tearen]
  )
  # select EEG leads only, for late N2
  tf_eeg$late = data.frame(
    Time = eegmerge$signal$EEG_A1_A2$t[tlatst:tlaten],
    EEG_A1_A2 = eegmerge$signal$EEG_A1_A2$data[tlatst:tlaten],
    EEG_O2_A1 = eegmerge$signal$EEG_O2_A1$data[tlatst:tlaten],
    EEG_C4_A1 = eegmerge$signal$EEG_C4_A1$data[tlatst:tlaten],
    EEG_C3_A2 = eegmerge$signal$EEG_C3_A2$data[tlatst:tlaten],
    EEG_O1_A2 = eegmerge$signal$EEG_O1_A2$data[tlatst:tlaten]
  )
  # loop over early and late
  for (i in 1:2){
    
    tslice = tf_eeg[[i]]
    
    earlylatef = c('early', 'late')
    timestartv = c(earlystart, latestart)
    # loop over all 6 intervals
    for (n in 1:6){
      
      tslice2 = tslice %>% 
        slice(teststart[n]:testend[n])
      # basic bandwidtch filter
      earlyfilt = eegfilter(tslice2[,2:6], 200, .5, 70, method = 'fir1')
      
      # loop over all 5 EEG leads
      for (j in 1:5){
        
        # estimate PSD
        PSDest = pspectrum(earlyfilt[,j], 200, plot = FALSE, AR = FALSE, Nyquist.normalize = TRUE)
        
        psdframe = data.frame(PSDest$freq, 
                              PSDest$spec)
        # relative PSD for each bandwidtch
        psdtot = psdframe %>% 
          dplyr::filter(PSDest.freq < 40 & PSDest.freq > 0.5)
        psddelta = psdtot %>% 
          dplyr::filter(PSDest.freq > 0.5 & PSDest.freq < 4)
        psdtheta = psdtot %>% 
          dplyr::filter(PSDest.freq > 4 & PSDest.freq < 7)
        psdalpha1 = psdtot %>% 
          dplyr::filter(PSDest.freq > 7 & PSDest.freq < 10)
        psdalpha2 = psdtot %>% 
          dplyr::filter(PSDest.freq > 10 & PSDest.freq < 13)
        psdbeta = psdtot %>% 
          dplyr::filter(PSDest.freq > 13 & PSDest.freq < 30)
        psdgamma = psdtot %>% 
          dplyr::filter(PSDest.freq > 30 & PSDest.freq < 40)
        # divide by total AUC
        totalpsd = auc(psdtot$PSDest.freq, psdtot$PSDest.spec)
        deltar = auc(psddelta$PSDest.freq, psddelta$PSDest.spec)/totalpsd
        thetar = auc(psdtheta$PSDest.freq, psdtheta$PSDest.spec)/totalpsd
        alpha1r = auc(psdalpha1$PSDest.freq, psdalpha1$PSDest.spec)/totalpsd
        alpha2r = auc(psdalpha2$PSDest.freq, psdalpha2$PSDest.spec)/totalpsd
        betar = auc(psdbeta$PSDest.freq, psdbeta$PSDest.spec)/totalpsd
        gammar = auc(psdgamma$PSDest.freq, psdgamma$PSDest.spec)/totalpsd
        # normalize
        wavefrac = c(deltar, thetar, alpha1r, alpha2r, betar, gammar)/sum(deltar, thetar, alpha1r, alpha2r, betar, gammar)
        # save to output frame
        outline = data.frame(id = ins,
                             smoke = smoke,
                             time = earlylatef[i],
                             timestart = timestartv[i],
                             delta = wavefrac[1],
                             theta = wavefrac[2],
                             alpha1 = wavefrac[3],
                             alpha2 = wavefrac[4],
                             beta = wavefrac[5],
                             gamma = wavefrac[6],
                             tensecinv = n,
                             sensor = sensorv[j])
        # bind final output
        ff = rbind(ff, outline)
      }
    }
  }
}

# write out table
write.table(ff, file = "PSD_mmjstudy.txt")
# read in out table
PSD_mmjstudy <- read.csv("~/sleepfiles/PSD_mmjstudy.txt", sep="")

# save to seperate DF
ff = PSD_mmjstudy

# average leads by subject
avepsd = ff[,5:10]
n = 5
avepsd2 = aggregate(avepsd,list(rep(1:(nrow(avepsd)%/%n+1),each=n,len=nrow(avepsd))),mean)[-1];
extractrow = ff[seq(1,2160, by = n), c(1:4,11)]
ff_av = cbind(extractrow, avepsd2)
ff_av2 = cbind(extractrow, avepsd2) %>% 
  dplyr::filter(time == 'late')

# early and late plots
plotint = ff_av
plotint2 = ff_av2
# plot labels
smoke.labs <- c("No THC", "Yes THC")
names(smoke.labs) <- c("no", "yes")
time.labs <- c("Early N2 Sleep", "Late N2 Sleep")
names(time.labs) <- c("early", "late")

# testing plot
ggplot(data = plotint, aes(y = gamma)) +
  geom_boxplot() + 
  labs(y = 'Relative PSD', title = 'Gamma') +
  facet_grid(time ~ smoke, 
             labeller = labeller(time = time.labs, smoke = smoke.labs))

# boxplots and two way ANOVA
{boxplot(delta ~ time + smoke, data = plotint)
  res.aov = aov(delta ~ time + smoke, data = plotint)
  summary(res.aov)}
{boxplot(theta ~ time + smoke, data = plotint)
  res.aov = aov(theta ~ time + smoke, data = plotint)
  summary(res.aov)}
{boxplot(alpha1 ~ time + smoke, data = plotint)
  res.aov = aov(alpha1 ~ time + smoke, data = plotint)
  summary(res.aov)}
{boxplot(alpha2 ~ time + smoke, data = plotint)
  res.aov = aov(alpha2 ~ time + smoke, data = plotint)
  summary(res.aov)}
{boxplot(beta ~ time + smoke, data = plotint)
  res.aov = aov(beta ~ time + smoke, data = plotint)
  summary(res.aov)}
{boxplot(gamma ~ time + smoke, data = plotint)
  res.aov = aov(gamma ~ time + smoke, data = plotint)
  summary(res.aov)}

# boxplots and one way ANOVA
{boxplot(delta ~ smoke, data = plotint2)
  res.aov = aov(delta ~ smoke, data = plotint2)
  summary(res.aov)}
{boxplot(theta ~ smoke, data = plotint2)
  res.aov = aov(theta ~ smoke, data = plotint2)
  summary(res.aov)}
{boxplot(alpha1 ~ smoke, data = plotint2)
  res.aov = aov(alpha1 ~ smoke, data = plotint2)
  summary(res.aov)}
{boxplot(alpha2 ~ smoke, data = plotint2)
  res.aov = aov(alpha2 ~ smoke, data = plotint2)
  summary(res.aov)}
{boxplot(beta ~ smoke, data = plotint2)
  res.aov = aov(beta ~ smoke, data = plotint)
  summary(res.aov)}
{boxplot(gamma ~ smoke, data = plotint2)
  res.aov = aov(gamma ~ smoke, data = plotint2)
  summary(res.aov)}



# PCA and K-means Testing, NOT USED IN REPORT
# Slight interesting groupings when principle component analysis done, may be useful moving forward
clusterdf = ff_av[,6:11]
n = 6
clusterdf2 = aggregate(clusterdf,list(rep(1:(nrow(clusterdf)%/%n+1),each=n,len=nrow(clusterdf))),mean)[-1]

namestest = ff_av$id[seq(1, 432, by = 12)]
namestest1 = paste0(namestest[1:18], rep('_early_smoke', 18))
namestest2 = paste0(namestest[1:18], rep('_late_smoke', 18))
namestest3 = paste0(namestest[19:36], rep('_early_nosmoke', 18))
namestest4 = paste0(namestest[19:36], rep('_late_nosmoke', 18))

namestest5 = c(rbind(namestest1, namestest2))
namestest6 = c(rbind(namestest3, namestest4))
namestest7 = c(namestest5, namestest6)

groups = c(c(rbind(rep('early_smoke', 18), rep('late_smoke', 18))), c(rbind(rep('early_nosmoke', 18), rep('late_nosmoke', 18))))

row.names(clusterdf2) = namestest7

k2 <- kmeans(clusterdf2, centers = 4, nstart = 25)
fviz_cluster(k2, data = clusterdf2)

res.pca <- prcomp(clusterdf2, scale. = TRUE)

fviz_pca_ind(res.pca,
             col.ind = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")
)

fviz_pca_var(res.pca,
             col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE
)


groups = as.factor(c(c(rbind(rep('Early THC', 18), rep('Late THC', 18))), c(rbind(rep('Early NoTHC', 18), rep('Late NoTHC', 18)))))

fviz_pca_ind(res.pca,
             habillage = groups,
             geom = 'point',
             addEllipses = TRUE,
             ellipse.type = "confidence",
             legend.title = "Groups"
)

distance <- get_dist(clusterdf2)
fviz_dist(distance, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))
