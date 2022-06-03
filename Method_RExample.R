# Method PSD
# Ian Arriaga MacKenzie

library(tidyverse)
library(reshape2)
library(eegkit)
library(psd)
library(RColorBrewer)
library(edf)


# Simulate data, similar to sampled EEG sleep data

# Sampling Rate (Hz)
n = 200
# 10 seconds of data
s = 10
# time vector at 200 Hz
t = seq(0, s, length.out = s*n)

# sinusoidal wave with signal at 1 Hz
s1 = sin(2*pi*t*1)
# signal at 5 Hz
s2 = sin(2*pi*t*5)
# signal at 10 Hz
s3 = sin(2*pi*t*10)
# signal at 20 Hz
s4 = sin(2*pi*t*20)
# signal at 45 Hz
s5 = sin(2*pi*t*45)

# combine all signals into one
comb_sig = s1 + s2 + s3 + s4 + s5

# generate a gaussian error and add it to the combined signal
g_err = rnorm(s * n, sd = 1)
comb_sig_err = comb_sig + g_err

# Apply a 3 Hz - 30 Hz band-pass FIR filter to the signal with error
fir_filt = eegfilter(comb_sig_err, Fs = n, lower = 3, upper = 30, method = "fir1", order = 300)

# combine all signals into a data frame for plotting
signalframe = data.frame(
  time = t,
  signal1 = s1,
  signal2 = s2,
  signal3 = s3,
  signal4 = s4,
  signal5 = s5,
  combsig = comb_sig,
  combsig_err = comb_sig_err,
  combsig_FIR = fir_filt
)

# select individual signals
signalframe_solo = signalframe %>% 
  select(time, signal1, signal2, signal2, signal3, signal4, signal5) %>% 
  rename("Signal 1" = signal1, "Signal 2" = signal2, "Signal 3" = signal3, "Signal 4" = signal4, "Signal 5" = signal5) %>% 
  melt(id.vars = "time")
# plot individual signals
ggplot(signalframe_solo, aes(time, value)) +
  geom_line() +
  facet_grid(variable ~ .) +
  labs(x = "Time", y = NULL) +
  theme_bw()

# select combined signals, error, and band pass
combframe = signalframe %>% 
  select(time, combsig, combsig_err, combsig_FIR) %>% 
  melt(id.vars = "time")
levels(combframe$variable) = c("Combined Signal", "Comb. Signal + Error", "Comb. Signal + FIR Filter")
# plot combined signals
ggplot(combframe, aes(time, value)) +
  geom_line() +
  facet_grid(variable ~ .) +
  theme_bw()


# Power Spectral Density - shows signals present in combined data
# decibals
# PSD of combined signal
eegpsd(signalframe$combsig_err, Fs = n, upper = 50, t = "b")
# PSD of band-pass signal
eegpsd(signalframe$combsig_FIR, Fs = n, upper = 50, t = "b")

# squared microvolts
# PSD of combined signal
eegpsd(signalframe$combsig_err, Fs = n, upper = 50, unit = "mV^2", t = "b")
# PSD of band-pass signal
eegpsd(signalframe$combsig_FIR, Fs = n, upper = 50, unit = "mV^2", t = "b")



# Use of multitaper power spectral density estimation

# PSD of combined signal
PSDest = pspectrum(signalframe$combsig_err, 200)
# Extract values for plotting
psdframe = data.frame(freq = PSDest$freq, 
                      spec = PSDest$spec)
# Trim data
psdframe_trim = psdframe[1:500,]
# Plot, vertical lines at simulated signal
linecol = brewer.pal(6, "Paired")[6]
ggplot(psdframe_trim, aes(freq, spec)) +
  geom_vline(xintercept = 1, color = linecol, size = 1)+
  geom_vline(xintercept = 5, color = linecol, size = 1)+
  geom_vline(xintercept = 10, color = linecol, size = 1)+
  geom_vline(xintercept = 20, color = linecol, size = 1)+
  geom_vline(xintercept = 45, color = linecol, size = 1) +
  geom_point(size = 2, alpha = 0.6) +
  labs(x = "Frequency", y = "Power") +
  theme_bw()


# PSD of band-pass signal
PSDest_fir = pspectrum(signalframe$combsig_FIR, 200)
# Extract values for plotting
psdframe_fir = data.frame(freq = PSDest_fir$freq, 
                      spec = PSDest_fir$spec)
# Trim data
psdframe_trim_fir = psdframe_fir[1:500,]
# Plot, vertical lines at simulated signal
ggplot(psdframe_trim_fir, aes(freq, spec)) +
  geom_vline(xintercept = 5, color = linecol, size = 1)+
  geom_vline(xintercept = 10, color = linecol, size = 1)+
  geom_vline(xintercept = 20, color = linecol, size = 1) +
  geom_point(size = 2, alpha = 0.5) +
  labs(x = "Frequency", y = "Power") +
  theme_bw()






# EEG sampling rate
samprate = 200

# load EDF file
eegmerge = read.edf(paste('~/NJT0884.edf.edf', 
                          sep = ''))

# early time
tstart = eegmerge[[1]]$timestamp.start
earlystart = c('2017-03-13 23:35:00 MDT')


# starting sample point
tearst = abs(as.numeric(difftime(tstart, earlystart, units = 'secs'))) * 200
tearen = tearst + 199

# select EEG lead, also FIR filter EEG lead
tf_eeg = data.frame(
  Time = 0:1999,
  EEG_A1_A2 = eegmerge$signal$EEG_A1_A2$data[tearst:tearen],
  FIR_A1_A2 = eegfilter(eegmerge$signal$EEG_A1_A2$data[tearst:tearen], Fs = 200, lower = 0.5, upper = 40, method = "fir1", order = 300)
)

# Plot raw/filtered EEG lead
eeg_plot = tf_eeg %>% 
  melt(id.vars = "Time")
levels(eeg_plot$variable) = c("A1 - A2", "A1 - A2, FIR Filter")
ggplot(eeg_plot, aes(Time, value)) +
  geom_line() +
  facet_grid(variable ~ .) +
  theme_bw() +
  labs(x = "Time (200 Hz Sampling Rate)", y = NULL)

# PSD estimation on raw/filtered signal
tf_psd = data.frame(
  PSD = pspectrum(tf_eeg$EEG_A1_A2, 200),
  PSD_FIR = pspectrum(tf_eeg$FIR_A1_A2, 200)
)

# Plot the PSD results, for raw and filtered signal
ggplot(tf_psd, aes(PSD.freq, PSD.spec)) +
  geom_line(size = 1) +
  xlim(0, 40) +
  labs(x = "Frequency", y = "Spectrum") +
  theme_bw()

ggplot(tf_psd, aes(PSD_FIR.freq, PSD_FIR.spec)) +
  geom_line(size = 1) +
  xlim(0, 40) +
  labs(x = "Frequency", y = "Spectrum") +
  theme_bw()


