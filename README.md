# EEG_OPMMEG_HealthyAdultStudy_YoungEpilepsy
Data and code for simultaneous EEG and OPM-MEG experiments conducted in healthy adults at Young Epilepsy, UK

Data are stored on the Open Science Framework. Project link: https://osf.io/z6w2t/
These data are 'pre-processed'. OPM-MEG data have had homogeneous field correction applied (as in Tierney et al. 2021, https://www.sciencedirect.com/science/article/pii/S1053811921007576), are demeaned, detrended, bandstop filtered between 49.5 and 50.5 (and 99.5 and 100.5), and bandpass filtered between 1 and 300Hz using FieldTrip. EEG data have not had homogeneous field correction applied but are referenced to the average as well as being demeaned, detrended, bandstop filtered between 49.5 and 50.5 (and 99.5 and 100.5), and bandpass filtered between 1 and 300Hz using FieldTrip.

The code changed throughout the course of the project, and the final versions uploaded here are unlikley to be 'cut and paste' and immediately useable on your PC. However, they are well commented and give you an understanding of the steps taken during analysis, allowing the results to be reproduced. If you have any questions at all, please email me at zseedat@youngepilepsy.org.uk and I will be very happy to help. 

For further methodlogical details, please read our paper: https://direct.mit.edu/imag/article/doi/10.1162/imag_a_00179/120899

Please reference this paper when using the data or code.

