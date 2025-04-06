FigureX.m is the MATLAB script to generate Figure X.
Please add 'tools' to the workspace.

yeEX: EOG data in Experiment X (time x channels x trials x subjects x conditions)
Channel 1 is vertical EOG and channel 2 is horizontal EOG.

ybEX: EOG data in Experiment X (time x channels x trials x subjects x conditions)
Only channel Cz is nonzero and the index of channel Cz is specified in each script.

etkEX: eyetracking data in Experiment X (time x channels x trials x subjects x conditions)
Channel 1 is blink; channels 2 is saccade; channel 3 is pupil size.


The function to do FDR correction can be downloaded from 
https://www.mathworks.com/matlabcentral/fileexchange/27418-fdr_bh

[h, crit_p, adj_p]=fdr_bh(pvals,q,method,report);


Please contact Dr. Nai Ding (ding_nai@zju.edu.cn) for any question.
