load('SPM.mat');

% Example for F-contrast
con = 1;  % contrast number
df1 = SPM.xCon(con).eidf;
df2 = SPM.xX.erdf;

% Load the SPM F-map
V = spm_vol('spmF_0001.nii');
Y = spm_read_vols(V);

% Compute uncorrected two-tailed p-values
p_uncorrected = 1 - fcdf(Y, df1, df2);

% Save the p-map
V.fname = 'pval_0001.nii';
spm_write_vol(V, p_uncorrected);
