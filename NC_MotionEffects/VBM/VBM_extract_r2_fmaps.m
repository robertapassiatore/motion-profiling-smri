% Load the SPM.mat file
load('SPM.mat');

% Set which contrast to use
con = 1;

% Get df1 and df2
df1 = SPM.xCon(con).eidf;  % Numerator degrees of freedom
df2 = SPM.xX.erdf;         % Denominator degrees of freedom

% Load the F-map
V = spm_vol('spmF_0001.nii');   % Replace with your actual F-map name
F_data = spm_read_vols(V);

% Compute R² voxel-wise
R2 = (F_data * df1) ./ (F_data * df1 + df2);

% Optional: Set invalid values (e.g. negative F) to NaN
R2(F_data < 0) = NaN;

% Save R² map as NIfTI
V.fname = 'r2_0001.nii';
spm_write_vol(V, R2);
