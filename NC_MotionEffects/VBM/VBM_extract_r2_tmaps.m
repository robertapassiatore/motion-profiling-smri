% --- Load the spmT map ---
spmT_path = 'spmT_0002.nii';
V = spm_vol(spmT_path);
T = spm_read_vols(V);

% --- Get residual degrees of freedom from the SPM.mat file ---
load('SPM.mat');  % loads SPM structure
df = SPM.xX.erdf;  % error degrees of freedom

% --- Convert T map to R² map ---
R2 = (T.^2) ./ (T.^2 + df);

% --- Save R² map ---
V_R2 = V;  % copy header
V_R2.fname = 'r2_0002.nii';
spm_write_vol(V_R2, R2);

% --- Load the spmT map ---
spmT_path = 'spmT_0003.nii';
V = spm_vol(spmT_path);
T = spm_read_vols(V);

% --- Get residual degrees of freedom from the SPM.mat file ---
load('SPM.mat');  % loads SPM structure
df = SPM.xX.erdf;  % error degrees of freedom

% --- Convert T map to R² map ---
R2 = (T.^2) ./ (T.^2 + df);

% --- Save R² map ---
V_R2 = V;  % copy header
V_R2.fname = 'r2_0003.nii';
spm_write_vol(V_R2, R2);
