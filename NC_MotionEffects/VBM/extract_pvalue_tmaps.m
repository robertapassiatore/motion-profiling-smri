% --- Load spmT map ---
spmT_path = 'spmT_0002.nii';
V = spm_vol(spmT_path);
T = spm_read_vols(V);

% --- Load degrees of freedom from SPM.mat ---
load('SPM.mat');
df = SPM.xX.erdf;

% --- Compute uncorrected two-tailed p-values ---
P_unc = 2 * (1 - tcdf(abs(T), df));

% --- Flatten and filter p-values ---
p_vec = P_unc(:);
valid_idx = ~isnan(p_vec) & p_vec > 0;
p_valid = p_vec(valid_idx);

% --- Benjamini-Hochberg FDR correction ---
[sorted_p, sort_idx] = sort(p_valid);
m = length(sorted_p);  % number of valid p-values
ranks = (1:m)';

% BH threshold: p(i) <= (i/m) * q
q = 0.05;  % FDR level
thresh_line = (ranks / m) * q;
below_thresh = sorted_p <= thresh_line;

% Find maximum significant p-value
if any(below_thresh)
    max_rank = find(below_thresh, 1, 'last');
    p_thresh = sorted_p(max_rank);
else
    p_thresh = 0;  % No significant voxels
end

% --- Create binary FDR significance mask ---
fdr_mask = zeros(size(p_vec));
fdr_mask(valid_idx) = p_valid <= p_thresh;
fdr_mask = reshape(fdr_mask, size(T));

% --- Save binary mask ---
V_mask = V;
V_mask.fname = 'spmT_0002_masked_FDR.nii';
spm_write_vol(V_mask, fdr_mask);
