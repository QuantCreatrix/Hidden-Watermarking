% =========================================================================
% LGAD-W: LBP-Guided Adaptive DCT Watermarking – DIFFERENTIAL EMBEDDING
% High PSNR & High Robustness
% =========================================================================
clc; clear; close all;

%% CONFIGURATION
IMAGE_FILE = 'cameraman.tif';
ALPHA = 0.15;                   
VARIANCE_PERCENTILE = 40;       
LBP_PERCENTILE      = 40;       
ARNOLD_ITER = 5;
OUTPUT_DIR = 'LGADW_Results';

%% LOAD IMAGE & DETERMINE WATERMARK SIZE
fprintf('\n========================================\n');
fprintf(' LGAD-W (Differential Embedding)\n');
fprintf('========================================\n\n');

if exist(IMAGE_FILE, 'file')
    I = imread(IMAGE_FILE);
else
    I = imread('cameraman.tif');
end
if size(I,3)==3, I = rgb2gray(I); end
I = im2double(I);
[rows, cols] = size(I);
fprintf('Host image: %d x %d\n', rows, cols);

[~, block_info] = compute_two_factor_map(I, VARIANCE_PERCENTILE, LBP_PERCENTILE);
max_bits = block_info.eligible_count;
max_wm_size = floor(sqrt(max_bits));
fprintf('Maximum watermark size: %d x %d\n', max_wm_size, max_wm_size);

WM_SIZE = max_wm_size;
if WM_SIZE < 8
    error('Image too small. Use a larger image or lower percentiles.');
end

W_orig = generate_robust_watermark(WM_SIZE);
fprintf('Using watermark size: %d x %d (%d bits)\n', WM_SIZE, WM_SIZE, WM_SIZE*WM_SIZE);

if ~exist(OUTPUT_DIR, 'dir'), mkdir(OUTPUT_DIR); end

%% EMBED (Differential)
fprintf('\n[STEP] Embedding watermark (differential DCT)...\n');
W_scrambled = arnold_transform(W_orig, ARNOLD_ITER, WM_SIZE);
[block_map, block_info] = compute_two_factor_map(I, VARIANCE_PERCENTILE, LBP_PERCENTILE);
fprintf('Eligible blocks: %d (out of %d)\n', block_info.eligible_count, block_info.total_blocks);

[I_watermarked, ~] = embed_watermark_diff(I, W_scrambled, block_map, ALPHA);

psnr_val = compute_psnr(I, I_watermarked);
ssim_val = compute_ssim_simple(I, I_watermarked);
fprintf('PSNR = %.2f dB, SSIM = %.4f\n', psnr_val, ssim_val);

%% ROUND-TRIP VERIFICATION
fprintf('\n[STEP] Verifying round-trip extraction...\n');
W_extracted = extract_watermark_diff(I_watermarked, block_map, ARNOLD_ITER, WM_SIZE, ALPHA);
nc_clean = compute_nc(W_orig, W_extracted);
bcr_clean = compute_bcr(W_orig, W_extracted);
fprintf('No‑attack NC = %.4f, BCR = %.4f\n', nc_clean, bcr_clean);

if nc_clean < 0.99
    warning('Round‑trip NC low (%f). Check embedding/extraction.', nc_clean);
end