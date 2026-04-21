% 
% LGAD-W: LBP-Guided Adaptive DCT Watermarking – DIFFERENTIAL EMBEDDING
% High PSNR & High Robustness
% 
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

%% ROBUSTNESS TESTS
fprintf('\n[STEP] Running robustness tests...\n');
results = struct();
sanitize = @(str) strrep(str, '.', 'p');

% JPEG
jpeg_qualities = [90, 80, 70, 60, 50, 40, 30];
fprintf('%-35s %8s %8s\n', 'Attack', 'NC', 'BCR');
fprintf('%s\n', repmat('-',1,55));
for q = jpeg_qualities
    I_att = jpeg_compress(I_watermarked, q);
    W_ext = extract_watermark_diff(I_att, block_map, ARNOLD_ITER, WM_SIZE, ALPHA);
    nc = compute_nc(W_orig, W_ext);
    bcr = compute_bcr(W_orig, W_ext);
    tag = sprintf('jpeg_q%d', q);
    results.(tag).nc = nc; results.(tag).bcr = bcr;
    results.(tag).image = I_att; results.(tag).wm = W_ext;
    fprintf('%-35s %8.4f %8.4f\n', sprintf('JPEG Q=%d',q), nc, bcr);
end

% Gaussian noise
for v = [0.001, 0.005, 0.01]
    I_att = imnoise(I_watermarked, 'gaussian', 0, v);
    W_ext = extract_watermark_diff(I_att, block_map, ARNOLD_ITER, WM_SIZE, ALPHA);
    nc = compute_nc(W_orig, W_ext);
    bcr = compute_bcr(W_orig, W_ext);
    tag = sanitize(sprintf('gauss_%.3f', v));
    results.(tag).nc = nc; results.(tag).bcr = bcr;
    results.(tag).image = I_att; results.(tag).wm = W_ext;
    fprintf('%-35s %8.4f %8.4f\n', sprintf('Gaussian var=%.3f',v), nc, bcr);
end

% Salt & Pepper
for d = [0.01, 0.02, 0.05]
    I_att = imnoise(I_watermarked, 'salt & pepper', d);
    W_ext = extract_watermark_diff(I_att, block_map, ARNOLD_ITER, WM_SIZE, ALPHA);
    nc = compute_nc(W_orig, W_ext);
    bcr = compute_bcr(W_orig, W_ext);
    tag = sanitize(sprintf('sp_%.2f', d));
    results.(tag).nc = nc; results.(tag).bcr = bcr;
    results.(tag).image = I_att; results.(tag).wm = W_ext;
    fprintf('%-35s %8.4f %8.4f\n', sprintf('Salt&Pepper d=%.2f',d), nc, bcr);
end

% Median filter
for f = [3,5]
    I_att = medfilt2(I_watermarked, [f f]);
    W_ext = extract_watermark_diff(I_att, block_map, ARNOLD_ITER, WM_SIZE, ALPHA);
    nc = compute_nc(W_orig, W_ext);
    bcr = compute_bcr(W_orig, W_ext);
    tag = sprintf('median_%d', f);
    results.(tag).nc = nc; results.(tag).bcr = bcr;
    results.(tag).image = I_att; results.(tag).wm = W_ext;
    fprintf('%-35s %8.4f %8.4f\n', sprintf('Median %dx%d',f,f), nc, bcr);
end

% Cropping
for r = [0.10, 0.25, 0.50]
    I_att = crop_image(I_watermarked, r);
    W_ext = extract_watermark_diff(I_att, block_map, ARNOLD_ITER, WM_SIZE, ALPHA);
    nc = compute_nc(W_orig, W_ext);
    bcr = compute_bcr(W_orig, W_ext);
    tag = sprintf('crop_%.0f', r*100);
    results.(tag).nc = nc; results.(tag).bcr = bcr;
    results.(tag).image = I_att; results.(tag).wm = W_ext;
    fprintf('%-35s %8.4f %8.4f\n', sprintf('Crop %d%%', round(r*100)), nc, bcr);
end

% Histogram equalization
I_att = histeq(I_watermarked);
W_ext = extract_watermark_diff(I_att, block_map, ARNOLD_ITER, WM_SIZE, ALPHA);
nc = compute_nc(W_orig, W_ext);
bcr = compute_bcr(W_orig, W_ext);
results.histeq.nc = nc; results.histeq.bcr = bcr;
results.histeq.image = I_att; results.histeq.wm = W_ext;
fprintf('%-35s %8.4f %8.4f\n', 'Histogram Equalization', nc, bcr);
fprintf('%s\n', repmat('-',1,55));

%% SAVE FIGURES
fprintf('\n[STEP] Saving figures...\n');
fig1 = figure('Visible','off','Position',[100 100 900 220]);
subplot(1,4,1); imshow(I); title('Host Image');
subplot(1,4,2); imshow(W_orig); title('Original WM');
subplot(1,4,3); imshow(W_scrambled); title('Scrambled WM');
subplot(1,4,4); imshow(I_watermarked); title(sprintf('Watermarked (PSNR=%.1fdB)',psnr_val));
saveas(fig1, fullfile(OUTPUT_DIR,'Fig1_overview.png'));

fig2 = figure('Visible','off','Position',[100 100 800 220]);
subplot(1,3,1); imshow(block_info.var_map_norm); title('Variance map');
subplot(1,3,2); imshow(block_info.lbp_map_norm); title('LBP texture map');
subplot(1,3,3); imshow(block_map); title('Two‑factor map (union)');
saveas(fig2, fullfile(OUTPUT_DIR,'Fig2_blockmaps.png'));

fig3 = figure('Visible','off','Position',[100 100 400 220]);
subplot(1,2,1); imshow(W_orig); title('Original');
subplot(1,2,2); imshow(W_extracted); title(sprintf('Extracted (NC=%.4f)',nc_clean));
saveas(fig3, fullfile(OUTPUT_DIR,'Fig3_noattack.png'));

fig4 = figure('Visible','off');
nc_vals = arrayfun(@(q) results.(sprintf('jpeg_q%d',q)).nc, jpeg_qualities);
plot(jpeg_qualities, nc_vals, 'b-o','LineWidth',1.5); grid on;
xlabel('JPEG Quality'); ylabel('NC'); title('Robustness vs JPEG');
saveas(fig4, fullfile(OUTPUT_DIR,'Fig4_jpeg.png'));

fig5 = figure('Visible','off','Position',[100 100 900 400]);
attacks = {'jpeg_q50', sanitize('gauss_0.010'), 'crop_25', sanitize('sp_0.05')};
labels = {'JPEG Q=50','Gaussian 0.01','Crop 25%','Salt&Pepper 0.05'};
for k = 1:4
    subplot(2,4,k); imshow(results.(attacks{k}).image); title(labels{k});
    subplot(2,4,k+4); imshow(results.(attacks{k}).wm); 
    title(sprintf('NC=%.3f',results.(attacks{k}).nc));
end
sgtitle('Attacked images (top) and extracted watermarks (bottom)');
saveas(fig5, fullfile(OUTPUT_DIR,'Fig5_attacks.png'));

fprintf('\n[DONE] Results saved in "%s"\n', OUTPUT_DIR);
fprintf('========================================\n');
fprintf(' SUMMARY\n');
fprintf('========================================\n');
fprintf('PSNR = %.2f dB\n', psnr_val);
fprintf('SSIM = %.4f\n', ssim_val);
fprintf('NC (no attack) = %.4f\n', nc_clean);
fprintf('BCR (no attack) = %.4f\n', bcr_clean);
fprintf('========================================\n\n');

%% =========================================================================
%  LOCAL FUNCTIONS
% =========================================================================

function [I_wm, log] = embed_watermark_diff(I, W_scrambled, block_map, alpha)
    % Differential embedding: Y(p1)-Y(p2) = +alpha (bit=1) or -alpha (bit=0)
    [rows,cols] = size(I);
    brows = floor(rows/8);
    bcols = floor(cols/8);
    
    mid_pos = [2,5; 3,4; 4,3; 5,2; 6,1; 1,6; 2,6; 3,5; 4,4; 5,3; 6,2; 7,1; 1,7; 2,7; 3,6; 4,5];
    pairs = [1 2; 3 4; 5 6; 7 8; 9 10; 11 12; 13 14; 15 16]; % indices into mid_pos
    nPairs = size(pairs,1);
    
    eligible = find(block_map);
    bits = W_scrambled(:);
    I_out = I;
    
    for k = 1:min(numel(eligible), numel(bits))
        [r,c] = ind2sub(size(block_map), eligible(k));
        blk = I_out((r-1)*8+1:r*8, (c-1)*8+1:c*8);
        Y = dct2(blk);
        
        % Choose pair for this block (cyclic)
        idx = mod(k-1, nPairs) + 1;
        p1 = mid_pos(pairs(idx,1), :);
        p2 = mid_pos(pairs(idx,2), :);
        d_target = alpha * (2*bits(k) - 1);  % +alpha for bit=1, -alpha for bit=0
        d_curr = Y(p1(1), p1(2)) - Y(p2(1), p2(2));
        delta = (d_target - d_curr) / 2;
        
        % Modify both coefficients minimally
        Y(p1(1), p1(2)) = Y(p1(1), p1(2)) + delta;
        Y(p2(1), p2(2)) = Y(p2(1), p2(2)) - delta;
        
        blk_new = idct2(Y);
        blk_new = max(0, min(1, blk_new));   % clip to valid range
        I_out((r-1)*8+1:r*8, (c-1)*8+1:c*8) = blk_new;
    end
    I_wm = I_out;
    log.embedded = min(numel(eligible), numel(bits));
end

function W_out = extract_watermark_diff(I_wm, block_map, arnold_iter, wm_size, alpha)
    [rows,cols] = size(I_wm);
    brows = floor(rows/8);
    bcols = floor(cols/8);
    
    mid_pos = [2,5; 3,4; 4,3; 5,2; 6,1; 1,6; 2,6; 3,5; 4,4; 5,3; 6,2; 7,1; 1,7; 2,7; 3,6; 4,5];
    pairs = [1 2; 3 4; 5 6; 7 8; 9 10; 11 12; 13 14; 15 16];
    nPairs = size(pairs,1);
    
    eligible = find(block_map);
    bits = zeros(wm_size*wm_size, 1);
    
    for k = 1:min(numel(eligible), wm_size*wm_size)
        [r,c] = ind2sub(size(block_map), eligible(k));
        blk = I_wm((r-1)*8+1:r*8, (c-1)*8+1:c*8);
        Y = dct2(blk);
        idx = mod(k-1, nPairs) + 1;
        p1 = mid_pos(pairs(idx,1), :);
        p2 = mid_pos(pairs(idx,2), :);
        d = Y(p1(1), p1(2)) - Y(p2(1), p2(2));
        bits(k) = double(d > 0);   % threshold at zero
    end
    
    W_scrambled = reshape(bits, wm_size, wm_size);
    W_out = arnold_inverse(W_scrambled, arnold_iter, wm_size);
end

function W = generate_robust_watermark(sz)
    W = zeros(sz, sz);
    border = round(sz/12);
    W(border+1:end-border, border+1:end-border) = 1;
    for i = 1:sz
        for j = 1:sz
            if abs(i-j) <= border || abs(i+j-sz-1) <= border
                W(i,j) = 1;
            end
        end
    end
    W = 1 - W;
end

function W_out = arnold_transform(W_in, iter, sz)
    W_out = W_in;
    for n = 1:iter
        W_new = zeros(sz,sz);
        for i = 0:sz-1
            for j = 0:sz-1
                ni = mod(i+j, sz) + 1;
                nj = mod(i+2*j, sz) + 1;
                W_new(ni,nj) = W_out(i+1,j+1);
            end
        end
        W_out = W_new;
    end
end

function W_out = arnold_inverse(W_in, iter, sz)
    W_out = W_in;
    for n = 1:iter
        W_new = zeros(sz,sz);
        for i = 0:sz-1
            for j = 0:sz-1
                ni = mod(2*i - j, sz) + 1;
                nj = mod(-i + j, sz) + 1;
                W_new(ni,nj) = W_out(i+1,j+1);
            end
        end
        W_out = W_new;
    end
end

function complexity = compute_lbp_complexity(block)
    [br,bc] = size(block);
    total = 0;
    for i = 2:br-1
        for j = 2:bc-1
            center = block(i,j);
            neigh = [block(i-1,j-1), block(i-1,j), block(i-1,j+1), ...
                     block(i,j+1), block(i+1,j+1), block(i+1,j), ...
                     block(i+1,j-1), block(i,j-1)];
            bits = double(neigh >= center);
            trans = sum(abs(diff([bits, bits(1)])));
            total = total + trans;
        end
    end
    complexity = total;
end

function [block_map, info] = compute_two_factor_map(I, var_pct, lbp_pct)
    [rows,cols] = size(I);
    brows = floor(rows/8);
    bcols = floor(cols/8);
    var_scores = zeros(brows,bcols);
    lbp_scores = zeros(brows,bcols);
    
    for r = 1:brows
        for c = 1:bcols
            blk = I((r-1)*8+1:r*8, (c-1)*8+1:c*8);
            var_scores(r,c) = var(blk(:));
            lbp_scores(r,c) = compute_lbp_complexity(blk);
        end
    end
    
    var_thresh = prctile(var_scores(:), 100-var_pct);
    lbp_thresh = prctile(lbp_scores(:), 100-lbp_pct);
    eligible = (var_scores >= var_thresh) | (lbp_scores >= lbp_thresh);
    block_map = eligible;
    
    var_map = zeros(rows,cols);
    lbp_map = zeros(rows,cols);
    v_norm = (var_scores - min(var_scores(:))) / (max(var_scores(:))-min(var_scores(:))+eps);
    l_norm = (lbp_scores - min(lbp_scores(:))) / (max(lbp_scores(:))-min(lbp_scores(:))+eps);
    for r = 1:brows
        for c = 1:bcols
            var_map((r-1)*8+1:r*8, (c-1)*8+1:c*8) = v_norm(r,c);
            lbp_map((r-1)*8+1:r*8, (c-1)*8+1:c*8) = l_norm(r,c);
        end
    end
    info.total_blocks = brows*bcols;
    info.eligible_count = sum(eligible(:));
    info.var_map_norm = var_map;
    info.lbp_map_norm = lbp_map;
end

function I_out = jpeg_compress(I, quality)
    tmp = [tempname '.jpg'];
    imwrite(I, tmp, 'jpg', 'Quality', quality);
    I_out = im2double(imread(tmp));
    if size(I_out,3)==3, I_out = rgb2gray(I_out); end
    delete(tmp);
end

function I_out = crop_image(I, ratio)
    [r,c] = size(I);
    I_out = I;
    cr = floor(r*ratio); cc = floor(c*ratio);
    I_out(1:cr, end-cc+1:end) = 0;
end

function val = compute_psnr(I1, I2)
    mse = mean((I1(:)-I2(:)).^2);
    val = 10*log10(1/(mse+eps));
end

function val = compute_ssim_simple(I1, I2)
    try
        val = ssim(I1,I2);
    catch
        C1=0.01^2; C2=0.03^2;
        mu1=mean(I1(:)); mu2=mean(I2(:));
        s1=std(I1(:)); s2=std(I2(:));
        s12=mean((I1(:)-mu1).*(I2(:)-mu2));
        val = (2*mu1*mu2+C1)*(2*s12+C2)/((mu1^2+mu2^2+C1)*(s1^2+s2^2+C2));
    end
end

function nc = compute_nc(W_orig, W_ext)
    w1 = double(W_orig(:)); w2 = double(W_ext(:));
    nc = sum(w1.*w2) / (sqrt(sum(w1.^2))*sqrt(sum(w2.^2))+eps);
    nc = max(0, min(1, nc));
end

function bcr = compute_bcr(W_orig, W_ext)
    bcr = mean(double(W_orig(:) > 0.5) == double(W_ext(:) > 0.5));
end
