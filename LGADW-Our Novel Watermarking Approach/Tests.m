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