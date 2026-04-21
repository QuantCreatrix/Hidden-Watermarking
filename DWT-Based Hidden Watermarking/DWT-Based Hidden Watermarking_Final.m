
% BLIND MEDICAL IMAGE WATERMARKING BASED ON LBP-DWT

clear; clc; close all;


%  PARAMETERS
% -------------------------------------------------------------------------
ARNOLD_KEY  = 5;    % Arnold transform iterations (secret key T)
REDUNDANCY  = 3;    % R: each watermark bit embedded R times (majority vote)
DWT_LEVEL   = 2;    % Wavelet decomposition level (1 = original, 2 = improved)
BLOCK_SIZE  = 3;    % LBP block size (fixed at 3x3 as in paper)


%  STEP 1: LOAD HOST IMAGE
% -------------------------------------------------------------------------
img_path = 'medical1.png';

if exist(img_path, 'file')
    host_raw = imread(img_path);
    if size(host_raw, 3) == 3
        host_raw = rgb2gray(host_raw);
    end
    host_raw = imresize(double(host_raw), [512 512]);
else
    warning('Image not found — using synthetic Shepp-Logan phantom.');
    host_raw = double(phantom('Modified Shepp-Logan', 512)) * 255;
end
host_raw = double(host_raw);

%  STEP 2: CREATE BINARY TEXT WATERMARK
%  Encodes patient information as a binary image (100x100 pixels).
%  You can replace this block with your own watermark image.
% -------------------------------------------------------------------------
wm_size    = [100, 100];
wm_img     = zeros(wm_size);
wm_img(10:30, 5:95) = 1;   % Simulated text line 1
wm_img(40:55, 5:70) = 1;   % Simulated text line 2
wm_img(65:80, 5:85) = 1;   % Simulated text line 3

W_original = wm_img;   % binary [0,1] watermark, size 100x100


%  STEP 3: EMBED WATERMARK
% -------------------------------------------------------------------------
fprintf('=== EMBEDDING ===\n');
[WMI, embed_info] = embed_watermark(host_raw, W_original, ...
                                    ARNOLD_KEY, REDUNDANCY, DWT_LEVEL);

PSNR_val = compute_psnr(host_raw, WMI);
SSIM_val = compute_ssim(host_raw, WMI);
fprintf('  PSNR = %.4f dB\n', PSNR_val);
fprintf('  SSIM = %.6f\n',   SSIM_val);

%  STEP 4: VISUALISE EMBEDDING RESULT
% -------------------------------------------------------------------------
figure('Name','Watermarking Result','NumberTitle','off','Position',[50 50 1100 380]);

subplot(1,3,1);
imshow(uint8(host_raw));
title('Original Medical Image','FontSize',11);

subplot(1,3,2);
imshow(uint8(WMI));
title(sprintf('Watermarked Image\nPSNR=%.2f dB  SSIM=%.4f', PSNR_val, SSIM_val),'FontSize',11);

subplot(1,3,3);
imshow(W_original);
title('Embedded Watermark','FontSize',11);

%  STEP 5: TEST AGAINST ATTACKS & EXTRACT WATERMARKS
% -------------------------------------------------------------------------
fprintf('\n=== ROBUSTNESS TESTS ===\n');

% Define attacks
attacks = define_attacks();
n_attacks = length(attacks);

NC_vals  = zeros(1, n_attacks);
BCR_vals = zeros(1, n_attacks);
EW_cells = cell(1, n_attacks);

for k = 1:n_attacks
    att = attacks(k);
    fprintf('  Attack [%d/%d]: %s ... ', k, n_attacks, att.name);

    % Apply attack
    attacked = apply_attack(WMI, att);

    % Extract watermark
    EW = extract_watermark(attacked, embed_info, ARNOLD_KEY, REDUNDANCY, DWT_LEVEL);

    % Compute metrics
    NC_vals(k)  = compute_nc(W_original, EW);
    BCR_vals(k) = compute_bcr(W_original, EW);
    EW_cells{k} = EW;

    fprintf('NC=%.4f  BCR=%.4f\n', NC_vals(k), BCR_vals(k));
end

%  STEP 6: VISUALISE EXTRACTED WATERMARKS UNDER EACH ATTACK
% -------------------------------------------------------------------------
cols = 4;
rows = ceil(n_attacks / cols);
figure('Name','Extracted Watermarks Under Attacks','NumberTitle','off', ...
       'Position',[50 50 1200 300*rows]);

for k = 1:n_attacks
    subplot(rows, cols, k);
    imshow(EW_cells{k});
    title(sprintf('%s\nNC=%.3f BCR=%.3f', attacks(k).name, NC_vals(k), BCR_vals(k)), ...
          'FontSize', 9, 'Interpreter', 'none');
end
sgtitle('Extracted Watermarks After Attacks','FontSize',12);

%  STEP 7: BAR CHART SUMMARY
% -------------------------------------------------------------------------
attack_names = {attacks.name};

figure('Name','NC and BCR Summary','NumberTitle','off','Position',[100 100 1000 420]);
x = 1:n_attacks;
bar_data = [NC_vals; BCR_vals]';
b = bar(x, bar_data, 0.7);
b(1).FaceColor = [0.25 0.45 0.78];
b(2).FaceColor = [0.85 0.35 0.25];
set(gca, 'XTickLabel', attack_names, 'XTickLabelRotation', 35, 'FontSize', 9);
ylim([0 1.1]);
yline(0.9,'--k','LineWidth',1.2);
legend({'NC','BCR'}, 'Location','southwest','FontSize',10);
ylabel('Value','FontSize',11);
title('Watermark Robustness: NC and BCR Under Attacks','FontSize',12);
grid on;

fprintf('\n=== DONE ===\n');

% ==========================================================================
%  FUNCTIONS
% ==========================================================================

function [WMI, info] = embed_watermark(host, W, arnold_key, R, dwt_level)
% EMBED_WATERMARK  Main embedding function.
%
%  Inputs:
%    host        - double grayscale host image (512x512)
%    W           - binary watermark (any size, will be used as-is)
%    arnold_key  - number of Arnold iterations (secret key T)
%    R           - redundancy factor (each bit embedded R times)
%    dwt_level   - DWT decomposition level (1 or 2)
%
%  Outputs:
%    WMI         - watermarked image (double, same size as host)
%    info        - struct containing metadata needed for extraction

    % --- Normalise host to [0,1] ---
    host_max = max(host(:));
    host_min = min(host(:));
    I_norm   = (host - host_min) / (host_max - host_min + eps);

    % --- [I1] Multi-level DWT decomposition ---
    [LL, detail_cells] = multilevel_dwt(I_norm, dwt_level);

    % Multiply by 10 (as in original paper) to work with integer-like values
    LL_scaled = LL * 10;

    % --- [I2] Apply rotation-invariant LBP to LL sub-band integer part ---
    LL_int  = floor(LL_scaled);
    LL_frac = LL_scaled - LL_int;
    LBP_map = compute_lbp_riu2(LL_int, 3);  % LBP^riu2 on 3x3 blocks

    % --- Prepare watermark ---
    % Arnold scrambling
    WA = arnold_transform(W, arnold_key);

    % ZigZag vectorisation
    wm_vector = zigzag_scan(WA);              % 1D binary vector
    wm_len    = length(wm_vector);

    % Segment into 8-bit groups
    n_blocks_needed = R * ceil(wm_len / 8);   % [I3] R copies of each segment

    % --- Embed into LL blocks ---
    [rows_ll, cols_ll] = size(LL_int);
    n_blocks_r = floor(rows_ll / 3);
    n_blocks_c = floor(cols_ll / 3);
    n_blocks   = n_blocks_r * n_blocks_c;

    if n_blocks < n_blocks_needed
        warning('Not enough blocks for R=%d redundancy. Reducing R.', R);
        R = floor(n_blocks * 8 / wm_len);
        R = max(R, 1);
        n_blocks_needed = R * ceil(wm_len / 8);
    end

    LL_watermarked = LL_int;
    block_idx = 0;

    for r = 1 : 3 : (n_blocks_r * 3)
        for c = 1 : 3 : (n_blocks_c * 3)
            if r + 2 > rows_ll || c + 2 > cols_ll
                continue;
            end
            block_idx = block_idx + 1;

            % Determine which watermark segment this block carries
            % [I3] Redundancy: block i carries segment mod(i-1, n_segs)+1
            n_segs   = ceil(wm_len / 8);
            seg_idx  = mod(block_idx - 1, n_segs) + 1;   % cycles over segments

            % Only embed if within redundancy budget
            copy_num = ceil(block_idx / n_segs);
            if copy_num > R
                continue;
            end

            % Extract 8-bit segment from watermark vector
            seg_start = (seg_idx - 1) * 8 + 1;
            seg_end   = min(seg_start + 7, wm_len);
            wm_seg    = zeros(1, 8);
            wm_seg(1 : seg_end - seg_start + 1) = wm_vector(seg_start : seg_end);

            % Get LBP code for this block
            lbp_code = LBP_map(ceil(r/3), ceil(c/3));
            lbp_bits = de2bi(lbp_code, 8, 'left-msb');

            % XOR watermark segment with LBP code
            Xj = bitxor(uint8(wm_seg), uint8(lbp_bits));

            % Pair-swap shuffle (Eq. 9 in paper)
            Yj = pair_swap(Xj);

            % LSB embedding into block coefficients
            PI_block = LL_watermarked(r:r+2, c:c+2);
            PI_block = embed_lsb_block(PI_block, Yj);

            % [I2+I3] PreserveLBP: adjust so neighborhood relations unchanged
            PI_center = PI_block(2,2);
            PI_orig   = LL_int(r:r+2, c:c+2);
            PI_block  = preserve_lbp(PI_block, PI_orig, PI_center);

            LL_watermarked(r:r+2, c:c+2) = PI_block;
        end
    end

    % Recombine integer and fractional parts
    LL_final = (LL_watermarked + LL_frac) / 10;

    % [I1] Multi-level inverse DWT
    I_norm_wm = multilevel_idwt(LL_final, detail_cells, dwt_level);

    % Denormalise
    WMI = I_norm_wm * (host_max - host_min) + host_min;
    WMI = max(0, min(255, WMI));

    % Store metadata for extraction
    info.host_max   = host_max;
    info.host_min   = host_min;
    info.wm_size    = size(W);
    info.wm_len     = wm_len;
    info.R          = R;
    info.block_idx_max = block_idx;
end

function EW = extract_watermark(attacked_img, info, arnold_key, R, dwt_level)
% EXTRACT_WATERMARK  Blind extraction (no original image required).

    % Normalise
    host_max = info.host_max;
    host_min = info.host_min;
    I_norm   = (attacked_img - host_min) / (host_max - host_min + eps);

    % [I1] Multi-level DWT
    [LL, ~] = multilevel_dwt(I_norm, dwt_level);
    LL_scaled = LL * 10;
    LL_int    = floor(LL_scaled);

    % [I2] Rotation-invariant LBP
    LBP_map = compute_lbp_riu2(LL_int, 3);

    [rows_ll, cols_ll] = size(LL_int);
    n_blocks_r = floor(rows_ll / 3);
    n_blocks_c = floor(cols_ll / 3);

    wm_len = info.wm_len;
    n_segs = ceil(wm_len / 8);

    % Accumulate votes for each bit: votes(bit_idx, value) -> count of 0s and 1s
    bit_votes = zeros(wm_len, 2);   % column 1=votes for 0, column 2=votes for 1

    block_idx = 0;
    for r = 1 : 3 : (n_blocks_r * 3)
        for c = 1 : 3 : (n_blocks_c * 3)
            if r + 2 > rows_ll || c + 2 > cols_ll
                continue;
            end
            block_idx = block_idx + 1;

            seg_idx  = mod(block_idx - 1, n_segs) + 1;
            copy_num = ceil(block_idx / n_segs);
            if copy_num > R
                continue;
            end

            seg_start = (seg_idx - 1) * 8 + 1;
            seg_end   = min(seg_start + 7, wm_len);

            % Get LBP code
            lbp_code = LBP_map(ceil(r/3), ceil(c/3));
            lbp_bits = de2bi(lbp_code, 8, 'left-msb');

            % Extract LSBs from block
            PI_block = LL_int(r:r+2, c:c+2);
            Yj = extract_lsb_block(PI_block);

            % Unshuffle (reverse of pair-swap)
            Xj = pair_swap(Yj);   % pair_swap is its own inverse

            % XOR to recover watermark segment
            wm_seg_extracted = bitxor(uint8(Xj), uint8(lbp_bits));

            % Accumulate votes
            for bit_k = 1 : (seg_end - seg_start + 1)
                abs_bit = seg_start + bit_k - 1;
                if abs_bit <= wm_len
                    v = wm_seg_extracted(bit_k) + 1;   % v=1 if bit=0, v=2 if bit=1
                    bit_votes(abs_bit, v) = bit_votes(abs_bit, v) + 1;
                end
            end
        end
    end

    % [I3] Majority vote
    wm_vector_extracted = double(bit_votes(:,2) > bit_votes(:,1))';

    % Reconstruct watermark matrix via inverse ZigZag
    WA_extracted = inv_zigzag_scan(wm_vector_extracted, info.wm_size);

    % Inverse Arnold transform
    EW = inv_arnold_transform(WA_extracted, arnold_key);
    EW = double(EW > 0.5);
end

%  ATTACK DEFINITIONS
% --------------------------------------------------------------------------
function attacks = define_attacks()
% Returns a struct array of attack definitions.

    attacks = struct('name', {}, 'type', {}, 'params', {});
    k = 0;

    k=k+1; attacks(k).name='No attack';          attacks(k).type='none';         attacks(k).params=struct();
    k=k+1; attacks(k).name='JPEG Q=50';          attacks(k).type='jpeg';         attacks(k).params=struct('quality',50);
    k=k+1; attacks(k).name='JPEG Q=70';          attacks(k).type='jpeg';         attacks(k).params=struct('quality',70);
    k=k+1; attacks(k).name='Gaussian noise';     attacks(k).type='gaussnoise';   attacks(k).params=struct('var',0.001);
    k=k+1; attacks(k).name='Salt & pepper';      attacks(k).type='sap';          attacks(k).params=struct('density',0.01);
    k=k+1; attacks(k).name='Median filter 3x3';  attacks(k).type='median';       attacks(k).params=struct('size',3);
    k=k+1; attacks(k).name='Gaussian filter 3x3';attacks(k).type='gaussfilt';    attacks(k).params=struct('size',3,'sigma',0.5);
    k=k+1; attacks(k).name='Crop 10%';           attacks(k).type='crop';         attacks(k).params=struct('fraction',0.10);
    k=k+1; attacks(k).name='Resize 256->512';    attacks(k).type='resize';       attacks(k).params=struct('scale',0.5);
end

function out = apply_attack(img, att)
% APPLY_ATTACK  Applies one attack to a double [0..255] image.

    out = img;
    switch att.type
        case 'none'
            % no-op

        case 'jpeg'
            % Write to temp JPEG and reload
            tmp = tempname;
            imwrite(uint8(img), [tmp '.jpg'], 'Quality', att.params.quality);
            out = double(imread([tmp '.jpg']));
            delete([tmp '.jpg']);

        case 'gaussnoise'
            noise = sqrt(att.params.var * 255^2) .* randn(size(img));
            out   = img + noise;
            out   = max(0, min(255, out));

        case 'sap'
            tmp = imnoise(uint8(img), 'salt & pepper', att.params.density);
            out = double(tmp);

        case 'median'
            out = double(medfilt2(uint8(img), [att.params.size att.params.size]));

        case 'gaussfilt'
            h   = fspecial('gaussian', att.params.size, att.params.sigma);
            out = imfilter(img, h, 'replicate');

        case 'crop'
            [H, W] = size(img);
            fr = att.params.fraction;
            r1 = round(H*fr)+1; r2 = H-round(H*fr);
            c1 = round(W*fr)+1; c2 = W-round(W*fr);
            cropped = img(r1:r2, c1:c2);
            out = imresize(cropped, [H W]);

        case 'resize'
            [H, W] = size(img);
            s   = att.params.scale;
            tmp = imresize(img, s);
            out = imresize(tmp, [H W]);
    end
end

% ==========================================================================
%  CORE ALGORITHM FUNCTIONS
% ==========================================================================

function [LL, detail_cells] = multilevel_dwt(I, level)
% MULTILEVEL_DWT  Applies Haar DWT iteratively to produce LL at given level.
% Returns the LL sub-band and a cell array of detail coefficients (for reconstruction).

    detail_cells = cell(level, 1);
    current = I;

    for lv = 1 : level
        [LL, LH, HL, HH] = haar_dwt2(current);
        detail_cells{lv}.LH = LH;
        detail_cells{lv}.HL = HL;
        detail_cells{lv}.HH = HH;
        current = LL;
    end
    LL = current;
end

function I_rec = multilevel_idwt(LL, detail_cells, level)
% MULTILEVEL_IDWT  Reconstructs image from LL and detail_cells.

    current = LL;
    for lv = level : -1 : 1
        current = haar_idwt2(current, detail_cells{lv}.LH, ...
                             detail_cells{lv}.HL, detail_cells{lv}.HH);
    end
    I_rec = current;
end

function [LL, LH, HL, HH] = haar_dwt2(I)
% HAAR_DWT2  Single-level 2D Haar DWT using the paper's formula (Eqs. 3-6).

    % Row-wise Haar
    Ir = (I(:,1:2:end) + I(:,2:2:end)) / 2;   % Low
    Ic = (I(:,1:2:end) - I(:,2:2:end)) / 2;   % High

    % Column-wise Haar on row results
    LL = (Ir(1:2:end,:) + Ir(2:2:end,:)) / 2;
    LH = (Ir(1:2:end,:) - Ir(2:2:end,:)) / 2;
    HL = (Ic(1:2:end,:) + Ic(2:2:end,:)) / 2;
    HH = (Ic(1:2:end,:) - Ic(2:2:end,:)) / 2;
end

function I = haar_idwt2(LL, LH, HL, HH)
% HAAR_IDWT2  Single-level 2D inverse Haar DWT.

    [H, W] = size(LL);
    Ir = zeros(2*H, W);
    Ic = zeros(2*H, W);

    Ir(1:2:end,:) = LL + LH;
    Ir(2:2:end,:) = LL - LH;
    Ic(1:2:end,:) = HL + HH;
    Ic(2:2:end,:) = HL - HH;

    I = zeros(2*H, 2*W);
    I(:,1:2:end) = Ir + Ic;
    I(:,2:2:end) = Ir - Ic;
end

function LBP_map = compute_lbp_riu2(I, ~)
% COMPUTE_LBP_RIU2  [I2] Rotation-invariant uniform LBP (LBP^riu2).
%
%  Standard LBP: compare each of 8 neighbors to center, MSB-first.
%  LBP^riu2: count the number of 1-bits in the pattern (0..8).
%  This mapping is invariant to rotation of the neighborhood.
%
%  Returns LBP_map of size [floor(R/3) x floor(C/3)] where each cell
%  contains the riu2 code (0..8) of the 3x3 block.

    [rows, cols] = size(I);
    n_br = floor(rows / 3);
    n_bc = floor(cols / 3);
    LBP_map = zeros(n_br, n_bc);

    % Neighbor offsets in clockwise order starting from top-left
    neighbors_r = [-1,-1,-1, 0, 1, 1, 1, 0];
    neighbors_c = [-1, 0, 1, 1, 1, 0,-1,-1];

    for br = 1 : n_br
        for bc = 1 : n_bc
            r0 = (br-1)*3 + 2;   % center row of block (1-indexed)
            c0 = (bc-1)*3 + 2;   % center col
            Pc = I(r0, c0);

            bits = zeros(1,8);
            for n = 1:8
                rn = r0 + neighbors_r(n);
                cn = c0 + neighbors_c(n);
                % Clamp to image boundary
                rn = max(1, min(rows, rn));
                cn = max(1, min(cols, cn));
                bits(n) = double(I(rn, cn) <= Pc);  % 1 if neighbor <= center
            end

            % LBP^riu2: count 1-bits (rotation-invariant measure)
            % (Instead of the binary code, we use the popcount, which equals
            %  the original LBP value mod its circular rotations)
            LBP_map(br, bc) = sum(bits);
        end
    end
end

function WA = arnold_transform(W, iterations)
% ARNOLD_TRANSFORM  Applies Arnold cat map to binary image W.

    [N, M] = size(W);
    WA = W;
    for t = 1 : iterations
        WA_new = zeros(N, M);
        for i = 1 : N
            for j = 1 : M
                i2 = mod(i + j - 2, N) + 1;
                j2 = mod(i + 2*j - 3, M) + 1;
                WA_new(i2, j2) = WA(i, j);
            end
        end
        WA = WA_new;
    end
end

function W = inv_arnold_transform(WA, iterations)
% INV_ARNOLD_TRANSFORM  Inverse Arnold cat map (inverse permutation matrix).

    [N, M] = size(WA);
    W = WA;
    for t = 1 : iterations
        W_new = zeros(N, M);
        for i = 1 : N
            for j = 1 : M
                i2 = mod(2*i - j - 1, N) + 1;
                j2 = mod(-i + j, M) + 1;
                W_new(i2, j2) = W(i, j);
            end
        end
        W = W_new;
    end
end

function v = zigzag_scan(M)
% ZIGZAG_SCAN  Reads matrix in ZigZag order, returns 1D vector.

    [H, W] = size(M);
    v = zeros(1, H*W);
    idx = 1;
    for s = 0 : H + W - 2
        if mod(s, 2) == 0
            % going up-right
            r = min(s, H-1);
            c = s - r;
            while r >= 0 && c < W
                v(idx) = M(r+1, c+1);
                idx = idx + 1;
                r = r - 1;
                c = c + 1;
            end
        else
            % going down-left
            c = min(s, W-1);
            r = s - c;
            while c >= 0 && r < H
                v(idx) = M(r+1, c+1);
                idx = idx + 1;
                r = r + 1;
                c = c - 1;
            end
        end
    end
end

function M = inv_zigzag_scan(v, sz)
% INV_ZIGZAG_SCAN  Reconstructs matrix from ZigZag-ordered vector.

    H = sz(1); W = sz(2);
    M = zeros(H, W);
    idx = 1;
    len = length(v);
    for s = 0 : H + W - 2
        if mod(s, 2) == 0
            r = min(s, H-1);
            c = s - r;
            while r >= 0 && c < W
                if idx <= len
                    M(r+1, c+1) = v(idx);
                    idx = idx + 1;
                end
                r = r - 1;
                c = c + 1;
            end
        else
            c = min(s, W-1);
            r = s - c;
            while c >= 0 && r < H
                if idx <= len
                    M(r+1, c+1) = v(idx);
                    idx = idx + 1;
                end
                r = r + 1;
                c = c - 1;
            end
        end
    end
end

function Yj = pair_swap(Xj)
% PAIR_SWAP  Swaps adjacent bit pairs (Eq. 9 in paper).
% Yj = [x1,x0, x3,x2, x5,x4, x7,x6]

    Yj = Xj;
    Yj(1) = Xj(2);  Yj(2) = Xj(1);
    Yj(3) = Xj(4);  Yj(4) = Xj(3);
    Yj(5) = Xj(6);  Yj(6) = Xj(5);
    Yj(7) = Xj(8);  Yj(8) = Xj(7);
end

function PI_out = embed_lsb_block(PI_block, bits)
% EMBED_LSB_BLOCK  Embeds 8 bits into the 8 neighbors of the 3x3 block.
% The center pixel is skipped; bits go into the 8 surrounding pixels.

    PI_out = PI_block;
    positions = [1,1; 1,2; 1,3; 2,3; 3,3; 3,2; 3,1; 2,1];  % clockwise
    for k = 1:8
        r = positions(k,1); c = positions(k,2);
        val = floor(PI_block(r,c));
        % Set LSB to bit value
        val = floor(val/2)*2 + bits(k);
        PI_out(r,c) = val;
    end
end

function bits = extract_lsb_block(PI_block)
% EXTRACT_LSB_BLOCK  Extracts 8 LSBs from the 8 surrounding neighbors.

    positions = [1,1; 1,2; 1,3; 2,3; 3,3; 3,2; 3,1; 2,1];
    bits = zeros(1,8);
    for k = 1:8
        r = positions(k,1); c = positions(k,2);
        val = floor(PI_block(r,c));
        bits(k) = mod(val, 2);
    end
end

function PI_star = preserve_lbp(PI_prime, PI_orig, PI_center)
% PRESERVE_LBP  Adjusts neighbor pixel values to maintain LBP code (Eq. 10).
% If embedding changed a neighbor so that the comparison with center flipped,
% add or subtract 2 to fix it.

    [H, W] = size(PI_prime);
    PI_star = PI_prime;
    positions = [1,1; 1,2; 1,3; 2,3; 3,3; 3,2; 3,1; 2,1];
    for k = 1:8
        r = positions(k,1); c = positions(k,2);
        orig_rel  = (PI_center >= PI_orig(r,c));    % original comparison
        prime_rel = (PI_center >= PI_prime(r,c));   % after LSB embedding

        if orig_rel ~= prime_rel
            % Relationship flipped -> fix with +/-2
            if ~orig_rel && prime_rel
                % Was: center < neighbor. Now: center >= neighbor. Push neighbor up.
                PI_star(r,c) = PI_prime(r,c) + 2;
            else
                % Was: center >= neighbor. Now: center < neighbor. Pull neighbor down.
                PI_star(r,c) = PI_prime(r,c) - 2;
            end
        end
    end
end

% ==========================================================================
%  QUALITY METRICS
% ==========================================================================

function p = compute_psnr(I, J)
% COMPUTE_PSNR  Peak signal-to-noise ratio (dB).

    mse = mean((double(I(:)) - double(J(:))).^2);
    if mse == 0
        p = Inf;
    else
        p = 10 * log10(255^2 / mse);
    end
end

function s = compute_ssim(I, J)
% COMPUTE_SSIM  Structural similarity index (window-based).

    I = double(I); J = double(J);
    C1 = (0.01*255)^2;
    C2 = (0.03*255)^2;
    win = fspecial('gaussian', 11, 1.5);

    mu1 = imfilter(I, win, 'replicate');
    mu2 = imfilter(J, win, 'replicate');
    mu1_sq = mu1.^2; mu2_sq = mu2.^2; mu1_mu2 = mu1.*mu2;
    sig1_sq = imfilter(I.^2, win, 'replicate') - mu1_sq;
    sig2_sq = imfilter(J.^2, win, 'replicate') - mu2_sq;
    sig12   = imfilter(I.*J, win, 'replicate') - mu1_mu2;

    ssim_map = ((2*mu1_mu2 + C1).*(2*sig12 + C2)) ./ ...
               ((mu1_sq + mu2_sq + C1).*(sig1_sq + sig2_sq + C2));
    s = mean(ssim_map(:));
end

function nc = compute_nc(W, EW)
% COMPUTE_NC  Normalized correlation between original and extracted watermark.

    W  = double(W(:));
    EW = double(EW(:));
    W  = W  - mean(W);
    EW = EW - mean(EW);
    denom = sqrt(sum(W.^2) * sum(EW.^2));
    if denom < eps
        nc = 0;
    else
        nc = sum(W .* EW) / denom;
    end
    nc = max(0, min(1, nc));
end

function bcr = compute_bcr(W, EW)
% COMPUTE_BCR  Bit correct rate = fraction of correctly recovered bits.

    W  = double(W(:) > 0.5);
    EW = double(EW(:) > 0.5);
    bcr = sum(W == EW) / numel(W);
end