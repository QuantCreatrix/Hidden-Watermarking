
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

%% --------------------------------------------------------------------------
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