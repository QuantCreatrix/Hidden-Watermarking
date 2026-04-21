%% =========================================================================
%  LOCAL FUNCTIONS
% =========================================================================

function [I_wm, log] = embed_watermark_diff(I, W_scrambled, block_map, alpha)
    % Differential embedding: Y(p1)-Y(p2) = +alpha (bit=1) or -alpha (bit=0)
    [rows,cols] = size(I);
    brows = floor(rows/8);
    bcols = floor(cols/8);
    
    % Define 8 pairs from the original 16 mid-frequency positions
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