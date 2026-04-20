clc; 
clear;
close all;

% Read images
X = im2double(imread('kohli.jpeg'));
if size(X,3)==3
    X = rgb2gray(X);
end
W = im2bw(imread('Watermark.png'));
W = imresize(W, [size(X,1)/8 size(X,2)/8]);  % FIX: changed /2 to /8

% Embed watermark
[Xw, key] = embed_watermark(X, W);

% Extract watermark
W_ext = extract_watermark(X, Xw, key);

% Similarity (NC)
nc = sum(sum(W .* W_ext)) / sum(sum(W.^2));

% Display
figure;
subplot(1,3,1); imshow(X); title('Original');
subplot(1,3,2); imshow(Xw); title('Watermarked');
subplot(1,3,3); imshow(W_ext); title(['Extracted, NC=' num2str(nc)]);


function [Xw, key] = embed_watermark(X, W)
[M,N] = size(X);
block = 8;

% DCT of image
Y = blockproc(X, [block block], @(b) dct2(b.data));

% Select middle-frequency mask
mask = zeros(8);
mask(3:6,3:6) = 1;

% Watermark already at [M/8, N/8] - no resize needed
wm = W;

% Pseudorandom permutation
rng(123);
perm = randperm(numel(wm));
wm_perm = reshape(wm(perm), size(wm));
key.perm = perm;

% Embedding
Yw = Y;
k = 1;
for i = 1:8:M
    for j = 1:8:N
        block_dct = Y(i:i+7, j:j+7);
        
        [r,c] = ind2sub(size(wm_perm), k);
        bit = wm_perm(r,c);
        
        if mask(4,4) == 1
            coeff = block_dct(4,4);
            if bit == 1
                block_dct(4,4) = coeff + 0.05;
            else
                block_dct(4,4) = coeff - 0.05;
            end
        end
        
        Yw(i:i+7, j:j+7) = block_dct;
        k = k + 1;
    end
end

% Inverse DCT
Xw = blockproc(Yw, [block block], @(b) idct2(b.data));
end


function W_ext = extract_watermark(X, Xw, key)
[M,N] = size(X);
block = 8;

Y  = blockproc(X,  [block block], @(b) dct2(b.data));
Yw = blockproc(Xw, [block block], @(b) dct2(b.data));

wm_size = [M/8 N/8];
W_ext_perm = zeros(wm_size);

k = 1;
for i = 1:8:M
    for j = 1:8:N
        block1 = Y(i:i+7, j:j+7);
        block2 = Yw(i:i+7, j:j+7);
        
        diff = block2(4,4) - block1(4,4);
        
        if diff > 0
            W_ext_perm(k) = 1;
        else
            W_ext_perm(k) = 0;
        end
        
        k = k + 1;
    end
end

% Reverse permutation
W_ext = zeros(wm_size);
W_ext(key.perm) = W_ext_perm(:);
W_ext = reshape(W_ext, wm_size);
end