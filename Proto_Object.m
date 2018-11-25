function [block,flag,S,E] = Proto_Object(smap,blocksize,thres)

% [RETURNS]
% block : block image (1: salient 0: not salient)
% flag : flag for saliency
% S: starting points of block (top-left)
% E: ending points of block (bottom-right)
%
% [PARAMETERS]
% smap   : saliency map
% blocksize : size of block
% thres : threshold for saliency
 

[M N] = size(smap);
% Make integral image of saliency map
smap_integral = cumsum(cumsum(smap),2);
x_block = blocksize(1);
y_block = blocksize(2);

flag = zeros(y_block,x_block);
block_width = floor(N/x_block);
block_height = floor(M/y_block);

block = zeros(M,N);

% Compute average saliency in the blocks and declare the region if the
% average saliency is over the threshold

for i = 1:x_block
    for j = 1:y_block
        S{i,j}.x = (i-1) * block_width + 1;
        S{i,j}.y = (j-1) * block_height + 1;
        E{i,j}.x = i * block_width;
        E{i,j}.y = j * block_height;
        
        f(i,j) = smap_integral(E{i,j}.y,E{i,j}.x) + smap_integral(S{i,j}.y,S{i,j}.x) ...
            - smap_integral(S{i,j}.y,E{i,j}.x) - smap_integral(E{i,j}.y,S{i,j}.x);
        f(i,j) = f(i,j)/(block_width*block_height);
        if f(i,j) > thres
            block(S{i,j}.y:E{i,j}.y-2, S{i,j}.x:E{i,j}.x-2) = 1;
            flag(i,j) = 1;
        end
    end
end


