function [img,LARK] = CompLARK(Inimg,conf)

% Compute Locally Adaptive Regression Kernels
%
% [USAGE]
% [img,LSK] = Computation_LARK_Mar23(img,conf) ;
%
% [RETURNS]
% LSK   : A collection of vectorized LSKs computed from every three pixels
% img   : Downscaled version of input image
%
% [PARAMETERS]
% Inimg   : the input image
% conf    : parameters
%


% Common for all the methods
wsize = conf.Wsize;
colormode = conf.colormode;
interval = conf.interval;

h = conf.h;
alpha = conf.alpha;

if size(Inimg,3) > 1
    Inimg = im2double(colorspace('Lab<-RGB',Inimg));
    if colormode == 1
        for i = 1:3
            LARK{i} =  subfunction(Inimg(:,:,i),wsize,h,interval,alpha) ;
            img{i} = imresize(Inimg(:,:,i),1/interval,'bicubic');
        end
    else
        Inimg = Inimg(:,:,1);
        Inimg = Inimg/std(Inimg(:));
        LARK =  subfunction(Inimg,wsize,h,interval,alpha) ;
        img = imresize(Inimg(:,:,1),1/interval,'bicubic');
    end

else
    Inimg = im2double(Inimg);
    Inimg = Inimg/std(Inimg(:));
    LARK =  subfunction(Inimg,wsize,h,interval,alpha) ;
    img = imresize(Inimg,1/interval,'bicubic');
end

function LARK = subfunction(img,wsize,h,interval,alpha)
[M,N] = size(img);
win = (wsize-1)/2;
%% Pilot estimation of graidents zx,zy

[zx,zy]=gradient(img);

%% Covariance compuatation

zx = padarray(zx,[win,win],'symmetric','both');
zy = padarray(zy,[win,win],'symmetric','both');
C11 = zeros(M,N);
C12 = zeros(M,N);
C22 = zeros(M,N);
tmp = zeros(2,2);
G = zeros(wsize^2,2);
gx = zeros(wsize,wsize);
gy = zeros(wsize,wsize);
K = fspecial('disk', win);
K = K ./ K(win+1, win+1);
len = sum(K(:));
for i = 1 : interval: M
    for j = 1 :interval: N
        gx = zx(i:i+wsize-1, j:j+wsize-1).*K;
        gy = zy(i:i+wsize-1, j:j+wsize-1).*K;
        G = [gx(:), gy(:)];
        
            [u s v] = svd(G,'econ');
            S1 = (s(1,1) + 1) / (s(2,2) + 1);
            S2 = 1/S1;
           tmp = (S1 * v(:,1) * v(:,1).' + S2 * v(:,2) * v(:,2).')*((s(1,1) * s(2,2) + 0.0000001)/len)^alpha;
            C11(i,j) = tmp(1,1);
            C12(i,j) = tmp(1,2);
            C22(i,j) = tmp(2,2);
        
    end
end

C11 = C11(1:interval:end,1:interval:end);
C12 = C12(1:interval:end,1:interval:end);
C22 = C22(1:interval:end,1:interval:end);
% 

[M,N] = size(C11);

C11 = padarray(C11,[win,win],'symmetric','both');
C12 = padarray(C12,[win,win],'symmetric','both');
C22 = padarray(C22,[win,win],'symmetric','both');


%% Polynonimal basis function
[x2,x1] = meshgrid(-win:win,-win:win);
x2 = x2;
x1 = x1;

x12 = 2*x1.*x2;
x11 = x1.^2;
x22 = x2.^2;

x1x1 = reshape(repmat(reshape(x11,[1 wsize^2]),M*N,1),[M,N,wsize wsize]);
x1x2 = reshape(repmat(reshape(x12,[1 wsize^2]),M*N,1),[M,N,wsize wsize]);
x2x2 = reshape(repmat(reshape(x22,[1 wsize^2]),M*N,1),[M,N,wsize wsize]);

%% LSK feature computation
LARK = zeros(M,N,wsize,wsize);
for i = 1:wsize
    for j = 1:wsize
        LARK(:,:,i,j) = C11(i:i+M-1,j:j+N-1).*x1x1(:,:,i,j)+ C12(i:i+M-1,j:j+N-1).*x1x2(:,:,i,j)+ C22(i:i+M-1,j:j+N-1).*x2x2(:,:,i,j);
    end
end

    LARK = exp(-(LARK)/h);


%% Normalization
LARK = reshape(LARK,[M N wsize^2]);
LARK = LARK./repmat(sum(LARK,3),[1 1 wsize^2]);



function DrawNonoverlappingLSKs(LSK,M,N,wsize_skr)
LSK = reshape(LSK,[M N wsize_skr^2]);
win_skr = (wsize_skr-1)/2;
for i = 1:wsize_skr:M-wsize_skr+1
    for j = 1:wsize_skr:N-wsize_skr+1
        Output_LSK(i:i+wsize_skr-1,j:j+wsize_skr-1) = reshape(LSK(i+win_skr,j+win_skr,:),[wsize_skr,wsize_skr]);
    end
end
figure(100),
subplot(1,1,1),imagesc(Output_LSK), axis image, axis off, colorbar;
