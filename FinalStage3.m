function [RM2,T] = FinalStage3(RM,Q,T,alpha,s_ind,SC,searchmode)

% Significance Test and Non-maxima Suppression
%
% [USAGE]
% RM2 = FinalStage(RM,Q,lpha,s_ind,SC,searchmode)

%
% [RETURNS]
% RM2   : Resemblance map with bounding boxes at detected objects
%
% [PARAMETERS]
% RM   : Resemblance map
% Q : query image
% alpha     : confidence level
% sind : scale indexes where objects are detected
% SC: scale factors



RM1 = RM;
RM2 = zeros(size(RM));


if max(RM(:)) > 0%0.7
    f_rho = RM(floor(size(Q,1)/2):end-floor(size(Q,1)/2),floor(size(Q,2)/2):end-floor(size(Q,2)/2));
    [E_pdf, ind] = ksdensity(f_rho(:));
    E_cdf = cumsum(E_pdf/sum(E_pdf));
    detection = find(E_cdf > alpha);
    %     figure, plot(ind,E_pdf/sum(E_pdf));
    T_n = ind(detection(1));

    if T_n < 0
        T_n = 0;
    end
    x = size(Q,1);
    y = size(Q,2);

    half_x1 = fix(x/2);
    half_y1 = fix(y/2);

    cnt = 0;
    while max(RM1(:)) > T_n
        %max(RM1(:))
        [x_ind,y_ind] = find(RM==max(RM1( :)));
        cnt = cnt+1;
        if numel(x_ind) > 1
            x_ind = x_ind(1);
            y_ind = y_ind(1);
        end
        %         SC(s_ind(x_ind,y_ind));
        %         [x_ind y_ind];
        if searchmode ~= 0
            scaled_x1 = floor(half_x1/SC(s_ind(x_ind,y_ind)));
            scaled_y1 = floor(half_y1/SC(s_ind(x_ind,y_ind)));
        else
            scaled_x1 = half_x1;
            scaled_y1 = half_y1;
        end
        x_b = x_ind-scaled_x1;
        y_b = y_ind-scaled_y1;
        x_e = x_ind+scaled_x1;
        y_e = y_ind+scaled_y1;


        [x1,x2] = meshgrid(-scaled_y1:scaled_y1,-scaled_x1:scaled_x1);


        kkk = exp(-(0.5/(1)^2) *(x1.^2 + x2.^2));

        if x_b <= 0
            x_b = 1;
            x_e = x_b+size(kkk,1)-1;
        end
        if x_e > floor(size(RM,1))
            x_e = floor(size(RM,1));
            x_b = x_e-size(kkk,1)+1;
        end

        if y_b <= 0
            y_b = 1;
            y_e = y_b+size(kkk,2)-1;
        end
        if y_e > floor(size(RM,2))
            y_e = floor(size(RM,2));
            y_b = y_e-size(kkk,2)+1;
        end
        %         RM2(x_b:x_e,y_b:y_e) = RM1(x_b:x_e,y_b:y_e).*kkk(1:size(RM1(x_b:x_e,y_b:y_e),1),1:size(RM1(x_b:x_e,y_b:y_e),2));



        Overlap{cnt} = zeros(size(RM));
        Overlap{cnt}(max(1,x_b):min(x_e,end),max(1,y_b):min(y_e,end)) = 1;

        if cnt > 1
            for cnt1 = 1:cnt-1
                intersection = Overlap{cnt}.*Overlap{cnt1};
                union = Overlap{cnt}+Overlap{cnt1};
                union = union>0;
                ol(cnt1) = sum(intersection(:))/sum(union(:));
                
%                                 figure(10000),
%                                 subplot(1,2,1),sc(intersection),
%                                 subplot(1,2,2),sc(union),pause;
            end
%             max(ol)
            if max(ol) < 0.01
                RM2 = DrawBox_bold(RM2,max(1,x_b),max(1,y_b),min(x_e,size(RM2,1)),min(y_e,size(RM2,2)),max(RM1(:)));
               T(max(1,x_b):min(x_e,end),max(1,y_b):min(y_e,end)) = imresize(Q,[size(RM2(max(1,x_b):min(x_e,end),max(1,y_b):min(y_e,end)),1),size(RM2(max(1,x_b):min(x_e,end),max(1,y_b):min(y_e,end)),2)]);
            RM2(max(1,x_b):min(x_e,end),max(1,y_b):min(y_e,end)) = RM2(max(1,x_b):min(x_e,end),max(1,y_b):min(y_e,end));
            RM1(max(1,x_b):min(x_e,end),max(1,y_b):min(y_e,end)) = 0;

            else
                RM1(max(1,x_b):min(x_e,end),max(1,y_b):min(y_e,end)) = 0;
            end
            clear ol;
        else
            RM2 = DrawBox_bold(RM2,max(1,x_b),max(1,y_b),min(x_e,size(RM2,1)),min(y_e,size(RM2,2)),max(RM1(:)));
            T(max(1,x_b):min(x_e,end),max(1,y_b):min(y_e,end)) = imresize(Q,[size(RM2(max(1,x_b):min(x_e,end),max(1,y_b):min(y_e,end)),1),size(RM2(max(1,x_b):min(x_e,end),max(1,y_b):min(y_e,end)),2)]);
            RM2(max(1,x_b):min(x_e,end),max(1,y_b):min(y_e,end)) = RM2(max(1,x_b):min(x_e,end),max(1,y_b):min(y_e,end));
            RM1(max(1,x_b):min(x_e,end),max(1,y_b):min(y_e,end)) = 0;

        end


        clear x1 x2 kkk;

    end
end


function img = DrawBox_bold(img, x0,y0,x1,y1,value)

img = func_DrawLine(img, x0, y0, x0, y1, value);
img = func_DrawLine(img, x0+1, y0+1, x0+1, y1-1, value);
img = func_DrawLine(img, x0+2, y0+2, x0+2, y1-2, value);
img = func_DrawLine(img, x0, y1, x1, y1, value);
img = func_DrawLine(img, x0+1, y1-1, x1-1, y1-1, value);
img = func_DrawLine(img, x0+2, y1-2, x1-2, y1-2, value);
img = func_DrawLine(img, x1, y1 ,x1, y0, value);
img = func_DrawLine(img, x1-1, y1-1 ,x1-1, y0+1, value);
img = func_DrawLine(img, x1-2, y1-2 ,x1-2, y0+2, value);
img = func_DrawLine(img, x1, y0, x0, y0, value);
img = func_DrawLine(img, x1-1, y0+1, x0+1, y0+1, value);
img = func_DrawLine(img, x1-2, y0+2, x0+2, y0+2, value);


function Img = func_DrawLine(Img, X0, Y0, X1, Y1, nG)
% Connect two pixels in an image with the desired graylevel
%
% Command line
% ------------
% result = func_DrawLine(Img, X1, Y1, X2, Y2)
% input:    Img : the original image.
%           (X1, Y1), (X2, Y2) : points to connect.
%           nG : the gray level of the line.
% output:   result
%
% Note
% ----
%   Img can be anything
%   (X1, Y1), (X2, Y2) should be NOT be OUT of the Img
%
%   The computation cost of this program is around half as Cubas's [1]
%   [1] As for Cubas's code, please refer
%   http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=4177
%
% Example
% -------
% result = func_DrawLine(zeros(5, 10), 2, 1, 5, 10, 1)
% result =
%      0     0     0     0     0     0     0     0     0     0
%      1     1     1     0     0     0     0     0     0     0
%      0     0     0     1     1     1     0     0     0     0
%      0     0     0     0     0     0     1     1     1     0
%      0     0     0     0     0     0     0     0     0     1
%
%


Img(X0, Y0) = nG;
Img(X1, Y1) = nG;
if abs(X1 - X0) <= abs(Y1 - Y0)
    if Y1 < Y0
        k = X1; X1 = X0; X0 = k;
        k = Y1; Y1 = Y0; Y0 = k;
    end
    if (X1 >= X0) & (Y1 >= Y0)
        dy = Y1-Y0; dx = X1-X0;
        p = 2*dx; n = 2*dy - 2*dx; tn = dy;
        while (Y0 < Y1)
            if tn >= 0
                tn = tn - p;
            else
                tn = tn + n; X0 = X0 + 1;
            end
            Y0 = Y0 + 1; Img(X0, Y0) = nG;
        end
    else
        dy = Y1 - Y0; dx = X1 - X0;
        p = -2*dx; n = 2*dy + 2*dx; tn = dy;
        while (Y0 <= Y1)
            if tn >= 0
                tn = tn - p;
            else
                tn = tn + n; X0 = X0 - 1;
            end
            Y0 = Y0 + 1; Img(X0, Y0) = nG;
        end
    end
else if X1 < X0
        k = X1; X1 = X0; X0 = k;
        k = Y1; Y1 = Y0; Y0 = k;
    end
    if (X1 >= X0) & (Y1 >= Y0)
        dy = Y1 - Y0; dx = X1 - X0;
        p = 2*dy; n = 2*dx-2*dy; tn = dx;
        while (X0 < X1)
            if tn >= 0
                tn = tn - p;
            else
                tn = tn + n; Y0 = Y0 + 1;
            end
            X0 = X0 + 1; Img(X0, Y0) = nG;
        end
    else
        dy = Y1 - Y0; dx = X1 - X0;
        p = -2*dy; n = 2*dy + 2*dx; tn = dx;
        while (X0 < X1)
            if tn >= 0
                tn = tn - p;
            else
                tn = tn + n; Y0 = Y0 - 1;
            end
            X0 = X0 + 1; Img(X0, Y0) = nG;
        end
    end
end
