function [F_Q,F_T] = PCAfeature(Q,T,wsize,W_Q,W_T,flag)
% [PCA]
% Dimension reduction step
%
% [USAGE]
% [F_T,F_Q,D] = PCA(Q,T,wsize,W_Q,W_T)
%
% [RETURNS]
% F_T   : Target Feature Matrix
% F_Q   : Query Feature Matrix

%
% [PARAMETERS]
% Q     : Query image
% T     : Target image
% wsize : a kernel size
% W_Q   : Query LSK matrix
% W_T   : Target LSK matrix
% flag : show features if 1 or do not show if 0




%% Principal Component Analysis
rsize = (wsize-1)/2;
LSK = W_Q(rsize+1:4:size(Q,1)-rsize,rsize+1:4:size(Q,2)-rsize,:);
WWW = reshape(LSK,[size(LSK,1)*size(LSK,2),size(LSK,3)*size(LSK,4)]);
Temp = sum(WWW,1);
cnt = size(WWW,1);

meanTemp = repmat(Temp',[1 cnt])./cnt;
Temp = Temp';

% compose covariance matrix of query LARK features and perform PCA
[v,d] = eig((WWW'-meanTemp)*(WWW'-meanTemp)');
[tmp l] = sort(diag(d),'descend');
v = v(:,l);
d = diag(tmp);
% remain top eigenvalues so that we can preserve 80% of energy
for i = 1:size(d,1)
    trace(d(1:i,1:i))./trace(d);
    if trace(d(1:i,1:i))./trace(d) > .8
        energy(i) = trace(d(1:i,1:i))./trace(d);
        valid = i;
        break;
    end
end

if flag == 1
    figure(10000);
    for i = 1:valid
        v1(:,:,i) = reshape(v(:,i),[wsize,wsize]);
        subplot(1,valid,i), imagesc(imresize(squeeze(v1(:,:,i)),5,'bicubic')), axis off, axis image,label( imresize(squeeze(v1(:,:,i)),5,'bicubic'),['eigen vector:' num2str(i)]);
    end
end

meanMat = repmat((Temp./cnt)',size(T,1)*size(T,2),1);

%% Obtain F_T and F_Q by projecting W_Q and W_T onto Projection Space V
%% (eveigenvector space)
F_T = (reshape(W_T,[size(T,1)*size(T,2) size(W_T,3)])-meanMat)*v(:,1:valid);
F_Q = (reshape(W_Q,[size(Q,1)*size(Q,2) size(W_Q,3)])-meanMat(1:size(Q,1)*size(Q,2),:))*v(:,1:valid);

F_T = reshape(F_T,[size(T,1) size(T,2) valid]);
F_Q = reshape(F_Q,[size(Q,1) size(Q,2) valid]);

if flag == 1
    figure(100);
    for i = 1:valid
        subplot(2,valid,i), imagesc(F_T(:,:,i)),colormap(gray), axis off, axis image; label(F_T(:,:,i),['F_T' num2str(i)]);
        subplot(2,valid,i+valid), imagesc(F_Q(:,:,i)),colormap(gray), axis off, axis image; label(F_Q(:,:,i),['F_Q' num2str(i)]);
    end
end

function label(im, str)
text(size(im, 2)/2, size(im, 1)+20, str, ...
    'Interpreter', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
return

