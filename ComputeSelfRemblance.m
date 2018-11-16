function SR = ComputeSelfRemblance(varargin)

% Compute Self-Resemblance

% [RETURNS]
% SR   : Saliency Map
%
% [PARAMETERS]
% img   : Input image
% LARK  : A collection of LARK descriptors
% param : parameters


% Nov 2, 2018 : Siddhant Katyan

img = varargin{1};
LARK = varargin{2};
param = varargin{3};


[M,N] = size(img);
win_N = (param.N-1)/2;
win_L = (param.L-1)/2;

% To avoid edge effect, we use mirror padding.
for k = 1:1
    for i = 1:size(LARK{1},3)
        LARK1{k}(:,:,i) = padarray(LARK{k}(:,:,i),[win_L,win_L],'symmetric','both');
        LARK2{k}(:,:,i) = padarray(LARK1{k}(:,:,i),[win_L,win_L],'symmetric','both');
    end
end

SR = zeros(M*N,1);

% Precompute Norm of center matrices
Norm_C = zeros(M*N,1);
for m = 1:2:param.L
    for n = 1:2:param.L
        Center{m,n} = [reshape(LARK1{1}(m:m+M-1,n:n+N-1,:),[M*N size(LARK{1},3)])];
        Norm_C = Norm_C + sum(Center{m,n}.^2,2);
    end
end

Norm_C = sqrt(Norm_C);

for i = 1:param.N
    for j = 1:param.N
        % compute Norm of surrounding matrices
        Norm_S = zeros(M*N,1);
        temp = zeros(M*N,1);
        for m = 1:2:param.L
            for n = 1:2:param.L
                Surround = [reshape(LARK2{1}(i+m-1:i+m-1+M-1,j+n-1:j+n-1+N-1,:),[M*N size(LARK{1},3)])];
                Norm_S = Norm_S + sum(Surround.^2,2);  
                temp = temp + sum(Center{m,n}.*Surround,2); % compute inner product between a center and surrounding matrices
            end
        end
        Norm_S = sqrt(Norm_S);
        SR(:,1) = SR(:,1) + exp( (-1+temp./(Norm_C.*Norm_S))/param.sigma^2); % compute self-resemblance using matrix cosine similarity
    end
end
SR = reshape(1./SR,[M,N]); %Final saliency map values
end
