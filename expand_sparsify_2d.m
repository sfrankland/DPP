%clear
%load 1x100_T.mat
%load 5x5_2d_space_laplace.mat
% # of active nodes in sparse layer.
pct_act = 0.2;

%size of expanded layer.
layer2=10000;

%# of observations
n=2000;

%# of nodes per variable.
d=15;

% for n=2:100;
% [spect F T P,Fc,Tc]=spectrogram(U(:,n), 50);
% [a,b]=max(max(real(spect)'));
% idx(n)=a;
% end


corner = ...
[0.6 .2 .2 0 0 0 0 0 0 0;
.2 .4 .2 .2 0 0 0 0 0 0;
0.2 .2 .4 .2 0.2 0 0 0 0 0;
0 0.2 .2 .4 .2  0 0 0 0;
0 0 0 .2 .4 .2 0 0 0 0;
0 0 0 0 .2 .4 .2 0 0 0;
0 0 0 0 0 0.2 .4 .2 0 0;
0 0 0 0 0 0 .2 .4 .2 0;
0 0 0 0 0 0 0 .2 .4 .2;
0 0 0 0 0 0 0 0 .2 .6];

side = diag([.2 .2 .2 .2 .2 .2 .2 .2 .2 .2]);
mid = ...
[.4 .2 0 0 0 0 0 0 0 0;
.2 .2 .2 0 0 0 0 0 0 0;
0 .2 .2 .2 0 0 0 0 0 0;
0 0 .2 .2 .2 0 0 0 0 0;
0 0 0 .2 .2 .2 0 0 0 0;
0 0 0 0 .2 .2 .2 0 0 0;
0 0 0 0 0 .2 .2 .2 0 0;
0 0 0 0 0 0 .2 .2 .2 0;
0 0 0 0 0 0 0 .2 .2 .2;
0 0 0 0 0 0 0 0 .2 .4];
z = zeros(10)

L = [
corner,side,z,z,z,z,z,z,z,z;
side,mid,side,z,z,z,z,z,z,z;
z,side,mid,side,z,z,z,z,z,z;
z,z,side,mid,side,z,z,z,z,z;
z,z,z,side,mid,side,z,z,z,z;
z,z,z,z,side,mid,side,z,z,z;
z,z,z,z,z,side,mid,side,z,z;
z,z,z,z,z,z,side,mid,side,z;
z,z,z,z,z,z,z,side,mid,side;
z,z,z,z,z,z,z,z,side,corner];


[U,V]=eig(L);
%if it's kwta, what's k?
kwta_options=[1 5 10 25 50 100 250 500];
kwta=0.005*layer2;
%for it=1:20;
%%%%%random walk context vector.%%%%%%%%%%%%
init = 0+randn(1,50);
for j=2:100;
   init(j,:) = init(j-1,:)+randn(1,50);
end
init=init';
%[U,V]=eig(cov(init));
%[U,V,P] = qr(cov(init));

%MEC
%load('/users/stevenmf/Desktop/transitions/1x100_T.mat')

%sensory LEC.
%%%%%random walk context vector.%%%%%%%%%%%%
% sensory = 0+randn(1,150);
% for j=2:150;
%    sensory(j,:) = sensory(j-1,:)+rand(1,50);
% end
% sensory=sensory';

%data = [zscore(U);zscore(init)]

[U,V]=eig(L);

Ushift1 = zscore(circshift(U,1));
Ushift2 = zscore(circshift(U,2));
Ushift3 = zscore(circshift(U,3));
Ushift4 = zscore(circshift(U,4));
Ushift5 = zscore(circshift(U,5));
Ushift6 = zscore(circshift(U,6));
Ushift7 = zscore(circshift(U,7));
Ushift8 = zscore(circshift(U,8));

%data = [zscore(U(:,2:10)) Ushift1(:, 2:10) Ushift2(:, 2:10) Ushift3(:, 2:10) Ushift4(:, 2:10) Ushift5(:, 2:10) Ushift6(:, 2:10) Ushift7(:, 2:10) Ushift8(:, 2:10)]';
data = [zscore(U(:,20:60)) Ushift1(:, 20:60) Ushift2(:, 20:60) Ushift3(:, 20:60) Ushift4(:, 20:60) Ushift5(:, 20:60) Ushift6(:, 20:60) Ushift7(:, 20:60) Ushift8(:, 20:60)]';

%data = [zscore(U(:,10:20)) Ushift1(:, 10:20) Ushift2(:, 10:20) Ushift3(:, 10:20) Ushift4(:, 10:20) Ushift5(:, 10:20)  Ushift6(:, 10:20) Ushift7(:, 10:20) Ushift8(:, 10:20)]';
%data = [zscore(U(:,1:20);zscore(init)]

%data = init;
%data = [zscore(U(:, 50:60)) Ushift1(:,50:60) Ushift2(:, 50:60) Ushift3(:, 50:60)  Ushift4(:, 50:60)  Ushift5(:, 50:60)  Ushift6(:, 50:60)  Ushift7(:, 50:60)  Ushift8(:,50:60)]';
%data = U;

%data = init;
%data=zscore(init);
%data = [zscore(randn(110,100))];
data=L;
[nfeat nobs] =size(data);
%data=zscore(sensory')';
n_states = length(L);
gamma = [0.95];
%for kloop=1:length(kwta_options);
 %   kwta=kwta_options(kloop);
% for p=1:length(gamma);
% cur_gamma=gamma(p);
% %SR
% init = inv(eye(n_states) - cur_gamma *T);
% [U,V,W]= eig(init);
% 

% %expand and sparsify w/ binary random projection.
for j=1:nfeat;
sprse(:,j)=randsample(layer2,layer2*pct_act);
weights_layer2(j,sprse(:, j))=1;
end
%for eigenvector case.
%data = L;
%output = data'*weights_layer2;

output = data'*weights_layer2;

% %kwta
for j=1:nobs;
    [a b] = sort(output(j,:), 'descend');
    output(j,b(kwta:end))=0;
    output(j,b(1:kwta))=1;
end
% 


%transpose random weights and reconstruct.
%init_recon = output*weights_layer2';

% for k=1:100;
%     for j=1:100;
% hdist_recon(k,j) = pdist([init_recon(k,:); init(j,:)], 'euclidean');
%     end
% end
% 
% [a,b]=min(hdist_recon);
% recon_perf(kloop, it) = mean(b==1:100);
% 


for k=1:nobs;
    for j=1:nobs;
hdist_vectors(k,j) = pdist([output(k,:); output(j,:)], 'hamming');
    end
end

% ctr=1;
% for j1=1:46;
%     for j2=j1+1:47;
%       for j3=j2+1:48;
%           for j4=j3+1:49;
%               for j5=j4+1:50;
% matrix(ctr, j1)=1;
% matrix(ctr, j2)=1;
% matrix(ctr, j3)=1;
% matrix(ctr, j4)=1;
% matrix(ctr, j5)=1;
% ctr=ctr+1;
%               end
%           end
%       end
%     end
% end

%repel within high-dimensional DG code.
%randomly sample a point in high dimensional space.
%what is the learned kernel for eigendecomposition?


% % sample a set Y from a dpp.  L is a decomposed kernel, and k is (optionally)
% % the size of the set to return.
% if ~exist('k2','var')  
%   % choose eigenvectors randomly
%   D = U ./ (1+U);
%   v = find(rand(length(D),1) <= D);
% else
%   % k-DPP
%   v = sample_k(U,k2);
% end
% 
% k = length(v);
% V = U(:,v);
% 
% % iterate
% Y = zeros(k,1);
% for i = k:-1:1
%   % compute probabilities for each item
%   P = sum(V.^2,2);
%   P = P / sum(P);
%   % choose a new item to include
%   Y(i) = find(rand <= cumsum(P),1);
%   % choose a vector to eliminate
%   j = find(V(Y(i),:),1);
%   Vj = V(:,j);
%   V = V(:,[1:j-1 j+1:end]);
%   % update V
%   V = V - bsxfun(@times,Vj,V(Y(i),:)/Vj(Y(i)));
% 
%   % orthogonalize
%   for a = 1:i-1
%     for b = 1:a-1
%       V(:,a) = V(:,a) - V(:,a)'*V(:,b)*V(:,b);
%     end
%     V(:,a) = V(:,a) / norm(V(:,a));
%     size(V)
%   end  
%  end
% 
% Y = sort(Y);






for k=2:100;
hdist_vectors(k) = pdist([output(k,:); output(k-1,:)], 'hamming');
actives = find(output(k,:)==1);
hdist_vectors_activeonly(k) = pdist([output(k,actives); output(k-1,actives)], 'hamming');
end

% %%%%similarity kernel and eigendecomposition of that data%%%%
% r = corr(init);
% [U,V]=eig(r);
% eig_output = U*weights_layer2;
% %kwta
% for j=1:100;
%     [a b] = sort(eig_output(j,:), 'descend');
%     eig_output(j,b(kwta:end))=0;
%     eig_output(j,b(1:kwta))=1;
% end
% for k=2:100;
% hdist_eigs(k) = pdist([eig_output(k,:); eig_output(k-1,:)], 'hamming');
% actives = find(eig_output(k,:)==1);
% hdist_eigs_activeonly(k) = pdist([eig_output(k,actives); eig_output(k-1,actives)], 'hamming');
% end
% 
% %for p1=1:100;
% %randvar = randn(100,10);
% randvar = U(randsample(100,100),:);
% rand_output = randvar*weights_layer2;
% %kwta
% for j=1:100;
%     [a b] = sort(rand_output(j,:), 'descend');
%     rand_output(j,b(kwta:end))=0;
%     rand_output(j,b(1:kwta))=1;
% end
% for k=2:100;
% hdist_rand(k) = pdist([rand_output(k,:);rand_output(k-1,:)], 'hamming');
% actives = find(rand_output(k,:)==1);
% hdist_rand_activeonly(k) = pdist([rand_output(k,actives); rand_output(k-1,actives)], 'hamming');
% end
% %end
% 
% 
% loop_results_eigs(kloop,p) = mean(hdist_eigs);
% loop_results_SR(kloop,p) = mean(hdist_vectors);
% loop_results_eigs_active(kloop,p) = mean(hdist_eigs);
% loop_results_SR_active(kloop,p) = mean(hdist_vectors);
% loop_results_rand(kloop,p) = mean(hdist_rand);
% loop_results_rand_active(kloop,p) = mean(hdist_rand_activeonly);
% end
% end
% 
% 
% lambda_lower = [1 10 20 30 40 50 60 70 80 90];
% 
% for p1=1:length(lambda_lower);
%   cur_lambda_lower = lambda_lower(p1);  
%   
% %%%%low freq eigenvectors%%%%
% r = corr(init);
% [U,V]=eig(r);
% eig_output = U(:,cur_lambda_lower:cur_lambda_lower+10)*weights_layer2(cur_lambda_lower:cur_lambda_lower+10,:);
% %kwta'
% for j=1:100;
%     [a b] = sort(eig_output(j,:), 'descend');
%     eig_output(j,b(kwta:end))=0;
%     eig_output(j,b(1:kwta))=1;
% end
% for k=2:100;
% hdist_inputs(k,p1) = pdist([U(k,cur_lambda_lower:cur_lambda_lower+10); U(k-1,cur_lambda_lower:cur_lambda_lower+10)], 'euclidean');
% hdist_eigs(k, p1) = pdist([eig_output(k,:); eig_output(k-1,:)], 'hamming');
% actives = find(eig_output(k,:)==1);
% hdist_eigs_activeonly(k,p1) = pdist([eig_output(k,actives); eig_output(k-1,actives)], 'hamming');
% end
% end
% 
% 
% lambda_lower = [1 10 20 30 40 50 60 70 80 90];
% 
% for p1=1:length(lambda_lower);
%   cur_lambda_lower = lambda_lower(p1);  
%   
% %%%%low freq eigenvectors%%%%
% r = corr(init);
% [U,V]=eig(r);
% eig_output = U(:,cur_lambda_lower:cur_lambda_lower+10)*weights_layer2(cur_lambda_lower:cur_lambda_lower+10,:);
% %kwta'
% for j=1:100;
%     [a b] = sort(eig_output(j,:), 'descend');
%     eig_output(j,b(kwta:end))=0;
%     eig_output(j,b(1:kwta))=1;
% end
% for k=2:100;
% hdist_inputs(k,p1) = pdist([U(k,cur_lambda_lower:cur_lambda_lower+10); U(k-1,cur_lambda_lower:cur_lambda_lower+10)], 'euclidean');
% hdist_eigs(k, p1) = pdist([eig_output(k,:); eig_output(k-1,:)], 'hamming');
% actives = find(eig_output(k,:)==1);
% hdist_eigs_activeonly(k,p1) = pdist([eig_output(k,actives); eig_output(k-1,actives)], 'hamming');
% end
% end



%for different distances.
%lambda_lower = [1 10 20 30 40 50 60 70 80 90];
lambda_lower = [2 11 21 31 41 51 61 71 81 91];
for p1=1:length(lambda_lower);
  cur_lambda_lower = lambda_lower(p1);

%%%%low freq eigenvectors%%%%
r = corr(init);
[U,V]=eig(r);
eig_output = U(:,cur_lambda_lower:cur_lambda_lower+9)*weights_layer2(cur_lambda_lower:cur_lambda_lower+9,:);
%kwta'
for j=1:100;
    [a b] = sort(eig_output(j,:), 'descend');
    eig_output(j,b(kwta:end))=0;
    eig_output(j,b(1:kwta))=1;
end
%sample different distances.
v=1:100;
% Absolute pairwise diifferences
dv = abs(bsxfun(@minus,v,v'));
%hdist_inputs(k,p1) = pdist([U(k,cur_lambda_lower:cur_lambda_lower+10); U(k-n,cur_lambda_lower:cur_lambda_lower+10)], 'euclidean');

%loop through and evaluate hamming for different distances.
for n=1:90;
   ctr=1;
   [a,b]=find(dv==n);
   for j=1:length(a); 
    hdist_inputs(ctr)=pdist([U(a(j),cur_lambda_lower:cur_lambda_lower+9); U(b(j),cur_lambda_lower:cur_lambda_lower+9)], 'euclidean');
    hdist_eigs(ctr) = pdist([eig_output(a(j),:); eig_output(b(j),:)], 'hamming');
    ctr=ctr+1;
   end
    avg_hdist_DG(p1,n) = mean(hdist_eigs);
    clear hdist_eigs
    clear hdist_inputs
end
end







for p1=1:length(lambda_lower);
  cur_lambda_lower = lambda_lower(p1);  
  
%%%%low freq eigenvectors%%%%
r = corr(init);
[U,V]=eig(r);
eig_output = U(:,cur_lambda_lower:cur_lambda_lower+10)*weights_layer2(cur_lambda_lower:cur_lambda_lower+10,:);
%kwta'
for j=1:100;
    [a b] = sort(eig_output(j,:), 'descend');
    eig_output(j,b(kwta:end))=0;
    eig_output(j,b(1:kwta))=1;
end

for k=2:100;
hdist_inputs(k,p1) = pdist([U(k,cur_lambda_lower:cur_lambda_lower+10); U(k-n,cur_lambda_lower:cur_lambda_lower+10)], 'euclidean');
hdist_eigs(k, p1) = pdist([eig_output(k,:); eig_output(k-1,:)], 'hamming');
actives = find(eig_output(k,:)==1);
hdist_eigs_activeonly(k,p1) = pdist([eig_output(k,actives); eig_output(k-1,actives)], 'hamming');
end
end



ivar = repmat(eye(10),10,1);
i_output = ivar(:, 1:10)*weights_layer2(1:10,:);
%kwta
for j=1:100;
    [a b] = sort(i_output(j,:), 'descend');
    i_output(j,b(kwta:end))=0;
    i_output(j,b(1:kwta))=1;
end
for k=2:100;
hdist_i(p1) = pdist([i_output(k,:);i_output(k-1,:)], 'hamming');
end




