% ns = [10000];
% for it=1:length(ns);
% clearvars -except ns it
m=4;
num_items=m^2;

act = [0.001 0.005 0.01 0.05 0.1 0.5];
kwta_options = [0.001 0.005 0.01 0.05 0.1];

act = [0.1];
kwta_options = [0.01];

%for it=1:length(kwta_options);
 %   for in_it=1:length(act);
 %clearvars -except act m num_items it binary_res order_res binary_sim_res order_sim_res binary_rand order_rand kwta_options c in_it
for sim_val=1:1000;
 clearvars catted
names = randn(m,100);
c=cov(names');

objects = randn(m,100);
c2=cov(objects');


% ctr=1;
% for k=1:4;
%     for j=1:4;
%        catted(ctr,:)=[names(k,:) objects(j,:)];
%        idx(ctr,:) = [k,j];
%        ctr=ctr+1;
%     end
% end

ctr=1;
for k=1:4;
    for j=1:4;
       catted(ctr,:)=names(k,:) + objects(j,:);
       idx(ctr,:) = [k,j];
       ctr=ctr+1;
    end
end

catted = [catted; names; objects];

%c   = c-min(c)/(max(c)-min(c));
    these_objs = randsample(m,m-1);
    these_words = randsample(m,m-1);
    o_complement = setdiff([1:m], these_objs);
    w_complement = setdiff([1:m], these_words);
 
%[U,V]=eig(c);
% # of active nodes in sparse layer.
pct_act = act;

%size of expanded layer.
layer2=10000;

%# of observations
n=16;

%# of nodes per variable.
d=15;

%if it's kwta, what's k?
%kwta_options=[1 5 10 25 50 100 250 500];

weights_layer2 = zeros(200, layer2);

%expand and sparsify w/ binary random projection.
% for j=1:200;
% sprse(:,j)=randsample(layer2,layer2*pct_act);
% weights_layer2(j,sprse(:, j))=1;
% end
% 
% [U,V]=eig(cov(catted'));
% %for eigenvector case.
% %output = U*weights_layer2;
% 
% output = catted*weights_layer2;
% 
% kwta=kwta_options(it)*layer2;
% % %kwta
% for j=1:num_items;
%     [a b] = sort(abs(output(j,:)), 'descend');
%     output(j,b(kwta:end))=0;
%     output(j,b(1:kwta))=1;
% end

% %kwta
% for j=1:num_items;
%     [a b] = sort(output(j,:), 'descend');
%     output(j,b(kwta:end))=0;
%     output(j,b(1:kwta))=1;
% end

%kernels
%c = cov(output');
%c = corr(output');
%c= (output*output')/norm(output);

%distance measures
% %hamming distance
% out=1-pdist(output,  'hamming');
% c=squareform(out);
% c=c+eye(m^2);
% c = (c - min(c)) ./ (max(c) - min(c)) ;
% %euclidean distance
% out=pdist(output,  'euclidean');
% c=squareform(out);
% c=c+eye(m^2);

%jaccard distance
%  out=1-pdist(output,  'jaccard');
%  c=squareform(out);
%  c=c+eye(m^2);
%  c = (c - min(c)) ./ (max(c) - min(c)) ;
%gaussian.
% nsq=sum(output.^2,2);
% K=bsxfun(@minus,nsq,(2*output)*output.');
% K=bsxfun(@plus,nsq.',K);
% K=exp(-K);
% c=K;


%c=normr(cov(catted'));

c=cov(catted');

%c=corr(catted');

% output=output';
% sigma=1;
% l=5;
% for j=1:num_items;
%      for k=1:num_items;     
% %c(j,k)=1-pdist([output(j,:); output(k,:)], 'hamming');
% c(j,k)=exp(-norm(output(j,:)-output(k,:)' ).^2)
% %c(j,k) = sigma^2*(exp( (output(j,:)-output(k,:)' ) / (2*l^2) ));
%      end
%  end
% % 

%exp(-norm( X(i,:) - X(j,:) ))^2);


%c=zscore(c);

%find(idx(:, 1)==these_words(p,:))
    for p=1:m-1;
    sub_pre(p) = intersect(find(idx(:, 1)==these_words(p,:)), find(idx(:, 2)==these_objs(p,:)));
    end
    target = [w_complement, o_complement];
    target_idx = find(idx==[w_complement o_complement]);
    target_idx = intersect(find(idx(:, 1)==w_complement), find(idx(:, 2)==o_complement));
    
%     I = eye(m-1);
%     for p=1:m-2;
%         for p2=p+1:m-1;
%             csub(ctr,1) = c(sub_pre(p), sub_pre(p2));
%             I(p,p2)=c(sub_pre(p), sub_pre(p2));
%             ctr=ctr+1;
%         end
%     end
%    Iout= triu(I)+triu(I,1)'  
   
 % I = eye(m);
 % I(1:m-1,1:m-1) = c(sub_pre, sub_pre);
 % I(:, m)=0;
 % I(m,:)=0;
 
%  [u,v]=eig(c);
%  [a,b]=max((u'*eye(m^2)));
%  eig_items(j) = sum(b);
 candidates = setdiff([1:num_items], sub_pre);
 %candidates = 1:24;
  %search for alternatives.
  for p=1:length(candidates);    
      %Itemp=I;
      cur_list = [sub_pre candidates(p)];
      present_absent = zeros(m^2+(m*2), m^2+(m*2));
      present_absent(cur_list, cur_list)=1;
     % for j=1:m-1;
       Itemp = c(cur_list, cur_list);
       %det_results(p) = det(Itemp);
       
       Iout= triu(Itemp);
       %det_results(p) =prod(diag(zscore(Iout)));
       %det_results(p) =prod(diag(normr(Iout)));
       %det_results(p) =sum(diag(normr(Iout)));
       
       out =  diag(c.*present_absent) ./ sum(c .* present_absent')+0.0001;
       Iout= out(cur_list, cur_list);
       
       %det_results(p) =sum(diag(Iout));
       det_results(p) =prod(diag(Iout));
       df
       
       %det_results(p) =det(Iout);
     % end
   
%    [V,D]=eig(Iout);
%    X = (V*diag(exp(diag(D)))/V);
%    det_results(p)=trace(expm(X));
   %sim_results(p) = mean(mean(corr((Iout))));
   
   
  end
  
  %[ares,bres]= min(det_results);
  [ares,bres]= max(det_results);
  [trsh srtd_res]= sort(det_results);
  binary_res(sim_val) = candidates(bres)==target_idx;
  %order_res(it, in_it,sim_val) = find(candidates(srtd_res)==target_idx);
%   [ares,bres]= max(sim_results);
%   [trsh srtd_res]= sort(sim_results);
%   binary_sim_res(sim_val) = candidates(bres)==target_idx;
%   %order_sim_res(it, in_it,sim_val) = find(candidates(srtd_res)==target_idx);
%   
%   binary_rand(sim_val) = randsample(length(candidates),1)==target_idx;
  %order_rand(it, in_it, sim_val)=find(randsample(length(candidates), 13)==target_idx);
end

%    end 
%end

%
% 
% m=4;
% 
% 
% %no severing eigenvectors.
% for j=1:1000;
% names = randn(m,100);
% c=cov(names');
% [L.U,L.D]=eig(c);
% [a,b]=max((L.U'*eye(m)));
% items = randsample(4, 3);
% %[a,b]=max((L.U*eye(m)));
% eig_items = sum(b);
% names = randn(m,1000);
% c=corr(names');
% [L.U,L.D]=eig(c);
% [a,b]=max((L.U'*eye(m)));
% %[a,b]=max((L.U*eye(m)));
% eig_items(j) = sum(b);
% % if b(1,1) == (1) && b(1,2) == 1; eig_items(j) = 1;
% % elseif b(1,1) == 1 && b(1,2) ==2; eig_items(j) = 2;
% % elseif b(1,1) == 2 && b(1,2) == 1; eig_items(j)=3;
% % elseif b(1,1) == 2 && b(1,2) == 2; eig_items(j)=4;   
% % end
% eig_samp(m,j)=length(unique(b));
% eig_samp_collisions(m,j)=(m-eig_samp(m,j))/m;
% end

