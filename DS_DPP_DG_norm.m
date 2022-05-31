% ns = [10000];
% for it=1:length(ns);
% clearvars -except ns it
m=2;
num_items=m^2;

act = [0.1];
%kwta_options = [0.01 0.05 0.1];

neuron_options = [10000];

%act = [1];
kwta_options = [0.1];
for it=1:length(neuron_options);
    for in_it=1:length(act);
 clearvars -except act m num_items it binary_res order_res binary_sim_res order_sim_res binary_rand order_rand kwta_options c in_it neuron_options binary_res_den
for sim_val=1:100;

names = randn(m,100);
c=cov(names');

objects = randn(m,100);
c2=cov(objects');



ctr=1;
for k=1:m;
    for j=1:m;
       catted(ctr,:)=[names(k,:) objects(j,:)];
       idx(ctr,:) = [k,j];
       ctr=ctr+1;
    end
end

%c   = c-min(c)/(max(c)-min(c));
    these_objs = randsample(m,m-1);
    these_words = randsample(m,m-1);
    o_complement = setdiff([1:m], these_objs);
    w_complement = setdiff([1:m], these_words);
 
[U,V]=eig(c);
% # of active nodes in sparse layer.
pct_act = act(in_it);

%size of expanded layer.
layer2=neuron_options(it);
%layer2=10000;
%# of observations
n=16;

%# of nodes per variable.
d=15;

%if it's kwta, what's k?
%kwta_options=[1 5 10 25 50 100 250 500];

weights_layer2 = zeros(200, layer2);

%expand and sparsify w/ binary random projection.
for j=1:200;
sprse(:,j)=randsample(layer2,layer2*pct_act);
weights_layer2(j,sprse(:, j))=1;
end

[U,V]=eig(cov(catted'));
%for eigenvector case.
%output = U*weights_layer2;

output = catted*weights_layer2;

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



%sdf
%kernels
%c = cov(output');
%c = corr(output');
%c= (output*output')/norm(output);



%distance measures

% %hamming distance
out=1-pdist(output,  'hamming');
c=squareform(out);
c=c+eye(m^2);
%sdf
c = (c - min(c)) ./ (max(c) - min(c)) ;

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
   
   %search for alternatives.
   
%   I = eye(m);
   candidates = setdiff([1:num_items], sub_pre);
%   
%   for p=1:length(candidates);    
%       Itemp=I;
%       for j=1:m-1;
%           Itemp(j,m) = c(candidates(p), sub_pre(j));
%       end
%       
%      for p1=1:m-2;
%         for p2=p+1:m-1;
%             csub(ctr,1) = c(sub_pre(p1), sub_pre(p2));
%             Itemp(p1,p2)=c(sub_pre(p1), sub_pre(p2));
%             ctr=ctr+1;
%         end
%      end
%     
%    Iout= triu(Itemp)+triu(Itemp,1)';
%    det_results(p) = sum(diag(normr(Iout)));
%    %det_results(p) = det(Iout);
%    sim_results(p) = mean(mean(corr((Iout))));
   
norm_1_output = (output - min(output) )./ (max(output) - min(output));
output = (output - min(output) )./ (max(output) - min(output));
%lateral_connections = zeros(layer2, layer2); 
%c=corr(catted');
DG = output(sub_pre,:)';
%DG=c(sub_pre,:)';
%outer product for normalization.
lateral_connections = (DG*DG');
s_lat = sum(lateral_connections);
%lateral_connections = (lat - min(lat) )./ (max(lat) - min(lat));
%norm_out = layer2 ./ (lateral_connections + 0.001) ;
%norm_out = output./s_lat;
norm_out = (catted*weights_layer2)./(s_lat+0.001);
%norm_out = (output)./(s_lat+0.001);
%norm_out = output+s_lat;
%max of max
[ares bres]=max(sum(norm_out'));

%mean of max
%[ares bres]=max(mean(norm_out'));

%norm_out = output ./ (output + (sum(lateral_connections) + 0.001));
%[ares bres]=max(sum(cov(output')));
%[ares bres]=max(sum(abs(norm_out')));
%[ares bres]=max(sum((output/norm_out)'));
%norm_out = layer2 ./ (lateral_connections + 0.001) ;

%norm_out = eye(m^2) ./ (lateral_connections + 0.001) ;
%select max normalized value.
%[ares,bres]=max(diag(norm_out));
% % %kwta
% for j=1:num_items;
%     [a b] = sort(abs(norm_out(j,:)), 'descend');
%     norm_out(j,b(kwta:end))=0;
%     norm_out(j,b(1:kwta))=1;
% end
  
%   [ares,bres]= max(det_results);
%   [trsh srtd_res]= sort(det_results);
  binary_res(it,in_it, sim_val) = bres==target_idx;
  %binary_res(it,  in_it, sim_val) = candidates(bres)==target_idx;
  %order_res(it, in_it,sim_val) = find(candidates(srtd_res)==target_idx);
  
%   [ares,bres]= max(sim_results);
%   [trsh srtd_res]= sort(sim_results);
%   binary_sim_res(it, in_it, sim_val) = candidates(bres)==target_idx;
%   %order_sim_res(it, in_it,sim_val) = find(candidates(srtd_res)==target_idx);
%   
   binary_rand(it, in_it, sim_val) = randsample(length(candidates),1)==target_idx;
  %order_rand(it, in_it, sim_val)=find(randsample(length(candidates), 13)==target_idx);
  
end

    end 
end


%project.
%before binarizing, normalize in high d.

%c, k,v



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

