
m=2;
num_items=m*m;

act = [0.001 0.005 0.01 0.05 0.1 0.5];

act = [0.1];
kwta_options = [0.001 0.005 0.01 0.05 0.1];

for it=1:length(kwta_options);
    for in_it=1:length(act);
 clearvars -except act m num_items it binary_res order_res binary_sim_res order_sim_res binary_rand order_rand kwta_options c in_it
for sim_val=1:100;


locations = randn(m,100);
c=cov(locations');

objects = randn(m,100);
c2=cov(objects');


ctr=1;
for k=1:m;
    for j=1:m;
       catted(ctr,:)=[locations(k,:) objects(j,:)];
       idx(ctr,:) = [k,j];
       ctr=ctr+1;
    end
end

c  = cov(catted');
%c   = c-min(c)/(max(c)-min(c));
    these_objs = randsample(m,m-1);
    these_words = randsample(m,m-1);
    o_complement = setdiff([1:m], these_objs);
    w_complement = setdiff([1:m], these_words);
 
pct_act = act(in_it);

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
for j=1:200;
sprse(:,j)=randsample(layer2,layer2*pct_act);
weights_layer2(j,sprse(:, j))=1;
end

[U,V]=eig(cov(catted'));
%for eigenvector case.
%output = U*weights_layer2;

output = catted*weights_layer2;

kwta=kwta_options(it)*layer2;
% %kwta
for j=1:num_items;
    [a b] = sort(abs(output(j,:)), 'descend');
    output(j,b(kwta:end))=0;
    output(j,b(1:kwta))=1;
end


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
out=1-pdist(output,  'hamming');
c=squareform(out);
c=c+eye(m^2);
%normalize hamming
c = (c - min(c)) ./ (max(c) - min(c)) ;


%     
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
   
  I = eye(m);
  candidates = setdiff([1:num_items], sub_pre);
  
  for p=1:length(candidates);    
      Itemp=I;
      for j=1:m-1;
          Itemp(j,m) = c(candidates(p), sub_pre(j));
      end
      
     for p1=1:m-2;
        for p2=p+1:m-1;
            csub(ctr,1) = c(sub_pre(p1), sub_pre(p2));
            Itemp(p1,p2)=c(sub_pre(p1), sub_pre(p2));
            ctr=ctr+1;
        end
     end
    
   Iout= triu(Itemp)+triu(Itemp,1)';
   det_results(p) = det(Iout);
   sim_results(p) = mean(mean(cov((Iout))));
  end

  [ares,bres]= max(det_results);
  [trsh srtd_res]= sort(det_results);
  binary_res(it, sim_val) = candidates(bres)==target_idx;
  %order_res(it, in_it,sim_val) = find(candidates(srtd_res)==target_idx);
  
  [ares,bres_sim]= max(sim_results);
  [trsh srtd_res]= sort(sim_results);
  binary_sim_res(it, sim_val) = candidates(bres_sim)==target_idx;
  %order_sim_res(it, in_it,sim_val) = find(candidates(srtd_res)==target_idx);
  
  binary_rand(it, in_it, sim_val) = randsample(candidates,1)==target_idx;
  
  %prediction error
  
  options = find(idx(:, 2)==o_complement);
  RU_object_choice = randsample(options,1);

  
  %ME
  
  %expected
  ME=options(options~=target_idx);
  MSE_ME(sim_val) = sum(catted(candidates(bres), :)-catted(target_idx, :)).^2;
  
  %unexpected.
  other=options(options~=target_idx);
  MSE_other(sim_val) = sum(catted(other, :)-catted(target_idx, :)).^2;
  %binary_rand(sim_val) = randsample(length(candidates),1)==target_idx;
  %order_rand(sim_val)=find(randsample(length(candidates), 13)==target_idx);
  
  
  %sim_expected
  MSE_sim_ME(sim_val) = sum(catted(candidates(bres_sim), :)-catted(target_idx, :)).^2;
  
  %sim unexpected.
  other=options(options~=target_idx);
  MSE_sim_other(sim_val) = sum(catted(candidates(bres_sim)-catted(other, :), :)).^2;
  
  
  
  %%sim random select
  MSE_random_ME(sim_val) = sum(catted(RU_object_choice, :)-catted(target_idx, :)).^2;
  MSE_random_other(sim_val) = sum(catted(RU_object_choice, :)-catted(other, :)).^2;
end

    end
end

    
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

