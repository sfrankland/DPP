
m=4;
num_items=m^2;

%for sim_val=1:1000;

names = randn(m,100);
c=cov(names');

objects = randn(m,100);
c2=cov(objects');


ctr=1;
for k=1:4;
    for j=1:4;
       catted(ctr,:)=[names(k,:) objects(j,:)];
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
 
[U,V]=eig(c);
% # of active nodes in sparse layer.
pct_act = 0.05;

%size of expanded layer.
layer2=1000;

%# of observations
n=16;

%# of nodes per variable.
d=15;

%if it's kwta, what's k?
kwta_options=[1 5 10 25 50 100 250 500];

weights_layer2 = zeros(200, layer2);

 % %expand and sparsify w/ binary random projection.
for j=1:num_items;
sprse(:,j)=randsample(layer2,layer2*pct_act);
weights_layer2(j,sprse(:, j))=1;
end
%for eigenvector case.
output = catted*weights_layer2;

kwta=0.005*layer2;
% %kwta
for j=1:num_items;
    [a b] = sort(output(j,:), 'descend');
    output(j,b(kwta:end))=0;
    output(j,b(1:kwta))=1;
end
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
  end
  
  [ares,bres]= max(det_results);
  [trsh srtd_res]= sort(det_results);
  binary_res(sim_val) = candidates(bres)==target_idx;
  order_res(sim_val) = find(candidates(srtd_res)==target_idx);
  
  %binary_rand(sim_val) = randsample(length(candidates),1)==target_idx;
  %order_rand(sim_val)=find(randsample(length(candidates), 13)==target_idx);
  
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

