% ns = [10000];
% for it=1:length(ns);
% clearvars -except ns it
m=4;
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
    clearvars -except act m num_items it binary_res order_res binary_sim_res order_sim_res binary_rand order_rand avg_overlap overall_mean sim_val
   component1 = randn(1,100);
   component2 = randn(1,100);
   component3 = randn(1,100);
   component4 = randn(1,100);
   
    route1 = [component1  randn(1,100)];
    route2 = [component1  randn(1,100)];
    route3 = [component2  randn(1,100)];
    route4 = [component2  randn(1,100)];
    route5 = [component3  randn(1,100)];
    route6 = [component3  randn(1,100)];
    route7 = [component4  randn(1,100)];
    route8 = [component4  randn(1,100)];
    
    %episodes
    route1_a = [route1 randn(1,50)];
    route2_a = [route2 randn(1,50)];
    route3_a = [route3 randn(1,50)];
    route4_a = [route4 randn(1,50)];
    route5_a = [route5 randn(1,50)];
    route6_a = [route6 randn(1,50)];
    route7_a = [route7 randn(1,50)];
    route8_a = [route8 randn(1,50)];

    route1_b = [route1 randn(1,50)];
    route2_b = [route2 randn(1,50)];
    route3_b = [route3 randn(1,50)];
    route4_b = [route4 randn(1,50)];
    route5_b = [route5 randn(1,50)];
    route6_b = [route6 randn(1,50)];
    route7_b = [route7 randn(1,50)];
    route8_b = [route8 randn(1,50)];
    
%     route1 = [component1  randn(1,100)];
%     route2 = [component1  randn(1,100)];
%     route3 = [component2  randn(1,100)];
%     route4 = [component2  randn(1,100)];
%     route5 = [component3  randn(1,100)];
%     route6 = [component3  randn(1,100)];
%     route7 = [component4  randn(1,100)];
%     route8 = [component4  randn(1,100)];
%     
%     
%     route1 = [component1  randn(1,100)];
%     route2 = [component1  randn(1,100)];
%     route3 = [component2  randn(1,100)];
%     route4 = [component2  randn(1,100)];
%     route5 = [component3  randn(1,100)];
%     route6 = [component3  randn(1,100)];
%     route7 = [component4  randn(1,100)];
%     route8 = [component4  randn(1,100)];
%     
    
    
    all_routes = [route1; route2; route3; route4; route5; route6; route7; route8];
    
        
    all_routes = [route1_a; route2_a; route3_a; route4_a; route5_a; route6_a; route7_a; route8_a; route1_b; route2_b; route3_b; route4_b; route5_b; route6_b; route7_b; route8_b];

    
    keys = randn(8, 100);
    %keys(1,:) = randn(1,100);
%     keys(2,:) = keys(1,:)+randn(1,100)+2;
%     keys(3,:) = keys(2,:)+randn(1,100)+3;
%     keys(4,:) = keys(3,:)+randn(1,100)+4;
%     keys(5,:) = keys(4,:)+randn(1,100)+5;
%     keys(6,:) = keys(5,:)+randn(1,100)+6;
%     keys(7,:) = keys(6,:)+randn(1,100)+7;
%     keys(8,:) = keys(7,:)+randn(1,100)+8;
    
    ctr=1;
    for j=1:16;
        for k=1:8;
            catted(ctr,:)=[keys(k,:) all_routes(j,:)];
            if j<9;
            idx(ctr,:) = [k,j];
            else 
                idx(ctr,:) = [k,j-8];
            end
            ctr=ctr+1;
        end
    end

  c=cov(catted');
 %randomly select an association to start.
 select = randsample(128,1);
 starter = catted(select, :);
 candidates = 1:8;
 candidates(idx(select, 2))=[];
 sub_pre = select;
 %('now here')
 for j=1:15;
     list = find(idx(:, 2)==candidates(j));
     for l=1:length(list);
        cur_list = [sub_pre list(l)];
        track_dets(l) = det(c(cur_list, cur_list));
     end
     [ares bres]=max(track_dets);
     sub_pre = [sub_pre list(bres)];
     clear track_dets
 end
 
  cmat = corr(catted(sub_pre, 1:100)', catted(sub_pre, 1:100)');
  %cmat = corr(catted(sub_pre, 1:100)', catted(sub_pre, 1:100)');
  idx = find(~eye(size(cmat)));
  avg_overlap(sim_val) = mean([cmat(1,2), cmat(3,4), cmat(5,6), cmat(7,8)]);
  overall_mean(sim_val) = mean(cmat(idx));
 
end
    end

 
 
%[U,V]=eig(c);
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

weights_layer2 = zeros(300, layer2);

%expand and sparsify w/ binary random projection.
for j=1:300;
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


%assign each route a key.






%     for p=1:m-1;
%     sub_pre(p) = intersect(find(idx(:, 1)==these_words(p,:)), find(idx(:, 2)==these_objs(p,:)));
%     end
%     
%     target = [w_complement, o_complement];
%     target_idx = find(idx==[w_complement o_complement]);
%     
%     target_idx = intersect(find(idx(:, 1)==w_complement), find(idx(:, 2)==o_complement));
%     
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
   %candidates = setdiff([1:num_items], sub_pre);
   candidates = 1:length(catted);
%   
  for p=1:length(candidates);    
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
norm_out = (catted*weights_layer2)./(s_lat+0.001);
%lateral_connections = (lat - min(lat) )./ (max(lat) - min(lat));
%norm_out = layer2 ./ (lateral_connections + 0.001) ;
%norm_out = output./s_lat;
%norm_out = output+s_lat;
%max of max
[ares bres]=max(max(norm_out'));

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

