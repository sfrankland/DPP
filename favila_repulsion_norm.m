% ns = [10000];
% for it=1:length(ns);
% clearvars -except ns it
m=4;
num_items=m^2;

act = [0.5];
%kwta_options = [0.01 0.05 0.1];

neuron_options = [10000];

%act = [1];
kwta_options = [0.1];
for it=1:length(neuron_options);
    for in_it=1:length(act);
 clearvars -except act m num_items it binary_res order_res binary_sim_res order_sim_res binary_rand order_rand kwta_options c in_it neuron_options binary_res_den no_overlap diff_overlap same_overlap m
for sim_val=1:100;
    disp(sim_val)
    clearvars -except act m num_items it binary_res order_res binary_sim_res order_sim_res binary_rand order_rand avg_overlap overall_mean sim_val no_overlap diff_overlap same_overlap neuron_options in_it
   component1 = randn(1,100);
   component2 = randn(1,100);
   component3 = randn(1,100);
   component4 = randn(1,100);

   nkeys=100;
   %generate similar scene pairs.
   pair1_a = randn(1,100);
   pair1_b = pair1_a+randn(1,100)*0.5;
   
   pair2_a = randn(1,100);
   pair2_b = pair2_a+randn(1,100)*0.5;
   
   pair3_a = randn(1,100);
   pair3_b = pair3_a+randn(1,100)*0.5;
   
   pair4_a = randn(1,100);
   pair4_b = pair4_a+randn(1,100)*0.5;
   
   pair5_a = randn(1,100);
   pair5_b = pair5_a+randn(1,100)*0.5;
   
   pair6_a = randn(1,100);
   pair6_b = pair6_a+randn(1,100)*0.5;
   
  
   %no association
   catted(1,:) = [pair1_a zeros(1,100)];
   catted(2,:) = [pair1_b zeros(1,100)];
   catted(3,:) = [pair2_a zeros(1,100)];
   catted(4,:) = [pair2_b zeros(1,100)];
   
   %similar scenes to different faces
   faces_diff_a = randn(1,100);
   faces_diff_b = faces_diff_a + randn(1,100)*3;
   
   catted(5,:) = [faces_diff_a pair3_a];
   catted(6,:) = [faces_diff_b pair3_b];
   
   faces_diff_a = randn(1,100);
   faces_diff_b = faces_diff_a + randn(1,100)*3; 
   catted(7,:) = [faces_diff_a pair4_a];
   catted(8,:) = [faces_diff_b pair4_b];
   
   
   %similar scenes to the same face.
   faces_diff_a = randn(1,100);
   catted(9, :) = [pair5_a faces_diff_a];
   catted(10,:) = [pair5_b faces_diff_a];
   
   
   faces_diff_a = randn(1,100);
   catted(11, :) = [pair6_a faces_diff_a];
   catted(12,:) = [pair6_b faces_diff_a];
   
   
   target_nontarget  = zeros(12, 12);
   target_nontarget(1:2,1:2)=1;
   target_nontarget(3:4,3:4)=1;
   target_nontarget(5:6,5:6)=1;
   target_nontarget(7:8,7:8)=1;
   target_nontarget(9:10,9:10)=1;
   target_nontarget(11:12,11:12)=1;
   
   
   
  no_overlap_mat = zeros(12,12);
  same_overlap_mat = zeros(12,12);
  diff_overlap_mat = zeros(12,12);
  
  no_overlap_mat(1,2)=1;
  no_overlap_mat(3,4)=1;
  
  diff_overlap_mat(5,6)=1;
  diff_overlap_mat(7,8)=1;

  same_overlap_mat(9,10)=1;
  same_overlap_mat(11,12)=1;
   
    keys = randn(nkeys, 500);
    %keys = [repmat(1:nkeys)'];
    %keys = zscore(repmat(1:nkeys, 100,1)');
    %keys(1,:) = randn(1,100);
%     keys(2,:) = keys(1,:)+randn(1,100)+2;
%     keys(3,:) = keys(2,:)+randn(1,100)+3;
%     keys(4,:) = keys(3,:)+randn(1,100)+4;
%     keys(5,:) = keys(4,:)+randn(1,100)+5;
%     keys(6,:) = keys(5,:)+randn(1,100)+6;
%     keys(7,:) = keys(6,:)+randn(1,100)+7;
%     keys(8,:) = keys(7,:)+randn(1,100)+8;
    




    ctr=1;
    for j=1:12;
        for k=1:nkeys;
            key_catted(ctr,:)=[keys(k,:) catted(j,:)];
            idx(ctr,:) = [k,j];
            ctr=ctr+1;
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
weights_layer2 = zeros(500, layer2);
%expand and sparsify w/ binary random projection.
for j=1:500;
sprse(:,j)=randsample(layer2,layer2*pct_act);
weights_layer2(j,sprse(:, j))=1;
end

%[U,V]=eig(cov(catted'));
%for eigenvector case.
%output = U*weights_layer2;

% output = key_catted*weights_layer2;
% %norm_1_output = (output - min(output) )./ (max(output) - min(output));
% output = (output - min(output) )./ (max(output) - min(output));

  c=corr(key_catted');
 %randomly select an association to start.
 select = randsample(nkeys*12,1);
 starter = key_catted(select, :);
 candidates = 1:nkeys;
 candidates(idx(select, 2))=[];
 sub_pre = select;
 %('now here')
 for j=1:11;
     list = find(idx(:, 2)==candidates(j));
     %[track_dets]=det(cov(output(sub_pre, :)'));
     for l=1:length(list);
        cur_list = [sub_pre list(l)];
        track_dets(l) = det(c(cur_list, cur_list));
        %hippocampal code for subset.
%         DG = output(sub_pre,:)';
%         %normalization
%         lateral_connections = (DG*DG');
%         s_lat = sum(lateral_connections);
%         %norm_out = (output(list, :)*weights_layer2')./(s_lat+0.001);
%         norm_out = (output(list, :)./(s_lat+0.001));
%         [ares bres]=max(sum(norm_out'));
%         sub_pre = [sub_pre list(bres)];
        %track_dets=max(max(norm_out')));
        %[track_dets(l,1)]=det(cov(norm_out'));
        % [track_dets(l,1)]=det(cov(output(sub_pre, :)'));
        %[track_dets(l,1)]=max(max(norm_out'));

%lateral_connections = (lat - min(lat) )./ (max(lat) - min(lat));
%norm_out = layer2 ./ (lateral_connections + 0.001) ;
%norm_out = output./s_lat;
%norm_out = output+s_lat;
%max of max

          end
     %[ares bres]=max(max(norm_out'));
      [ares bres]=max(track_dets);
       sub_pre = [sub_pre list(bres)];
     clear track_dets
 end
 
  cmat = corr(key_catted(sub_pre, 1:500)', key_catted(sub_pre, 1:500)');
  %cmat = corr(key_catted(sub_pre, :)', key_catted(sub_pre,:)');
  r_idx = find(~eye(size(cmat)));
  %no_overlap(sim_val) = mean(mean([cmat(1:4, 1:4)]));
  %diff_overlap(sim_val) = mean(mean([cmat(5:8, 5:8)]));
  %same_overlap(sim_val) = mean(mean([cmat(9:12, 9:12)]));
  
 no_overlap(sim_val) = mean(cmat(no_overlap_mat==1)) - mean(cmat(target_nontarget==0));
 diff_overlap(sim_val) = mean(cmat(diff_overlap_mat==1)) - mean(cmat(target_nontarget==0));
 same_overlap(sim_val) = mean(cmat(same_overlap_mat==1)) - mean(cmat(target_nontarget==0));
  %overall_mean(sim_val) = mean(cmat(target_nontarget==0));
 
    end
end
end
