
for m=1:20;
    %for keyn=1:19
    %disp(keyn)
    nkeys=m;
    nvalues=m;
for i=1:100;
%data.
values = randn(nvalues, 100);
keys = randn(nkeys, 100);

%conjunctions.
clear idx
clear state_mat
ctr=1;
for k=1:nvalues;
    for j=1:nkeys;
      state_mat(ctr,:) = [values(k,:) keys(j,:)];
      idx(ctr,1)=k; idx(ctr, 2)=j;
ctr=ctr+1;
    end
end
    
c = corr(state_mat');
 
%Neural algorithm for DPP-based Key Selection.
act = [0.1];
%kwta_options = [0.01 0.05 0.1];
neuron_options = [10000];
  % # of active nodes in sparse layer.
pct_act = act;
%size of expanded layer.
layer2=neuron_options(1);
%layer2=10000;
%if it's kwta, what's k?
%kwta_options=[1 5 10 25 50 100 250 500];
weights_layer2 = zeros(nkeys^2, layer2);
%expand and sparsify w/ binary random projection.
for j=1:200;
sprse(:,j)=randsample(layer2,layer2*pct_act);
weights_layer2(j,sprse(:, j))=1;
end
output = state_mat*weights_layer2;
    %random starting index.
%    for i=1:10;
    clear dets
    %p=randsample(m,1);
    for j=1:m;
        if j==1;
            p=randsample(nkeys*nvalues,1);
            %p=randsample([15; 19; 2; 6; 23], 1);
            key_select(j) = (idx(p,2));
            idx_select = p;
            %candidate_sub_matrix=  this_state(key_select(j), :);
        else
%             weights_layer2 = zeros(nkeys^2, layer2);
%             %expand and sparsify w/ binary random projection.
%             for j=1:500;
%                 sprse(:,j)=randsample(layer2,layer2*pct_act);
%                 weights_layer2(j,sprse(:, j))=1;
%             end
            this_state= intersect(find(idx(:,1) == j), find(idx(:,2) == p));
            list = 1:(nvalues*nkeys);
            %make the sub-matrix;
            %c(idx_select, idx_select)
            for h=1:length(list);
                temp_test = [idx_select list(h)];
                 %hippocampal code for subset.
                 DG = output(temp_test,:)';
%                %normalization
                lateral_connections = (DG*DG');
                s_lat = sum(lateral_connections);
                norm_out = sum(sum(output(temp_test, :))./(s_lat+0.001));
%         norm_out = (output(list, :)./(s_lat+0.001));
              [ares(h) bres(h)]=max(norm_out);
%         sub_pre = [sub_pre list(bres)];
                %dets(h) = det(c(temp_test, temp_test));
            end
            [a,b]=max(ares);
            idx_select(j)=list(b);
            key_select(j)= idx(list(b),2);
        end
    end
    
 dpp_key_select = key_select;
 [trash dpp_rep(m,i)]=size(unique(dpp_key_select));
 dpp_rep_collisions(m,i)=(m-dpp_rep(m,i))/m;
 
end
end
