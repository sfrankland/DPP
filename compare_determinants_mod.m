
for L=1:100;
% ns = [100 500 1000 5000 10000];
 
% for it=1:length(ns);
% clearvars -except ns it

m=4;
names = randn(m,100);
objects = randn(m,100);
num_items=m^2;
ctr=1;
for k=1:4;
    for j=1:4;
       catted(ctr,:)=names(k,:) + objects(j,:);
       idx(ctr,:) = [k,j];
       ctr=ctr+1;
    end
end


ctr=1;
for k=1:4;
    for j=1:4;
        if ctr==1;
       catted_keys(ctr,:)=names(k,:);
        else
              catted_keys(ctr,:)=names(k,:) + objects(j,:);
       idx(ctr,:) = [k,j];
        end
    ctr=ctr+1;
    end
end




%size of expanded layer.
layer2=10000;
%layer2=10000;
%# of observations
n=16;
pct_act=0.5;
%# of nodes per variable.
d=15;
kwta=0.01*layer2;
%if it's kwta, what's k?
%kwta_options=[1 5 10 25 50 100 250 500];
weights_layer2 = zeros(100, layer2);
%expand and sparsify w/ binary random projection.
for j=1:100;
sprse(:,j)=randsample(layer2,layer2*pct_act);
weights_layer2(j,sprse(:, j))=1;
end


output = catted_keys*weights_layer2;
%output = (output - min(output) )./ (max(output) - min(output));

% %kwta
% for j=1:num_items;
%     [a b] = sort(output(j,:), 'descend');
%     output(j,b(kwta:end))=0;
%     output(j,b(1:kwta))=1;
% end


    these_objs = randsample(m,m-1);
    these_words = randsample(m,m-1);
    o_complement = setdiff([1:m], these_objs);
    w_complement = setdiff([1:m], these_words);
    
%c=1-squareform(pdist(output, 'hamming'));
c=cov(output');
c1=cov(catted');
candidates = 1:16;

sub_pre = randsample(16,1);
        DG = output(sub_pre,:)';
        lateral_connections = (DG*DG');
        s_lat = sum(lateral_connections);
%     for p=1:m-1;
%     sub_pre(p) = intersect(find(idx(:, 1)==these_words(p,:)), find(idx(:, 2)==these_objs(p,:)));
%     end
  for p=1:length(candidates);    
      %Itemp=I;
      cur_list = [sub_pre candidates(p)];
      present_absent = zeros(m^2, m^2);
      present_absent(cur_list, cur_list)=1;
      Itemp_orig = c1(cur_list, cur_list);
      Itemp_hipp = c(cur_list, cur_list);
       %hippocampal code for subset.
%      %normalization
         %DG = output(sub_pre,:)';
         %DG = output(sub_pre,:)';
         %lateral_connections = (DG.*DG);
         %s_lat = sum(lateral_connections);
         
         %     norm_out = (output(candidates(p), :)./ sum(lateral_connections)) + 0.001;
        
%         %norm_out = (output(list, :)*weights_layer2')./(s_lat+0.001);
        %out =  diag(triu(Itemp_orig)) ./ (normr(triu(Itemp_orig))+0.0001);
       out =  diag(triu(Itemp_hipp)) ./ (normr(triu(Itemp_hipp))+0.0001);

       % out=triu(Itemp_hipp)/(sum(triu(Itemp_hipp)')+0.001);
        %out =  diag(triu(Itemp_hipp)) ./ (sum(triu(Itemp_hipp)./length(Itemp_hipp))+0.0001);
        %out =  diag(triu(Itemp_hipp)) ./ (normr(Itemp_hipp)+0.0001);
        %norm_out = (catted_keys*weights_layer2)./(s_lat+0.001);
       % norm_out = (output(candidates(p), :)./ (sum((s_lat.*output(candidates(p),:))+0.001)));
        
        
     [det_norm(L,p)] = sum(sum(out'));
%         sub_pre = [sub_pre list(bres)];
      
      det_hipp_prod(L,p) = prod(diag(Itemp_hipp));
      det_hipp_sum(L,p) = sum(diag(Itemp_hipp));
      det_mat(L,p) = det(triu(Itemp_orig));   
      det_sum(L,p) = sum(diag(Itemp_orig));

%out =  diag(c.*present_absent) ./ sum(c .* present_absent')+0.0001;
%out =  diag(c.*eye(m^2)) ./ sum(c .* eye(m^2)')+0.0001;
%out =  diag(triu(Itemp)) ./ normr(triu(Itemp))+0.0001;
%det_norm_prod(L,p) = prod(diag(triu(out)));
  end
  avg_corr_mat(L) = corr(det_norm(L,:)', det_mat(L,:)', 'type', 'Spearman');
  %avg_corr_mat(L) = corr(det_hipp_prod(L,:)', det_mat(L,:)', 'type', 'pearson');
   [a1,b1]=max(det_hipp_prod(L,:));
   [a2,b2]=max(det_mat(L,:));
   max_equal(L) = b1==b2;
  avg_corr_sum(L) = corr(det_hipp_sum(L,:)', det_sum(L,:)');
end
