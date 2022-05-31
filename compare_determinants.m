
for L=1:1000;
    m=4;
names = randn(m,100);
objects = randn(m,100);

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
layer2=1000;
%layer2=10000;
%# of observations
n=16;
pct_act=1;
%# of nodes per variable.
d=15;

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

    these_objs = randsample(m,m-1);
    these_words = randsample(m,m-1);
    o_complement = setdiff([1:m], these_objs);
    w_complement = setdiff([1:m], these_words);
c=cov(catted');
candidates = 1:16;
    for p=1:m-1;
    sub_pre(p) = intersect(find(idx(:, 1)==these_words(p,:)), find(idx(:, 2)==these_objs(p,:)));
    end
  for p=1:length(candidates);    
      %Itemp=I;
      cur_list = [sub_pre candidates(p)];
      present_absent = zeros(m^2, m^2);
      present_absent(cur_list, cur_list)=1;
      Itemp = c(cur_list, cur_list);

       %hippocampal code for subset.
%      %normalization
        DG = output(sub_pre,:)';
        lateral_connections = (DG*DG');
        s_lat = sum(lateral_connections);
%         %norm_out = (output(list, :)*weights_layer2')./(s_lat+0.001);
        
        out =  diag(triu(Itemp)) ./ normr(triu(Itemp))+0.0001;
        norm_out = (output(cur_list, :)./(s_lat+0.001));
        %norm_out = (output(candidates(p), :)./ (sum((s_lat.*output(candidates(p),:))+0.001)));
        
        
        [det_norm(L,p)] = mean(sum(norm_out'));
%         sub_pre = [sub_pre list(bres)];
      
      det_prod(L,p) = prod(diag(Itemp));
      det_mat(L,p) = det(triu(Itemp));   
      det_sum(L,p) = sum(diag(Itemp));

%out =  diag(c.*present_absent) ./ sum(c .* present_absent')+0.0001;
%out =  diag(c.*eye(m^2)) ./ sum(c .* eye(m^2)')+0.0001;
%out =  diag(triu(Itemp)) ./ normr(triu(Itemp))+0.0001;
%det_norm_prod(L,p) = prod(diag(triu(out)));
  end
  %avg_corr_mat(L) = corr(det_norm_prod(L,:)', det_mat(L,:)');
  avg_corr_mat(L) = corr(det_norm(L,:)', det_mat(L,:)', 'type', 'Spearman');
   [a1,b1]=max(det_norm(L,:));
   [a2,b2]=max(det_mat(L,:));
   max_equal(L) = b1==b2;
  avg_corr_sum(L) = corr(det_norm(L,:)', det_sum(L,:)');
end
