



x=[1:100];
y=[1:100];
ctr=1;
    for j=1:100;
        for k=1:100;
            data_mat(ctr,1)=j;
            data_mat(ctr,2)=k;
            ctr=ctr+1;
        end
    end
    
    ctr=1;
for p=1:length(data_mat);
    for q=1:length(data_mat);
        d(ctr)=sqrt((data_mat(p, 1) - data_mat(q, 1) + data_mat(p, 2) - data_mat(q, 2))^2); 
ctr=ctr+1;
    end
end


norm_d = 1 - (d - min(d)) / (max(d) - min(d));


    %norm_d=18-d;
%euclidean distance.
% ctr=1;
%     for j=1:10;
%         for k=1:10;
%             for m=1:10;
%                 for n=1:10;
%             d(ctr) = sqrt(((x(j)-x(k)) + (y(m)-y(n)))^2);
%             ctr=ctr+1;
%         end
%             end
%         end
%     end
% 
k = reshape(norm_d, 10000, 10000);
%k=A;
sub_pre = randsample(10000,1);
%sub_pre=1;

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
   
  I = eye(100);
  candidates = 1:10000;
  m=100;
  
  num_samples=1000;
%   
%    for o1=1:length(candidates);
%        for o2=1:length(candidates);
%            for o3=1:length(candidates);
%                for o4=1:length(candidates);
%                     for o5=1:length(candidates);
%                         for o6=1:length(candidates);
%                             for o7=1:length(candidates);
%                                 for o8=1:length(candidates);
%                                     k(candidates([o1, o2, o3, o4, o5, o6, o7, o8]
%            
             


test_list = randsample(10000,num_samples);



     init_mat= k(test_list, test_list);
     init_det = det(current_mat);
     obj=0;
     
     while obj ~= met;
     item_to_adjust=test_list(randsample(num_samples));
     test_list(item_to_adjust) = d(test_list(item_to_adjust));
     
     current_mat= k(test_list, test_list);
     this_det = det(current_mat);
     %update items, and state.
     if this_det > det_state;
         test_list = cur_list;
         det_state = this_det;
     else item_to_adjust = test_list(randsample(num_samples));
     end
     
     end
     
     
     
     [a,choice]=max(these_dets);






  for l=1:num_samples;
      if l==1;
            cur_list = sub_pre;
      end
      for p=1:length(candidates);
        test_list = [candidates(p) cur_list];
          current_mat= k(test_list, test_list);
          these_dets(p) = det(current_mat);
     end
      [a,choice]=max(these_dets);
      %add to the list.
      cur_list = [cur_list choice];
      clear these_dets current_mat test_list
  end
  
      
%      for p1=1:100-2;
%         for p2=p+1:100-1;
%             csub(ctr,1) = k(sub_pre(p1), sub_pre(p2));
%             Itemp(p1,p2)=k(sub_pre(p1), sub_pre(p2));
%             ctr=ctr+1;
%         end
%      end
%     
   Iout= triu(Itemp)+triu(Itemp,1)';
   det_results(p) = det(Iout);
   sim_results(p) = mean(mean(corr((Iout))));
   
   
  end

%10000