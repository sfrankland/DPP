warning('off')
for m=1:20;
    %for keyn=1:19
        %disp(keyn)
nkeys=m;
nvalues=m;
% 
% for j=1:1000;
% a=randsample(m,m, true);
% [r_rep(m,j) trash]=size(unique(a));
% rand_items(j) = sum(a);
% r_rep_collisions(m,j)=(m-r_rep(m,j))/m;
% end
% 

for i=1:500;
  
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
    
c = cov(state_mat');


%Random uniform
a=randsample(nkeys,nvalues, true);
[rand_key_select(m,i) trash]=size(unique(a));
r_rep_collisions(m,i)=(m-rand_key_select(m,i))/m;
%end
%avg_ran_key_select(keyn,m) = mean(rand_key_select(m,:)');
%    end
%end

% %Similarity-Based Key Selection.
%     %random starting index.
    p=randsample(m,1);
    for j=1:m;
        if j==1;
            key_select(j) = (idx(p,2));
        else
            this_state= intersect(find(idx(:,1) == j), find(idx(:,2) == p));
            list = 1:(nkeys*nvalues);
            list(this_state)=[];
            [a,b]=max(c(this_state, :));
            key_select(j)= idx(b,2);
        end
    end
    
sim_key_select=key_select;
[trash sim_rep(m,i)]=size(unique(sim_key_select));
sim_rep_collisions(m,i)=(m-sim_rep(m,i))/m;
% key_select=[];
     
%DPP-based Key Selection.
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
            this_state= intersect(find(idx(:,1) == j), find(idx(:,2) == p));
            list = 1:(nvalues*nkeys);
            %make the sub-matrix;
            %c(idx_select, idx_select)
            for h=1:length(list);
                temp_test = [idx_select list(h)];
                
                %[V,D]=eig(c(temp_test, temp_test));
                %dets(h)=trace(V*diag(exp(diag(D)))/V);
                dets(h) = det(c(temp_test, temp_test));
            end
            [a,b]=max(dets);
            %[a,b]=min(dets);
            idx_select(j)=list(b);
            key_select(j)= idx(list(b),2);
        end
    end
    
 dpp_key_select = key_select;
 [trash dpp_rep(m,i)]=size(unique(dpp_key_select));
 dpp_rep_collisions(m,i)=(m-dpp_rep(m,i))/m;
end
  %avg_dpp_key_select(keyn,m) = mean(dpp_rep(keyn,:));  
  %  end
    
end

 
%  
% 
% for k=1:m^2;
%     for j=1:m^2;
% gauss_kern(k,j) = exp(-norm( state_mat(k,:) - state_mat(j,:) ))^2;
%     end
% end
% 
% 
% 
% 
% 
% 
% 
% 
% %no severing eigenvectors.
% for j=1:1000;
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
% 
% 
% %name case
% for j=1:1000;
% names = eye(m);
% c=cov(names');
% [L.U,L.D]=eig(c);
% [a,b]=max((real(L.U)'*eye(m)));
% %[a,b]=max((L.U*eye(m)));
% name_items(m,j) = sum(b);
% name_samp(m,j)=length(unique(b));
% name_samp_collisions(m,j)=(m-eig_samp(m,j))/m;
% end
% 
% for j=1:1000;
% names = randn(m,1000);
% c=corr(names');
% [a,b]=max((names'*eye(m)));
% %[a,b]=max((L.U*eye(m)));
% corr_items(j) = sum(b);
% corr_samp(m,j)=length(unique(b));
% corr_samp_collisions(m,j)=(m-eig_samp(m,j))/m;
% end
% 
% 
% %simple severing.
% % for j=1:1000;
% % names = randn(m,1000);
% % c=corr(names');
% % [L.U,L.D]=eig(c);
% % [a,b(k)]=max((L.U(k,:)*eye(m-k+1)));
% % L.U(k,:)=[];
% % end
%  
% %[a,b]=max((L.U*eye(m)));
% % eig_items(j) = sum(b);
% % eig_samp(m,j)=length(unique(b));
% % eig_samp_collisions(m,j)=(m-eig_samp(m,j))/m;
% end
% 
% 
% 
% %bar([length(find(eig_items==2)); length(find(eig_items==3)); length(find(eig_items==4)); 0 ;length(find(rand_items==2)); length(find(rand_items==3)); length(find(rand_items==4)); 0; length(find(name_items==2)); length(find(name_items==3)); length(find(name_items==4))])
% 
% 
% 
% % %fourier
% for j=1:1000;
% names = randn(m,1000);
% c=corr(names');
% [L.U]=fft(c);
% [a,b]=max((real(L.U)'*eye(m)));
% fourier_items(j) = sum(a);
% fourier_samp(m,j)=length(unique(b));
% fourier_samp_collisions(m,j)=(m-fourier_samp(m,j))/m;
% end
% 
% %QR
% for j=1:1000;
% names = randn(m,1000);
% c=corr(names');
% [q,r]=qr(c);
% [a,b]=max((real(q)'*eye(m)));
% q_samp(m,j)=length(unique(b));
% q_samp_collisions(m,j)=(m-fourier_samp(m,j))/m;
% end
% 
% 
% end
% 
% 
% 
% names = randn(4,1000);
% c=corr(names');
% [u,v]=eig(c);
% [a,b]=max((u'*eye(m)));
% 
% %no severing eigenvectors.
% for j=1:1000;
% 
% 
% names = randn(m,1000);
% c=corr(names');
% [L.U,L.D]=eig(c);
% 
% for s=1:100;
%     
% [a,b]= max(sum(L.D));
% [a_eig,b_eig]= max(L.U(:,b)'*eye(m+1-s));
% 
% %cut
% L.U(b,:)=[];
% L.D(b,:)=[];
% idx_list(s) = b;
% end
% 
% [a,b]=max((u'*eye(m)));
% eig_samp(m,j)=length(unique(b));
% eig_samp_collisions(m,j)=(m-eig_samp(m,j))/m;
% end
% 
% 
%   % choose eigenvectors randomly
%   D = L.D/ (1+L.D);
%   v = find(rand(length(D),1) <= D);
% L.U=L.U(:)
% k = length(v);    
% V = L.U(:,v);
% 
% %severing eigenvectors.
% % iterate
% Y = zeros(k,1);
% for i = k:-1:1
%   
%   % compute probabilities for each item
%   P = sum(V.^2,2);
%   P = P / sum(P);
%   
%   % choose a new item to include
%   Y(i) = find(rand <= cumsum(P),1);
% 
%   % choose a vector to eliminate
%   j = find(V(Y(i),:),1);
%   Vj = V(:,j);
%   V = V(:,[1:j-1 j+1:end]);
% 
%   % update V
%   V = V - bsxfun(@times,Vj,V(Y(i),:)/Vj(Y(i)));
% 
%   % orthogonalize
%   for a = 1:i-1
%     for b = 1:a-1
%       V(:,a) = V(:,a) - V(:,a)'*V(:,b)*V(:,b);
%     end
%     V(:,a) = V(:,a) / norm(V(:,a));
%     size(V)
%   end
%   
%  end
% 
% Y = sort(Y);
% 
% disp(m)
% 
% end
% 
% %basic collision probability with random uniform is 1/m. r_rep should track
% %that.
% 
% %plot the eigenvector improvement relative to 1/m.
