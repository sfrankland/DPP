
for it=1:1000;
    
objects = randn(2,100);
relations= randn(2,100);
relations(3,:)= relations(2,:)*-1;

state1 = [objects(1,:) relations(1,:)];
state2 = [objects(1,:) relations(2,:)];
state3 = [objects(1,:) relations(3,:)];
state4 = [objects(2,:) relations(1,:)];
state5 = [objects(2,:) relations(2,:)];
state6 = [objects(2,:) relations(3,:)];


state_mat = [state1; state2; state3; state4; state5; state6];

gauss_kern = eye(6,6);
%gaussian kernel.
for k=1:6;
for j=1:6;
  % if j~=k;
gauss_kern(k,j) = exp(-norm( state_mat(k,:) - state_mat(j,:) ))^2;
end
end
%end
%gauss_kern_z = zscore(gauss_kern, 0, 2);

det_sub_mats= eye(6,6)*-0.000000000000000001;
for k=1:6;
    for j=1:6;
        if k~=j;
        %det_sub_mats(k,j) = det([gauss_kern(k,k) gauss_kern(k,j); gauss_kern(j,k) gauss_kern(j,j)]);
        det_sub_mats(k,j)= gauss_kern(k,k)*gauss_kern(j,j) - gauss_kern(j,k)*gauss_kern(k,j);
    end
end
end

%[a,b]=max(det_sub_mats*-1);

[a,b]= min(gauss_kern);

%predicted model.
% 1 should be randomly varying between 5,6
%2 should be 6.
%3 should be 5.
%4 should be randomly varying 2 and 3.
%5 should be 3.
%6 should be 2.

repel=0;
if b(2)==6; repel=repel+1;
   if b(3)==5; repel=repel+1;
       if b(5)==3; repel=repel+1;
           if b(6)==2; repel=repel+1;
           end
       end
   end
end

perf(it) = repel/4;

end
