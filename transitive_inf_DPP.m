
for it=1:1000;
    clearvars -except it perf TI
        nRelations=6;
 
a = randn(1,50)+1;
 b = randn(1,50);
 c = randn(1,50)+2;

a=a+randn(1,50)*0.5;
b=b+randn(1,50)*0.5;
c=c+randn(1,50)*0.5;
objects = [a;b;c];
relation = randn(nRelations,50);


%first a>b.
%c>a.


ctr=1;
for j=1:3;
for k=1:nRelations;
catted(ctr,:) = [objects(j,:) relation(k,:)];
idx(ctr,1)=j;
idx(ctr,2)=k;
ctr=ctr+1;
end
end
c=corr(catted');
[a,b_max]=max(mean(relation'));
[a,b_min]=min(mean(relation'));

itemC = intersect(find(idx(:,1)==3), find(idx(:,2)==b_max));
itemB = intersect(find(idx(:,1)==2), find(idx(:,2)==b_min));

subset=[itemC itemB];

for j=1:3*nRelations;
    current_subset=[subset j];
    data= [catted(current_subset(1), :);catted(current_subset(2), :);catted(current_subset(3),:)];
    vol_data(j) = det(corr(data'));
end

%a>b, && c>a.
itemA = find(idx(:,1)==1);
itemB = find(idx(:,1)==2);
itemC = find(idx(:,1)==3);


[a,b]=max(vol_data);
t=idx(b,:);
perf(it) = t(1)==1 && t(2)~=b_max && t(2)~=b_min;

%sign(mean(relation(t(2),:)')-mean(relation(b_max,:))') . == sign(mean(relation(t(b_min),:)')-mean(relation(b_max,:))');
TI(it) = sign(mean(relation(t(2),:)')-mean(relation(b_max,:))') == sign(mean(relation(b_min,:)')-mean(relation(b_max,:))');




end

