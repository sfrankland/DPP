
m=4;
locations = randn(m,100);
c=cov(locations');

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


% assume symmetric matrix.
k=triu(catted);
det=prod(diag(k));




