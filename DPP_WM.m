

center = randn(30,1)
keymat = randn(20,30);
for k=1:20;
small_region(k,:) = center+randn(30,1);
mid_region(k,:) = center+randn(30,1)*2;
large_region(k,:) = center+randn(30,1)*6;
end

ctr=1;
for k=1:20;
    for j=1:20;
        sum_of_close_codes(ctr, :) = small_region(k,:) + keymat(j,:);
        sum_of_mid_codes(ctr, :) = mid_region(k,:) + keymat(j,:);
        sum_of_far_codes(ctr, :) = large_region(k,:) + keymat(j,:);
        ctr=ctr+1;
    end
end

c_close = cov(sum_of_close_codes');
c_mid = cov(sum_of_mid_codes');
c_far = cov(sum_of_far_codes');

s = squareform(pdist(small_region, 'euclidean'));
m = squareform(pdist(mid_region, 'euclidean'));
l = squareform(pdist(large_region, 'euclidean'));


% items that are close together (e.g., s) will have a small determinant.
% far apart start with more dissimilar codes, larger determinant.


%map percept to word. 
%
s = normr(1-normr(squareform(pdist(small_region, 'euclidean'))));
m = normr(1-normr(squareform(pdist(mid_region, 'euclidean'))));
l = normr(1-normr(squareform(pdist(large_region, 'euclidean'))));