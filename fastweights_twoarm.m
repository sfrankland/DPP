
rate_ctr=1;
decay_ctr=1;
lr=0.5;

for lr=0.5
for e=0.5
for it=1:1000;

sequence = rand(9,100);
Wx = rand(100,50); 
Wh = eye(50);
%Wh=randn(50,50);
h = rand(1,50);

%carry through (MPP). standard fast weights.
for step=2:9;
    h_in = sequence(step, :) * Wx;
    hout = e*(h(step-1,:)'*h(step-1,:));
    %h(step,:) = zscore(h_in*hout);
    %similarity to previous state, and add input/hidden).
    h(step,:) = zscore(  lr*(hout  * (h(step-1, :)*Wh' + h_in)') + (h(step-1, :)*Wh' + h_in)');
    %h(step,:) = (hout  * (h(step-1, :)*Wh' + h_in)') + (h(step-1, :)*Wh' + h_in)';
    %h(step,:) = zscore(h_in*hout) +  zscore(h(step-1, :)*Wh);
    %h(step,:) = h_in*hout;
    adj(step-1)=corr(h(step, :)', h(step-1,:)');
end 

mpp_compare(it, rate_ctr) = mean(mean(abs(corr(h'))));
mpp_compare_noabs(it, rate_ctr) = mean(mean(corr(h')));
mpp_adj(it, rate_ctr) = mean(abs(adj));
mpp_det(it) = det(cov(h'));


%individuate (DPP).
for step=2:9;
    h_in = sequence(step, :) * Wx;
    hout = e*(h(step-1,:)'*h(step-1,:));
    %h(step,:) = zscore(h_in/hout) + zscore(h(step-1, :)*Wh);
    %h(step,:) = zscore((h(step-1, :)*Wh+h_in / (hout)+0.001));
    h(step,:) = zscore((lr*(h(step-1, :)*Wh' + h_in) / hout) + (h(step-1, :)*Wh' + h_in));
    %h(step,:) = zscore((h(step-1, :)*Wh+h_in) / (randn(50,50)+0.001));
    %h(step,:) = h_in/hout;
    adj(step-1)=corr(h(step, :)', h(step-1,:)');
end 

dpp_compare(it, rate_ctr) = mean(mean(abs(corr(h'))));
dpp_compare_noabs(it, rate_ctr) = mean(mean(corr(h')));
dpp_adj(it, rate_ctr) = mean(abs(adj));
dpp_det(it) = det(cov(h'));

%rand numerator
for step=2:9;
    h_in = sequence(step, :) * Wx;
    %hout = h(step-1,:)'*h(step-1,:);
    hout=e*(randn(50,1)*randn(1,50));
    %hout = hout*0.5;
    %h(step,:) = zeros(1,50) + zscore(h(step-1, :)*Wh); 
    h(step,:) = zscore( lr* (hout  * (h(step-1, :)*Wh' + h_in)') + (h(step-1, :)*Wh' + h_in)');
    h(step,:) = zscore(rand(1,50));
    %h(step,:) = zscore(hout*randn(1,50)');
    adj(step-1)=corr(h(step, :)', h(step-1,:)');
   
    
end 

rand_compare(it, rate_ctr) = mean(mean(abs(corr(h'))));
rand_compare_noabs(it, rate_ctr) = mean(mean(corr(h')));
rand_adj(it, rate_ctr) = mean(abs(adj));
rand_det(it) = det(cov(h'));

%rand denominator
for step=2:9;
    h_in = sequence(step, :) * Wx;
    %hout = h(step-1,:)'*h(step-1,:);
    hout=e*(randn(50,1)*randn(1,50));
    %hout = hout*0.5;
    %h(step,:) = zeros(1,50) + zscore(h(step-1, :)*Wh); 
    h(step,:) = zscore((lr*(h(step-1, :)*Wh' + h_in) / hout) + (h(step-1, :)*Wh' + h_in));
    %h(step,:) = zscore(hout*randn(1,50)');
    adj(step-1)=corr(h(step, :)', h(step-1,:)');
end 

rand_denom_compare(it, rate_ctr, decay_ctr) = mean(mean(abs(corr(h'))));
rand_denom_compare_noabs(it, rate_ctr, decay_ctr) = mean(mean(corr(h')));
rand_denom_adj(it, rate_ctr, decay_ctr) = mean(abs(adj));

end
rate_ctr=rate_ctr+1;
end
decay_ctr=decay_ctr+1;
end

