clc;clear all;
syms x

M = 1e4;
T_vec = [2e4 3e4 4e4];

u_1 = 0;
K_vec = [16 32 64];

res_exp5_uni = [];
for ii = 1:length(K_vec)
    K = K_vec(ii);
    T = T_vec(ii);
    count=0;
    L = (sqrt(2*log(floor(T/(K*(log2(K)))))/(floor(T/(K*(log2(K)))))))/8;
    eqn = ((((x+u_1)/2).^2)./(L+((x-u_1).^2)/12)) == 1;
    roots1 = double(solve(eqn,x));

    a = zeros(1, K);
    b = zeros(1, K);
    b(1) = roots1(2);  

    arms_p = (0.98).^(2:K);
    roots2 = cell2mat(arrayfun(@(i) double(solve((((x/2).^2)./(L+(x.^2)/12)) == arms_p(i),x)),1:length(arms_p),'UniformOutput',false));

    if K==16
        b(2:K) = roots2(2,:);
    elseif K==32
        b(2:K) = roots2(2,:);
    elseif K==64
        b(2:K) = roots2(2,:);
    end
    mean_ = (b+a)/2;
    variance = (b-a).^2/12;
    sharpe = (mean_.^2)./(L+variance);
    [~, ind3] = sort(sharpe,'descend');
%     for i=1:M
    parfor (i=1:M,4)
        uu = [];
        for j=1:K
            xx = a(j) + (b(j) - a(j))*rand(1,floor(T/K));
            sr = (mean(xx).^2)/(L + var(xx));
            uu = [uu sr];
        end
        [~, p] = max(uu);
        if ind3(1) == p
            count = count+1;
        end
    end
    
    res_exp5_uni = [res_exp5_uni (1 - (count/M))];
end
