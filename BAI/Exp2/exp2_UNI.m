clc;clear all;
syms x

M = 1e4;
T = 3e4;
u_1 = 0;
K_vec = [16 32 64];

res_exp2_uni = [];
for K = K_vec
    count=0;
    L = (sqrt(2*log(floor(T/(K*(log2(K)))))/(floor(T/(K*(log2(K)))))))/8;
    eqn = ((((x+u_1)/2).^2)./(L+((x-u_1).^2)/12)) == 1;
    roots1 = double(solve(eqn,x));

    roots2 = double(solve((((x/2).^2)./(L+(x.^2)/12)) == 0.925,x));
    roots3 = double(solve((((x/2).^2)./(L+(x.^2)/12)) == 0.875,x));
    a = zeros(1, K);
    b = zeros(1, K);
    b(1) = roots1(2);
    if K==16
        rr1 = 2:6;
        rr2 = 7:K;
        b(rr1) = roots2(2)*ones(1,length(rr1));
        b(rr2) = roots3(2)*ones(1,length(rr2));
    elseif K==32
        rr1 = 2:14;
        rr2 = 15:K;
        b(rr1) = roots2(2)*ones(1,length(rr1));
        b(rr2) = roots3(2)*ones(1,length(rr2));
    elseif K==64
        rr1 = 2:30;
        rr2 = 31:K;
        b(rr1) = roots2(2)*ones(1,length(rr1));
        b(rr2) = roots3(2)*ones(1,length(rr2));
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
    
    res_exp2_uni = [res_exp2_uni (1 - (count/M))];
end
