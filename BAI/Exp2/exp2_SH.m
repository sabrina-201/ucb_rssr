%% Sequential Halving for Sharpe ratio
clc;clear all;
syms x

M = 1e4;
T = 3e4;
u_1 = 0;
K_vec = [16 32 64];

res_exp2_sh = [];
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
    
%     for m=1:M
    parfor (m=1:M,4)
        true_mean = (a+b)/2;
        true_var = ((b-a).^2)/12;
        true_sr = (true_mean.^2)./(L + true_var);
        [true_best_var, true_best_arm] = sort(true_sr, 'descend');
        arm_set = 1:K;
        dd=[];
        for r=0:ceil(log2(K))-1
            X_bar_t = zeros(1,K);
            X_squared_t = zeros(1,K);
            V_t = zeros(1,K);
            S_r = ceil(K/2^r);
            t_r = floor(T/(S_r*ceil(log2(K))));
            dd=[dd t_r];
            for k=arm_set
                for t=1:t_r
                    X = a(k) + (b(k) - a(k))*rand;
                    X_bar_t(k) = ((t - 1)*X_bar_t(k) + X)/t;
                    X_squared_t(k) = ((t - 1)*X_squared_t(k) + X^2)/t;
                    V_t(k) = X_squared_t(k) - ((X_bar_t(k))^2);
                end
            end
            [~, sorted_arms] = sort((X_bar_t(arm_set).^2)./(L + V_t(arm_set)), 'descend');
            arm_set = arm_set(sorted_arms);
            arm_set = arm_set(1:ceil(length(arm_set)/2));
         end
         
        if arm_set == true_best_arm(1)
            count = count+1;
        end
    end
    
    res_exp2_sh = [res_exp2_sh (1 - (count/M))];
end
