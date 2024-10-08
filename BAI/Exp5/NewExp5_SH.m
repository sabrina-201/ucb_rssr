%% Sequential Halving for Sharpe ratio
clc;clear all;
syms x

M = 1e4;
T_vec = [2e4 3e4 4e4];
u_1 = 0;
K = 32;

res_NewExp5_sh = [];
for ii = 1:length(T_vec)
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
    
    res_NewExp5_sh = [res_NewExp5_sh (1 - (count/M))];
end
