%Sequential Halving for Sharpe ratio
clc;clear all;
M = 1e5; %number of simulations

T = 3e4;
K_vec = [16 32 64]; %number of arms

res_exp6_sh = [];
res_exp6_sr = [];
res_exp6_uni = [];
for K = K_vec

    count1=0;
    count2=0;
    count3=0;
    
    L = (sqrt(2*log(floor(T/(K*(log2(K)))))/(floor(T/(K*(log2(K)))))))/8;
    
%     for m=1:M
    parfor (m=1:M,4)
        b = rand(1,K);
        a = b.*rand(1,K);
    
        true_mean = (a+b)/2;
        true_var = ((b-a).^2)/12;
        true_sr = (true_mean.^2)./(L + true_var);
        [true_best_var, true_best_arm] = sort(true_sr, 'descend');
        arm_set = 1:K;
    
        %% SH
        
        for r=0:ceil(log2(K))-1
            X_bar_t = zeros(1,K);
            X_squared_t = zeros(1,K);
            V_t = zeros(1,K);
            S_r = ceil(K/2^r);
            t_r = floor(T/(S_r*ceil(log2(K))));%floor(T/(S_r*(log2(K))));
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
            count1 = count1+1;
        end
    
        %% SR
        B=b;
        A=a;
        mean_ = (B+A)/2;
        variance = (B-A).^2/12;
        sharpe = (mean_.^2)./(L+variance);
    
        [~, ind3] = sort(sharpe,'descend');
        
        arm_set = 1:K;
        X_bar_t = zeros(1,K);
        X_squared_t= zeros(1,K);
        count_arm = zeros(1,K);
        s_mean= zeros(1,K);
        V_t = zeros(1,K);
        for i=1:K-1
            if i==1
                t_n = ceil((1/log_bar(K)) * (T-K)/(K+1-i));
            else
                t_n = ceil((1/log_bar(K)) * (T-K)/(K+1-i)) - ceil((1/log_bar(K)) * (T-K)/(K+1-(i-1)));
            end
            for k=arm_set
                for t=1:1:t_n
                    X = A(k) + (B(k) - A(k))*rand;
                    count_arm(k) = count_arm(k) + 1;
                    X_bar_t(k) = ((count_arm(k) - 1)*X_bar_t(k) + X)/count_arm(k);
                    X_squared_t(k) = ((count_arm(k) - 1)*X_squared_t(k) + X^2)/count_arm(k);
                    V_t(k) = X_squared_t(k) - ((X_bar_t(k))^2);
                end
            end
            [~, temp_ind] = sort((X_bar_t(arm_set).^2)./((V_t(arm_set))+L), 'descend'); 
            arm_set(temp_ind(end))=[];
        end
        if arm_set==ind3(1)
            count2 =count2+1;
        end
    
        %% Unifrom
        uu = [];
        for j=1:K
            xx = a(j) + (b(j) - a(j))*rand(1,floor(T/K));
            sr = mean(xx).^2/(L + var(xx));
            uu = [uu sr];
        end
        [~, p] = max(uu);
        if ind3(1) == p
            count3 = count3+1;
        end
    
    end
    
    res_exp6_sh = [res_exp6_sh (1 - (count1/M))];
    res_exp6_sr = [res_exp6_sr (1 - (count2/M))];
    res_exp6_uni = [res_exp6_uni (1 - (count3/M))];
end
