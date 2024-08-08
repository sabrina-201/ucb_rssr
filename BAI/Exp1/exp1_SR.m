%% Successive Rejects
clc;clear;
syms x

M = 1e4; 
T = 3e4;

K_vec = [16 32 64];
res_exp1_sr = [];

for K = K_vec
    count = 0;
    L = (sqrt(2*log(floor(T/(K*(log2(K)))))/(floor(T/(K*(log2(K)))))))/8;

    eqn = (((x/2).^2)./(L+(x.^2)/12)) == 1;
    roots1 = double(solve(eqn,x));
    
    eqn = (((x/2).^2)./(L+(x.^2)/12)) == 0.9;
    roots2 = double(solve(eqn,x));

    B = zeros(1,K);
    A = zeros(1,K);
    B(1) = roots1(2);
    B(2:end) = roots2(2)*ones(1,K-1);
    for m=1:M
%     parfor (m=1:M,4)
        mean_ = (B+A)/2;
        variance = (B-A).^2/12;
        sharpe = (mean_.^2)./(L+variance);
        [~, ind3] = sort(sharpe,'descend');
        
        arm_set = 1:K;
        X_bar_t = zeros(1,K);
        X_squared_t= zeros(1,K);
        s_var = zeros(1,K);
        s_mean= zeros(1,K);
        V_t = zeros(1,K);
        count_arm = zeros(1,K);
        dd=[];
        for i=1:K-1
            if i==1
                t_n = ceil((1/log_bar(K)) * (T-K)/(K+1-i));
            else
                t_n = ceil((1/log_bar(K)) * (T-K)/(K+1-i)) - ceil((1/log_bar(K)) * (T-K)/(K+1-(i-1)));
            end
            dd=[dd t_n];
            for k=arm_set
                for t=1:1:t_n
                    X = A(k) + (B(k) - A(k))*rand;
                    count_arm(k) = count_arm(k)+1;
                    X_bar_t(k) = ((count_arm(k) - 1)*X_bar_t(k) + X)/count_arm(k);
                    X_squared_t(k) = ((count_arm(k) - 1)*X_squared_t(k) + X^2)/count_arm(k);
                    V_t(k) = X_squared_t(k) - (X_bar_t(k)^2);
                end
            end
            [~, temp_ind] = sort((X_bar_t(arm_set).^2)./((V_t(arm_set))+L), 'descend'); 
            arm_set(temp_ind(end))=[];
        end
        if arm_set==ind3(1)
            count =count+1;
        end
    end
    res_exp1_sr = [res_exp1_sr (1 - (count/M))];
end

   