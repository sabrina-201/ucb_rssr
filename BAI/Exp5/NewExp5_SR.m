%% Successive Rejects
clc;clear;
syms x

M = 1e4; 
T_vec = [2e4 3e4 4e4];

K = 32;
res_NewExp5_sr = [];

for ii = 1:length(T_vec)
    T = T_vec(ii);
    count = 0;
    L = (sqrt(2*log(floor(T/(K*(log2(K)))))/(floor(T/(K*(log2(K)))))))/8;

    eqn = (((x/2).^2)./(L+(x.^2)/12)) == 1;
    roots1 = double(solve(eqn,x));

    B = zeros(1,K);
    A = zeros(1,K);
    B(1) = roots1(2);
    
    arms_p = (0.98).^(2:K);
    roots2 = cell2mat(arrayfun(@(i) double(solve((((x/2).^2)./(L+(x.^2)/12)) == arms_p(i),x)),1:length(arms_p),'UniformOutput',false));

    if K==16
        B(2:K) = roots2(2,:);
    elseif K==32
        B(2:K) = roots2(2,:);
    elseif K==64
        B(2:K) = roots2(2,:);
    end
%     for m=1:M
    parfor (m=1:M,4)
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
    res_NewExp5_sr = [res_NewExp5_sr (1 - (count/M))];
end

   