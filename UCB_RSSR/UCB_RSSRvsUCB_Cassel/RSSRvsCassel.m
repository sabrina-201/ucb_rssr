%% Performance comparison of UCB-RSSR v/s USB-Cassel

% For the given parameters of A, B and L only Fig. 2(a) is generated.
% Change the values of A, B, K, and \Delta to get desired Figures 2(b) and 2(c). 

clc;clear all;
T = 1e4; %number of pulls of arms

M = 1e3; %number of simulations
K = 2; %number of arms

r_U_UCB = 0; % threshold for average reward
epsi_U_UCB_vec = [0.5 1 15 20]; % regularization paramter
gamma_U_UCB = 2.5;
nu_U_UCB = 0.5;
q_U_UCB = 2;
L_vec = epsi_U_UCB_vec;

rssr_res = [];
cassel_res = [];
delta_vec = [];

for ii = 1:length(L_vec)
    L = L_vec(ii);
    epsi_U_UCB = epsi_U_UCB_vec(ii);
    b_U_UCB = 4*max([epsi_U_UCB^(-1) 2*epsi_U_UCB^(-1.5)*(abs(r_U_UCB)+1)]);
    s_loser_Sharpe = zeros(M, T);
    s_loser_U_UCB = zeros(M, T);
    A = [0.5 0.6];
    B = [13.5 13.55];
    sigma_square = (1/12)*(B - A).^2;
    mu = (B + A)/2;
    mu_by_sigma = mu.^2./(epsi_U_UCB+ sigma_square);
    delta_vec = [delta_vec abs(diff(mu_by_sigma))];
    [~, true_winner_arm] = max(mu_by_sigma);
    for m=1:M
        sigma_square = zeros(1, K); 
        mu = zeros(1, K);
        mu_by_sigma = zeros(1, K);
    
        X_bar_t_Sharpe = zeros(1,K);
        X_squared_t_Sharpe = zeros(1,K);
        s_Sharpe = zeros(1,K);
        V_t_Sharpe = zeros(1,K);
        B_t_Sharpe = zeros(1,K);
    
        X_bar_t_U_UCB = zeros(1,K);
        X_squared_t_U_UCB = zeros(1,K);
        s_U_UCB = zeros(1,K);
        V_t_U_UCB = zeros(1,K);
        B_t_U_UCB = zeros(1,K);
        exp_phase = ceil(8*log(T)/L^2);
        for k = 1:1:K
            for t=1:exp_phase
                X = A(k) + rand()*(B(k) - A(k));
                X_bar_t_Sharpe(k) = ((t - 1)*X_bar_t_Sharpe(k) + X)/t;
                X_squared_t_Sharpe(k) = ((t - 1)*X_squared_t_Sharpe(k) + X^2)/t;
                V_t_Sharpe(k) = X_squared_t_Sharpe(k) - ((X_bar_t_Sharpe(k))^2);
                s_Sharpe(k) = s_Sharpe(k) + 1;
                B_t_Sharpe(k)  = (X_bar_t_Sharpe(k)^2)/(L+V_t_Sharpe(k));
    
                X_bar_t_U_UCB(k) = ((t - 1)*X_bar_t_U_UCB(k) + X)/t;
                X_squared_t_U_UCB(k) = ((t - 1)*X_squared_t_U_UCB(k) + X^2)/t;
                V_t_U_UCB(k) = X_squared_t_U_UCB(k) - ((X_bar_t_U_UCB(k))^2);
                s_U_UCB(k) = s_U_UCB(k) + 1;
                B_t_U_UCB(k)  = (X_bar_t_U_UCB(k)^2-r_U_UCB)/(epsi_U_UCB + V_t_U_UCB(k));
            end
        end
        
        [~, winner_arm_Sharpe] = max(B_t_Sharpe);
    
        [~, winner_arm_U_UCB] = max(B_t_U_UCB);
     
        for t = exp_phase*K + 1:1:T
            %% SR-UCB
            X = A(winner_arm_Sharpe) + rand()*(B(winner_arm_Sharpe) - A(winner_arm_Sharpe));
            X_bar_t_Sharpe(winner_arm_Sharpe) = (s_Sharpe(winner_arm_Sharpe)*X_bar_t_Sharpe(winner_arm_Sharpe) + X)/(s_Sharpe(winner_arm_Sharpe) + 1);
            X_squared_t_Sharpe(winner_arm_Sharpe) = (s_Sharpe(winner_arm_Sharpe)*X_squared_t_Sharpe(winner_arm_Sharpe) + X^2)/(s_Sharpe(winner_arm_Sharpe) + 1);
            V_t_Sharpe(winner_arm_Sharpe) = X_squared_t_Sharpe(winner_arm_Sharpe) - ((X_bar_t_Sharpe(winner_arm_Sharpe))^2);
            s_Sharpe(winner_arm_Sharpe) = s_Sharpe(winner_arm_Sharpe) + 1;
            for k = 1:1:K
               epsilon = sqrt(2.0*log(t)/s_Sharpe(k));
               B_t_Sharpe(k) = ((X_bar_t_Sharpe(k))^2/(L+V_t_Sharpe(k))) + ((X_squared_t_Sharpe(k) + V_t_Sharpe(k) + 2*epsilon+L)*epsilon/((V_t_Sharpe(k)+ L)*(V_t_Sharpe(k) + L - 3*epsilon) + 0*2*epsilon^2));
            end
            if ((X_squared_t_Sharpe(k) + V_t_Sharpe(k) + 2*epsilon+L)*epsilon/((V_t_Sharpe(k)+ L)*(V_t_Sharpe(k) + L - 3*epsilon) + 0*2*epsilon^2))<0
                disp('a');
            end
            [~, winner_arm_Sharpe] = max(B_t_Sharpe);
            s_loser_Sharpe(m, t) = (t - s_Sharpe(true_winner_arm));
    
            %% U-UCB
            X = A(winner_arm_U_UCB) + rand()*(B(winner_arm_U_UCB) - A(winner_arm_U_UCB));
            X_bar_t_U_UCB(winner_arm_U_UCB) = (s_U_UCB(winner_arm_U_UCB)*X_bar_t_U_UCB(winner_arm_U_UCB) + X)/(s_U_UCB(winner_arm_U_UCB) + 1);
            X_squared_t_U_UCB(winner_arm_U_UCB) = (s_U_UCB(winner_arm_U_UCB)*X_squared_t_U_UCB(winner_arm_U_UCB) + X^2)/(s_U_UCB(winner_arm_U_UCB) + 1);
            V_t_U_UCB(winner_arm_U_UCB) = X_squared_t_U_UCB(winner_arm_U_UCB) - ((X_bar_t_U_UCB(winner_arm_U_UCB))^2);
            s_U_UCB(winner_arm_U_UCB) = s_U_UCB(winner_arm_U_UCB) + 1;
            for k = 1:1:K
                B_t_U_UCB(k) = (((X_bar_t_U_UCB(k))-r_U_UCB)/sqrt(epsi_U_UCB + V_t_U_UCB(k))) + 2*b_U_UCB*max([sqrt(gamma_U_UCB*log(t)/(s_U_UCB(k)*nu_U_UCB))...
                    (gamma_U_UCB*log(t)/(s_U_UCB(k)*nu_U_UCB))^(q_U_UCB/2)]);
            end
    
            [~, winner_arm_U_UCB] = max(B_t_U_UCB);
            s_loser_U_UCB(m, t) = (t - s_U_UCB(true_winner_arm));  
        end
    end
    rssr_res = [rssr_res; mean(s_loser_Sharpe)];
    cassel_res = [cassel_res; mean(s_loser_U_UCB)];
end

%% Plotting

delta = delta_vec(1);
rssr_res(1,1:590) = rssr_res(1,591);
plot(1:200:t,delta*rssr_res(1,1:200:end),'-d','DisplayName','$\texttt{UCB-RSSR}: L = 0.5$','LineWidth',2);
grid on;
hold on;
cassel_res(1,1:590) = cassel_res(1,591);
plot(1:200:t,delta*cassel_res(1,1:200:end),'-o','DisplayName','$\texttt{U-UCB}: \epsilon_0 = 0.5$','LineWidth',2);

delta = delta_vec(2);
rssr_res(1,1:148) = rssr_res(1,149);
plot(1:200:t,delta*rssr_res(2,1:200:end),'-d','DisplayName','$\texttt{UCB-RSSR}: L = 1$','LineWidth',2);
cassel_res(1,1:148) = cassel_res(1,149);
plot(1:200:t,delta*cassel_res(2,1:200:end),'-o','DisplayName','$\texttt{U-UCB}: \epsilon_0 = 1$','LineWidth',2);

delta = delta_vec(3);
rssr_res(1,1:2) = rssr_res(1,3);
plot(1:200:t,delta*rssr_res(3,1:200:end),'-d','DisplayName','$\texttt{UCB-RSSR}: L = 15$','LineWidth',2);
cassel_res(1,1:2) = cassel_res(1,3);
plot(1:200:t,delta*cassel_res(3,1:200:end),'-o','DisplayName','$\texttt{U-UCB}: \epsilon_0 = 15$','LineWidth',2);

delta = delta_vec(4);
rssr_res(1,1:2) = rssr_res(1,3);
plot(1:200:t,delta*rssr_res(4,1:200:end),'-d','DisplayName','$\texttt{UCB-RSSR}: L = 20$','LineWidth',2);
cassel_res(1,1:2) = cassel_res(1,3);
plot(1:200:t,delta*cassel_res(4,1:200:end),'-o','DisplayName','$\texttt{U-UCB}: \epsilon_0 = 20$','LineWidth',2);

extraInputs = {'interpreter','latex','fontsize',15};
xlabel('Time steps',extraInputs{:});
ylabel('Regret',extraInputs{:});
lgd = legend('Interpreter','latex','location','northwest');
ax=gca;
ax.FontSize = 25;
lgd.FontSize = 19;
pbaspect([1 0.8 1]);
ax.TickLabelInterpreter = "latex";
