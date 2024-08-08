%% Performance of UCB-RSSR vs Modified-GRA-UCB and Modified-MVTS

% Output: tx, ty and tz are the Regret for UCB-RSSR, Modified-GRA-UCB, and Modified-MVTS respectively.
% Change value of K to complete Fig 3(a).

clc;clear all;
T = 1e4; %number of pulls of arms
M = 1e1; %number of simulations
K = 2; %number of arms
alpha=0.05;
ro=1;

tx=zeros(M, T);
ty=zeros(M, T);
tz=zeros(M, T);
L = 1;%sqrt(2*log(T)/(PILOT_FRACTION*T/K))/8;
L_G=L;
L_M=L;
load('arms_set2.mat')
A = arms_set2.A(1:K);
B = arms_set2.B(1:K);

% parfor (m=1:M,4)
for m=1:M
    s_loser_Sharpe = zeros(1, T);
    s_loser_GRA = zeros(1, T);
    s_loser_MVTS = zeros(1, T);

    X_bar_t_Sharpe = zeros(1,K);
    X_squared_t_Sharpe = zeros(1,K);
    s_Sharpe = zeros(1,K);
    V_t_Sharpe = zeros(1,K);
    B_t_Sharpe = zeros(1,K);

    X_bar_t_GRA = zeros(1,K);
    X_squared_t_GRA = zeros(1,K);
    s_GRA = zeros(1,K);
    V_t_GRA = zeros(1,K);
    B_t_GRA = zeros(1,K);

    X_bar_t_MVTS = zeros(1,K);
    X_squared_t_MVTS = zeros(1,K);
    s_MVTS = zeros(1,K);
    alph_a = 0.5*ones(1,K);
    bet_a = 0.5*ones(1,K);
    sigma_square = (1/12)*(B - A).^2;
    mu = (B + A)/2;
    mu_by_sigma = (mu.^2)./(L+sigma_square);
    [~, true_winner_arm] = max(mu_by_sigma);
    exp_phase = ceil(8*log(T)/(L^2));
    for k = 1:1:K
        for t=1:exp_phase
            X = A(k) + rand()*(B(k) - A(k));
            X_bar_t_Sharpe(k) = ((t - 1)*X_bar_t_Sharpe(k) + X)/t;
            X_squared_t_Sharpe(k) = ((t - 1)*X_squared_t_Sharpe(k) + X^2)/t;
            V_t_Sharpe(k) = X_squared_t_Sharpe(k) - ((X_bar_t_Sharpe(k))^2);
            s_Sharpe(k) = s_Sharpe(k) + 1;
            B_t_Sharpe(k)  = (X_bar_t_Sharpe(k)^2)/(L+V_t_Sharpe(k));

            X_bar_t_GRA(k) = ((t - 1)*X_bar_t_GRA(k) + X)/t;
            X_squared_t_GRA(k) = ((t - 1)*X_squared_t_GRA(k) + X^2)/t;
            V_t_GRA(k) = X_squared_t_GRA(k) - ((X_bar_t_GRA(k))^2);
            s_GRA(k) = s_GRA(k) + 1;
            B_t_GRA(k)  = (X_bar_t_GRA(k)^2+ sqrt(log(t)/s_GRA(k)))...
                /(L_G+(((V_t_GRA(k)*(s_GRA(k)-1))/(chi2inv(1-alpha,s_GRA(k)-1)))));
            alph_a(k) = alph_a(k)+1;
            bet_a(k) = max(bet_a(k), X);
            s_MVTS(k) = s_MVTS(k) + 1;
        end
    end
    
    [~, winner_arm_Sharpe] = max(B_t_Sharpe);

    [~, winner_arm_GRA] = max(B_t_GRA);

    draw = zeros(1,K);
    for k = 1:1:K
       draw(k) = gprnd(1./alph_a(k),bet_a(k)./alph_a(k),bet_a(k));
    end
    [~, winner_arm_MVTS] = max(((draw/2).^2)./(L_M + (draw.^2)/12));
    
    for t = exp_phase*K + 1:1:T
        X = A(winner_arm_Sharpe) + rand()*(B(winner_arm_Sharpe) - A(winner_arm_Sharpe));
        X_bar_t_Sharpe(winner_arm_Sharpe) = (s_Sharpe(winner_arm_Sharpe)*X_bar_t_Sharpe(winner_arm_Sharpe) + X)/(s_Sharpe(winner_arm_Sharpe) + 1);
        X_squared_t_Sharpe(winner_arm_Sharpe) = (s_Sharpe(winner_arm_Sharpe)*X_squared_t_Sharpe(winner_arm_Sharpe) + X^2)/(s_Sharpe(winner_arm_Sharpe) + 1);
        V_t_Sharpe(winner_arm_Sharpe) = X_squared_t_Sharpe(winner_arm_Sharpe) - ((X_bar_t_Sharpe(winner_arm_Sharpe))^2);
        s_Sharpe(winner_arm_Sharpe) = s_Sharpe(winner_arm_Sharpe) + 1;
        for k = 1:1:K
            epsilon = sqrt(2*log(t)/s_Sharpe(k));
            B_t_Sharpe(k) = ((X_bar_t_Sharpe(k)))^2/(L+V_t_Sharpe(k)) + ((X_squared_t_Sharpe(k) + V_t_Sharpe(k) + 2*epsilon+L)*epsilon/((V_t_Sharpe(k)+L)*(V_t_Sharpe(k) + L - 3*epsilon)+ 0*2*epsilon^2));
        end

        [~, winner_arm_Sharpe] = max(B_t_Sharpe);
        s_loser_Sharpe(t) = (t - s_Sharpe(true_winner_arm));

        X = A(winner_arm_GRA) + rand()*(B(winner_arm_GRA) - A(winner_arm_GRA));
        X_bar_t_GRA(winner_arm_GRA) = (s_GRA(winner_arm_GRA)*X_bar_t_GRA(winner_arm_GRA) + X)/(s_GRA(winner_arm_GRA) + 1);
        X_squared_t_GRA(winner_arm_GRA) = (s_GRA(winner_arm_GRA)*X_squared_t_GRA(winner_arm_GRA) + X^2)/(s_GRA(winner_arm_GRA) + 1);
        V_t_GRA(winner_arm_GRA) = X_squared_t_GRA(winner_arm_GRA) - ((X_bar_t_GRA(winner_arm_GRA))^2);
        s_GRA(winner_arm_GRA) = s_GRA(winner_arm_GRA) + 1;
        B_t_GRA(winner_arm_GRA) =  (ro*(X_bar_t_GRA(winner_arm_GRA)^2+sqrt(log(t)/s_GRA(winner_arm_GRA))))...
            /(L_G+(((V_t_GRA(winner_arm_GRA)*(s_GRA(winner_arm_GRA)-1))/(chi2inv(1-alpha,s_GRA(winner_arm_GRA)-1)))));
        [~, winner_arm_GRA] = max(B_t_GRA);
        s_loser_GRA(t) = (t - s_GRA(true_winner_arm));

        X = A(winner_arm_MVTS) + rand()*(B(winner_arm_MVTS) - A(winner_arm_MVTS));
        X_bar_t_MVTS(winner_arm_MVTS) = (s_MVTS(winner_arm_MVTS)*X_bar_t_MVTS(winner_arm_MVTS) + X)/(s_MVTS(winner_arm_MVTS) + 1);
        alph_a(winner_arm_MVTS) = alph_a(winner_arm_MVTS)+1;
        bet_a(winner_arm_MVTS) = max(bet_a(winner_arm_MVTS), X);
        s_MVTS(winner_arm_MVTS) = s_MVTS(winner_arm_MVTS) + 1;
        for k = 1:1:K
           draw(k) = gprnd(1./alph_a(k),bet_a(k)./alph_a(k),bet_a(k));
        end
        [~, winner_arm_MVTS] = max(((draw/2).^2)./(L_M + (draw.^2)/12));
        s_loser_MVTS(t) = (t - s_MVTS(true_winner_arm));   
    end
    tx(m,:) = s_loser_Sharpe;
    ty(m,:) = s_loser_GRA;
    tz(m,:) = s_loser_MVTS;
end

%%

fig1=figure;
A = arms_set2.A;
B = arms_set2.B;

sigma_square = (1/12)*(B - A).^2;
mu = (B + A)/2;

L = 1;
mu_by_sigma = mu.^2./(L+ sigma_square);
delta = abs(mu_by_sigma(1) - mu_by_sigma(2));

tem = mean(tx);
tem(1:370) = tem(371);
plot(1:200:t,delta*tem(1:200:end),'-d','DisplayName','$\texttt{UCB-RSSR}: K = 5$','LineWidth',2);
hold on;
grid on;
tem = mean(ty);
tem(1:370) = tem(371);
plot(1:200:t,delta*tem(1:200:end),'-*','DisplayName','$\texttt{Modified-GRA-UCB}: K = 5$','LineWidth',2);
tem = mean(tz);
tem(1:370) = tem(371);
plot(1:200:t,delta*tem(1:200:end),'-x','DisplayName','$\texttt{Modified-MVTS}: K = 5$','LineWidth',2);

extraInputs = {'interpreter','latex','fontsize',15};
xlabel('Time steps',extraInputs{:});
ylabel('Regret',extraInputs{:});
lgd = legend('Interpreter','latex','Position',[0.505730644861857 0.436882546652031 0.23495945930481 0.247804610318331]);
ax=gca;
ax.FontSize = 25;
lgd.FontSize = 15;
pbaspect([1 0.8 1]);
ax.TickLabelInterpreter = "latex";
set(gca, 'YScale', 'log')
