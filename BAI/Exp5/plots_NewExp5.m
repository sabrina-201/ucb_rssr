%% Comp with cassel New
clc;clear all; close all


%% Exp 5

fig1=figure;

load('res_NewExp5_sh.mat');
load('res_NewExp5_sr.mat');
load('res_NewExp5_uni.mat');

sh = res_NewExp5_sh;
sr = res_NewExp5_sr;
us = res_NewExp5_uni;

bh = bar([sh; sr; us],1);

bh(1).FaceColor = "#7E2F8E";
bh(2).FaceColor = "#77AC30";
bh(3).FaceColor = "#4DBEEE";
grid;

Exp = {'$\texttt{SHSR}$', '$\texttt{SuRSR}$', 'Uniform'};
xticks(1:numel(Exp)); % Setting x-axis ticks
xticklabels(Exp); % Setting x-axis labels
 

extraInputs = {'interpreter','latex','fontsize',15};
ylabel('Probability of Error - $e_n$',extraInputs{:});
lgd = legend('$n = 2e4$','$n = 3e4$', '$n = 4e4$','Interpreter','latex','location','northwest');
ax=gca;
ax.FontSize = 25;
lgd.FontSize = 22;
pbaspect([1 0.8 1]);
ax.TickLabelInterpreter = "latex";
title('Experiment 5','interpreter','latex')
% ylim([0 0.02])
