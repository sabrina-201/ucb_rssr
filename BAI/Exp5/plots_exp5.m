%% Comp with cassel New
clc;clear all; close all


%% Exp 5

fig1=figure;

load('res_exp5_sh.mat');
load('res_exp5_sr.mat');
load('res_exp5_uni.mat');

sh = res_exp5_sh;
sr = res_exp5_sr;
us = res_exp5_uni;

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
lgd = legend('$K=16, \; n = 2e4$','$K=32, \; n = 3e4$', '$K=64, \; n = 4e4$','Interpreter','latex','location','northwest');
ax=gca;
ax.FontSize = 25;
lgd.FontSize = 22;
pbaspect([1 0.8 1]);
ax.TickLabelInterpreter = "latex";
title('Experiment 5','interpreter','latex')
% ylim([0 0.02])
