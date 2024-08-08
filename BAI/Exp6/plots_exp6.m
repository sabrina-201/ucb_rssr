%% Comp with cassel New
clc;clear all; close all


%% Exp 4

fig1=figure;

load('res_exp6_sh.mat');
load('res_exp6_sr.mat');
load('res_exp6_uni.mat');

sh = res_exp6_sh;
sr = res_exp6_sr;
us = res_exp6_uni;

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
lgd = legend('$K=16$','$K=32$', '$K=64$','Interpreter','latex','location','northwest');
ax=gca;
ax.FontSize = 25;
lgd.FontSize = 22;
pbaspect([1 0.8 1]);
ax.TickLabelInterpreter = "latex";
title('Experiment 6','interpreter','latex')
% ylim([0 0.6])
