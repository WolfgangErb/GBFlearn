% GBFlearn: a toolbox for graph signal interpolation
% and classification with graph basis functions (GBFs)
% (C) W. Erb 15.01.2023

% Name: GBF_example_ITP_sensor2_comparison.m
% This example script compares GBF interpolation with
% bandlimited interpolation on a sensor graph. For bandlimited
% interpolation strong Runge-type artifacts are visible.

% Test scenario:
% graph: small sensor graph
% kernel: diffusion and bandlimited 
% number of sampling nodes: m = 40;
% problem: compare GBF interpolation with bandlimited interpolation

clear all
close all

% Paths
addpath(genpath('./core/'))
addpath(genpath('./data/'))

%Choose graph
G.type = 'sensor2';

%Generate graph
[G.nodes,G.edges,G.A] = GBF_gengraph(G.type);

%Calculate the graph Laplacian
G.N = length(G.nodes(:,1));
G.deg = sum(G.A,1);
isD = diag(1./sqrt(G.deg));
G.L = eye(G.N) - isD*G.A*isD;

%Choose plotting parameter
plotpar.MM = 2;                 %size of dots
plotpar.ub = 0.02;              %upper boundary
plotpar.lb = 0.02;              %left boundary
plotpar.fontsize = 14;          %fontsize
plotpar.colorbar = 'y';         %set 'y' if you want a colorbar
plotpar.edge = 1;               %set 1 if you want to plot edges
plotpar.edgewidth = 1;          %width of edges

%Create phantom and interpolation data on graph
M = 40;
idxW = (1:M)';

findP = find(G.nodes(:,1) >= 0.5);
fOriginal = zeros(G.N,1);
fOriginal(findP) = 1;
yW = fOriginal(idxW);

%Kernel parameters
type1 = 'diffusion';
alpha1 = 4;

type2 = 'bandlimited';
alpha2 = [0.93,200];

%Calculate GBF interpolants
bf1 = GBF_genGBFeff(G.L, idxW,type1,alpha1);
bf2 = GBF_genGBFeff(G.L, idxW,type2,alpha2);

s1 = GBF_itpGBF(bf1, idxW, yW);
s2 = GBF_itpGBF(bf2, idxW, yW);

%Plot GBF and bandlimited interpolation
figure('Units', 'pixels', ...
'Position', [0 50 1200 400]);

subplot(2,6,[1,2,7,8]),
plotpar.uaxis = 1.1;
plotpar.laxis = -0.1;
GBF_drawsignal(G.nodes,G.edges,fOriginal,plotpar);
title('Original');
set(gca,'XTick',[], 'YTick', [])
hold off;

subplot(2,6,[3,4,9,10]), 
plotpar.uaxis = max(s1);
plotpar.laxis = min(s1);
GBF_drawsignal(G.nodes,G.edges,s1,plotpar);
title('GBF interpolation');
hold on;
plot( G.nodes(idxW,1),G.nodes(idxW,2),'o','color',[0, 0, 0]/255,'LineWidth',1,'MarkerSize',8)
set(gca,'XTick',[], 'YTick', [])
hold off;

subplot(2,6,[5,6,11,12]),
plotpar.uaxis = max(s2);
plotpar.laxis = min(s2);
GBF_drawsignal(G.nodes,G.edges,s2,plotpar);
title('Bandlimited interpolation');
hold on;
plot( G.nodes(idxW,1),G.nodes(idxW,2),'o','color',[0, 0, 0]/255,'LineWidth',1,'MarkerSize',8)
set(gca,'XTick',[], 'YTick', [])
hold off;