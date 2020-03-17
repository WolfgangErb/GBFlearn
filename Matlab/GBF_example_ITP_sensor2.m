% GBFlearn: a toolbox for graph signal interpolation
% and classification with graph basis functions (GBFs)
% (C) W. Erb 01.03.2020

% GBF_example_ITP_sensor2 is a very simple script to show how a GBF interpolant
% can be calculated and plotted by the toolbox GBFlearn

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

% calculate GFT
[G.U,G.Lambda] = GBF_spectrum(G.L,'ascend');

%Choose plotting parameter
plotpar.MM = 2;                 %size of dots
plotpar.ub = 0.02;              %upper boundary
plotpar.lb = 0.02;              %left boundary
plotpar.uaxis = 1.1;            %upper colorbar boundary
plotpar.laxis = -0.1;           %lower colorbar boundary
plotpar.fontsize = 14;          %fontsize
plotpar.colorbar = 'y';         %set 'y' if you want a colorbar
plotpar.edge = 1;               %set 1 if you want to plot edges
plotpar.edgewidth = 1;          %width of edges

%Choose interpolation nodes and values on the graph
M = 12;
IntIndex = (1:M)';
y = zeros(M,1); y(1:M/2) = ones(M/2,1);

%Kernel parameter
type = 'diffusion';
alpha = -10;

%Calculate GBF interpolant
bf = GBF_genGBF(G.U,G.Lambda, IntIndex,type,alpha);
s = GBF_itpGBF(bf, IntIndex, y);

%Plot GBF interpolant
figure
GBF_drawsignal(G.nodes,G.edges,s,plotpar);
title('GBF interpolant with diffusion kernel');
hold on;
plot( G.nodes(IntIndex,1),G.nodes(IntIndex,2),'o','color',[0, 0, 0]/255,'LineWidth',1,'MarkerSize',8)
set(gca,'XTick',[], 'YTick', [])
hold off;
