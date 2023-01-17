% GBFlearn: a toolbox for graph signal interpolation
% and classification with graph basis functions (GBFs)
% (C) W. Erb 15.01.2023

% Name: GBF_example_EFF_verylargesensor.m
% Example script for the usage of efficient block Krylov space methods on a
% very large sensor network

% Test scenario:
% graph: very large sensor network (n = 25000)
% kernel: diffusion with parameter alpha = 80
% number of sampling nodes: M = 200
% problem: calculate kernel classificator using a RLS solution for
% labeled data on M = 200 sampling nodes 
% Regularization parameter: gamma = 0.001

clear all, close all

%Paths
addpath(genpath('./core/'))
addpath(genpath('./data/'))

%Choose graph
G.type = 'sensorverylarge';

%Generate graph
tic;
[G.nodes,G.edges,G.A] = GBF_gengraph(G.type);
                             
fprintf(1, 'Time to generate very large graph: '); fprintf(1,'%5f \n', toc); 
%Calculate the graph Laplacian
tic;
G.N = length(G.nodes(:,1));
G.deg = sum(G.A,1)';
isD = spdiags(1./sqrt(G.deg),0,G.N,G.N);
G.L = speye(G.N) - isD*G.A*isD;
fprintf(1, 'Time to generate graph Laplacian:  '); fprintf(1,'%5f \n', toc); 

%Choose plotting parameters
plotpar.MM = 1;                 %size of dots
plotpar.ub = 0.1;               %upper boundary
plotpar.lb = 0.1;               %left boundary
plotpar.fontsize = 14;          %fontsize
plotpar.uaxis = 1;              %upper boundary for colormap
plotpar.laxis = 0;              %lower boundary for colormap
plotpar.colorbar = 'n';         %put 'y' if you want a colorbar
plotpar.edge = 0;               %put 1 if you want to plot edges
plotpar.edgewidth = 0.5;        %width of edges

%Choose labeled nodes of the graph
M = 200;
idxW = randsample(G.N,M);
f = 2*(G.nodes(:,1) > 0.5)-1;
y = f(idxW);

%Kernel parameters for interpolation
type = 'diffusion';
alpha = [80,20,1e-8];
method = 'sbl';
lambda = 0.001;
   
% Calculate GBF approximant
tic; bf = GBF_genGBFeff(G.L,idxW,type,alpha,method);
s  = GBF_RLSGBF(bf, idxW, y,lambda);
fprintf(1, 'Time to calculate classfication:   '); fprintf(1,'%5f \n', toc);
    
% Calculate classification based on GBF approximation
idxsup = find(s>=0); idxsdown = find(s<0);
sclass(idxsup) = 1; sclass(idxsdown) = -1;

% Plots
plotpar.uaxis = max(s);               %upper boundary for colormap
plotpar.laxis = min(s);              %lower boundary for colormap

figure('Units', 'pixels', 'Position', [0 50 600 600]);
GBF_drawsignal(G.nodes,G.edges,s,plotpar);
title('Kernel approximation based on GBFRLS')
hold on
plot( G.nodes(idxW,1),G.nodes(idxW,2),'o','color',[0, 0, 0]/255,'LineWidth',1,'MarkerSize',5)
% set(gca,'XTick',[], 'YTick', [])
hold off

plotpar.uaxis = 1;               %upper boundary for colormap
plotpar.laxis = -1;              %lower boundary for colormap

figure('Units', 'pixels', 'Position', [0 50 600 600]);
GBF_drawsignal(G.nodes,G.edges,sclass,plotpar);
title('Kernel classification based on GBFRLS')
hold on
plot( G.nodes(idxW,1),G.nodes(idxW,2),'o','color',[0, 0, 0]/255,'LineWidth',1,'MarkerSize',5)
% set(gca,'XTick',[], 'YTick', [])
hold off
