% GBFlearn: a toolbox for graph signal interpolation
% and classification with graph basis functions (GBFs)
% (C) W. Erb 15.01.2023

% Name: GBF_example_EFF_minnesota.m
% Test script to compare kernel calculation using an efficent block Krylov solver and 
% an inefficient eigenvalue decomposition of the graph Laplacian 

% Test scenario:
% graph: minnesota
% kernel: diffusion
% number of sampling nodes: M = 200 (random)
% problem: calculate kernel predictor using regularized least squares
% solution for a data set on M = 200 sampling nodes
% regularization parameter: gamma = 0.001
% used methods: global block Krylov vs eigenvalue decomposition of L

clear all, close all;

%Paths
addpath(genpath('./core/'))
addpath(genpath('./data/'))

%Choose graph
G.type = 'minnesota';

%Generate graph
tic; [G.nodes,G.edges,G.A] = GBF_gengraph(G.type);
fprintf(1, 'Time to generate Minnesota graph:  '); fprintf(1,'%5f \n', toc); 

%Calculate the graph Laplacian

tic; G.N = length(G.nodes(:,1));
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
plotpar.colorbar = 'y';         %put 'y' if you want a colorbar
plotpar.edge = 0;               %put 1 if you want to plot edges
plotpar.edgewidth = 0.5;        %width of edges

%Choose labeled nodes of the graph
M = 200;
idxW = randsample(G.N,M);
f = 2*((G.nodes(:,1)+92).^2+ (G.nodes(:,2)-44).^2 < 8)-1;
y = f(idxW);

%Kernel parameters for interpolation
type = 'diffusion';
alpha = [20,50,1e-14];
lambda = 0.001;

% Calculate GBF approximant (slow method for diffusion)
tic; bf1 = GBF_genGBF(G.L,idxW,type,alpha);
fprintf(1, 'Time to calculate basis functions (slow): '); fprintf(1,'%5f \n', toc);

% Calculate GBF approximant (fast iterative method for diffusion)
tic; bf2 = GBF_genGBFeff(G.L,idxW,type,alpha);
fprintf(1, 'Time to calculate basis functions (fast): '); fprintf(1,'%5f \n', toc);

tic; s1  = GBF_RLSGBF(bf1, idxW, y,lambda);
s2  = GBF_RLSGBF(bf2, idxW, y,lambda);

% Calculate classification based on GBF-RLS solution
sclass1 = sign(s1);
sclass2 = sign(s2);

fprintf(1, 'Time to calculate classficiation:  '); fprintf(1,'%5f \n', toc);

% Plots
plotpar.uaxis = max(s1);              %upper boundary for colormap
plotpar.laxis = min(s1);              %lower boundary for colormap

figure('Units', 'pixels', 'Position', [0 50 600 600]);
GBF_drawsignal(G.nodes,G.edges,s1,plotpar);
title('RLS solution using diffusion GBF')
hold on
plot( G.nodes(idxW,1),G.nodes(idxW,2),'o','color',[0, 0, 0]/255,'LineWidth',1,'MarkerSize',5)
% set(gca,'XTick',[], 'YTick', [])
hold off

plotpar.uaxis = 1;               %upper boundary for colormap
plotpar.laxis = -1;              %lower boundary for colormap

figure('Units', 'pixels', 'Position', [0 50 600 600]);
GBF_drawsignal(G.nodes,G.edges,sclass1,plotpar);
title('Classification based on RLS solution')
hold on
plot( G.nodes(idxW,1),G.nodes(idxW,2),'o','color',[0, 0, 0]/255,'LineWidth',1,'MarkerSize',5)
% set(gca,'XTick',[], 'YTick', [])
hold off

% Compute errors between the two solutions
rmaerr = norm(s1-s2,inf)/norm(s1,inf);            % relative max absolute error
rrmserr = norm(s1-s2)/norm(s1)/sqrt(length(s1));  % relative rms error

% Print errors 
fprintf('no. points\t rmaerr \t rrmserr \n')
fprintf('\t%4d\t %.3e\t %.3e\n ',M,rmaerr,rrmserr);