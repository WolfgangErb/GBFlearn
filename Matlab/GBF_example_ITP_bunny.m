% GBFlearn: a toolbox for graph signal interpolation
% and classification with graph basis functions (GBFs)
% (C) W. Erb 15.01.2023

% Name: GBF_example_ITP_bunny.m
% In this example script we plot different GBFs and show how GBF interpolation performs 
% in comparison to bandlimited interpolation. 

% Test scenario:
% graph: bunny (2D projection of the well-known Stanford bunny)
% kernel: variational spline, diffusion GBF and bandlimited kernel
% number of sampling nodes: N = 10
% problem: calculate and plot GBF interpolant for given sampling data

clear all
close all

% Paths
addpath(genpath('./core/'))
addpath(genpath('./data/'))

%Choose graph
G.type = 'bunny';

%Generate graph
[G.nodes,G.edges,G.A] = GBF_gengraph(G.type);

%Calculate the normalized graph Laplacian
G.N = length(G.nodes(:,1));
G.deg = sum(G.A,1);
isD = diag(1./sqrt(G.deg));
G.L = eye(G.N) - isD*G.A*isD;

%Choose plotting parameter
plotpar.MM = 1;                 %size of dots
plotpar.ub = 0.02;              %upper boundary
plotpar.lb = 0.02;              %left boundary
plotpar.uaxis = 1.25;           %upper colorbar boundary
plotpar.laxis = -0.25;          %lower colorbar boundary
plotpar.fontsize = 14;          %fontsize
plotpar.colorbar = 'y';         %set 'y' if you want a colorbar
plotpar.edge = 0;               %set 1 if you want to plot edges
plotpar.edgewidth = 1;          %width of edges

%Choose interpolation nodes and values on the graph
idxW = 17;
yW = ones(size(idxW));

%Kernel parameters for interpolation
type1 = 'diffusion';
alpha1 = 10;

type2 = 'bandlimited';
alpha2 = [0.1, 100];

type3 = 'varspline';
alpha3 = [2, 0.02];

%Calculate GBF interpolation
bf1 = GBF_genGBFeff(G.L,idxW,type1,alpha1);
bf2 = GBF_genGBFeff(G.L,idxW,type2,alpha2);
bf3 = GBF_genGBFeff(G.L,idxW,type3,alpha3);

s1 = GBF_itpGBF(bf1, idxW, yW);
s2 = GBF_itpGBF(bf2, idxW, yW);
s3 = GBF_itpGBF(bf3, idxW, yW);

% Plot the basis functions for the first three kernels

figure('Units', 'pixels', ...
'Position', [0 50 1200 400]);

subplot(2,6,[1,2,7,8]),
GBF_drawsignal(G.nodes,G.edges,s1,plotpar);
title('Diffusion GBF');
hold on;
plot( G.nodes(idxW,1),G.nodes(idxW,2),'o','color',[0, 0, 0]/255,'LineWidth',1,'MarkerSize',5)
set(gca,'XTick',[], 'YTick', [])
hold off;

subplot(2,6,[3,4,9,10]), 
GBF_drawsignal(G.nodes,G.edges,s2,plotpar);
title('Bandlimited GBF');
hold on;
plot( G.nodes(idxW,1),G.nodes(idxW,2),'o','color',[0, 0, 0]/255,'LineWidth',1,'MarkerSize',5)
set(gca,'XTick',[], 'YTick', [])
hold off;

subplot(2,6,[5,6,11,12]),
GBF_drawsignal(G.nodes,G.edges,s3,plotpar);
title('Variational spline GBF');
hold on;
plot( G.nodes(idxW,1),G.nodes(idxW,2),'o','color',[0, 0, 0]/255,'LineWidth',1,'MarkerSize',5)
set(gca,'XTick',[], 'YTick', [])
hold off;


%Choose second set of interpolation nodes and values on the graph
idxW = (1:10)';
yW = zeros(size(idxW));
yW(1:5) = ones(5,1);

%Calculate GBF interpolation of this second set with third and fourth kernel
bfvs = GBF_genGBFeff(G.L,idxW,type3,alpha3);
bfband = GBF_genGBFeff(G.L,idxW,type2,alpha2);
[svs,~] = GBF_itpGBF(bfvs,idxW,yW);
[sband,~] = GBF_itpGBF(bfband,idxW,yW);

%Plot the resulting GBF interpolant
figure('Units', 'pixels', ...
'Position', [0 50 1000 500]);

subplot(2,4,[1,2,5,6]),
plotpar.MM = 2;                     %size of dots
plotpar.uaxis = max(svs);           %upper boundary for colormap
plotpar.laxis = min(svs);           %lower boundary for colormap
GBF_drawsignal(G.nodes,G.edges,svs,plotpar);
title('Variational spline GBF interpolant');
hold on;
plot( G.nodes(idxW,1),G.nodes(idxW,2),'o','color',[0, 0, 0]/255,'LineWidth',1,'MarkerSize',8)
set(gca,'XTick',[], 'YTick', [])
hold off;

subplot(2,4,[3,4,7,8]),
plotpar.MM = 2;                       %size of dots
plotpar.uaxis = max(sband);           %upper boundary for colormap
plotpar.laxis = min(sband);           %lower boundary for colormap
GBF_drawsignal(G.nodes,G.edges,sband,plotpar);
title('Bandlimited interpolant');
hold on;
plot( G.nodes(idxW,1),G.nodes(idxW,2),'o','color',[0, 0, 0]/255,'LineWidth',1,'MarkerSize',8)
set(gca,'XTick',[], 'YTick', [])
hold off;