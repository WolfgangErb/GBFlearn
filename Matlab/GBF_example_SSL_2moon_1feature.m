% GBFlearn: a toolbox for graph signal interpolation
% and classification with graph basis functions (GBFs)
% (C) W. Erb 01.03.2020

% GBF_example_SSL_2moon_1feature shows how SSL based on a
% feature-augmented GBF-RLS solution can be used to improve
% classification accuracy of a supervised scheme.
% The included feature is a unsupervised classification of the graph 
% based on spectral clustering. The test data set consists of two 
% point clouds in the form of two nested moons.  

clear all
close all

%Paths
addpath(genpath('./core/'))
addpath(genpath('./data/'))

%Choose graph
G.type = '2moon';
load data_2moon.mat;

%Generate graph
[G.nodes,G.edges,G.A] = GBF_gengraph(G.type);

%Calculate the normalized graph Laplacian
G.N = length(G.nodes(:,1));
G.deg = sum(G.A,1);
isD = diag(1./sqrt(G.deg));
G.L = eye(G.N) - isD*G.A*isD;

%Calculate spectrum of graph
[G.U,G.Lambda] = GBF_spectrum(G.L,'ascend');

%Choose plotting parameters
plotpar.MM = 0.5;               %size of dots
plotpar.ub = 0.1;               %upper boundary
plotpar.lb = 0.1;               %left boundary
plotpar.fontsize = 14;          %fontsize
plotpar.uaxis = 1;              %upper boundary for colormap
plotpar.laxis = -1;             %lower boundary for colormap
plotpar.colorbar = 'n';         %put 'y' if you want a colorbar
plotpar.edge = 0;               %put 1 if you want to plot edges
plotpar.edgewidth = 0.5;        %width of edges

%Choose labeled nodes of the graph
idxW = [5,6]';
yW = label(idxW);

%Kernel parameter
type = 'diffusion';             %Type of GBF
alpha = -50;                    %Shape parameter of kernel
lambda = 0.0001;                %Regularization parameter
gamma = -1;                     %Shape parameter of feature kernel

%Calculate standard GBF-RLS solution
bf = GBF_genGBF(G.U,G.Lambda, idxW,type,alpha);
s = GBF_RLSGBF(bf, idxW, yW,lambda);

%Calculate (unsupervised) classification based on spectral clustering
%Corresponds to the normalized cut of Shi and Malik
%Instead of the median you can also take the mean or simply cut = 0. 
scut = G.U(:,2); 
cut = median(scut);

idxcutup = find(scut>=cut); idxcutdown = find(scut<cut);
scut(idxcutup)=1; scut(idxcutdown) = -1;

%Calculate PSI-GBF-RLS solution
binK = GBF_genbinK(idxW,scut,gamma);
[sPSI,~] = GBF_RLSGBF(bf.*binK, idxW, yW, lambda);

%Calculate classifications
idxsup = find(s>=0); idxsdown = find(s<0);
idxsPSIup = find(sPSI>=0); idxsPSIdown = find(sPSI<0);

sclass(idxsup) = 1; sclass(idxsdown) = -1;
sPSIclass(idxsPSIup) = 1; sPSIclass(idxsPSIdown) = -1;

%Plot 1: supervised classification

figure('Units', 'pixels', ...
'Position', [0 50 1200 300]);

subplot(2,8,[1,2,9,10]),
GBF_drawsignal(G.nodes,G.edges,label,plotpar);
title('Original');
set(gca,'XTick',[], 'YTick', [])
hold off;

subplot(2,8,[3,4,11,12]),
GBF_drawsignal(G.nodes,G.edges,scut,plotpar);
title('Spectral classification');
set(gca,'XTick',[], 'YTick', [])
hold off;

subplot(2,8,[5,6,13,14]),
GBF_drawsignal(G.nodes,G.edges,s,plotpar);
title('GBF-RLS solution');
hold on;
plot( G.nodes(idxW,1),G.nodes(idxW,2),'o','color',[0, 0, 0]/255,'LineWidth',1,'MarkerSize',5)
set(gca,'XTick',[], 'YTick', [])
hold off;

subplot(2,8,[7,8,15,16]),
GBF_drawsignal(G.nodes,G.edges,sclass,plotpar);
title('GBF-RLS classification');
hold on;
plot( G.nodes(idxW,1),G.nodes(idxW,2),'o','color',[0, 0, 0]/255,'LineWidth',1,'MarkerSize',5)
set(gca,'XTick',[], 'YTick', [])
hold off;

%Plot 2: semi-supervised classification with feature-augmented GBF

figure('Units', 'pixels', ...
'Position', [0 50 1200 300]);

subplot(2,8,[1,2,9,10]),
GBF_drawsignal(G.nodes,G.edges,label,plotpar);
title('Original');
set(gca,'XTick',[], 'YTick', [])
hold off;

subplot(2,8,[3,4,11,12]),
GBF_drawsignal(G.nodes,G.edges,scut,plotpar);
title('Spectral classification');
set(gca,'XTick',[], 'YTick', [])
hold off;

subplot(2,8,[5,6,13,14]),
GBF_drawsignal(G.nodes,G.edges,sPSI,plotpar);
title('\psi-GBF-RLS solution');
hold on;
plot( G.nodes(idxW,1),G.nodes(idxW,2),'o','color',[0, 0, 0]/255,'LineWidth',1,'MarkerSize',5)
set(gca,'XTick',[], 'YTick', [])
hold off;

subplot(2,8,[7,8,15,16]),
GBF_drawsignal(G.nodes,G.edges,sPSIclass,plotpar);
title('\psi-GBF-RLS classification');
hold on;
plot( G.nodes(idxW,1),G.nodes(idxW,2),'o','color',[0, 0, 0]/255,'LineWidth',1,'MarkerSize',5)
set(gca,'XTick',[], 'YTick', [])
hold off;

