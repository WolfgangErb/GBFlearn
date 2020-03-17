% GBFlearn: a toolbox for graph signal interpolation
% and classification with graph basis functions (GBFs)
% (C) W. Erb 01.03.2020

% GBF_example_SSL_2moon_2feature shows how SSL based on a
% feature-augmented GBF-RLS solution can be used to improve
% classification accuracy of a supervised scheme.
% The two included features are an unsupervised classification result  
% based on spectral clustering and a given left-right separation 
% of the data. The test data set consists of two 
% point clouds in the form of nested moons.  

clear all
close all

% Paths
addpath(genpath('./core/'))
addpath(genpath('./data/'))

%Choose graph
G.type = '2moon';
load data_2moon.mat;

%Generate graph
[G.nodes,G.edges,G.A] = GBF_gengraph(G.type);

%Calculate the graph Laplacian
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
plotpar.uaxis = 3;              %upper boundary for colormap
plotpar.laxis = -3;             %lower boundary for colormap
plotpar.colorbar = 'n';         %put 'y' if you want a colorbar
plotpar.edge = 0;               %put 1 if you want to plot edges
plotpar.edgewidth = 0.5;        %width of edges

%Choose labeled nodes of the graph
IntIndex = [6,7,11,12]';
y = [label(IntIndex),sign(G.nodes(IntIndex,1))];

%Kernel parameter
type = 'diffusion';             %Type of GBF
alpha = -50;                    %Shape parameter of kernel
lambda = 0.0001;                %Regularization parameter
gamma = -1;                     %Shape parameter of 1. feature kernel
gamma2 = 0.1;                   %Shape parameter of 2. feature kernel

%Calculate standard GBF-RLS solution
bf = GBF_genGBF(G.U,G.Lambda,IntIndex,type,alpha);
s = GBF_RLSGBF(bf,IntIndex,y,lambda);

%Initiate 1. feature based on spectral clustering
%(Shi-Malik normalized cut using the median)
scut = G.U(:,2);
cut = median(scut);

%Initiate 2. feature based on given left-right-prior
scut2 = G.nodes(:,1);
cut2 = mean(scut2);

%Calculate binary classification features
idxcutup = find(scut>=cut); idxcutup2 = find(scut2>=cut2);
idxcutdown = find(scut<cut); idxcutdown2 = find(scut2<cut2);

scut(idxcutup)=1; scut(idxcutdown) = -1;
scut2(idxcutup2)=1; scut2(idxcutdown2) = -1;

%Calculate PSI-GBF-RLS solution
binK  = GBF_genbinK(IntIndex,scut,gamma);
binK2 = GBF_genbinK(IntIndex,scut2,gamma2);

sPSI = GBF_RLSGBF(bf.*binK.*binK2, IntIndex, y, lambda);

%Calculate corresponding classifications
idxlab1 = find((label(:,1)>=0).*(G.nodes(:,1)<0));
idxlab2 = find((label(:,1)>=0).*(G.nodes(:,1)>=0));
idxlab3 = find((label(:,1)<0).*(G.nodes(:,1)<0));
idxlab4 = find((label(:,1)<0).*(G.nodes(:,1)>=0));

idxs1 = find((s(:,1)>=0).*(s(:,2)<0));
idxs2 = find((s(:,1)>=0).*(s(:,2)>=0));
idxs3 = find((s(:,1)<0).*(s(:,2)<0));
idxs4 = find((s(:,1)<0).*(s(:,2)>=0));

idxsPSI1 = find((sPSI(:,1)>=0).*(sPSI(:,2)<0));
idxsPSI2 = find((sPSI(:,1)>=0).*(sPSI(:,2)>=0));
idxsPSI3 = find((sPSI(:,1)<0).*(sPSI(:,2)<0));
idxsPSI4 = find((sPSI(:,1)<0).*(sPSI(:,2)>=0));

labclass(idxlab1) = 3;
labclass(idxlab2) = 1;
labclass(idxlab3) = -1;
labclass(idxlab4) = -3;
sclass(idxs1) = 3;
sclass(idxs2) = 1;
sclass(idxs3) = -1;
sclass(idxs4) = -3;
sPSIclass(idxsPSI1) = 3;
sPSIclass(idxsPSI2) = 1;
sPSIclass(idxsPSI3) = -1;
sPSIclass(idxsPSI4) = -3;

%Plot: comparison of supervised and semi-supervised classification

figure('Units', 'pixels', ...
'Position', [0 50 900 300]);

colormap(flipud(copper(4)));
labels = {'(-1,-1)','(-1,1)','(1,-1)','(1,1)'};
h = lcolorbar(labels,'fontweight','bold');
set(h,'Position',[0.92 0.365 0.02 0.3])

subplot(2,6,[1,2,7,8]),
GBF_drawsignal(G.nodes,G.edges,labclass,plotpar);
title('Original');
set(gca,'XTick',[], 'YTick', [])
hold off;

subplot(2,6,[3,4,9,10]),
GBF_drawsignal(G.nodes,G.edges,sclass,plotpar);
title('GBF-RLS classification');
hold on;
plot( G.nodes(IntIndex,1),G.nodes(IntIndex,2),'o','color',[0, 0, 0]/255,'LineWidth',1,'MarkerSize',5)
set(gca,'XTick',[], 'YTick', [])
hold off;

subplot(2,6,[5,6,11,12]),
GBF_drawsignal(G.nodes,G.edges,sPSIclass,plotpar);
title('\psi-GBF-RLS classification');
hold on;
plot( G.nodes(IntIndex,1),G.nodes(IntIndex,2),'o','color',[0, 0, 0]/255,'LineWidth',1,'MarkerSize',5)
set(gca,'XTick',[], 'YTick', [])
hold off;

