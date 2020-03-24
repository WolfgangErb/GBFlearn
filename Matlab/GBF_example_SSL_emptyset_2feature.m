% GBFlearn: a toolbox for graph signal interpolation
% and classification with graph basis functions (GBFs)
% (C) W. Erb 01.03.2020

% GBF_example_SSL_emptyset_2feature shows how SSL based on a
% feature-augmented GBF-RLS solver can be used to improve
% classification accuracy in an example where geometric features are
% available. The test data set hast the form of a slashed O. 

clear all
close all

% Paths
addpath(genpath('./core/'))
addpath(genpath('./data/'))

%Choose graph
G.type = 'emptyset';
load data_emptyset.mat;

%Generate graph
[G.nodes,G.edges,G.A] = GBF_gengraph(G.type);

%Calculate the graph Laplacian
G.N = length(G.nodes(:,1));
G.deg = sum(G.A,1);
isD = diag(1./sqrt(G.deg));
G.L = eye(G.N) - isD*G.A*isD;

%Calculate Spectrum of graph
[G.U,G.Lambda] = GBF_spectrum(G.L,'ascend');

%Choose plotting parameter
plotpar.MM = 1;                 %size of dots
plotpar.ub = 0.05;              %upper boundary
plotpar.lb = 0.05;              %left boundary
plotpar.uaxis = 1;              %upper colormap boundary
plotpar.laxis = -1;             %lower colormap boundary
plotpar.fontsize = 14;          %fontsize
plotpar.colorbar = 'n';         %set 'y' if you want a colorbar
plotpar.edge = 0;               %set 1 if you want to plot edges
plotpar.edgewidth = 1;          %width of edges

%Choose interpolation nodes and values on the graph
M = 18;
idxW = (1:M)';
yW = label(1:M);

%Kernel parameters
type = 'diffusion';             %Type of GBF
alpha = -5;                     %Shape parameter
lambda = 0.001;                 %Regularization parameter
gamma1 = 10;                    %Similarity parameter feature 1
gamma2 = 10;                    %Similarity parameter feature 2

%Calculate GBF-RLS solution
bf = GBF_genGBF(G.U,G.Lambda, idxW, type, alpha);
s = GBF_RLSGBF(bf, idxW, yW,lambda);

% Extract geometric features
idxplus = find(yW==1); idxminus = find(yW==-1);

X1 = nodes(idxplus,:);
X2 = [nodes(idxminus,1),ones(numel(nodes(idxminus,1)),1)];

z1 = [X1, ones(size(X1,1), 1)]\sum(X1.*X1, 2);
ccircle = 0.5*z1(1:2);
rcircle = sqrt(z1(3) + ccircle'*ccircle);

z2 = (X2'*X2)\(X2'*nodes(idxminus,2));
aline = z2(1);
bline = z2(2);

%1. geometric feature: distance to reference circle
dist1 = nodes - ones(G.N,1)*ccircle';
dist1 = abs(sqrt(dist1(:,1).^2 + dist1(:,2).^2)-rcircle);

%2. geometric feature: distance to reference line
dist2 = abs(nodes(:,2) - aline*nodes(:,1) - bline);

%Joint feature information in similarity graph
scut = [2*exp(-20.*dist1.^2)-1,2*exp(-5.*dist2.^2)-1];

%Calculate PSI-GBF-RLS solution for 1 and 2 additional features
simK1 = GBF_genexpK(idxW,scut(:,1),gamma1);
simK2 = GBF_genexpK(idxW,scut(:,2),gamma2);

sAdd1 = GBF_RLSGBF(bf.*simK1, idxW, yW, lambda);
sAdd2 = GBF_RLSGBF(bf.*simK1.*simK2, idxW, yW, lambda);

%Calculate RLS classifications

idxsup = find(s>=0);
idxsdown = find(s<0);
idxsAdd2up = find(sAdd2>=0);
idxsAdd2down = find(sAdd2<0);
idxsAdd1up = find(sAdd1>=0);
idxsAdd1down = find(sAdd1<0);

sclass(idxsup) = 1;
sclass(idxsdown) = -1;
sAdd1class(idxsAdd1up) = 1;
sAdd1class(idxsAdd1down) = -1;
sAdd2class(idxsAdd2up) = 1;
sAdd2class(idxsAdd2down) = -1;

%Plot original classification and the two applied geometric features

figure('Units', 'pixels', ...
'Position', [0 50 900 300]);

subplot(2,6,[1,2,7,8]),
GBF_drawsignal(G.nodes,G.edges,label,plotpar);
title('Original');
set(gca,'XTick',[], 'YTick', [])
hold off;

subplot(2,6,[3,4,9,10]),
GBF_drawsignal(G.nodes,G.edges,scut(:,1),plotpar);
title('Geometric feature 1');
hold on;
plot( G.nodes(idxW,1),G.nodes(idxW,2),'o','color',[0, 0, 0]/255,'LineWidth',1,'MarkerSize',5)
set(gca,'XTick',[], 'YTick', [])
hold off;

subplot(2,6,[5,6,11,12]),
GBF_drawsignal(G.nodes,G.edges,scut(:,2),plotpar);
title('Geometric feature 2');
hold on;
plot( G.nodes(idxW,1),G.nodes(idxW,2),'o','color',[0, 0, 0]/255,'LineWidth',1,'MarkerSize',5)
set(gca,'XTick',[], 'YTick', [])
hold off;

%Plot supervised classification and two semi-supervised classifications
%based on one and two additional features

figure('Units', 'pixels', ...
'Position', [0 50 900 300]);

subplot(2,6,[1,2,7,8]),
GBF_drawsignal(G.nodes,G.edges,sclass,plotpar);
title('GBF-RLS classification');
hold on;
plot( G.nodes(idxW,1),G.nodes(idxW,2),'o','color',[0, 0, 0]/255,'LineWidth',1,'MarkerSize',5)
set(gca,'XTick',[], 'YTick', [])
hold off;

subplot(2,6,[3,4,9,10]),
GBF_drawsignal(G.nodes,G.edges,sAdd1class,plotpar);
title('\psi-GBF-RLS class. (1 feat.)');
hold on;
plot( G.nodes(idxW,1),G.nodes(idxW,2),'o','color',[0, 0, 0]/255,'LineWidth',1,'MarkerSize',5)
set(gca,'XTick',[], 'YTick', [])
hold off;

subplot(2,6,[5,6,11,12]),
GBF_drawsignal(G.nodes,G.edges,sAdd2class,plotpar);
title('\psi-GBF-RLS class. (2 feat.)');
hold on;
plot( G.nodes(idxW,1),G.nodes(idxW,2),'o','color',[0, 0, 0]/255,'LineWidth',1,'MarkerSize',5)
set(gca,'XTick',[], 'YTick', [])
hold off;