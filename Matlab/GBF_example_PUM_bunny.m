% GBFlearn: a toolbox for graph signal interpolation
% and classification with graph basis functions (GBFs)
% (C) W. Erb 15.01.2023

% Name: GBF_example_PUM_bunny.m
% Example script to see how a partition of unity method (PUM) can be combined
% with a GBF interpolation. 

% Test scenario:
% graph: bunny
% kernel: variational spline with s = 2 and epsilon = 0.01
% number of sampling nodes: N = 100
% PUM parameters: J=6 subdomains and overlapping radius R = 2 
% problem: calculate global GBFPUM interpolant for N = 100 sampling data

clear all; close all; format short e; 

% Paths
addpath(genpath('./core/'))
addpath(genpath('./data/'))
addpath(genpath('./PUM/'))

%Choose graph
G.type = 'bunny';
load data_bunny.mat        % loads nodeselect

%Generate graph
[G.nodes,G.edges,G.A] = GBF_gengraph(G.type);

%Calculate the normalized graph Laplacian
G.N = length(G.nodes(:,1));
G.deg = sum(G.A,1);
isD = diag(1./sqrt(G.deg));
G.L = eye(G.N) - isD*G.A*isD;

%Calculate Spectrum of Laplacian
[G.U,G.Lambda] = GBF_spectrum(G.L,'ascend');

%Choose plotting parameter
plotpar.MM = 1;                 %size of dots
plotpar.ub = 0.02;              %upper boundary
plotpar.lb = 0.02;              %left boundary
plotpar.uaxis = 1.15;           %upper colorbar boundary
plotpar.laxis = -0.15;          %lower colorbar boundary
plotpar.fontsize = 14;          %fontsize
plotpar.colorbar = 'n';         %set 'y' if you want a colorbar
plotpar.edge = 0;               %set 1 if you want to plot edges
plotpar.edgewidth = 1;          %width of edges

%Extract bandlimited function u4
f = G.U(:,4);

% Choose number of interpolation nodes and sampling data
N = 100;                        % Number of interpolation nodes
idxW = nodeselect(1:N)'; y = f(idxW);

% PUM parameters
JJ = 6;               % Number of clusters
RR = 2;               % Increasing parameter for subdomains

% Kernel parameters for interpolation (variational splines)
type = 'varspline';
alpha = [2,0.01];

% Generate Partition of Unity

% 1a) Calculate J-center clustering of graph in JJ components
[idxcluster,idxQ] = GBF_Jcenters_greedy(G.A,JJ,idxW,idxW(1));

% 1b) Augment the clusters to subdomains for PU
[idxdomain] = GBF_domainaugment(G.edges,idxcluster,RR);

% 1c) Generate PU out of the single domains
[phi,idxWdomain,ydomain] = GBF_genPUM(G.edges,idxcluster,idxdomain,idxW,y);

% 2)3)4) Use PU to calculate global kernel-PUM interpolant
s = GBF_RLSGBFPUM(G.L,idxdomain,phi,idxWdomain,ydomain,type,alpha,0);

% Compute errors
rmaerr = norm(s(:)-f,inf)/norm(f,inf);            % relative max absolute error
rrmserr = norm((s(:)-f))/norm(f)/sqrt(length(f));  % relative rms error

% Print errors
fprintf('no. points\t rmaerr \t rrmserr \n')
fprintf('\t%4d\t %.3e\t %.3e\n ',N,rmaerr,rrmserr);

% Plot partition of unity
for k = 1:JJ
  figure('Units', 'pixels', 'Position', [0 50 600 600]);
  GBF_drawsignal(G.nodes,G.edges,phi(:,k),plotpar);
  title([num2str(k),'. Subdomain'])
  hold on
  inddiff = setdiff((1:G.N)',idxdomain{k}); 
  plot( G.nodes(inddiff,1),G.nodes(inddiff,2),'.','color',[208,209,230]/255,'LineWidth',1,'MarkerSize',10)
  hold on
  plot( G.nodes(idxdomain{k}(idxWdomain{k}),1),G.nodes(idxdomain{k}(idxWdomain{k}),2),'o','color',[0, 0, 0]/255,'LineWidth',1,'MarkerSize',5)
  hold on
  plot( G.nodes(idxQ(k),1),G.nodes(idxQ(k),2),'o','color','r','LineWidth',2,'MarkerSize',8)
  set(gca,'XTick',[], 'YTick', [])
  hold off
end

plotpar.colorbar = 'y';           %put 'y' if you want a colorbar
plotpar.uaxis = max(s);           %upper colorbar boundary
plotpar.laxis = min(s);           %lower colorbar boundary

figure('Units', 'pixels', 'Position', [0 50 600 600]);
GBF_drawsignal(G.nodes,G.edges,s,plotpar);
title('Global GBF-PUM Interpolation')
hold on
plot( G.nodes(idxW,1),G.nodes(idxW,2),'o','color',[0, 0, 0]/255,'LineWidth',1,'MarkerSize',5)
hold on
plot( G.nodes(idxQ,1),G.nodes(idxQ,2),'o','color','r','LineWidth',2,'MarkerSize',8)
set(gca,'XTick',[], 'YTick', [])
hold off