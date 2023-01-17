% GBFlearn: a toolbox for graph signal interpolation
% and classification with graph basis functions (GBFs)
% (C) W. Erb 15.01.2023

% Name: GBF_example_EFF_path_varspline_cheb.m
% Test script to compare kernel approximations using Chebyshev polynomials

% Test scenario:
% graph: path
% kernel: variational spline
% number of sampling nodes: M = 1 (center of paht)
% problem: calculate kernel column corresponding to center node of path
% used methods: Chebyshev methods with different degree m

clear all

%Paths
addpath(genpath('./core/'))
addpath(genpath('./data/'))

%Choose graph
G.type = 'path';

%Generate graph
tic;
[G.nodes,G.edges,G.A] = GBF_gengraph(G.type);
fprintf(1, 'Time to generate Path graph:       '); fprintf(1,'%5f \n', toc); 
%Calculate the graph Laplacian
tic;
G.N = length(G.nodes(:,1));
G.deg = sum(G.A,1)';
isD = spdiags(1./sqrt(G.deg),0,G.N,G.N);
G.L = speye(G.N) - isD*G.A*isD;
fprintf(1, 'Time to generate graph Laplacian:  '); fprintf(1,'%5f \n', toc); 

%Choose labeled nodes of the graph
M = 101;               %center node of path graph
idxW = M;              
y = ones(size(idxW));

% Iteration numbers 
KK2 = 120;
KK3 = 80;
KK4 = 40;
% Kernel parameters (variational spline)
s = 2;
epsilon = 0.001;
tol = 1e-14;

type1 = 'varspline';
alpha1 = [s,epsilon];
method1 = 'direct';
type2 = 'varspline';
alpha2 = [s,epsilon,KK2];
method2 = 'cheb';
type3 = 'varspline';
alpha3 = [s,epsilon,KK3];
method3 = 'cheb';
type4 = 'varspline';
alpha4 = [s,epsilon,KK4];
method4 = 'cheb';
lambda = 0;

% Calculate GBF approximant (fast Chebyshev iterative method for diffusion)

bf1 = GBF_genGBFeff(G.L,idxW,type1,alpha1);
bf2 = GBF_genGBFeff(G.L,idxW,type2,alpha2,method2);
bf3 = GBF_genGBFeff(G.L,idxW,type3,alpha3,method3);
bf4 = GBF_genGBFeff(G.L,idxW,type4,alpha4,method4);

s1  = GBF_RLSGBF(bf1, idxW, y,lambda);
s2  = GBF_RLSGBF(bf2, idxW, y,lambda);
s3  = GBF_RLSGBF(bf3, idxW, y,lambda);
s4  = GBF_RLSGBF(bf4, idxW, y,lambda);

error2 = norm(s2-s1,inf);
error3 = norm(s3-s1,inf);
error4 = norm(s4-s1,inf);

% Plots
figure('Units', 'pixels', 'Position', [0 50 1200 300]);
subplot(1,4,1),
h1 = stem(G.nodes(:,1),s1,'k-.','linewidth',0.2);
set(h1,'MarkerFaceColor','red')
set(h1,'markersize',2)
axis square;
axis([0,1,-0.1,1.1]);
set(gca,'XTick',[0, 0.5, 1])
set(gca,'XTickLabel',{'1','101','201'})
xlabel('Exact variational spline');
subplot(1,4,2),
h2 = stem(G.nodes(:,1),s2,'k-.','linewidth',0.2);
set(h2,'MarkerFaceColor','red')
set(h2,'markersize',2)
hold on
line([0.5-KK2/200 0.5+KK2/200], [-0.05 -0.05],'color','b','linewidth',2);
hold off
axis square;
axis([0,1,-0.1,1.1]);
set(gca,'XTick',[0, 0.5, 1])
set(gca,'XTickLabel',{'1','101','201'})
xlabel(['m = ',num2str(KK2), ', error = ', num2str(error2)]);
subplot(1,4,3),
h3 = stem(G.nodes(:,1),s3,'k-.','linewidth',0.2);
set(h3,'MarkerFaceColor','red')
set(h3,'markersize',2)
hold on
line([0.5-KK3/200 0.5+KK3/200], [-0.05 -0.05],'color','b','linewidth',2);
hold off
axis square;
axis([0,1,-0.1,1.1]);
set(gca,'XTick',[0, 0.5, 1])
set(gca,'XTickLabel',{'1','101','201'})
xlabel(['m = ',num2str(KK3), ', error = ', num2str(error3)]);
subplot(1,4,4),
h4 = stem(G.nodes(:,1),s4,'k-.','linewidth',0.2);
set(h4,'MarkerFaceColor','red')
set(h4,'markersize',2)
hold on
line([0.5-KK4/200 0.5+KK4/200], [-0.05 -0.05],'color','b','linewidth',2);
hold off
axis square;
axis([0,1,-0.1,1.1]);
set(gca,'XTick',[0, 0.5, 1])
set(gca,'XTickLabel',{'1','101','201'})
xlabel(['m = ',num2str(KK4), ', error = ', num2str(error4)]);
% title('Different GBFs on path graph')
% set(gca,'XTick',[], 'YTick', [])
hold off