% GBFlearn: a toolbox for graph signal interpolation
% and classification with graph basis functions (GBFs)
% (C) W. Erb 15.01.2023

% Name: GBF_example_EFF_bunny_krylov_symmetry.m
% Test script to check symmetry and positive definiteness of the
% collocation matrices when using Krylov subspace methods for the
% approximation of the kernels. 

% Test scenario:
% graph: bunny
% kernel: diffusion with parameter alpha = 20
% number of sampling nodes: N = 40
% number of Krylov space iterations: m = 6

clear all, close all

% Paths
addpath(genpath('./core/'))
addpath(genpath('./data/'))

%Choose graph
G.type = 'bunny';
load data_bunny.mat 

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
plotpar.edge = 0;               %set 1 if you want to plot edges
plotpar.edgewidth = 1;          %width of edges

%Initialize
MM = 6;             %number of iterations 
NN = 40;            %number of nodes
idxW = nodeselect(1:NN)';
yW = ones(size(idxW));

%Kernel parameters for interpolation
alpha = 20;         %diffusion parameter
tol = 1e-20;

type1 = 'diffusion';
method1 = 'cbl';
type2 = 'diffusion';
method2 = 'gbl';
type3 = 'diffusion';
method3 = 'sbl';
type4 = 'diffusion';
method4 = 'cheb';
type5 = 'diffusion';
method5 = 'cheb2';

%Start evaluation   
bf1 = GBF_genGBFeff(G.L,idxW,type1,[alpha,MM,tol],method1);
bf2 = GBF_genGBFeff(G.L,idxW,type2,[alpha,MM,tol],method2);
bf3 = GBF_genGBFeff(G.L,idxW,type3,[alpha,MM,tol],method3);
bf4 = GBF_genGBFeff(G.L,idxW,type4,[alpha,MM,tol],method4);
bf5 = GBF_genGBFeff(G.L,idxW,type5,[alpha,MM,tol],method5);

[~,~,KW1] = GBF_itpGBF(bf1, idxW,yW);
[~,~,KW2] = GBF_itpGBF(bf2, idxW,yW);
[~,~,KW3] = GBF_itpGBF(bf3, idxW,yW);
[~,~,KW4] = GBF_itpGBF(bf4, idxW,yW);
[~,~,KW5] = GBF_itpGBF(bf5, idxW,yW);

error1 = norm(KW1-KW1')/norm(KW1);
error2 = norm(KW2-KW2')/norm(KW2);
error3 = norm(KW3-KW3')/norm(KW3);
error4 = norm(KW4-KW4')/norm(KW4);
error5 = norm(KW5-KW5')/norm(KW5);

mineig1 = eig(KW1);
mineig2 = eig(KW2);
mineig3 = eig(KW3);
mineig4 = eig(KW4);
mineig5 = eig(KW5);

fprintf('Number of Krylov iterations m = %2d \n',MM);
fprintf('Number of nodes (block size) N = %2d \n',NN);
fprintf('Symmetry error for sbl in iteration step %2d: %20.16f \n',MM,error3);

%%
figure('Units', 'pixels', ...
    'Position', [50 50 600 600]);
hold on;

hker1 = plot(real(mineig1),imag(mineig1),'-k');
hker2 = plot(real(mineig2),imag(mineig2),'-r');
hker3 = plot(real(mineig3),imag(mineig3),'-g');
hker4 = plot(real(mineig4),imag(mineig4),'-b');
hker5 = plot(real(mineig5),imag(mineig5),'-m');

set(hker1                         , ...
  'Linestyle'          , 'none'      , ...
  'Marker'          , 'd'      , ...
  'Markersize'       , 12           , ...
  'MarkerEdgeColor'           , [0,0,0]/255 , ...
  'MarkerFaceColor'           , [0,0,0]/255    );
set(hker2                         , ...
  'Linestyle'          , 'none'      , ...
  'Marker'          , 'o'      , ...
  'Markersize'       , 11           , ...
  'MarkerEdgeColor'           , [153,52,4]/255 , ...
  'MarkerFaceColor'           , [153,52,4]/255    );
set(hker3                         , ...
  'Linestyle'          , 'none'      , ...
  'Marker'          , 'd'      , ...
  'Markersize'       , 6           , ...
  'MarkerEdgeColor'           , [217,95,14]/255 , ...
  'MarkerFaceColor'           , [217,95,14]/255    );
set(hker4                         , ...
  'Linestyle'          , 'none'      , ...
  'Marker'          , 'o'      , ...
  'Markersize'       , 5           , ...
  'MarkerEdgeColor'           , [254,153,41]/255 , ...
  'MarkerFaceColor'           , [254,153,41]/255    );
set(hker5                         , ...
  'Linestyle'          , 'none'      , ...
  'Marker'          , 's'      , ...
  'Markersize'       , 5           , ...
  'MarkerEdgeColor'           , [254,217,142]/255 , ...
  'MarkerFaceColor'           , [254,217,142]/255    );

hTitle  = title ('');
hXLabel = xlabel('Real part of eigenvalue','Interpreter','latex' );
hYLabel = ylabel('Imaginary part of eigenvalue of Interpolation matrix','Interpreter','latex' );

hLegend = legend( ...
  [hker1, hker2, hker3, hker4, hker5], ...
  'diffusion - cbl'   , ...
  'diffusion - gbl'   , ...
  'diffusion - sbl'   , ...
  'diffusion - cheb'   , ...
  'diffusion - cheb^2'   , ...
  'location', 'NorthEast');

set( gca                    , ...
  'FontName'   , 'Helvetica', ...
  'FontSize'   ,   14       , ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'on'      , ...
  'XGrid'       , 'on'      , ...
  'xlim'        , [-0.02,0.07]    , ...
  'ylim'        , [-0.0002,0.0002]    , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1         );

set([hTitle, hXLabel, hYLabel], ...
    'FontName'   , 'AvantGarde' ,...
    'FontSize'   ,   14  );
set( hLegend             , ...
    'FontSize'   , 14    );
set( hTitle                    , ...
    'FontWeight' , 'bold'      );

TT = get(gca,'TightInset');
TT(1) = TT(1)-0.001;
set(gca,'Position', [TT(1) TT(2), 1 - TT(1)-TT(3),1-TT(2)-TT(4)]);
hold off;