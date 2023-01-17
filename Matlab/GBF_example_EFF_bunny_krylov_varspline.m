% GBFlearn: a toolbox for graph signal interpolation
% and classification with graph basis functions (GBFs)
% (C) W. Erb 15.01.2023

% Name: GBF_example_EFF_bunny_krylov_varspline.m
% Test script to test convergence of block Krylov space methods for 
% increasing iteration number m. For five Krylov methods and increasing m, approximate
% kernels and kernel interpolants are calculated and compared. 

% Test scenario:
% graph: bunny
% kernel: variational spline with s =  2, epsilon = 0.05
% number of sampling nodes: N = 20
% problem: calculate kernel interpolant for data given on sampling nodes. 

clear all, close all

% Paths
addpath(genpath('./core/'))
addpath(genpath('./data/'))

%Choose graph
G.type = 'bunny';
load data_bunny.mat        % loads nodeselect

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

%Choose interpolation nodes and values on the graph
NN = 20;
idxW = (1:NN)';
yW = double(G.nodes(idxW,1)<0); %binary labels

%Kernel parameters for interpolation (variational spline)
alpha = [2,0.05];
tol = 1e-20;

% Calculate exact kernel and interpolant
bff = GBF_genGBFeff(G.L,idxW,'varspline',alpha,'direct');
f = GBF_itpGBF(bff, idxW,yW);

type1 = 'varspline';
method1 = 'cbl';
type2 = 'varspline';
method2 = 'gbl';
type3 = 'varspline';
method3 = 'sbl';
type4 = 'varspline';
method4 = 'cheb';
type5 = 'varspline';
method5 = 'cheb2';


% Initialize tests
MM = 2:2:50;
nMM = length(MM);

s1 = zeros(G.N,nMM); s2 = zeros(G.N,nMM); 
s3 = zeros(G.N,nMM); s4 = zeros(G.N,nMM); s5 = zeros(G.N,nMM);

error1 = zeros(1,nMM); error2 = zeros(1,nMM);
error3 = zeros(1,nMM); error4 = zeros(1,nMM); error5 = zeros(1,nMM);

% Start evaluation loop
for i = 1:nMM
    
bf1 = GBF_genGBFeff(G.L,idxW,type1,[alpha,MM(i),tol],method1);
bf2 = GBF_genGBFeff(G.L,idxW,type2,[alpha,MM(i),tol],method2);
bf3 = GBF_genGBFeff(G.L,idxW,type3,[alpha,MM(i),tol],method3);
bf4 = GBF_genGBFeff(G.L,idxW,type4,[alpha,MM(i),tol],method4);
bf5 = GBF_genGBFeff(G.L,idxW,type5,[alpha,MM(i),tol],method5);

s1(:,i) = GBF_itpGBF(bf1, idxW,yW);
s2(:,i) = GBF_itpGBF(bf2, idxW,yW);
s3(:,i) = GBF_itpGBF(bf3, idxW,yW);
s4(:,i) = GBF_itpGBF(bf4, idxW,yW);
s5(:,i) = GBF_itpGBF(bf5, idxW,yW);

error1(i) = norm(s1(:,i)-f,inf)/norm(f,inf);
error2(i) = norm(s2(:,i)-f,inf)/norm(f,inf);
error3(i) = norm(s3(:,i)-f,inf)/norm(f,inf);
error4(i) = norm(s4(:,i)-f,inf)/norm(f,inf);
error5(i) = norm(s5(:,i)-f,inf)/norm(f,inf);

fprintf('Error for %3d interpolation points   : %20.16f \n',MM(i),error1(i));
end

%%
figure('Units', 'pixels', ...
    'Position', [50 50 600 600]);
hold on;

hker1 = plot(MM,log10(error1),'k');
hker2 = plot(MM,log10(error2),'r');
hker3 = plot(MM,log10(error3),'g');
hker4 = plot(MM,log10(error4),'b');
hker5 = plot(MM,log10(error5),'m');

set(hker1                         , ...
  'LineStyle'       , '-'      , ...
  'LineWidth'       , 2           , ...
  'Color'           , [0,0,0]/255    );
set(hker2                         , ...
  'LineStyle'       , '-'      , ...
  'LineWidth'       , 2           , ...
  'Color'           , [153,52,4]/255    );
set(hker3                         , ...
  'LineStyle'       , '-'      , ...
  'LineWidth'       , 2           , ...
  'Color'           , [217,95,14]/255    );
set(hker4                         , ...
  'LineStyle'       , '-'      , ...
  'LineWidth'       , 2           , ...
  'Color'           , [254,153,41]/255    );
set(hker5                         , ...
  'LineStyle'       , '-'      , ...
  'LineWidth'       , 2           , ...
  'Color'           , [254,217,142]/255    );

hTitle  = title ('');
hXLabel = xlabel('Number $m$ of iterations','Interpreter','latex' );
hYLabel = ylabel('Log approximation error','Interpreter','latex' );

hLegend = legend( ...
  [hker1, hker2, hker3, hker4, hker5], ...
  'var. spline - cbl'   , ...
  'var. spline - gbl'   , ...
  'var. spline - sbl'   , ...
  'var. spline - cheb'   , ...
  'var. spline - cheb^2'   , ...
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
  'ylim'        , [-12,2]    , ...
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

%%
figure('Units', 'pixels', ...
'Position', [0 50 800 600]);

K = 25;
plotpar.uaxis = max(s1(:,K));
plotpar.laxis = min(s1(:,K));
GBF_drawsignal(G.nodes,G.edges,s1(:,K),plotpar);
title('Kernel predictor (variational spline, s = 2)');
hold on;
plot( G.nodes(idxW,1),G.nodes(idxW,2),'o','color',[0, 0, 0]/255,'LineWidth',1,'MarkerSize',8)
set(gca,'XTick',[], 'YTick', [])
hold off;
