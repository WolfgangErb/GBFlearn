% GBFlearn: a toolbox for graph signal interpolation
% and classification with graph basis functions (GBFs)
% (C) W. Erb 01.03.2020

% GBF_example_ITP_bunny_approximation_u4 is a script to show how GBF 
% interpolation performs in terms of approximation power when the number
% of interpolation nodes is increased. The approximated 
% signal is the eigenvector u_4 of the graph Laplacian corresponding to the
% fourth smallest eigenvalue. 

clear all
close all

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

%Calculate Spectrum of graph
[G.U,G.Lambda] = GBF_spectrum(G.L,'ascend');

%Choose plotting parameter
plotpar.MM = 2;                 %size of dots
plotpar.ub = 0.02;              %upper boundary
plotpar.lb = 0.02;              %left boundary
plotpar.fontsize = 14;          %fontsize
plotpar.colorbar = 'y';         %set 'y' if you want a colorbar
plotpar.edge = 0;               %set 1 if you want to plot edges
plotpar.edgewidth = 1;          %width of edges

%Extract bandlimited function u4
f = G.U(:,4);

%Kernel parameters for interpolation
type1 = 'diffusion';
alpha1 = -20;

type2 = 'polydecay';
alpha2 = -2;

type3 = 'polydecay';
alpha3 = -4;

type4 = 'varspline';
alpha4 = [-2,0.01];

% Initialize
MM = 2:2:100;
nMM = length(MM);

s1 = zeros(G.N,nMM); s2 = zeros(G.N,nMM); 
s3 = zeros(G.N,nMM); s4 = zeros(G.N,nMM);

error1 = zeros(1,nMM); error2 = zeros(1,nMM);
error3 = zeros(1,nMM); error4 = zeros(1,nMM);

% Start evaluation loop
for i = 1:nMM
    
idxW = nodeselect(1:MM(i));
yW = f(idxW);

bf1 = GBF_genGBF(G.U,G.Lambda,idxW,type1,alpha1);
bf2 = GBF_genGBF(G.U,G.Lambda,idxW,type2,alpha2);
bf3 = GBF_genGBF(G.U,G.Lambda,idxW,type3,alpha3);
bf4 = GBF_genGBF(G.U,G.Lambda,idxW,type4,alpha4);

s1(:,i) = GBF_itpGBF(bf1, idxW,yW);
s2(:,i) = GBF_itpGBF(bf2, idxW,yW);
s3(:,i) = GBF_itpGBF(bf3, idxW,yW);
s4(:,i) = GBF_itpGBF(bf4, idxW,yW);

error1(i) = norm(s1(:,i)-f,inf)/norm(f,inf);
error2(i) = norm(s2(:,i)-f,inf)/norm(f,inf);
error3(i) = norm(s3(:,i)-f,inf)/norm(f,inf);
error4(i) = norm(s4(:,i)-f,inf)/norm(f,inf);

fprintf('Error for %3d interpolation points   : %20.16f \n',MM(i),error1(i));
end

figure('Units', 'pixels', ...
    'Position', [50 50 600 600]);
hold on;

hker1 = plot(MM,log10(error1),'k');
hker2 = plot(MM,log10(error2),'r');
hker3 = plot(MM,log10(error3),'g');
hker4 = plot(MM,log10(error4),'b');

set(hker1                         , ...
  'LineStyle'       , '-'      , ...
  'LineWidth'       , 2           , ...
  'Color'           , [0/255,0/255,0/255]    );
set(hker2                         , ...
  'LineStyle'       , '-'      , ...
  'LineWidth'       , 2           , ...
  'Color'           , [139/255,69/255,19/255]    );
set(hker3                         , ...
  'LineStyle'       , '-'      , ...
  'LineWidth'       , 2           , ...
  'Color'           , [230/255,120/255,15/255]    );
set(hker4                         , ...
  'LineStyle'       , '-'      , ...
  'LineWidth'       , 2           , ...
  'Color'           , [255/255,190/255,0/255]    );


% hTitle  = title ('Approximation error of GBF-interpolants');
hTitle  = title ('');
hXLabel = xlabel('Number of interpolation nodes','Interpreter','latex' );
hYLabel = ylabel('Log approximation error','Interpreter','latex' );

hLegend = legend( ...
  [hker1, hker2, hker3, hker4], ...
  'Diffusion kernel t = 20'   , ...
  'Polynomial decay s = 2'   , ...
  'Polynomial decay s = 4'   , ...
  'Variation spline s = 2'   , ...
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
  'ylim'        , [-4,1]    , ...
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
TT(1) = TT(1)-0.1;
%TT(2) = TT(2)-0.01;
set(gca,'Position', [TT(1) TT(2), 1 - TT(1)-TT(3),1-TT(2)-TT(4)]);
hold off;

figure('Units', 'pixels', ...
'Position', [0 50 1000 500]);

subplot(2,4,[1,2,5,6]),
KK = 20;
plotpar.uaxis = max(s3(:,KK));
plotpar.laxis = min(s3(:,KK));
GBF_drawsignal(G.nodes,G.edges,s3(:,KK),plotpar);
title('GBF interpolant (polynomial, s = 4)');
hold on;
plot( G.nodes(idxW(1:MM(KK)),1),G.nodes(idxW(1:MM(KK)),2),'o','color',[0, 0, 0]/255,'LineWidth',1,'MarkerSize',8)
set(gca,'XTick',[], 'YTick', [])
hold off;

subplot(2,4,[3,4,7,8]), 
plotpar.uaxis = max(abs(f-s3(:,KK)))/max(f);
plotpar.laxis = 0;
GBF_drawsignal(G.nodes,G.edges,abs(f-s3(:,KK))/max(f),plotpar);
title('Interpolation error (polynomial, s = 4)');
hold on;
plot( G.nodes(idxW(1:MM(KK)),1),G.nodes(idxW(1:MM(KK)),2),'o','color',[0, 0, 0]/255,'LineWidth',1,'MarkerSize',8)
set(gca,'XTick',[], 'YTick', [])
hold off;