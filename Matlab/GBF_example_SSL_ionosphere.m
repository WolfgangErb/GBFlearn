% GBFlearn: a toolbox for graph signal interpolation
% and classification with graph basis functions (GBFs)
% (C) W. Erb 01.03.2020

% GBF_example_SSL_ionosphere shows how SSL based on a
% feature-augmented GBF-RLS solution can be used to improve
% classification accuracy of a supervised scheme.
% The included feature is a unsupervised classification of the graph 
% based on spectral clustering. The test data set consists of the 
% ionosphere data set from the UCI machine learning repository.  

clear all
close all

% Paths
addpath(genpath('./core/'))
addpath(genpath('./data/'))

%Choose graph
G.type = 'ionosphere';
load data_ionosphere.mat;

%Generate graph
[G.nodes,G.edges,G.A] = GBF_gengraph(G.type);

%Calculate the graph Laplacian
G.N = length(G.nodes(:,1));
G.deg = sum(G.A,1);
isD = diag(1./sqrt(G.deg));
G.L = eye(G.N) - isD*G.A*isD;

%Calculate Spectrum of graph
[G.U,G.Lambda] = GBF_spectrum(G.L,'ascend');

%Calculate (unsupervised) classification based on spectral clustering
%Corresponds to the normalized cut of Shi and Malik
%Instead of median you can also take the mean or simply cut = 0. 
sCut = G.U(:,2);
cut = median(sCut);

idxcutup = find(sCut>=cut);
idxcutdown = find(sCut<cut);

sCutclass = zeros(size(sCut));
sCutclass(idxcutup)=1;
sCutclass(idxcutdown) = -1;

%Kernel parameters
type = 'diffusion';             % Type of GBF
alpha = -5;                     % Shape parameter
lambda = 0.001;                 % Regularization parameter
gamma1 = 0.5;                   % Binary kernel parameter 1. test
gamma2 = -0.5;                  % Binary kernel parameter 2. test

%Start randomized supervised training
M = [10,20,30,40,50,60];
R = 100;

%Initialize
sRLSclass = zeros(G.N,R); 
sTP1class = zeros(G.N,R);
sTP2class = zeros(G.N,R);

meanclasserror = zeros(4,length(M));

%Calculate the mean classification errors of R random experiments
for j = 1:length(M)
  for r = 1:R
    IntIndex = randsample(G.N,M(j));
    y = label(IntIndex);

    %Calculate standard RLS solution
    bf = GBF_genGBF(G.U,G.Lambda,IntIndex,type,alpha);
    sRLS = GBF_RLSGBF(bf,IntIndex,y,lambda);

    %Calculate feature adapted RLS solution
    bin1 = GBF_genbinK(IntIndex,sCutclass,gamma1);
    bin2 = GBF_genbinK(IntIndex,sCutclass,gamma2);
    sTP1 = GBF_RLSGBF(bf.*bin1, IntIndex, y, lambda);
    sTP2 = GBF_RLSGBF(bf.*bin2, IntIndex, y, lambda); 
  
    %Calculate classifications
    idxsRLSup = find(sRLS>=0); idxsRLSdown = find(sRLS<0);
    idxsTP1up = find(sTP1>=0); idxsTP1down = find(sTP1<0);
    idxsTP2up = find(sTP2>=0); idxsTP2down = find(sTP2<0);

    sRLSclass(idxsRLSup,r) = 1; sRLSclass(idxsRLSdown,r) = -1;
    sTP1class(idxsTP1up,r) = 1; sTP1class(idxsTP1down,r) = -1;
    sTP2class(idxsTP2up,r) = 1; sTP2class(idxsTP2down,r) = -1;
  end

  classerror = zeros(3,R);
  classerror(1,:) = 1-sum(abs(sRLSclass - repmat(label,1,R)),1)/(G.N)/2;
  classerror(2,:) = 1-sum(abs(sTP1class - repmat(label,1,R)),1)/(G.N)/2;
  classerror(3,:) = 1-sum(abs(sTP2class - repmat(label,1,R)),1)/(G.N)/2;

  meanclasserror(1,j) = 1-sum(abs(sCutclass - label))/(G.N)/2;
  meanclasserror(2,j) = mean(classerror(1,:));
  meanclasserror(3,j) = mean(classerror(2,:));
  meanclasserror(4,j) = mean(classerror(3,:));
end

fprintf('Mean classification error for SSL on the ionosphere data set\n');
fprintf('------------------------------------------------------------\n');
fprintf(1, 'Number of labeled nodes:       '); fprintf(1,'%5d \t', [10,20,30,40,50,60]);   fprintf(1, '\n');
fprintf(1, 'Spectral clustering:           '); fprintf(1,'%1.4f \t', meanclasserror(1,:)); fprintf(1, '\n');
fprintf(1, 'GBF-RLS classification:        '); fprintf(1,'%1.4f \t', meanclasserror(2,:)); fprintf(1, '\n');
fprintf(1, 'Augm. GBF-RLS class., Test 1:  '); fprintf(1,'%1.4f \t', meanclasserror(3,:)); fprintf(1, '\n');
fprintf(1, 'Augm. GBF-RLS class., Test 2:  '); fprintf(1,'%1.4f \t', meanclasserror(4,:)); fprintf(1, '\n');
