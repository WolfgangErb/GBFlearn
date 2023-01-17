% GBFlearn: a toolbox for graph signal interpolation
% and classification with graph basis functions (GBFs)
% (C) W. Erb 01.03.2020

function [nodes,edges,A] = GBF_gengraph(type)

% Provides the data for different graph examples

% INPUT:    
% type         : '2moon'           : The two-moon dataset
%              : 'bunny'           : The stanford bunny
%              : 'emptyset'        : A dataset in form of a slashed O
%              : 'ionosphere'      : From UCI repository
%              : 'minnesota'       : The Minnesota graph
%              : 'path'            : A simple path graph
%              : 'rand'            : A random sensor network
%              : 'sensor1'         : A simple sensor network
%              : 'sensor2'         : A larger sensor network
%              : 'sensorverylarge' : A very large sensor network
%              : 'star'            : A simple star graph
%              : 'wbc'             : Wisconsin breast cancer (UCI repository)
%
% OUTPUT:    
% nodes        : Nodes of the graph G
% edges        : Edges of the graph G
% A            : Adjacency matrix of the graph G

switch type
      
case '2moon'
    
  load data_2moon.mat
  r = 0.5;
  [edges,A] = GBF_sim_rballs(nodes,r);

case 'bunny'
    
  load data_bunny.mat
  nodes = [bunny(:,1),bunny(:,2)];
  thresh = 0.0025;
  N = size(nodes,1);

  stp = 0;
  for i = 1 : N
    for j = i+1 : N
      if norm( nodes(i,:) - nodes(j,:), 2 ) <= thresh
        stp = stp + 1;
        idx(stp) = j;
      end
    end
  end

  nodes(idx,:)=[];
  r = 0.01;
  [edges,A] = GBF_sim_rballs(nodes,r);

case 'emptyset'
    
  load data_emptyset.mat
  r = 0.2;
  [edges,A] = GBF_sim_rballs(nodes,r);

case 'ionosphere'
    
  load data_ionosphere
  alpha = 1;
  [edges,A] = GBF_sim_gauss(nodes,alpha);
  
case 'minnesota'
    
  load data_minnesota.mat
  nodes = xy;
  % A = full(A);
  A(349,355)=1;
  A(355,349)=1;
  [row,col] = find(A~=0);
  idxx = find(row < col);
  edges = [row(idxx),col(idxx)];
  
case 'path'
    
  N = 201;
  nodes = zeros(N,2);
  edges = zeros(N-1,2);
  nodes(:,1) = linspace(0,1,N)';
  edges(:,1) = (1:N-1)';
  edges(:,2) = (2:N)';
  A = GBF_adjmat(edges,N);  
  
case 'rand'
    
  N = 300;
  thresh = 0.02;
  shift = [0,0];
  D = [1,1];
  nodes = ones(N,1)*D.*rand(N,2)+ones(N,1)*shift;
  stp = 0;
  for i = 1 : N
    for j = i+1 : N
      if norm( nodes(i,:) - nodes(j,:), 2 ) <= thresh
        stp = stp + 1;
        idx(stp) = j;
      end
    end
  end

  nodes(idx,:)=[];
  r = 1/6;
  [edges,A] = GBF_sim_rballs(nodes,r);
  
case 'sensor1'
    
  load data_sensor1.mat
  r = 1/6;
  [edges,A] = GBF_sim_rballs(nodes,r);

case 'sensor2'
    
  load data_sensor2.mat
  r = 1/6;
  [edges,A] = GBF_sim_rballs(nodes,r);
  
case 'sensorverylarge'
    
  load data_sensorverylarge.mat
  
case 'star'
    
  N = 40;
  nodes = zeros(N,2);
  edges = ones(N-1,2);
  nodes(2:N,:) = [cos((0:N-2)'*2*pi/(N-1)),sin((0:N-2)'*2*pi/(N-1))];
  edges(:,2) = (2:N)';
  A = GBF_adjmat(edges,N);
  
case 'wbc'
    
  load data_WBC
  alpha = 0.05;
  [edges,A] = GBF_sim_gauss(nodes,alpha);

end

return