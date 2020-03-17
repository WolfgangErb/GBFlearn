% GBFlearn: a toolbox for graph signal interpolation
% and classification with graph basis functions (GBFs)
% (C) W. Erb 01.03.2020

function GBF_drawsignal(nodes,edges,signal,plotpar)

% This function draws a signal on the nodes of the graph G

% INPUT:    
% nodes        : Nodes of the graph G
% edges        : Edges of the graph G
% signal       : The signal on the Nodes of G
% plotpar      : The following parameters are relevant:
%                MM          : size of nodes
%                lb          : additional space for left and right boundary
%                ub          : additional space for upper and lower boundary
%                uaxis       : upper bound for representation of values
%                laxis       : lower bound for representation of values
%                fontsize    : fontsize
%                colorbar    : 'y' or 'n'
%                edge        : 0 (no) or 1 (yes)

 if ~exist('plotpar','var')
      plotpar.flag = 'plotpar created automatically';
 end
 
 if ~isfield(plotpar,'MM')
      plotpar.MM = 4;
 end
 
 if ~isfield(plotpar,'lb')
     plotpar.lb = 0.01;
 end
 
 if ~isfield(plotpar,'ub')
     plotpar.ub = 0.01;
 end
 
 if ~isfield(plotpar,'uaxis')
     plotpar.uaxis = 1;
 end
 
 if ~isfield(plotpar,'laxis')
     plotpar.laxis = 0;
 end
  
 if ~isfield(plotpar,'fontsize')
     plotpar.fontsize = 18;
 end
  
 if ~isfield(plotpar,'colorbar')
    plotpar.colorbar = 'y';
 end
    
 if ~isfield(plotpar,'edge')
      plotpar.edge = 1;
 end
 
 if ~isfield(plotpar,'edgewidth')
      plotpar.edgewidth = 1;
 end
 
% plot the graph with nodes and edges

% draw the edges
if (plotpar.edge == 1)
    for k = 1 : size(edges,1)
      i = edges(k,1);
      j = edges(k,2);    
      plot( [nodes(i,1) nodes(j,1)], [nodes(i,2) nodes(j,2)], 'color',[0.6,0.6,0.6], 'LineWidth',plotpar.edgewidth);
      hold on
    end
end

% draw the signal on the nodes
scatter(nodes(:,1),nodes(:,2), plotpar.MM*100, signal, '.');
hold off;
caxis([plotpar.laxis plotpar.uaxis] ) ;
axis square;
axis([min(nodes(:,1))-plotpar.lb max(nodes(:,1))+plotpar.lb min(nodes(:,2))-plotpar.ub max(nodes(:,2))+plotpar.ub]);
h = get(gcf,'CurrentAxes');
set(h, 'FontName', 'cmr10', 'FontSize', plotpar.fontsize)
colormap(flipud(copper));
if plotpar.colorbar == 'y'
   colorbar;
end

return