
function plotSelectedPoints(X,alreadySelected,totClusters)
% colors={'b','g','k','w','c','r','y'};
colors={'w','k','r','w','c','r','y'};

hold on;
%alreadySelected
%totClusters

for i=1:totClusters
    
scatter3(X(alreadySelected==i,1),X(alreadySelected==i,2),X(alreadySelected==i,3),'LineWidth',3,'Marker','p','MarkerEdgeColor',colors{i});
end
hold off;
