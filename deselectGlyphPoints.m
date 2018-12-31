
function SS=deselectGlyphPoints(X,alreadySelected,totClusters)


% %axis on;
XX=X(:,1);
YY=X(:,2);
ZZ=X(:,3);
[SS firstpoint secondpoint]=rbb3select(XX,YY,ZZ);
SS=uint16(SS);
for i=1:size(SS,1)
    if SS(i,1)
        SS(i,1)=totClusters;
        if ~alreadySelected(i,1)
            SS(i,1)=0;
        elseif min([norm([XX(i) YY(i) ZZ(i)]-firstpoint) norm([XX(i) YY(i) ZZ(i)]-secondpoint)])>1
            SS(i,1)=0;
        end
    end
    
end



