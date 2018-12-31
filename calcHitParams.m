function colvals=calcHitParams(w,assoc)
colvals=zeros(size(w,1),1);
for i=1:size(w,1)
    colvals(i)=size(assoc(assoc==i),1);
end