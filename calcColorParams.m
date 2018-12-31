

function colvals=calcColorParams(w,C,r)
colvals=zeros(size(w,1),1);
for i=1:size(w,1)
    for j=1:r
        colvals(i)=colvals(i)+mean(sqrt(sum((bsxfun(@minus,w(C{i,j},:),w(i,:))).^2,2)));
    end
    
end
colvals=colvals./r;
end