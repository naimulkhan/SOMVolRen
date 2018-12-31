
% load c3-12.mat;
% glyphHard(X);
% % %axis on;
% XX=X(:,1);
% YY=X(:,2);
% ZZ=X(:,3);
% [SS firstpoint secondpoint]=rbb3select(XX,YY,ZZ);
% hold on
% Points=[XX(SS) YY(SS) ZZ(SS)];
% size(SS)
% newPoints=[];
% for i=1:size(SS,1)
%     if SS(i,1) && min([norm([XX(i) YY(i) ZZ(i)]-firstpoint) norm([XX(i) YY(i) ZZ(i)]-secondpoint)])>1
%         SS(i,1)=false;
%     end
% end
% % size(Points)
% % size(newPoints)
% scatter3(XX(SS),YY(SS),ZZ(SS),'o','r')
% hold off


curdir='../MyVTKData/';
curdata='Foot';

load c3-12.mat;
load(strcat(curdir,curdata,'features.mat'));
load(strcat(curdir,curdata,'SOFM.mat'));

glyph(X,ones(size(X,1),1),calcColorParams(w,C,1));


