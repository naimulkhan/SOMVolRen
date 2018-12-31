load c3-12.mat;
 

glyphHard(X);
     %zoom(10);
     disp('Click anywhere on the surface, then hit return')
     pause
     [p v vi face facei] = select3d;
     
      disp(sprintf('\nYou clicked at\nX: %.2f\nY: %.2f\nZ: %.2f',p(1),p(2),p(3)'))
     disp(sprintf('\nThe nearest vertex is\nX: %.2f\nY: %.2f\nZ: %.2f',v(1),v(2),v(3)'))
     hold on
     scatter3(v(1),v(2),v(3),'o','r')