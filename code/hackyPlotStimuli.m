Spos = [0.731,0.3208,1];
Sneg = [0.2690,0.6792,0];
LminusMpos = [0.9780    0.4064    0.4978];
LminusMneg = [0.0220    0.5936    0.5022];

figure
subplot(2,2,1)
rectangle('Position',[2 4 2 2],'Curvature',[1 1],'FaceColor',Spos)
axis square
axis off
title('S+');
subplot(2,2,2)
rectangle('Position',[2 4 2 2],'Curvature',[1 1],'FaceColor',Sneg)
axis square
axis off
title('S-');
subplot(2,2,3)
rectangle('Position',[2 4 2 2],'Curvature',[1 1],'FaceColor',LminusMpos)
axis square
axis off
title('L-M');
subplot(2,2,4)
rectangle('Position',[2 4 2 2],'Curvature',[1 1],'FaceColor',LminusMneg)
axis square
axis off
title('M-L');
