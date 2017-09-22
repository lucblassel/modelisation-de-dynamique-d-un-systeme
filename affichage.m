global superp_C
green = [0 1-color 0];
red=[1-color 0 0];
greeninit=[0.2 0.6 0.2];
redinit=[0.8 0.2 0];
orange=[1 0.6 0];
%######### PREDATEURS ################
s(1)=subplot(2,2,4); %position de la fenêtre
plot(P0,'--','Color',redinit); %affichage courbe initiale
hold on
plot(P,'Color',red);
hold on
plot(V,'Color', orange);
legend('predateur t0','predateurs','vitesse predateurs','Location',...
 'southoutside','Orientation','horizontal')
title(s(1),'PREDATEURS','Color',red)
xlabel(s(1),'position longitudinale')
ylabel(s(1),'densite de population')
axis(s(1),[0 ndx -5 15]);
set(gca, 'Ticklength', [0 0])
if superp_C~=1 % voir fichier resolution
 hold off
end
%######### PROIES ################
s(2)=subplot(2,2,3);
plot(N0,'--','Color', greeninit);
hold on
plot(N,'Color',green);
hold on
plot(gradN,'m');
legend('proies t0','proies','gradient de proies','Location',...
 'southoutside','Orientation','horizontal')
title(s(2),'PROIES','Color',green)
xlabel(s(2),'position longitudinale')
ylabel(s(2),'densite de population')
set(gca, 'Ticklength', [0 0])
if superp_C~=1
 hold off
end
%######### PREDATEURS vs. PROIES ################
s(3)=subplot(2,2,1:2);
%suptitle(['ITERATION NUMERO : ' num2str(t) ]);%affichage n° d'itération
plot(N0,'--','Color',greeninit);
hold on
plot(P0,'--','Color',redinit);
hold on
plot(N,'Color',green);
hold on
plot(P,'Color',red);
axis(s(3),[0 ndx -1 15]);
legend('proies t0','predateur
t0','proies','predateurs','Location','eastoutside','Orientation','vertical')
title(s(3),'PROIES vs.PREDATEURS','Color','k')
xlabel(s(3),'position longitudinale')
ylabel(s(3),'Densite de population')
set(gca, 'Ticklength', [0 0])
20 
% ######### enregistrement des figures #########
if enregistrer_image==1
 T=num2str(t);
 base='figure';
 name=strcat(base,T,'.png');
 namex=strcat(base,T,'.emf');

 print -dmeta -r0 name
 saveas(gcf,name);
 saveas(gcf,namex);
end
if superp_C~=1
 hold off
end
pause(0.01)
