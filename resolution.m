
close
clear all
global h a e m K kappa r deltaN deltaP deltaV L ndt ndx deltax superp_C
enregistrer_image=0; % vaut 1 si on souhaite neregistrer.
change_color=0; %vaut 1 pour des courbes dont la couleur change
superp_C=0; %vaut 1 si on souhaite superposer les courbes
larg_marche=4; % si condition initiale en escalier
init;
ndt = 1000; % Nombre de pas de temps
ncourbe=100;%nombre de courbe à afficher


deltat=1;%le pas de temps vaut 1
ndx=1000; % Nombre d'éléments

deltax=L/ndx;%maille régulière
ntemps=floor(ndt/deltat);%nombre de calculs dans la boucle en temps
compteC=floor(ntemps/ncourbe);%tout les compteC courbe on affiche
k=1;%compteur de courbe
color=0; %pour changement de couleur des courbes
grey=[0.8 0.8 0.8];
fig=figure('units','normalized','outerposition',[0 0 1 1]);
set(fig, 'PaperPositionMode','auto');
set(gcf,'Color',grey)
set(gcf,'inverthardcopy','off')

%#################### Conditions Initiales #############################
 % N(i)= 15*exp(-10*(i.*i)); %densité initiale de ravageurs biblio
 % P(i)=6*exp(-100*(i.*i)); %densité initiale de prédateurs biblio
 N=zeros(1, ndx);
 Y=marche(larg_marche,round(0.3*ndx));
 %haut de marche
 i=1:(round(0.3*ndx)-3);
 N(i)=1;
 %descente de marche
 i=(round(0.3*ndx)-2):(round(0.3*ndx)+2);
 N(i)=Y(1)*i.*i.*i+Y(2)*i.*i+Y(3)*i+Y(4);
 %après la marche
 i=(round(0.3*ndx)+3):ndx;
 N(i)=0;
 P=zeros(1, ndx);
 Y=marche(-larg_marche,round(0.7*ndx));
 %bas de marche
 i=1:(round(0.7*ndx)-3);
 P(i)=0;
 %montée de marche
 i=(round(0.7*ndx)-2):(round(0.7*ndx)+2);
 P(i)=Y(1)*i.*i.*i+Y(2)*i.*i+Y(3)*i+Y(4);
 %haut de marche
 i=(round(0.7*ndx)+3):ndx;
 P(i)=1;

22 
 V=zeros(1, ndx);%vitesse nulle initialement
 N0=N;
 P0=P;

 FoN=(deltaN*deltat)/(deltax^2);%analogue du nbr de Fourrier proies.
 FoP=(deltaP*deltat)/(deltax^2);%analogue du nbr de Fourrier prédateurs.
 FoV=(deltaV*deltat)/(deltax^2);%analogue du nbr de Fourrier vitesse.



 % Boucle en temps
 %===========================================================
 for t=1:deltat:ndt
%########### Calcul champ de VITESSE ################
 %=======================================================

 gradN= gradiant(N);%gradient de proies

 % Condition aux limites à l'entrée
 V(1)=kappa*gradN(1);%pas d'acceleration
 % Condition en sortie
 V(ndx)=kappa*gradN(ndx);
 % Intégration de l'EDP

 %########### V #################
 % Condition aux limites à l'entrée
 aCV(1) = 1 ;
 aEV(1) = 0;
 bV(1) = V(1);
 % Intégration de l'EDP
 for i=2:ndx-1 ;
 aWV(i) = -FoV ;
 aCV(i) = 1+2*FoV ;
 aEV(i) = -FoV ;
 bV(i) = V(i)+deltat*(kappa*gradN(i)) ;
 end
 % Condition en sortie
 aWV(ndx) = 0 ;
 aCV(ndx) = 1 ;
 bV(ndx) = V(ndx);
% Résolution V et déduction NcP

 V=tridiag(aWV,aCV,aEV,bV);
i=1:ndx;
NcP(i)=V(i)*deltat/deltax;




23 
%########### Calcul champ des différentes populations P etN #############

 % Condition aux limites à l'entrée
 aCN(1) = 1+FoN;
 aEN(1) = -FoN ;
 bN(1) = N(1)+...
 (deltat)*(r*N(1)*(1-N(1)/K)-(a*N(1)*P(1))/(1+a*h*N(1)));

 aCP(1) =1+FoP;
 aEP(1) =-FoP;
 bP(1) = NcP(2)*P(2)+(1-NcP(1))*P(1)... convection
 +(deltat)*(e*(a*N(1)*P(1))/(1+a*h*N(1))-m*P(1)) ;

 % Intégration de l'EDP
 i=2:ndx-1 ;
 aWN(i) = -FoN*ones(1,ndx-2) ;
 aCN(i) = (1+2*FoN)*ones(1,ndx-2) ;
 aEN(i) = -FoN*ones(1,ndx-2) ;
 bN(i) =N(i) ...
 +(deltat)*(r*N(i).*(1-N(i)/K)-(a*N(i).*P(i))./(1+a*h*N(i)));

 aWP(i) = -FoP*ones(1,ndx-2) ;
 aCP(i) = (1+2*FoP)*ones(1,ndx-2) ;
 aEP(i) = -FoP*ones(1,ndx-2) ;
 bP(i) = NcP(i-1).*P(i-1)+(1-NcP(i)).*P(i)...
 +(deltat)*(e*(a*N(i).*P(i))./(1+a*h*N(i))-m*P(i));

 % Condition en sortie
 aWN(ndx) = -FoN ;
 aCN(ndx) = 1+FoN ;
 bN(ndx) = N(ndx)+...
 (deltat)*(r*N(ndx)*(1-N(ndx)/K)-(a*N(ndx)*P(ndx))/(1+a*h*N(ndx))) ;

 aWP(ndx) = -FoP ;
 aCP(ndx) = 1+FoP ;
 bP(ndx) = NcP(ndx-1)*P(ndx-1)+(1-NcP(ndx))*P(ndx)...
 +(deltat)*(e*(a*N(ndx)*P(ndx))/(1+a*h*N(ndx))-m*P(ndx));
 % Résolution
 N=tridiag(aWN,aCN,aEN,bN);
 P=tridiag(aWP,aCP,aEP,bP);

%################ Affichage des courbes ########################

 if k==compteC || t==ndt
 affichage;
 disp('ITERATION NUMERO : ')
 t
 k=1;
 if color<=0.95 && change_color==1%changement couleur des courbes
 color=color+0.050;
 end
 else
 k=k+1;
 end

 end
