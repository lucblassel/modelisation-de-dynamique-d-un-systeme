%fonction escalier continue sur une bande de largeur a et centrée en x
function Y=marche(a,x)
%y=A*x^3+B*x²+C*x+D
%Résolution d'un système d'équation de 4 équations à 4 inconnues
%f(x-a)=1
%f(x+a)=0
%f'(x-a)=0
%f'(x+a)=0
a=a/2;
A=[ (x-a)^3 (x-a)^2 (x-a) 1;...
 (x+a)^3 (x+a)^2 (x+a) 1;...
 3*(x-a)^2 2*(x-a) 1 0;...
 3*(x+a)^2 2*(x+a) 1 0];
B=[1;0;0;0];
Y=A\B;% coefficient A,B,C et D
