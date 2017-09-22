function gradN=gradiant(N)
global deltax
n=length(N);
gradN=zeros(1,n);
gradN(1)=(N(2)-N(1))/deltax; %condition de Neumann sur N
gradN(n)=(N(n)-N(n-1))/deltax;
i=2:n-1;
gradN(i)=(N(i+1)-N(i-1))/(2*deltax);
gradN;
19 
