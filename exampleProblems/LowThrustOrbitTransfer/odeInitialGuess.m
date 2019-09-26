function dx=odeInitialGuess(t,x,f,data)

theta=0;
phi=30;

u1=-sind(theta);
u2=cosd(theta)*cosd(phi);
u3=-cosd(theta)*sind(phi);
u=[u1 u2 u3];
p=[-8];


% Evaluate ODE right-hand side
dx=f(x',u,p,t,data)';

end