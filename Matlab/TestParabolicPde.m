clear; close all
global U D beta CL L0;

U=1;
D = 0.5;
beta = -0.1;
CL = 1.0;
L0 = 5.0;

x = linspace(0,1,21);
t = linspace(0,5,11);

m = 1;
sol = pdepe(m,@pdex1pde,@pdex1ic,@pdex1bc,x,t);

u = sol(:,:,1);
surf(x,t,u)
title('Numerical solution computed with 20 mesh points')
xlabel('Distance x')
ylabel('Time t')

figure

for i=1:length(u(:,1))
    L(i) = L0+U*t(i)
    plot(x*L(i),u(i,:))
    %plot(x,u(i,:))
    hold on 
end
title('Solution')
xlabel('Distance x')
ylabel('u(x,2)')

% Plot equilibrium soln (for fixed domain)
% R = max(L);
% r = linspace(0,R,50)
% u_exact = CL*besseli(0,r*sqrt(-beta/D))/besseli(0,R*sqrt(-beta/D));
% 
% plot(r,u_exact,'ro')


function [c,f,s] = pdex1pde(x,t,u,dudx) % Equation to solve
global U D beta L0;
L = L0 + U*t;
dLdt = U;
% 
% c = L^2/D;
% f = dudx;
% s = x*dLdt*L*dudx/D+beta/D*u;
c = L^2/D;
f = dudx;
s = L^2*beta*u/D + x*L*dLdt*dudx/D;

end

%----------------------------------------------
function u0 = pdex1ic(x) % Initial conditions
global CL;

u0 = CL;
end

%----------------------------------------------
function [pl,ql,pr,qr] = pdex1bc(xl,ul,xr,ur,t) % Boundary conditions
global CL;

pl = 0;
ql = 1;
pr = ur-CL;
qr = 0;
end
%----------------------------------------------