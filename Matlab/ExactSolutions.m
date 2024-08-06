% Compare numerical Pde solution to radial solution
clear, close all

global R;

%Poisson's Equation on Unit Disk
C0=1;
R=5; %% Can only do this currently
k=0.1;


%Create the PDE model and include the geometry.
model = createpde();
g = geometryFromEdges(model,@radius_circleg);

applyBoundaryCondition(model,"dirichlet", ...
                             "Edge",1:model.Geometry.NumEdges, ...
                             "u",C0);
specifyCoefficients(model,"m",0,"d",0,"c",1,"a",0,"f",-k);

hmax = 0.1;
generateMesh(model,"Hmax",hmax);


% figure
% pdemesh(model); 
% axis equal

results = solvepde(model);
u = results.NodalSolution;
figure
pdeplot(model,"XYData",u)
title("Uniform Uptake");
xlabel("x")
ylabel("y")

%Compare this result with the exact analytical solution and plot the error.
mesh = model.Mesh;
p = mesh.Nodes;
u = results.NodalSolution;
r = sqrt(p(1,:).^2 + p(2,:).^2);
exact = C0-k/4*(R^2 - r.^2);
figure
pdeplot(model,"XYData",u-exact')
title("Uniform Uptake Error");
xlabel("x")
ylabel("y")
 

%Constant uptake 

%Create the PDE model and include the geometry.
model = createpde();
geometryFromEdges(model,@radius_circleg);

applyBoundaryCondition(model,"dirichlet", ...
                             "Edge",1:model.Geometry.NumEdges, ...
                             "u",C0);
specifyCoefficients(model,"m",0,"d",0,"c",1,"a",k,"f",0);

hmax = 0.1;
generateMesh(model,"Hmax",hmax);


% figure
% pdemesh(model); 
% axis equal

results = solvepde(model);
u = results.NodalSolution;
figure
pdeplot(model,"XYData",u)
title("Constant Uptake");
xlabel("x")
ylabel("y")

%Compare this result with the exact analytical solution and plot the error.
mesh = model.Mesh;
p = mesh.Nodes;
u = results.NodalSolution;
r=sqrt(p(1,:).^2 + p(2,:).^2);
exact = C0*besseli(0,r*sqrt(k))/besseli(0,R*sqrt(k));

figure
pdeplot(model,"XYData",u - exact')
title("Constant Uptake Error");
xlabel("x")
ylabel("y")
 

