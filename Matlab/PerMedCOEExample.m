%%
clear, close all 

model = createpde();

width = 220; %240-20 as were tracking centre of cells.
t_end = 10;
diffusion = 2000;
cell_level_uptake = 20;
cell_width = 20;

gm = multicuboid(width,width,width);
model.Geometry = gm;

pdegplot(gm,"FaceLabels","on","FaceAlpha",0.5)


figure

applyBoundaryCondition(model,"dirichlet", ...
                             "Face",1:model.Geometry.NumFaces, ...
                             "u",10);


specifyCoefficients(model,"m",0,...
                          "d",1,...
                          "c",diffusion,...
                          "a",0,...
                          "f",-cell_level_uptake/cell_width/cell_width/cell_width);


setInitialConditions(model,0);
setInitialConditions(model,10,"Face",1);
setInitialConditions(model,10,"Face",2);
setInitialConditions(model,10,"Face",3);
setInitialConditions(model,10,"Face",4);
setInitialConditions(model,10,"Face",5);
setInitialConditions(model,10,"Face",6);


generateMesh(model,"Hmax",20);

tlist = 0:0.1:t_end;

results = solvepde(model,tlist);

figure
pdeplot3D(model,"ColorMapData",results.NodalSolution(:,2))


%% Plot Cross section 
xx = -width/2:width/20:width/2;
yy = xx;
[XX,YY] = meshgrid(xx,yy);
ZZ = width/2*ones(size(YY));
uintrp = interpolateSolution(results,XX,YY,ZZ,1:length(tlist));

figure
subplot(1,2,1)
usol = uintrp(:,1);
usol = reshape(usol,size(XX));
surf(XX,YY,usol) 
title("t = 0")
subplot(1,2,2)
usol = uintrp(:,end);
usol = reshape(usol,size(XX));
surf(XX,YY,usol) 
title("t = t_end")


%% Interpolate onto cells and create timeseries

xx = -width/2+20:20:width/2-20;
yy = xx;
zz = 20:20:width-20;
[XX,YY,ZZ] = meshgrid(xx,yy,zz);

uintrp = interpolateSolution(results,XX,YY,ZZ,1:length(tlist));

figure

plot(tlist,sum(uintrp)/(length(xx)^3)*cell_width^3);
axis([0,10,0,10]);


