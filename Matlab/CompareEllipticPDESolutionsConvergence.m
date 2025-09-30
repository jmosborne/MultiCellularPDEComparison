clear, close all

global U D beta CL L0;

% Match these to Chaste sims
U = 1;
D = 0.5;
beta = -0.1;
CL = 1.0;
L0 = 5.0;

plot_time_step = 21;

addpath ../../Chaste/anim/matlab/

% Dirs to load

BaseType = 'Elliptic';

TissueTypes = {'StaticDisc'};


DomainTypes = {'BoxDomain'};

PdeTypes = {'UniformPde','VolumeScaledAveragedPde'};

colors = {'r.';
          'b.'};

legends = {'Box Vol Ave';
           'Box Ad Vol Ave'};

BoxHs = {1.0/2.0,1.0/4.0,1.0/8.0,1.0/16};
colors = {'k.','bx','rd','go'}'


for TissueTypeIndex = 1:length(TissueTypes)
    TissueType = TissueTypes{TissueTypeIndex};

    next_index = 1; %To keep legends aligned
    
    for DomainTypeIndex = 1:length(DomainTypes)
        DomainType = DomainTypes{DomainTypeIndex};

        for PdeTypeIndex = 1:size(PdeTypes,2)
            PdeType = PdeTypes{DomainTypeIndex,PdeTypeIndex};
            
            figure

            for BoxHIndex = 1:length(BoxHs)
                BoxH = BoxHs{BoxHIndex};

                dir = [BaseType, '/', TissueType, '/', DomainType, '/', PdeType, '/Box_', num2str(BoxH)]

                data = LoadNonConstantLengthData(['../../../testoutput/',dir,'/results_from_time_0/celldata_oxygen.dat'])
    
                timestep = min([plot_time_step,length(data)]);
    
                end_time = data{timestep}(1);
                node_id = data{timestep}(2:5:end-4);
                cell_id = data{timestep}(3:5:end-3);
                x = data{timestep}(4:5:end-2);
                y = data{timestep}(5:5:end-1);
                c = data{timestep}(6:5:end);
            
                r = sqrt(x.*x+y.*y);
    
                %plot(r,c,'.','Color',[1.0-(BoxHIndex-1)/(length(BoxHs)-1),0, (BoxHIndex-1)/(length(BoxHs)-1),1.0])
                plot(r,c,colors{BoxHIndex});
                hold on
    
                % main_legend{next_index} = legends{DomainTypeIndex,PdeTypeIndex};
                % next_index = next_index + 1;
            end

            title([TissueType,DomainType,PdeType])
            SaveAsPngEpsAndFig(-1,['Figs/',BaseType,TissueType,PdeType,'_no_exact'], 12, 7/5, 12);
            
            R=sqrt(max(x.*x+y.*y));
            C0=max(c);
            r=linspace(0,R,1000);
            k=0.1;
            
            exact = C0*besseli(0,r*sqrt(k))/besseli(0,R*sqrt(k));
            
            plot(r,exact,'k--');
                    % SaveAsPngEpsAndFig(-1,['Figs/',BaseType,TissueType,PdeType], 12, 7/5, 12);
        end
    end
end




function [c,f,s] = pdepde(x,t,u,dudx) % Equation to solve
global U D beta L0;
L = L0 + U*t;
dLdt = U;

c = L^2/D;
f = dudx;
s = L^2*beta*u/D + x*L*dLdt*dudx/D;
end

function [c,f,s] = advectionpdepde(x,t,u,dudx) % Equation to solve
global U D beta L0;
L = L0 + U*t;
dLdt = U;

c = L^2/D;
f = dudx;
s = L^2*beta*u/D + x*L*dLdt*dudx/D - L*U*x*dudx/D;
end

%----------------------------------------------
function u0 = pdeic(x) % Initial conditions
global CL;

u0 = CL;
end

%----------------------------------------------
function [pl,ql,pr,qr] = pdebc(xl,ul,xr,ur,t) % Boundary conditions
global CL;

pl = 0;
ql = 1;
pr = ur-CL;
qr = 0;
end
%----------------------------------------------
