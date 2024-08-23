clear, close all

global U D beta CL L0;

% Match these to Chaste sims
U = 1;
D = 0.5;
beta = -0.1;
CL = 1.0;
L0 = 5.0;

plot_time_step = 31;

addpath ../../Chaste/anim/matlab/

% Dirs to load

BaseType = 'Parabolic';

TissueTypes = {'StaticDisc','GrowingDisc','ProliferatingDisc'};

DomainTypes = {'GrowingDomain', 'BoxDomain','BoxDomainAdvection'};

PdeTypes = {'UniformPde', 'CellwisePde', 'VolumeScaledCellwisePde';
            'UniformPde', 'AveragedPde', 'VolumeScaledAveragedPde';
            'UniformPde', 'AveragedPde', 'VolumeScaledAveragedPde'};

colors = {'kd','rs','g^';
          'k*','rx','g+';
          'k.','r.','g.'}

legends = {'Grow Unif','Grow Cell','Grow Vol Cell';
           'Box Unif','Box Ave','Box Vol Ave';
           'Box Ad Unif','Box Ad Ave','Box Ad Vol Ave'};

for TissueTypeIndex = 1:3 %3
    TissueType = TissueTypes{TissueTypeIndex};

    next_index = 1; %To keep legends aligned
    
    figure
    
    for DomainTypeIndex = 1:3 %3
        DomainType = DomainTypes{DomainTypeIndex};

        for PdeTypeIndex = 1:3 %3
            PdeType = PdeTypes{DomainTypeIndex,PdeTypeIndex};

            dir = [BaseType, '/', TissueType, '/', DomainType, '/', PdeType]

            data = LoadNonConstantLengthData(['../../../testoutput/',dir,'/results_from_time_0/celldata_morphogen.dat'])

            timestep = min([plot_time_step,length(data)]);

            end_time = data{timestep}(1);
            node_id = data{timestep}(2:5:end-4);
            cell_id = data{timestep}(3:5:end-3);
            x = data{timestep}(4:5:end-2);
            y = data{timestep}(5:5:end-1);
            c = data{timestep}(6:5:end);
        
            r = sqrt(x.*x+y.*y);

            plot(r,c,colors{DomainTypeIndex,PdeTypeIndex})
            hold on

            main_legend{next_index} = legends{DomainTypeIndex,PdeTypeIndex};
            next_index = next_index + 1;


            % % Plot velocities
            % velocity_data = LoadNonConstantLengthData(['../../../testoutput/',dir,'/results_from_time_0/cellvelocities.dat'])
            % 
            % timestep = length(velocity_data);
            % 
            % end_time = velocity_data{timestep}(1);
            % cell_id = velocity_data{timestep}(2:5:end-4);
            % x = velocity_data{timestep}(3:5:end-3);
            % y = velocity_data{timestep}(4:5:end-2);
            % u = velocity_data{timestep}(5:5:end-1);
            % v = velocity_data{timestep}(6:5:end);
            % 
            % r = sqrt(x.*x+y.*y);
            % ur = sqrt(u.*u+v.*v);
            % 
            % bin_edges = 0:10;
            % bin_centres = (1:10)-0.5;
            % bin_u = cell(10,1);
            % bin_count = zeros(length(bin_centres),1);
            % 
            % for i=1:length(r)
            %     bin = find(r(i)<bin_edges & r(i)+1>bin_edges)-1;
            %     assert(bin<=length(bin_centres));
            %     bin_count(bin) = bin_count(bin) + 1;
            %     bin_u{bin}(bin_count(bin)) = ur(i);
            % end
            % 
            % 
            % for i=1:length(bin_centres)
            %     ave_ur(i) = mean(bin_u{i});
            % end
            % 
            % plot(bin_centres,ave_ur,colors{DomainTypeIndex,PdeTypeIndex})
            % %plot(r,ur,colors{DomainTypeIndex,PdeTypeIndex})
            % hold on
            % plot([0,max(r)],[0,0.5],'k--')
        end
    end
    title(TissueType)
    legend(main_legend,'Location','northwest')
    SaveAsPngEpsAndFig(-1,['Figs/',BaseType,TissueType,'_no_exact'], 12, 7/5, 12);

    U = (max(r)-L0)/end_time;

    x = linspace(0,1,21);
    t = linspace(0,end_time,11);
    
    m = 1;
    sol = pdepe(m,@pdepde,@pdeic,@pdebc,x,t);
    advection_sol = pdepe(m,@advectionpdepde,@pdeic,@pdebc,x,t);
    
    u = sol(:,:,1);
    advection_u = advection_sol(:,:,1);
    
    L = L0+U*end_time
    plot(x*L,u(end,:),'m--')
    plot(x*L,advection_u(end,:),'b--')
   
    main_legend{next_index} = 'Exact';
    main_legend{next_index+1} = 'Advec Exact';
    legend(main_legend,'Location','westoutside')
    
   
    %xlim([0,10])
    %ylim([0.0,1.1])

    % SaveAsPngEpsAndFig(-1,['Figs/',BaseType,TissueType], 7, 7/5, 9);
    SaveAsPngEpsAndFig(-1,['Figs/',BaseType,TissueType], 12, 7/5, 12);
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
