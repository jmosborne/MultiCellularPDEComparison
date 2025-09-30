clear, close all

global U D beta CL L0;

% Match these to Chaste sims
U = 1;
D = 0.1;
beta = -0.1;
CL = 1.0;
L0 = 5.0;


% Change both of these
plot_time_step = 11;
max_radius = 10;

addpath ../../Chaste/anim/matlab/

% Dirs to load

BaseType = 'Parabolic';

%TissueTypes = {'StaticDisc','GrowingDisc','ProliferatingDisc'};
TissueTypes = {'GrowingDisc','ProliferatingDisc'};

%DomainTypes = {'GrowingDomain', 'BoxDomain','BoxDomainAdvection'};
DomainTypes = {'BoxDomain','BoxDomainAdvection'};

% PdeTypes = {'UniformPde', 'CellwisePde', 'VolumeScaledCellwisePde';
%             'UniformPde', 'AveragedPde', 'VolumeScaledAveragedPde';
%             'UniformPde', 'AveragedPde', 'VolumeScaledAveragedPde'};
% 
% colors = {'kd','rs','g^';
%           'k*','rx','g+';
%           'k.','r.','g.'}
% 
% legends = {'Grow Unif','Grow Cell','Grow Vol Cell';
%            'Box Unif','Box Ave','Box Vol Ave';
%            'Box Ad Unif','Box Ad Ave','Box Ad Vol Ave'};

PdeTypes = {'UniformPde','VolumeScaledAveragedPde';
            'UniformPde','VolumeScaledAveragedPde'};

BoxHs = {2,1,1.0/2.0,1.0/4.0};
BoxHs = {1.0/2.0,1.0/4.0,1.0/8.0};
colors = {'k.','b.','r.','g.'};
line_colors = {'k-','b-','r-','g-'};
boundary_line_colors = {'k:','b:','r:','g:'};


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

                pde_data = LoadNonConstantLengthData(['../../../testoutput/',dir,'/results_from_time_0/results.vizpdesolution'])
    
                timestep = min([plot_time_step,length(pde_data)]);
    
                end_time = pde_data{timestep}(1);
                node_id = pde_data{timestep}(2:4:end-3);
                x = pde_data{timestep}(3:4:end-2);
                y = pde_data{timestep}(4:4:end-1);
                c = pde_data{timestep}(5:4:end);
            
                r = sqrt(x.*x+y.*y);

                nodes_to_plot = find(r<=max_radius);


                r_keep = r(nodes_to_plot);
                c_keep = c(nodes_to_plot);

                step = 0.5    
                r_bins = 0:step:max_radius 
                r_bin_centres = r_bins(2:end) - 0.5*step;

                for i=1:length(r_bins)-1
                    c_set = c(all([r>=r_bins(i);r<r_bins(i+1)],1));
                    c_bin_mean(i) = mean(c_set);
                    c_bin_min(i) = min(c_set);
                    c_bin_max(i) = max(c_set);
                end



                %plot(r,c,'.','Color',[1.0-(BoxHIndex-1)/(length(BoxHs)-1),0, (BoxHIndex-1)/(length(BoxHs)-1),1.0])
                plot(r_bin_centres,c_bin_mean,line_colors{BoxHIndex});
                hold on
                plot(r_bin_centres,c_bin_min,boundary_line_colors{BoxHIndex});
                plot(r_bin_centres,c_bin_max,boundary_line_colors{BoxHIndex});

                clear c_bin_min c_bin_max c_bin_mean

                % main_legend{next_index} = legends{DomainTypeIndex,PdeTypeIndex};
                % next_index = next_index + 1;

                data = LoadNonConstantLengthData(['../../../testoutput/',dir,'/results_from_time_0/celldata_morphogen.dat'])

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
                alpha(0.5)
                hold on
                % 
                % % main_legend{next_index} = legends{DomainTypeIndex,PdeTypeIndex};
                % % next_index = next_index + 1;
            end

            title([TissueType,DomainType,PdeType])
            SaveAsPngEpsAndFig(-1,['Figs/',BaseType,TissueType,PdeType,'_no_exact'], 12, 7/5, 12);
            
            U = (max_radius-L0)/end_time;
        
            x = linspace(0,1,21);
            t = linspace(0,end_time,11);
            
            m = 1;
            if strcmp(DomainType,'BoxDomain')
                sol = pdepe(m,@pdepde,@pdeic,@pdebc,x,t);
            elseif strcmp(DomainType,'BoxDomainAdvection')
                sol = pdepe(m,@advectionpdepde,@pdeic,@pdebc,x,t);
            else
                assert(0);
            end
            
            u = sol(:,:,1);
            
            L = L0+U*end_time;
            plot(x*L,u(end,:),'k--')
            SaveAsPngEpsAndFig(-1,['Figs/',BaseType,TissueType,PdeType,'_no_legend'], 12, 7/5, 12);


            % main_legend{next_index} = 'Exact';
            % main_legend{next_index+1} = 'Advec Exact';
            % legend(main_legend,'Location','northwest')
            % 
            % 
            % xlim([0,L])
            % ylim([0.2,1.0])
            % 
            % % SaveAsPngEpsAndFig(-1,['Figs/',BaseType,TissueType,PdeType], 7, 7/5, 9);
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
