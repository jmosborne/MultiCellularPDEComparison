clear, close all

addpath ../../Chaste/anim/matlab/
addpath MinBoundSuite/MinBoundSuite/

% Dirs to load 

BaseType = 'EllipticConvergence';

TissueTypes = {'GrowingTissue'};

DomainTypes = {'BoxDomain'};

PdeTypes = {'VolumeScaledAveragedPde'};

ComparisonPdeDir =  'GrowingTissue/GrowingDomain/CellwisePde/MeshBc/8'

BcTypes = {'BoxBc'};
%BcTypes = {'TissueBc'};
%BcTypes = {'MeshBc'};


%BoxHs = [8,4,2,1.0,1.0/2.0,1.0/4.0,1.0/8.0];
%BoxHs = [1.0/16.0,1.0/8.0,1.0/4.0,1.0/2.0,1.0,2.0,4.0,8.0];
BoxHs = [1.0/16,1.0/8.0,1.0/4.0,1.0/2.0,1.0,2.0,4.0,8.0];
BoxHsString = {'1/16','1/8','1/4','1/2','1','2','4','8'}
assert(length(BoxHsString)==length(BoxHs));


Runs = [0:11, 20:31, 40:49, 60:71, 80:91];
%Runs = [0:20:99];

colors = {'kx-','rx-'};

legends = {'All','Prolif'};
    

for RunIndex = 1:1%00
    Run = RunIndex-1;

    dir = [BaseType, '/Run', num2str(Run), '/',ComparisonPdeDir]

    data = LoadNonConstantLengthData(['../../../testoutput/',dir,'/results_from_time_0/results.vizcelltypes'])

    timestep = length(data);


    cell_types = data{timestep}(2:end);
    
    num_all_cells_comparison(RunIndex) = length(cell_types);
    num_prolif_cells_comparison(RunIndex) = length(find(cell_types==1));
end

mean_num_all_cells_comparison = mean(num_all_cells_comparison);
mean_num_prolif_cells_comparison = mean(num_prolif_cells_comparison)

for TissueTypeIndex = 1:length(TissueTypes)
    TissueType = TissueTypes{TissueTypeIndex};
    
    for DomainTypeIndex = 1:length(DomainTypes)
        DomainType = DomainTypes{DomainTypeIndex};

        for BcTypeIndex = 1:length(BcTypes)
            BcType = BcTypes{BcTypeIndex};

            for PdeTypeIndex = 1:length(PdeTypes)
                PdeType = PdeTypes{DomainTypeIndex,PdeTypeIndex};
                    
                clear c_0;


                figure 
                hold on

                for BoxHIndex = 1:length(BoxHs)
                    BoxH = BoxHs(BoxHIndex);

                    for RunIndex = 1:length(Runs)
                        Run = Runs(RunIndex)

                        dir = [BaseType, '/Run', num2str(Run), '/', TissueType, '/', DomainType, '/', PdeType, '/', BcType, '/', num2str(BoxH)]
    
                        data = LoadNonConstantLengthData(['../../../testoutput/',dir,'/results_from_time_0/results.vizcelltypes']);
            
                        timestep = length(data);
            
            
                        cell_types = data{timestep}(2:end);
                        
                        num_all_cells(BoxHIndex,RunIndex) = length(cell_types);
                        num_prolif_cells(BoxHIndex,RunIndex) = length(find(cell_types==1));

                        % Now plot tissue Boundary
                        node_data = LoadNonConstantLengthData(['../../../testoutput/',dir,'/results_from_time_0/results.viznodes']);
                        boundary_data = LoadNonConstantLengthData(['../../../testoutput/',dir,'/results_from_time_0/results.vizboundarynodes']);
            
                        node_points = node_data{timestep};
                        boundary_points = boundary_data{timestep}(2:end);

                        node_x = node_points(2:2:end-1);
                        node_y = node_points(3:2:end);
                        boundary_nodes = find(boundary_points);

                        % plot(node_x,node_y,'ko')
                        % hold on 
                        % plot(node_x(boundary_nodes),node_y(boundary_nodes),'rx')
                        
                        [outer_center,outer_radius] = minboundcircle(node_x,node_y);

                        % Don't need boundary_nbodes as uses convex hull 
                        % [inner_center,inner_radius] = incircle(node_x(boundary_nodes),node_y(boundary_nodes));
                        [inner_center,inner_radius] = incircle(node_x,node_y);
                        

                        % theta = linspace(0,2*pi,100);
                        % xc = outer_center(1) + outer_radius*cos(theta);
                        % yc = outer_center(2) + outer_radius*sin(theta);
                        % plot(xc,yc,'b-')
                        % xc = inner_center(1) + inner_radius*cos(theta);
                        % yc = inner_center(2) + inner_radius*sin(theta);
                        % plot(xc,yc,'r-')
                        
                        bounding_sphere(BoxHIndex,RunIndex) = outer_radius;
                        bounding_sphere_ratio(BoxHIndex,RunIndex) = inner_radius/outer_radius;

                        K = convhull(node_x,node_y);

                        plot(node_x(K),node_y(K),'-','Color',[1.0-(BoxHIndex-1)/(length(BoxHs)-1),0, (BoxHIndex-1)/(length(BoxHs)-1),0.2])

                        hold on
                    end
                end

                axis equal
                axis off

                SaveAsPngEpsAndFig(-1,['Figs/Boundary_',BaseType,TissueType,DomainType,BcType,PdeType], 7, 7/5, 9);
                
  
                figure
                boxchart(num_all_cells', 'BoxFaceColor', "#0072BD")
                xticklabels(BoxHsString)
                hold on 
                plot(1:length(BoxHs), mean(num_all_cells,2),'bx-')
                hold on
                if (BcType(1) ~= 'B') 
                    plot([0,length(BoxHs)+1],mean_num_all_cells_comparison*[1,1],'k--')
                end
                title('Num Cells')
                SaveAsPngEpsAndFig(-1,['Figs/CellNumBox_',BaseType,TissueType,DomainType,BcType,PdeType], 7, 7/5, 9);
                
                figure 
                boxchart(num_prolif_cells', 'BoxFaceColor',"#D95319")
                xticklabels(BoxHsString)
                hold on 
                plot(1:length(BoxHs), mean(num_prolif_cells,2),'rx-')
                hold on 
                if (BcType(1) ~= 'B') 
                    plot([0,length(BoxHs)+1],mean_num_prolif_cells_comparison*[1,1],'k--')
                end
                title('Num Prolif Cells')
                SaveAsPngEpsAndFig(-1,['Figs/ProlifCellNumBox_',BaseType,TissueType,DomainType,BcType,PdeType], 7, 7/5, 9);

                figure 
                boxchart(bounding_sphere_ratio', 'BoxFaceColor',"#D95319")
                xticklabels(BoxHsString)
                hold on 
                plot(1:length(BoxHs), mean(bounding_sphere_ratio,2),'rx-')
                title('Bounding Sphere Ratio')
                SaveAsPngEpsAndFig(-1,['Figs/BoundingSphereRatioBox_',BaseType,TissueType,DomainType,BcType,PdeType], 7, 7/5, 9);
                
                figure 
                boxchart(bounding_sphere', 'BoxFaceColor',"#D95319")
                xticklabels(BoxHsString)
                hold on 
                plot(1:length(BoxHs), mean(bounding_sphere,2),'rx-')
                title('Bounding Sphere Radius')
                SaveAsPngEpsAndFig(-1,['Figs/BoundingSphereBox_',BaseType,TissueType,DomainType,BcType,PdeType], 7, 7/5, 9);

            end
        end
    end
    

    % R=sqrt(max(x.*x+y.*y));
    % C0=max(c);
    % r=linspace(0,R,1000);
    % k=0.1;
    % 
    % exact = C0*besseli(0,r*sqrt(k))/besseli(0,R*sqrt(k));
    % 
    % plot(r,exact,'k--');

    % SaveAsPngEpsAndFig(-1,['Figs/',BaseType,TissueType], 7, 7/5, 9);
    
end
