clear, close all

addpath ../../Chaste/anim/matlab/

% Dirs to load 

BaseType = 'Elliptic';

TissueTypes = {'StaticDisc','StaticDiscApoptotic'};

DomainTypes = {'GrowingDomain', 'BoxDomain'};

PdeTypes = {'UniformPde', 'CellwisePde', 'VolumeScaledCellwisePde';
            'UniformPde', 'AveragedPde', 'VolumeScaledAveragedPde'};

colors = {'kd','rs','g^';
    'k*','rx','g+';
    'k.','r.','g.'};

legends = {'Grow Unif','Grow Cell','Grow Vol Cell';
           'Box Unif','Box Ave','Box Vol Ave'};

for TissueTypeIndex = 1:length(TissueTypes)
    TissueType = TissueTypes{TissueTypeIndex};

    figure
    next_index = 1; %To keep legends aligned

    for DomainTypeIndex = 1:length(DomainTypes)
        DomainType = DomainTypes{DomainTypeIndex};

        for PdeTypeIndex = 1:length(PdeTypes)
            PdeType = PdeTypes{DomainTypeIndex,PdeTypeIndex};

            dir = [BaseType, '/', TissueType, '/', DomainType, '/', PdeType]

            data = LoadNonConstantLengthData(['../../../testoutput/',dir,'/results_from_time_0/celldata_oxygen.dat'])

            timestep = length(data);


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
            
        end
    end
    title(TissueType)
    legend(main_legend,'Location','northwest')


    R=sqrt(max(x.*x+y.*y));
    C0=max(c);
    r=linspace(0,R,1000);
    k=0.1;
    
    exact = C0*besseli(0,r*sqrt(k))/besseli(0,R*sqrt(k));
    
    plot(r,exact,'k--');

    % SaveAsPngEpsAndFig(-1,['Figs/',BaseType,TissueType], 7, 7/5, 9);
    SaveAsPngEpsAndFig(-1,['Figs/',BaseType,TissueType], 12, 7/5, 12);
end
