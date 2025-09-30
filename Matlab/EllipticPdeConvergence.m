clear, close all

addpath ../../Chaste/anim/matlab/

% Dirs to load 

BaseType = 'EllipticConvergence';

TissueTypes = {'StaticDisc'};%,'StaticDiscApoptotic'};

DomainTypes = {'BoxDomain'};

PdeTypes = {'UniformPde', 'VolumeScaledAveragedPde'};

BcTypes = {'Box', 'Tissue'};

BoxHs = [4,2,1,1.0/2.0,1.0/4.0,1.0/8.0,1.0/16.0];





colors = {'kx-','rx-';
           'bx-','gx-'};

legends = {'Unif Box','Unif Tiss';
           'Ave Box','Ave Tiss'};

for TissueTypeIndex = 1:length(TissueTypes)
    TissueType = TissueTypes{TissueTypeIndex};

    for DomainTypeIndex = 1:length(DomainTypes)
        DomainType = DomainTypes{DomainTypeIndex};

        for BcTypeIndex = 1:length(BcTypes)
            BcType = BcTypes{BcTypeIndex};
            
            figure(2*(BcTypeIndex-1)+1)
            figure(2*(BcTypeIndex-1)+2)
            next_index = 1; %To keep legends aligned

            for PdeTypeIndex = 1:length(PdeTypes)
                PdeType = PdeTypes{DomainTypeIndex,PdeTypeIndex};
                    
                clear c_0;

                for BoxHIndex = 1:length(BoxHs)
                    BoxH = BoxHs(BoxHIndex);

                    dir = [BaseType, '/', TissueType, '/', DomainType, '/', PdeType, '/', BcType, '/', num2str(BoxH)]

                    data = LoadNonConstantLengthData(['../../../testoutput/',dir,'/results_from_time_0/celldata_oxygen.dat'])
        
                    timestep = length(data);
        
        
                    node_id = data{timestep}(2:5:end-4);
                    cell_id = data{timestep}(3:5:end-3);
                    x = data{timestep}(4:5:end-2);
                    y = data{timestep}(5:5:end-1);
                    c = data{timestep}(6:5:end);
        
                    r = sqrt(x.*x+y.*y);
        
                    [r_0(BoxHIndex) index] = min(r);
                    c_0(BoxHIndex)  = c(index);
       
                end
                figure(2*(BcTypeIndex-1)+1)
                loglog(BoxHs(1:end-1),abs(c_0(1:end-1)-c_0(end)),colors{PdeTypeIndex,BcTypeIndex});
                hold on

                figure(2*(BcTypeIndex-1)+2)
                semilogx(BoxHs(1:end),abs(c_0(1:end)),colors{PdeTypeIndex,BcTypeIndex});
                hold on

                main_legend{next_index} = legends{PdeTypeIndex,BcTypeIndex};
                next_index = next_index + 1;
            end
            
            figure(2*(BcTypeIndex-1)+1)
            title(BcType)
            legend(main_legend,'Location','northwest')
            SaveAsPngEpsAndFig(-1,['Figs/Error_',BaseType,TissueType,DomainType,BcType], 7, 7/5, 9);

            figure(2*(BcTypeIndex-1)+2)
            title(BcType)
            legend(main_legend,'Location','northwest')
            SaveAsPngEpsAndFig(-1,['Figs/Concentration_',BaseType,TissueType,DomainType,BcType], 7, 7/5, 9);

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
