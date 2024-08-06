clear, close all

addpath ../../Chaste/anim/matlab/

% Dirs to load 



dir = {'TestEllipticBoxDomainAveragedPde'
       'TestEllipticBoxDomainVolumeDependentAveragedPde'
       'TestEllipticBoxDomainUniformPde'
       'TestEllipticGrowingDomainCellwisePde'
       'TestEllipticGrowingDomainVolumeDependentCellwisePde'
       'TestEllipticGrowingDomainUniformPde'}

colors = {'k.','r.','g.','ko','rx','g+'}
legends = {'Box Averaged',
           'Box Volume Averaged',
           'Box Uniform',
           'Growing Cellwise',
           'Growing Volume Cellwise',
           'Growing Uniform',
           'Exact Solution'}


for sim = 1:length(dir)
    data = LoadNonConstantLengthData(['../../../testoutput/',dir{sim},'/results_from_time_0/celldata_oxygen.dat'])
    
    timestep = 1
    
    
    node_id = data{1}(2:5:end-4);
    cell_id = data{1}(3:5:end-3);
    x = data{1}(4:5:end-2);
    y = data{1}(5:5:end-1);
    c = data{1}(6:5:end);
    
    
    r = sqrt(x.*x+y.*y);
    
    plot(r,c,colors{sim})
    hold on
end

R=sqrt(max(x.*x+y.*y));
C0=max(c);
r=linspace(0,R,1000);
k=0.1;

exact = C0*besseli(0,r*sqrt(k))/besseli(0,R*sqrt(k));

plot(r,exact,'k--');


legend(legends,'Location','northwest')