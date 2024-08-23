% Code to plot velocities within growing disk shaped tissue
clear; close all;



num_sims = 6
max_top_bin = 30

bin_u = cell(1,max_top_bin);

bin_edges = 0:max_top_bin;
bin_centres = (1:max_top_bin)-0.5;
bin_u = cell(1,max_top_bin);
bin_count = zeros(length(bin_centres),1);

for sim_index = 0:(num_sims-1)
    
    velocity_data = LoadNonConstantLengthData(['../../../testoutput/Parabolic/ExponentialVelocityTest/Run_',num2str(sim_index),'/results_from_time_0/cellvelocities.dat'])
    
    timestep = length(velocity_data);
    
    end_time = velocity_data{timestep}(1);
    cell_id = velocity_data{timestep}(2:5:end-4);
    x = velocity_data{timestep}(3:5:end-3);
    y = velocity_data{timestep}(4:5:end-2);
    u = velocity_data{timestep}(5:5:end-1);
    v = velocity_data{timestep}(6:5:end);
    
    r = sqrt(x.*x+y.*y);
    ur = sqrt(u.*u+v.*v);
    
    max_r = max(r);
    
    for i=1:length(r)
        bin = find(r(i)<bin_edges & r(i)+1>bin_edges)-1;
        assert(bin<=length(bin_centres));
        bin_count(bin) = bin_count(bin) + 1;
        bin_u{bin}(bin_count(bin),1) = ur(i);
    end
    
end

figure
    
    % pad with NaNs
    padded_bin_u = cellfun(@(x)[x(1:end); NaN(max(bin_count)-length(x),1)],bin_u,'UniformOutput',false);
     
    % make a matrix
    matrix_bin_u = cell2mat(padded_bin_u);
    
    boxplot(matrix_bin_u,bin_centres)
    
    
    % figure
    % for i=1:length(bin_centres)
    %   ave_ur(i) = mean(bin_u{i});
    % end
    % plot(bin_centres,ave_ur,'ko')
    % plot(r,ur,colors{DomainTypeIndex,PdeTypeIndex})
    hold on
    
    R = ceil(max_r)

    plot([0,R]+0.5,[0,0.5],'k--')
    x = linspace(0,R);
    y = 0.5*x.*x./R./R;
    plot(x,y,'k:')
    xlim([0,R]+0.5)
    ylim([-0.1,1.1])

    

