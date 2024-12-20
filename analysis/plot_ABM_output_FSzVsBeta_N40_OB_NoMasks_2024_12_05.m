% visualisation of ABM output summary stats as heatmap


clear all
close all

%% plot results from ABM 
% 
home_dirname = ['C:/Users/czachreson/Desktop/AIRBORNE/ward_model/ward_model_code/ward_model_repo/',...
                'TransmissionDynamics/',...
                'output/',...
                'shared_room_experiments/']

% home_dirname = ['/home/cameron/Desktop/compositions_in_progress/AIRBORNE/ward_model/',...
%                 'TransmissionDynamics/',...
%                 'output/',...
%                 'shared_room_experiments/']

experiment_type = 'OB'
%experiment_type = 'R0'

%index_case_type = 'nurse'
index_case_type = 'patient';

% index_case_config_vals = {'same', 'different'};
%index_case_config = 'different';


N = '1000'


% experiment_label = ['N' N '_full_run_pShared_N4_' index_case_config '_2024_11_25'];
% 
% experiment_dirname = [experiment_type '_p_vs_f_' index_case_type '_' experiment_label];

output_table = table(); 

% variables: 

index_case_config_vals = {'different'};%{'same', 'different'};



%ACH_p_vals = {'6p0', '12p0'};
ACH_p_vals = {'6p0'};%, '12p0'};

ACH_br_vals = {'2p0'};
%ACH_br_vals = {'2p0'};


Vol_fac_dbl_vals = {'1p5'};
%Vol_fac_dbl_vals = {'2p0'};

%pShared_vals = {'0p0', '0p5', '1p0'};
pShared_vals = {'1p0', '0p0'};%, '0p0'};%, '0p5', '1p0'};


beta_vals = {'100p0', '250p0' , '500p0', '1000p0', '2500p0', '5000p0', '10000p0', '25000p0', '50000p0', '100000p0'}

n_cols = 18;
n_rows = numel(ACH_p_vals) * numel(ACH_br_vals) * numel(Vol_fac_dbl_vals) * numel(pShared_vals) ;

vartypes = {'double', 'double', 'double', 'double','double','string'...
            'double', 'double', 'double', 'double',...
            'double', 'double', 'double', 'double',...
            'double', 'double', 'double', 'double'};
        
varnames = {'VDbl', 'ACHp', 'ACHb', 'pShared', 'beta', 'index_cfg',...
            'p1_AR', 'p10_AR', 'q90_AR', 'mean_AR',...
            'p1_ARpatient','p10_ARpatient', 'q90_ARpatient', 'mean_ARpatient',...
            'p1_ARnurse','p10_ARnurse', 'q90_ARnurse', 'mean_ARnurse'} ;


summary_output = table('Size', [n_rows, n_cols], 'VariableTypes', vartypes,...
                        'VariableNames', varnames); 

scenario_index = 0;


for index_case_cfg = index_case_config_vals
    index_case_config = index_case_cfg{:};


    
    experiment_label = ['N' N '_full_run_FSvsBeta_14d_sync_2024_12_03'];

    experiment_dirname = [experiment_type '_p_vs_f_' index_case_type '_' experiment_label];




for pShared = pShared_vals
    pShared_label = pShared{:};

    for beta = beta_vals
    beta_label = beta{:};

    
    for Vol_fac_dbl = Vol_fac_dbl_vals
        Vol_fac_dbl_label = Vol_fac_dbl{:};
    
        for ACH_p = ACH_p_vals 
            ACH_p_label = ACH_p{:};
            for ACH_br = ACH_br_vals
                ACH_br_label = ACH_br{:};


                
                scenario_index = scenario_index + 1; 

                data_flabel = ['pShared_' pShared_label,...
                    '_f_0p0_ACHpr_' ACH_p_label,...
                    '_ACHbr_' ACH_br_label,...
                    '_VfacDbl_' Vol_fac_dbl_label,...
                    '_beta_' beta_label '.csv'];

                data_fname = fullfile(home_dirname, experiment_dirname, 'detailed_output', data_flabel);
                data_table = readtable(data_fname);
                
                % compute summary statistics
                [q90_AR, p1_AR, p10_AR, avg_AR] = compute_summary_stats(data_table.AR) ;
                [q90_ARpatient, p1_ARpatient, p10_ARpatient, avg_ARpatient] = compute_summary_stats(data_table.AR_patient) ;
                [q90_ARnurse, p1_ARnurse, p10_ARnurse, avg_ARnurse] = compute_summary_stats(data_table.AR_nurse) ;
                
                % fill in scenario variables: 
                %{'VDbl', 'ACHp', 'ACHb', 'pShared', 'p10', 'q90', 'mean'} ;
                VDbl = str2double(strrep(Vol_fac_dbl_label, 'p', '.'));
                ACHp = str2double(strrep(ACH_p_label, 'p', '.'));
                ACHb = str2double(strrep(ACH_br_label, 'p', '.'));
                pSh = str2double(strrep(pShared_label, 'p', '.'));
                Bta = str2double(strrep(beta_label, 'p', '.'));
                
                summary_output.VDbl(scenario_index) = VDbl;
                summary_output.ACHp(scenario_index) = ACHp;
                summary_output.ACHb(scenario_index) = ACHb;
                summary_output.pShared(scenario_index) = pSh;
                summary_output.index_cfg(scenario_index) = index_case_config;
                summary_output.beta(scenario_index) = Bta;
                
                %p1 
                summary_output.p1_AR(scenario_index) = p1_AR;
                summary_output.p1_ARpatient(scenario_index) = p1_ARpatient;
                summary_output.p1_ARnurse(scenario_index) = p1_ARnurse;


                %p10
                summary_output.p10_AR(scenario_index) = p10_AR;
                summary_output.p10_ARpatient(scenario_index) = p10_ARpatient;
                summary_output.p10_ARnurse(scenario_index) = p10_ARnurse;
                
                %q90
                summary_output.q90_AR(scenario_index) = q90_AR;
                summary_output.q90_ARpatient(scenario_index) = q90_ARpatient;
                summary_output.q90_ARnurse(scenario_index) = q90_ARnurse;
                
                %mean
                summary_output.mean_AR(scenario_index) = avg_AR;
                summary_output.mean_ARpatient(scenario_index) = avg_ARpatient;
                summary_output.mean_ARnurse(scenario_index) = avg_ARnurse;
                
                % compile into full output table with extended column
                % labels: 
                
                % append unique label to dataframe columns
                var_label_prefix = ['pS' pShared_label,...
                    '_ACHp' ACH_p_label,...
                    '_ACHb' ACH_br_label,...
                    '_V' Vol_fac_dbl_label,...
                    '_Beta' beta_label];
                
                for i = 1:size(data_table.Properties.VariableNames, 2)
                       var_label_i = [var_label_prefix, '_', data_table.Properties.VariableNames{i}];
                       data_table.Properties.VariableNames{i} = var_label_i;
                end
                
                output_table = [output_table, data_table];
                
                
                
            end 
        end
    end
end

end

end


% visualise some comparisons: 

if strcmp(experiment_type, 'OB')
    xlabel_text = ['outbreak size (total cases)'];
elseif strcmp(experiment_type, 'R0')
    xlabel_text = ['secondary attack rate (1st generation only)'];
end

subtitle_text = ['index case type: ' index_case_type ', ' index_case_config, ' room'];

%% case 1: violin plots of outbreak size vs beta (single occupancy)
% compare to nurses

v1 = output_table.pS0p0_ACHp6p0_ACHb2p0_V1p5_Beta100p0_AR;
v2 = output_table.pS0p0_ACHp6p0_ACHb2p0_V1p5_Beta250p0_AR;
v3 = output_table.pS0p0_ACHp6p0_ACHb2p0_V1p5_Beta500p0_AR;
v4 = output_table.pS0p0_ACHp6p0_ACHb2p0_V1p5_Beta1000p0_AR;
v5 = output_table.pS0p0_ACHp6p0_ACHb2p0_V1p5_Beta2500p0_AR;
v6 = output_table.pS0p0_ACHp6p0_ACHb2p0_V1p5_Beta5000p0_AR;
v7 = output_table.pS0p0_ACHp6p0_ACHb2p0_V1p5_Beta10000p0_AR;
v8 = output_table.pS0p0_ACHp6p0_ACHb2p0_V1p5_Beta25000p0_AR;
v9 = output_table.pS0p0_ACHp6p0_ACHb2p0_V1p5_Beta50000p0_AR;
v10 = output_table.pS0p0_ACHp6p0_ACHb2p0_V1p5_Beta100000p0_AR;


input_mat = [v3, v4, v5, v6, v7, v8, v9, v10];

qs = 0.5%[q1, q2]

fid = 1;

x_vals = [1:3:24];
bandwidth = 1;
jitter_width = 2
c_code = 'b'

point_cloud_with_mean_and_quantile(input_mat, x_vals, qs, jitter_width, bandwidth, fid, c_code )

% case 2: violin plots of outbreak size vs beta (double-occupancy)
% compare to nurses

v1 = output_table.pS1p0_ACHp6p0_ACHb2p0_V1p5_Beta100p0_AR;
v2 = output_table.pS1p0_ACHp6p0_ACHb2p0_V1p5_Beta250p0_AR;
v3 = output_table.pS1p0_ACHp6p0_ACHb2p0_V1p5_Beta500p0_AR;
v4 = output_table.pS1p0_ACHp6p0_ACHb2p0_V1p5_Beta1000p0_AR;
v5 = output_table.pS1p0_ACHp6p0_ACHb2p0_V1p5_Beta2500p0_AR;
v6 = output_table.pS1p0_ACHp6p0_ACHb2p0_V1p5_Beta5000p0_AR;
v7 = output_table.pS1p0_ACHp6p0_ACHb2p0_V1p5_Beta10000p0_AR;
v8 = output_table.pS1p0_ACHp6p0_ACHb2p0_V1p5_Beta25000p0_AR;
v9 = output_table.pS1p0_ACHp6p0_ACHb2p0_V1p5_Beta50000p0_AR;
v10 = output_table.pS1p0_ACHp6p0_ACHb2p0_V1p5_Beta100000p0_AR;

input_mat = [v3, v4, v5, v6, v7, v8, v9, v10];

x_vals = [1:3:24]+1;
c_code = 'r'

point_cloud_with_mean_and_quantile(input_mat, x_vals, qs, jitter_width, bandwidth, fid, c_code )
xticks([1:3:24]+0.5)
xticklabels([500, 1000, 2500, 5000, 10000, 25000, 50000, 100000])
title(['outbrak size at day 14'])
xlabel('infectiousness')
ylabel(['outbreak size (total infections)'])
xlim([-1, 25])







%% local functions 


function violin_with_mean_and_quantile_nzv(input_matrix, q, figure_index)

    input_matrix_g0 = input_matrix;
    input_matrix_g0(input_matrix <= 0) = NaN;

    mean_vals = mean(input_matrix);
    q_vals = quantile(input_matrix, q);

    figure(figure_index)
    violinplot(input_matrix_g0);

    hold on
    plot(mean_vals, '-o')
    plot(q_vals, '-o')
    hold off


end

% fiddly. 


function point_cloud_plot(input_matrix, x_vals, jitter_width, bandwidth, figure_index)

    figure(figure_index)
    hold on

    % point clouds

    jitter = zeros(size(input_matrix));

    for col = 1:size(input_matrix, 2)
        [jitter(:, col), input_matrix(:, col)] = ksdensity(input_matrix(:, col), input_matrix(:, col), 'Bandwidth', bandwidth);
        jitter(:, col) = jitter(:, col) .* (1/max(jitter(:, col)));
    end

    x_jitter = repmat(x_vals, [size(input_matrix, 1), 1]);
    x_rnd = (0.5 - rand(size(x_jitter))) .* (jitter .* jitter_width); 

    x_jitter = x_jitter + x_rnd; 
    b = scatter(x_jitter, input_matrix, '.');

    f = [];
    xf = [];

    % violin plots: 
    for col = 1:size(input_matrix, 2)    
        [f(:, col), xf(:, col)] = kde(input_matrix(:, col), Bandwidth=bandwidth);

        %make my own damn violin plot.... 
        xvals = ones(size(f(:, col), 1), 1)*col - f(:, 1) .* (1/max(f(:, col)) .* (jitter_width/2) );
        plot(xvals, xf(:, 1), 'k-')
        
        xvals = ones(size(f(:, col), 1), 1)*col + f(:, 1) .* (1/max(f(:, col)) .* (jitter_width/2) );
        plot(xvals, xf(:, 1), 'k-')
        
        ylim([0, max(max(input_matrix))])
    end

    hold off


end

function point_cloud_with_mean_and_quantile(input_matrix, x_vals, q, jitter_width, bandwidth, figure_index, c_code)

    figure(figure_index)
    hold on

    % find maximum for scaling: 
    p_max_vals = [];
    for col = 1:size(input_matrix, 2)
    
        lower_bound = min(min(input_matrix));
        upper_bound = max(max(input_matrix));
        delta = 0.5;

        [p_kd, ~] = ksdensity(input_matrix(:, col), lower_bound:delta:upper_bound, 'Bandwidth', bandwidth);

        p_max_vals = [p_max_vals, max(p_kd)]; 

    end

    p_max = max(p_max_vals);



    % point clouds

   
    for col = 1:size(input_matrix, 2)
        
        [p_kde, ~] = ksdensity(input_matrix(:, col), input_matrix(:, col), 'Bandwidth', bandwidth);
        
        jitter_scale = p_kde .* (1/p_max) .* jitter_width;

        x_rnd = (0.5 - rand(size(p_kde))) .* jitter_scale; 

        x_jitter = x_vals(col) + x_rnd; 
        
        scatter(x_jitter, input_matrix(:, col), '.', c_code);

    end




    mean_vals = mean(input_matrix);
    q_vals = quantile(input_matrix, q);

    % violin plots: 
    for col = 1:size(input_matrix, 2)    

        upper_bound = max(input_matrix(:, col));

        [f, xf] = ksdensity(input_matrix(:, col), lower_bound:delta:upper_bound, 'Bandwidth', bandwidth);

        %make my own damn violin plot.... 
        xv = ones(size(f, 1), 1)*x_vals(col) - f .* (1/p_max) .* (jitter_width/2) ;
        plot(xv, xf, 'k-')
        xv1 = xv;

        xv = ones(size(f, 1), 1)*x_vals(col) + f .* (1/p_max) .* (jitter_width/2) ;
        plot(xv, xf, 'k-')
        xv2 = xv;
        
        % mean and quantile
        y_mean = mean_vals(col);
        y_quant = q_vals(col);
        [~, xid_mean] = min(abs(y_mean - xf)); 
        [~, xid_quant] = min(abs(y_quant - xf));

    



        x_mean = [xv1(xid_mean), xv2(xid_mean)] ;
        x_quant = [xv1(xid_quant), xv2(xid_quant)]; 

        min_line_length = jitter_width/4;

        if (x_mean(2) - x_mean(1)) < (min_line_length)
            x_mean(1) = -1 * min_line_length/2 + x_vals(col);
            x_mean(2) = min_line_length/2 + x_vals(col);
        end

        if (x_quant(2) - x_quant(1)) < (min_line_length)
            x_quant(1) = -1 * min_line_length/2 + x_vals(col);
            x_quant(2) = min_line_length/2 + x_vals(col);
        end



        plot(x_mean, [y_mean, y_mean], 'k-', 'Linewidth', 2)
        plot(x_quant, [y_quant, y_quant], 'k-')

    end

    % set x ticks. 
    xticks(x_vals)

    ylim([0, max(max(input_matrix))])

    hold off


end

function point_cloud_with_mean_and_quantile_NZ(input_matrix, x_vals, q, jitter_width, bandwidth, figure_index)


    input_matrix(input_matrix <= 0) = NaN;


    figure(figure_index)
    hold on

    % find maximum for scaling: 
    p_max_vals = [];
    for col = 1:size(input_matrix, 2)
    
        lower_bound = min(min(input_matrix));
        upper_bound = max(max(input_matrix));
        delta = 0.5;

        [p_kd, ~] = ksdensity(input_matrix(:, col), lower_bound:delta:upper_bound, 'Bandwidth', bandwidth);

        p_max_vals = [p_max_vals, max(p_kd)]; 

    end

    p_max = max(p_max_vals);



    % point clouds

   
    for col = 1:size(input_matrix, 2)
        
        [p_kde, ~] = ksdensity(input_matrix(:, col), input_matrix(:, col), 'Bandwidth', bandwidth);
        
        jitter_scale = p_kde .* (1/p_max) .* jitter_width;

        x_rnd = (0.5 - rand(size(p_kde))) .* jitter_scale; 

        x_jitter = x_vals(col) + x_rnd; 
        
        scatter(x_jitter, input_matrix(:, col), '.');

    end


    % violin plots: 
    for col = 1:size(input_matrix, 2)    

        upper_bound = max(input_matrix(:, col));

        [f, xf] = ksdensity(input_matrix(:, col), lower_bound:delta:upper_bound, 'Bandwidth', bandwidth);

        %make my own damn violin plot.... 
        xv = ones(size(f, 1), 1)*x_vals(col) - f .* (1/p_max) .* (jitter_width/2) ;
        plot(xv, xf, 'k-')
        xv1 = xv;

        xv = ones(size(f, 1), 1)*x_vals(col) + f .* (1/p_max) .* (jitter_width/2) ;
        plot(xv, xf, 'k-')
        xv2 = xv;
        
        % mean and quantile
        y_vals = input_matrix(:, col);
        y_vals = y_vals(~isnan(y_vals));

        y_mean = mean(y_vals);

        y_quant = quantile(y_vals, q);

        [~, xid_mean] = min(abs(y_mean - xf)); 
        [~, xid_quant] = min(abs(y_quant - xf));

        x_mean = [xv1(xid_mean), xv2(xid_mean)] ;
        x_quant = [xv1(xid_quant), xv2(xid_quant)]; 

        min_line_length = jitter_width/4;

        if (x_mean(2) - x_mean(1)) < (min_line_length)
            x_mean(1) = -1 * min_line_length/2 + x_vals(col);
            x_mean(2) = min_line_length/2 + x_vals(col);
        end

        if (x_quant(2) - x_quant(1)) < (min_line_length)
            x_quant(1) = -1 * min_line_length/2 + x_vals(col);
            x_quant(2) = min_line_length/2 + x_vals(col);
        end



        plot(x_mean, [y_mean, y_mean], 'k-', 'Linewidth', 2)
        plot(x_quant, [y_quant, y_quant], 'k-')

    end

    % set x ticks. 
    xticks(x_vals)

    ylim([0, max(max(input_matrix))])

    hold off


end

function [q90, p1, p10, avg] = compute_summary_stats(input_vec) 

        q90 = quantile(input_vec, 0.9);
        p10 = sum(input_vec > 10) / size(input_vec, 1);
        p1 = sum(input_vec >= 1) / size(input_vec, 1);
        avg = mean(input_vec);
end

function histogram_plot_with_quantiles_logy(input_matrix, qs, fid, leg_txt, nbins)
    % compute histogram
    max_x = max(max(input_matrix));
    dx = ceil(max_x/nbins);
    min_x = 0-dx/2;
    vmin = 1; 
    for i = 1:size(input_matrix, 2)
        figure(fid)
        h = histogram(input_matrix(:, i), 'BinEdges', min_x:dx:nbins*dx, 'Normalization','probability')
        vals = h.Values
        if vmin > min(vals(vals>0))
            vmin = min(vals(vals>0));
        end
        q_y = max(vals) + max(vals)/10
        % need a legend? 
     
        hold on 
        for j = 1:size(qs, 2)
            q_j  = quantile(input_matrix(:, i), qs(j));
            plot([q_j, q_j], [0.001, q_y])
        end
    end
    set(gca, 'YScale', 'log')
    legend(leg_txt)
    ylim([vmin/2, 1])
    hold off 
end

function histogram_plot_with_quantiles(input_matrix, qs, fid, leg_txt, nbins)
    % compute histogram
    max_x = max(max(input_matrix));
    dx = ceil(max_x/nbins);
    min_x = 0-dx/2;
    vmin = 1; 
    for i = 1:size(input_matrix, 2)
        figure(fid)
        h = histogram(input_matrix(:, i), 'BinEdges', min_x:dx:nbins*dx, 'Normalization','probability')
        vals = h.Values
        if vmin > min(vals(vals>0))
            vmin = min(vals(vals>0));
        end
        q_y = max(vals) + max(vals)/10
        % need a legend? 
     
        hold on 
        for j = 1:size(qs, 2)
            q_j  = quantile(input_matrix(:, i), qs(j));
            plot([q_j, q_j], [0.001, q_y])
        end
    end
    %set(gca, 'YScale', 'log')
    legend(leg_txt)
    %ylim([vmin/2, 1])
    hold off 
end

function histogram_bar_plot(input_matrix, qs, fid, leg_txt, nbins)
    % compute histogram
    max_x = max(max(input_matrix));
    dx = ceil(max_x/nbins);
    min_x = 0-dx/2;
    vmin = 1; 
    xtick_labels = 0:dx:(nbins-1)*dx;

    counts = [];

    for i = 1:size(input_matrix, 2)
        figure(fid)
        h = histcounts(input_matrix(:, i), 'BinEdges', min_x:dx:nbins*dx, 'Normalization','probability')
        counts = [counts, h']
        
        vals = h
        if vmin > min(vals(vals>0))
            vmin = min(vals(vals>0));
        end
        q_y = max(vals) + max(vals)/10
        % need a legend? 
     
        hold on 
        for j = 1:size(qs, 2)
            q_j  = quantile(input_matrix(:, i), qs(j));
            plot([q_j, q_j], [0.001, q_y])
        end
    end

    bar(xtick_labels', counts, 'grouped')
    
     xticks(xtick_labels)
    % 
     xticklabels(xtick_labels) 
     
    %set(gca, 'YScale', 'log')
    legend(leg_txt)
    %ylim([vmin/2, 1])
    hold off 
end