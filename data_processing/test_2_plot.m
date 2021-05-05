clear; close all; clc;

% Add functions path
addpath('utils/')
addpath(genpath('/home/hernan/Git/matlab_tools/'))

% Workspace folder
path2folder = 'outputs/noisy/';

% Load workspace
load([path2folder,'workspace.mat'],'mean_SinMod_mag','mean_HARP_mag',...
    'mean_DENSE_mag','mean_SinMod_ang','mean_HARP_ang','mean_DENSE_ang',...
    'mean_SinMod_CC','mean_HARP_CC','mean_DENSE_CC','mean_SinMod_RR',...
    'mean_HARP_RR','mean_DENSE_RR')
load([path2folder,'workspace.mat'],'std_SinMod_mag','std_HARP_mag',...
    'std_DENSE_mag','std_SinMod_ang','std_HARP_ang','std_DENSE_ang',...
    'std_SinMod_CC','std_HARP_CC','std_DENSE_CC','std_SinMod_RR',...
    'std_HARP_RR','std_DENSE_RR')

% Plot settings
api = struct(...
    'AxesFontSize',  20,...
    'AxesLineWidth', 2,...
    'LegendFontSize', 17,...
    'Axis', [],...
    'XLabel', true,...
    'YLabel', true,...
    'XLabelStr', 'Pixel size (mm)',...
    'YLabelStr', [],...
    'YAxisTickValues', []);
plot_line_width = 3.5;
plot_marker_size = 8.5;

% Plots visibility
visibility = 'off';

% Pixel size
pxsz = 0.001;
dr = 0.17;

% Cardiac phases
cp = [1 1.5 2.0 2.5 3.0];
labels = [true,false,false,false];

%% ERROR PLOTS
% Cardiac phase, resolution and tag spacing of the plots
t = 8;
% t = 4;
r = 5:-1:1;

% Number of different colors and markers
% co = linspecer(5);
% co = distinguishable_colors(5,'w');
co = [0.0000 0.4470 0.7410
      0.8500 0.3250 0.0980
      0.9290 0.6940 0.1250
      0.4940 0.1840 0.5560
      0.4660 0.6740 0.1880
      0.3010 0.7450 0.9330
      0.6350 0.0780 0.1840];
lt = {'s','-','-.'};

%% MAGNITUDE

% Legends holders
h1 = zeros([1 4]);
h2 = zeros([1 2]);
h3 = zeros([1 4]);

% Plot error for each tag frequency
for foo=1:2

    figure('Visible',visibility)
    for f=1:4

        % Plot results
        if foo==1
            h1(f) = plot(cp,squeeze(mean_HARP_mag(f,r,t)),lt{2},'Color',co(f,:),...
                'LineWidth',plot_line_width); hold on        
        else
            h1(f) = plot(cp,squeeze(mean_SinMod_mag(f,r,t)),lt{3},'Color',co(f,:),...
                'LineWidth',plot_line_width); hold on
        end
        h3(f) = plot(cp,squeeze(mean_DENSE_mag(f,r,t)),lt{1},'Color',co(f,:),...
                'LineWidth',0.65*plot_line_width,'MarkerFaceColor',co(f,:),...
                'MarkerSize',plot_marker_size); hold on;
        
        % Get current axis and figure
        if f == 1
            ref_ax = gca;
            ref_fig = gcf;
            fp = get(ref_fig,'position');
            ap = get(ref_ax,'position');
            set(gcf,'position',[fp(1) fp(2) 0.7*fp(3) fp(4)])
            set(gca,'position',[1.5*ap(1) 2*ap(2) 0.85*ap(3) 0.9*ap(4)])
        end

%         % Print results
%         fprintf('\nnRMSE results\n')
%         fprintf('--------------------------------------------------------\n')
%         fprintf('HARP mean nRMSE: [%.1f, %.1f, %.1f, %.1f, %.1f]\n',squeeze(mean_HARP_mag(f,2:end)));
%         fprintf('HARP std nRMSE:  [%.1f, %.1f, %.1f, %.1f, %.1f]\n',squeeze(std_HARP_mag(f,2:end)))
%         fprintf('--------------------------------------------------------\n')
%         fprintf('SinMod mean nRMSE: [%.1f, %.1f, %.1f, %.1f, %.1f]\n',squeeze(mean_SinMod_mag(f,2:end)));
%         fprintf('SinMod std nRMSE:  [%.1f, %.1f, %.1f, %.1f, %.1f]\n',squeeze(std_SinMod_mag(f,2:end)))
%         fprintf('--------------------------------------------------------\n')
%         fprintf('HARP-I mean nRMSE: [%.1f, %.1f, %.1f, %.1f, %.1f]\n',squeeze(mean_DENSE_mag(f,2:end)));
%         fprintf('HARP-I std nRMSE:  [%.1f, %.1f, %.1f, %.1f, %.1f]\n',squeeze(std_DENSE_mag(f,2:end)))

    end

    % Send DENSE plots to the front
    for i=1:4
      uistack(h3(i),'top')
    end    
    
    % Dummy plots
    h2(1) = plot(NaN,NaN,lt{foo+1},'Color','k','LineWidth',plot_line_width); hold on;
    h2(2) = plot(NaN,NaN,lt{1},'Color','k','LineWidth',plot_line_width); hold off;

    % Plot formatting
    if foo==1
        api.XLabel = false;
        api.YLabel = labels(1);
        api.Legend = {'SP-HR','DENSEanalysis'};
    else
        api.XLabel = false;
        api.YLabel = false;        
        api.Legend = {'SinMod','DENSEanalysis'};
    end
    api.YLabelStr = 'nRMSE (\%)';
%     api.Axis = [0.75 3.25 0 15];
    api.Axis = [0.75 3.25 0 20];
    api.YAxisTickValues = 0:5:20;
    api.Legend1 = h1;
    api.Legend2 = h2;
    nice_plot_1(api);

    % Print results
    print('-depsc','-r600',[path2folder,sprintf('MAG_%01d',foo)])
    
end


%% ANGULAR

% Legends holders
h1 = zeros([1 4]);
h2 = zeros([1 2]);
h3 = zeros([1 4]);

% Plot error for each tag frequency
for foo=1:2

    figure('Visible',visibility)
    for f=1:4

        % Plot results
        if foo==1
            h1(f) = plot(cp,squeeze(mean_HARP_ang(f,r,t)),lt{2},'Color',co(f,:),...
                'LineWidth',plot_line_width); hold on        
        else
            h1(f) = plot(cp,squeeze(mean_SinMod_ang(f,r,t)),lt{3},'Color',co(f,:),...
                'LineWidth',plot_line_width); hold on
        end
        h3(f) = plot(cp,squeeze(mean_DENSE_ang(f,r,t)),lt{1},'Color',co(f,:),...
                'LineWidth',0.65*plot_line_width,'MarkerFaceColor',co(f,:),...
                'MarkerSize',plot_marker_size); hold on;

        % Get current axis and figure
        if f == 1
            ref_ax = gca;
            ref_fig = gcf;
            fp = get(ref_fig,'position');
            ap = get(ref_ax,'position');
            set(gcf,'position',[fp(1) fp(2) 0.7*fp(3) fp(4)])
            set(gca,'position',[1.5*ap(1) 2*ap(2) 0.85*ap(3) 0.9*ap(4)])
        end

    end

    % Send DENSE plots to the front
    for i=1:4
      uistack(h3(i),'top')
    end    
    
    % Dummy plots
    h2(1) = plot(NaN,NaN,lt{foo+1},'Color','k','LineWidth',plot_line_width); hold on;
    h2(2) = plot(NaN,NaN,lt{1},'Color','k','LineWidth',plot_line_width); hold off;

    % Plot formatting
    if foo==1
        api.XLabel = false;
        api.YLabel = labels(1);
        api.Legend = {'SP-HR','DENSEanalysis'};
    else
        api.XLabel = false;
        api.YLabel = false;        
        api.Legend = {'SinMod','DENSEanalysis'};
    end
    api.YLabelStr = 'DE ($^o$)';
    api.Axis = [0.75 3.25 0 15];
    api.YAxisTickValues = 0:5:25;
    api.Legend1 = h1;
    api.Legend2 = h2;
    nice_plot_1(api);

    % Print results
    print('-depsc','-r600',[path2folder,sprintf('ANG_%01d',foo)])
    
end


%% CIRCUMFERENTIAL

% Legends holders
h1 = zeros([1 4]);
h2 = zeros([1 2]);
h3 = zeros([1 4]);

% Plot error for each tag frequency
for foo=1:2

    figure('Visible',visibility)
    for f=1:4

        % Plot results
        if foo==1
            h1(f) = plot(cp,squeeze(mean_HARP_CC(f,r,t)),lt{2},'Color',co(f,:),...
                'LineWidth',plot_line_width); hold on        
        else
            h1(f) = plot(cp,squeeze(mean_SinMod_CC(f,r,t)),lt{3},'Color',co(f,:),...
                'LineWidth',plot_line_width); hold on
        end
        h3(f) = plot(cp,squeeze(mean_DENSE_CC(f,r,t)),lt{1},'Color',co(f,:),...
                'LineWidth',0.65*plot_line_width,'MarkerFaceColor',co(f,:),...
                'MarkerSize',plot_marker_size); hold on;
        
        % Get current axis and figure
        if f == 1
            ref_ax = gca;
            ref_fig = gcf;
            fp = get(ref_fig,'position');
            ap = get(ref_ax,'position');
            set(gcf,'position',[fp(1) fp(2) 0.7*fp(3) fp(4)])
            set(gca,'position',[1.5*ap(1) 2*ap(2) 0.85*ap(3) 0.9*ap(4)])
        end

    end

    % Send DENSE plots to the front
    for i=1:4
      uistack(h3(i),'top')
    end
    
    % Dummy plots
    h2(1) = plot(NaN,NaN,lt{foo+1},'Color','k','LineWidth',plot_line_width); hold on;
    h2(2) = plot(NaN,NaN,lt{1},'Color','k','LineWidth',plot_line_width); hold off;

    % Plot formatting
    if foo==1
        api.XLabel = true;
        api.YLabel = true;
        api.Legend = {'SP-HR','DENSEanalysis'};
    else
        api.XLabel = true;
        api.YLabel = false;        
        api.Legend = {'SinMod','DENSEanalysis'};
    end
    api.YLabelStr = 'nRMSE CC (\%)';
    api.Axis = [0.75 3.25 0 50];
    api.YAxisTickValues = 0:10:50;
    api.Legend1 = h1;
    api.Legend2 = h2;
    nice_plot_1(api);

    % Print results
    print('-depsc','-r600',[path2folder,sprintf('CC_%01d',foo)])

end


%% RADIAL

% Legends holders
h1 = zeros([1 4]);
h2 = zeros([1 2]);
h3 = zeros([1 4]);

% Plot error for each tag frequency
for foo=1:2

    figure('Visible',visibility)
    for f=1:4

        % Plot results
        if foo==1
            h1(f) = plot(cp,squeeze(mean_HARP_RR(f,r,t)),lt{2},'Color',co(f,:),...
                'LineWidth',plot_line_width); hold on        
        else
            h1(f) = plot(cp,squeeze(mean_SinMod_RR(f,r,t)),lt{3},'Color',co(f,:),...
                'LineWidth',plot_line_width); hold on
        end
        h3(f) = plot(cp,squeeze(mean_DENSE_RR(f,r,t)),lt{1},'Color',co(f,:),...
                'LineWidth',0.65*plot_line_width,'MarkerFaceColor',co(f,:),...
                'MarkerSize',plot_marker_size); hold on;
        
        % Get current axis and figure
        if f == 1
            ref_ax = gca;
            ref_fig = gcf;
            fp = get(ref_fig,'position');
            ap = get(ref_ax,'position');
            set(gcf,'position',[fp(1) fp(2) 0.7*fp(3) fp(4)])
            set(gca,'position',[1.5*ap(1) 2*ap(2) 0.85*ap(3) 0.9*ap(4)])
        end

    end
    
    % Send DENSE plots to the front
    for i=1:4
      uistack(h3(i),'top')
    end

    % Dummy plots
    h2(1) = plot(NaN,NaN,lt{foo+1},'Color','k','LineWidth',plot_line_width); hold on;
    h2(2) = plot(NaN,NaN,lt{1},'Color','k','LineWidth',plot_line_width); hold off;

    % Plot formatting
    if foo==1
        api.XLabel = true;
        api.YLabel = true;
        api.Legend = {'SP-HR','DENSEanalysis'};
    else
        api.XLabel = true;
        api.YLabel = false;        
        api.Legend = {'SinMod','DENSEanalysis'};
    end
    api.YLabelStr = 'nRMSE RR (\%)';
    api.Axis = [0.75 3.25 0 65];
    api.YAxisTickValues = 0:20:90;
    api.Legend1 = h1;
    api.Legend2 = h2;
    nice_plot_1(api);

    % Print results
    print('-depsc','-r600',[path2folder,sprintf('RR_%01d',foo)])

end

% clear all