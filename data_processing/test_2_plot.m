clear; close all; clc;

% Add functions path
addpath('utils/')
addpath(genpath('/home/hernan/Git/matlab_tools/'))

% Workspace folder
path2folder = 'outputs/noisy/';

% Load workspace
load([path2folder,'workspace.mat'],...
         'nRMSE','MDE','nRMSE_CC','nRMSE_RR')

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
plot_line_width = 2.0;
plot_marker_size = 3.5;

% Plots visibility
visibility = 'off';

% Pixel size
% pxsz = 2*[1 1.5 2 2.5 3];
pxsz = 2*[3 2.5 2 1.5 1];
% dr = linspace(-2*0.17,2*0.17,5);
dr = linspace(-2*0.17,2*0.17,4);

% Cardiac phases
cp = [1 1.5 2.0 2.5 3.0];
labels = [true,false,false,false];
nlevel = 1:4;

%% ERROR PLOTS
% Cardiac phase of the plots
t = 10;

% Outlier marker
outlier = '';

% Number of different colors and markers
co = linspecer(5);
% co = distinguishable_colors(5,'w');

%% MAGNITUDE

% HARP
figure('Visible',visibility)
for r=1:5
    for n=nlevel
        if pxsz(r) == pxsz(r)+dr(n)
            boxplot(100*squeeze(nRMSE.HARP(:,n,r,t)),'Positions',[pxsz(r)+dr(n)],...
                    'Labels',{num2str(pxsz(r))},'symbol',outlier); hold on
        else
            boxplot(100*squeeze(nRMSE.HARP(:,n,r,t)),'Positions',[pxsz(r)+dr(n)],...
                    'symbol',outlier); hold on
        end
        set(findobj(gca,'type','line'),'linew',plot_line_width)
    end
end
hold off

% Plot formatting
api.XLabel = false;
api.YLabel = true;
api.YLabelStr = 'nRMSE (\%)';
api.Axis = [1.35 6.65 0 20];
api.YAxisTickValues = 0:5:100;
nice_boxplot(api);
print('-depsc','-r600',[path2folder,'1_MAG'])

% Get current axis and figure
ref_ax = gca;
ref_fig = gcf;
fig_pos = get(ref_fig,'position');
ax_pos = get(ref_ax,'position');
set(gcf,'position',[fig_pos(1) fig_pos(2) fig_pos(3) fig_pos(4)])
dx = 0.1;
set(gca,'position',[ax_pos(1) ax_pos(2)+dx ax_pos(3) ax_pos(4)-0.6*dx])


% SinMod
figure('Visible',visibility)
for r=1:5
    for n=nlevel
        if pxsz(r) == pxsz(r)+dr(n)
            boxplot(100*squeeze(nRMSE.SinMod(:,n,r,t)),'Positions',[pxsz(r)+dr(n)],...
                    'Labels',{num2str(pxsz(r))},'symbol',outlier); hold on
        else
            boxplot(100*squeeze(nRMSE.SinMod(:,n,r,t)),'Positions',[pxsz(r)+dr(n)],...
                    'symbol',outlier); hold on
        end
        set(findobj(gca,'type','line'),'linew',plot_line_width)        
    end
end
hold off
api.YLabel = false;
nice_boxplot(api);
print('-depsc','-r600',[path2folder,'2_MAG'])

% Format as the reference figure
set(gcf,'position',[fig_pos(1) fig_pos(2) fig_pos(3) fig_pos(4)])
set(gca,'position',[ax_pos(1) ax_pos(2)+dx ax_pos(3) ax_pos(4)-0.6*dx])

% DENSE
figure('Visible',visibility)
for r=1:5
    for n=nlevel
        if pxsz(r) == pxsz(r)+dr(n)
            boxplot(100*squeeze(nRMSE.DENSE(:,n,r,t)),'Positions',[pxsz(r)+dr(n)],...
                    'Labels',{num2str(pxsz(r))},'symbol',outlier); hold on
        else
            boxplot(100*squeeze(nRMSE.DENSE(:,n,r,t)),'Positions',[pxsz(r)+dr(n)],...
                    'symbol',outlier); hold on
        end
        set(findobj(gca,'type','line'),'linew',plot_line_width)        
    end
end
hold off
nice_boxplot(api);
print('-depsc','-r600',[path2folder,'3_MAG'])

% Format as the reference figure
set(gcf,'position',[fig_pos(1) fig_pos(2) fig_pos(3) fig_pos(4)])
set(gca,'position',[ax_pos(1) ax_pos(2)+dx ax_pos(3) ax_pos(4)-0.6*dx])


%% ANGULAR

% HARP
figure('Visible',visibility)
for r=1:5
    for n=nlevel
        if pxsz(r) == pxsz(r)+dr(n)
            boxplot(squeeze(MDE.HARP(:,n,r,t)),'Positions',[pxsz(r)+dr(n)],...
                    'Labels',{num2str(pxsz(r))},'symbol',outlier); hold on
        else
            boxplot(squeeze(MDE.HARP(:,n,r,t)),'Positions',[pxsz(r)+dr(n)],...
                    'symbol',outlier); hold on
        end
        set(findobj(gca,'type','line'),'linew',plot_line_width)    
    end
end
hold off

% Plot formatting
api.XLabel = true;
api.YLabel = true;
api.YLabelStr = 'DE ($^o$)';
api.Axis = [1.35 6.65 0 25];
api.YAxisTickValues = 0:5:100;
nice_boxplot(api);
print('-depsc','-r600',[path2folder,'1_DE'])

% Format as the reference figure
set(gcf,'position',[fig_pos(1) fig_pos(2) fig_pos(3) fig_pos(4)])
set(gca,'position',[ax_pos(1) ax_pos(2)+dx ax_pos(3) ax_pos(4)-0.6*dx])


% SinMod
figure('Visible',visibility)
for r=1:5
    for n=nlevel
        if pxsz(r) == pxsz(r)+dr(n)
            boxplot(squeeze(MDE.SinMod(:,n,r,t)),'Positions',[pxsz(r)+dr(n)],...
                    'Labels',{num2str(pxsz(r))},'symbol',outlier); hold on
        else
            boxplot(squeeze(MDE.SinMod(:,n,r,t)),'Positions',[pxsz(r)+dr(n)],...
                    'symbol',outlier); hold on
        end
        set(findobj(gca,'type','line'),'linew',plot_line_width)    
    end
end
hold off
api.YLabel = false;
nice_boxplot(api);
print('-depsc','-r600',[path2folder,'2_DE'])

% Format as the reference figure
set(gcf,'position',[fig_pos(1) fig_pos(2) fig_pos(3) fig_pos(4)])
set(gca,'position',[ax_pos(1) ax_pos(2)+dx ax_pos(3) ax_pos(4)-0.6*dx])

% DENSE
figure('Visible',visibility)
for r=1:5
    for n=nlevel
        if pxsz(r) == pxsz(r)+dr(n)
            boxplot(squeeze(MDE.DENSE(:,n,r,t)),'Positions',[pxsz(r)+dr(n)],...
                    'Labels',{num2str(pxsz(r))},'symbol',outlier); hold on
        else
            boxplot(squeeze(MDE.DENSE(:,n,r,t)),'Positions',[pxsz(r)+dr(n)],...
                    'symbol',outlier); hold on
        end
        set(findobj(gca,'type','line'),'linew',plot_line_width)    
    end
end
hold off
nice_boxplot(api);
print('-depsc','-r600',[path2folder,'3_DE'])

% Format as the reference figure
set(gcf,'position',[fig_pos(1) fig_pos(2) fig_pos(3) fig_pos(4)])
set(gca,'position',[ax_pos(1) ax_pos(2)+dx ax_pos(3) ax_pos(4)-0.6*dx])


%% CC

% HARP
figure('Visible',visibility)
for r=1:5
    for n=nlevel
        if pxsz(r) == pxsz(r)+dr(n)
            boxplot(squeeze(100*nRMSE_CC.HARP(:,n,r,t)),'Positions',[pxsz(r)+dr(n)],...
                    'Labels',{num2str(pxsz(r))},'symbol',outlier); hold on
        else
            boxplot(squeeze(100*nRMSE_CC.HARP(:,n,r,t)),'Positions',[pxsz(r)+dr(n)],...
                    'symbol',outlier); hold on
        end
        set(findobj(gca,'type','line'),'linew',plot_line_width)    
    end
end
hold off

% Plot formatting
api.XLabel = false;
api.YLabel = true;
api.YLabelStr = 'nRMSE CC (\%)';
api.Axis = [1.35 6.65 0 35];
api.YAxisTickValues = 0:10:100;
nice_boxplot(api);
print('-depsc','-r600',[path2folder,'1_CC'])

% Format as the reference figure
set(gcf,'position',[fig_pos(1) fig_pos(2) fig_pos(3) fig_pos(4)])
set(gca,'position',[ax_pos(1) ax_pos(2)+dx ax_pos(3) ax_pos(4)-0.6*dx])


% SinMod
figure('Visible',visibility)
for r=1:5
    for n=nlevel
        if pxsz(r) == pxsz(r)+dr(n)
            boxplot(squeeze(100*nRMSE_CC.SinMod(:,n,r,t)),'Positions',[pxsz(r)+dr(n)],...
                    'Labels',{num2str(pxsz(r))},'symbol',outlier); hold on
        else
            boxplot(squeeze(100*nRMSE_CC.SinMod(:,n,r,t)),'Positions',[pxsz(r)+dr(n)],...
                    'symbol',outlier); hold on
        end
        set(findobj(gca,'type','line'),'linew',plot_line_width)        
    end
end
hold off
api.YLabel = false;
nice_boxplot(api);
print('-depsc','-r600',[path2folder,'2_CC'])

% Format as the reference figure
set(gcf,'position',[fig_pos(1) fig_pos(2) fig_pos(3) fig_pos(4)])
set(gca,'position',[ax_pos(1) ax_pos(2)+dx ax_pos(3) ax_pos(4)-0.6*dx])

% DENSE
figure('Visible',visibility)
for r=1:5
    for n=nlevel
        if pxsz(r) == pxsz(r)+dr(n)
            boxplot(squeeze(100*nRMSE_CC.DENSE(:,n,r,t)),'Positions',[pxsz(r)+dr(n)],...
                    'Labels',{num2str(pxsz(r))},'symbol',outlier); hold on
        else
            boxplot(squeeze(100*nRMSE_CC.DENSE(:,n,r,t)),'Positions',[pxsz(r)+dr(n)],...
                    'symbol',outlier); hold on
        end
        set(findobj(gca,'type','line'),'linew',plot_line_width)        
    end
end
hold off
nice_boxplot(api);
print('-depsc','-r600',[path2folder,'3_CC'])

% Format as the reference figure
set(gcf,'position',[fig_pos(1) fig_pos(2) fig_pos(3) fig_pos(4)])
set(gca,'position',[ax_pos(1) ax_pos(2)+dx ax_pos(3) ax_pos(4)-0.6*dx])


%% RR

% HARP
figure('Visible',visibility)
for r=1:5
    for n=nlevel
        if pxsz(r) == pxsz(r)+dr(n)
            boxplot(squeeze(100*nRMSE_RR.HARP(:,n,r,t)),'Positions',[pxsz(r)+dr(n)],...
                    'Labels',{num2str(pxsz(r))},'symbol',outlier); hold on
        else
            boxplot(squeeze(100*nRMSE_RR.HARP(:,n,r,t)),'Positions',[pxsz(r)+dr(n)],...
                    'symbol',outlier); hold on
        end
        set(findobj(gca,'type','line'),'linew',plot_line_width)        
    end
end
hold off

% Plot formatting
api.XLabel = true;
api.YLabel = true;
api.YLabelStr = 'nRMSE RR (\%)';
api.Axis = [1.35 6.65 0 70];
api.YAxisTickValues = 0:10:150;
nice_boxplot(api);
print('-depsc','-r600',[path2folder,'1_RR'])

% Format as the reference figure
set(gcf,'position',[fig_pos(1) fig_pos(2) fig_pos(3) fig_pos(4)])
set(gca,'position',[ax_pos(1) ax_pos(2)+dx ax_pos(3) ax_pos(4)-0.6*dx])


% SinMod
figure('Visible',visibility)
for r=1:5
    for n=nlevel
        if pxsz(r) == pxsz(r)+dr(n)
            boxplot(squeeze(100*nRMSE_RR.SinMod(:,n,r,t)),'Positions',[pxsz(r)+dr(n)],...
                    'Labels',{num2str(pxsz(r))},'symbol',outlier); hold on
        else
            boxplot(squeeze(100*nRMSE_RR.SinMod(:,n,r,t)),'Positions',[pxsz(r)+dr(n)],...
                    'symbol',outlier); hold on
        end
        set(findobj(gca,'type','line'),'linew',plot_line_width)        
    end
end
hold off
api.YLabel = false;
nice_boxplot(api);
print('-depsc','-r600',[path2folder,'2_RR'])

% Format as the reference figure
set(gcf,'position',[fig_pos(1) fig_pos(2) fig_pos(3) fig_pos(4)])
set(gca,'position',[ax_pos(1) ax_pos(2)+dx ax_pos(3) ax_pos(4)-0.6*dx])

% DENSE
figure('Visible',visibility)
for r=1:5
    for n=nlevel
        if pxsz(r) == pxsz(r)+dr(n)
            boxplot(squeeze(100*nRMSE_RR.DENSE(:,n,r,t)),'Positions',[pxsz(r)+dr(n)],...
                    'Labels',{num2str(pxsz(r))},'symbol',outlier); hold on
        else
            boxplot(squeeze(100*nRMSE_RR.DENSE(:,n,r,t)),'Positions',[pxsz(r)+dr(n)],...
                    'symbol',outlier); hold on
        end
        set(findobj(gca,'type','line'),'linew',plot_line_width)        
    end
end
hold off
nice_boxplot(api);
print('-depsc','-r600',[path2folder,'3_RR'])

% Format as the reference figure
set(gcf,'position',[fig_pos(1) fig_pos(2) fig_pos(3) fig_pos(4)])
set(gca,'position',[ax_pos(1) ax_pos(2)+dx ax_pos(3) ax_pos(4)-0.6*dx])