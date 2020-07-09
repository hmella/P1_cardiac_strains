%% PATHS
% Inputs path
input_folder = '../data_generation/inputs/';
addpath(input_folder)

% functions paths
addpath(genpath('utils/'))
addpath(genpath('/home/hernan/cardiac_motion'))
addpath(genpath('/home/hernan/denseanalysis'))


% Output folder
outputs = {'outputs/3D/Exact/',...
           'outputs/3D/HARP/',...
           'outputs/3D/SinMod/',...
           'outputs/3D/DENSE/',};
for i=1:numel(outputs)
   mkdir(outputs{i});
end
   
%% INPUT DATA
% Analysis to be performed
RUN_EXACT   = false;
RUN_HARP    = false;
RUN_SinMod  = false;
RUN_DENSE   = false;
STRAIN_CURVES = true;

% Number of cardiac phases
fr  = 1:20;
Nfr = numel(fr);

% Temporal fitting parameters
tfit_args = struct('Mask',           [],...
                   'Frames',         1:Nfr,...
                   'TemporalOrder',  10,...
                   'Show',           false);

%% FILTERS SPECS (for image processing)
% Filter specs
KspaceFilter = 'Transmission';
BTW_cutoff = 1;
BTW_order  = [];
KspaceFollowing = false;

%% IMAGING PARAMETERS
% Resolutions and FOV
resD = [129 129];
resC = [257 257];
FOV = [0.35 0.35];
pxszD = FOV./resD;
pxszC = FOV./resC;

% Encoding frequencies
tag_spac = 0.0080;          % [m]
ke_spamm = 2*pi./tag_spac;  % [rad/m]
ke_dense = 0.12*1000;       % [rad/m]

%% MAIN CODE

% Plot kspaces
for view={'base','mid','apex'}
    % Load SPAMM image
    filename = sprintf('CI_%s.mat',view{1});
    IPath = [input_folder,'3D_experiments/',filename];
    load(IPath);
    M = squeeze(I.M(:,:,1,1,:));
    IC = (squeeze(I.I1 - I.I2));
    IC = squeeze(I.I1);

    % Re-format the images
    Nfr = size(IC,4);
    for i=1:Nfr
        M(:,:,i) = M(:,:,i)';
        IC(:,:,1,i) = IC(:,:,1,i)';
        IC(:,:,2,i) = IC(:,:,2,i)';
    end
    
    % Load DENSE image
    filename = sprintf('DI_%s.mat',view{1});
    IPath = [input_folder,'3D_experiments/',filename];
    load(IPath);
    M = squeeze(I.M(:,:,1,1,:));
    ID = squeeze(I.I1 - I.I2);
    ID = squeeze(I.I1);

    % Re-format the images
    Nfr = size(ID,4);
    for i=1:Nfr
        M(:,:,i) = M(:,:,i)';
        ID(:,:,1,i) = ID(:,:,1,i)';
        ID(:,:,2,i) = ID(:,:,2,i)';
    end

%     % Image filtering
%     h = ButterworthFilter(Isz(1:2),Isz(1:2)/2,40,5);
%     ID = ktoi(h.*itok(ID));

    % Plots
    h = figure('Visible','off');
    t = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
    nexttile
    imagesc(abs(itok(IC(:,:,1,8)))); colormap gray
    axis equal off; caxis([0 3.0e+04])
    nexttile
    imagesc(abs((IC(:,:,1,8)))); colormap gray
    axis equal off;
    nexttile
    imagesc(abs(itok(ID(:,:,1,8)))); colormap gray
    axis equal off; caxis([0 1.75e+04])
    nexttile
    imagesc(angle((ID(:,:,1,8)))); colormap gray
    axis equal off
    print('-depsc','-r1200',sprintf('kspaces_%s',view{1}))

end
    
%% EXACT ANALYSIS
if RUN_EXACT
    for view={'base','mid','apex'}
        for name={'RCI','RDI'}
            
            % Debug
            fprintf('\n  Processing reference data (%s view)',view{1})

            % Load data
            filename = sprintf('%s_%s.mat',name{1},view{1});
            IPath = [input_folder,'3D_experiments/',filename];
            load(IPath);
            M = squeeze(I.M(:,:,1,1,:));
            I = squeeze(I.I1);

            % Re-format the images
            Nfr = size(I,4);
            for i=1:Nfr
                M(:,:,i) = M(:,:,i)';
                I(:,:,1,i) = I(:,:,1,i)';
                I(:,:,2,i) = I(:,:,2,i)';
            end

            % Image size
            Isz = size(M);
            
            % Pixel size
            if name{1}=='RDI'
                pxsz = pxszD;
            else
                pxsz = pxszC;
            end
            
            % Get displacements in pixels
            ue = 0.01*angle(I)/pxsz(1);
            
            % Displacement
            ux_EXACT = squeeze(ue(:,:,1,:));
            uy_EXACT = squeeze(ue(:,:,2,:));                
            dxe = ux_EXACT(:,:,fr);    % millimeters
            dye = uy_EXACT(:,:,fr);    % millimeters
            dxe(~repmat(M(:,:,1),[1 1 Nfr])) = nan;
            dye(~repmat(M(:,:,1),[1 1 Nfr])) = nan;         
            
            % Strain
            [X, Y] = meshgrid(1:size(ux_EXACT,2), 1:size(ux_EXACT,1));
            options = struct(...
                'X', X,...
                'Y', Y,...
                'mask',M(:,:,1),...
                'times',1:Nfr,...
                'dx', dxe,...
                'dy', dye,...
                'Origin', [],...
                'Orientation', []);
            st = mypixelstrain(options);
            RR_EXACT = NaN([Isz(1) Isz(2) Nfr]);
            CC_EXACT = NaN([Isz(1) Isz(2) Nfr]);
            RR_EXACT(repmat(st.maskimage(:,:,1),[1 1 Nfr])) = st.RR(:);
            CC_EXACT(repmat(st.maskimage(:,:,1),[1 1 Nfr])) = st.CC(:);

%             figure(1)
%             sgtitle(sprintf('%s (%s)',name{1},view{1}))
%             subplot 121;
%             imagesc(CC_EXACT(:,:,8),'AlphaData',st.maskimage);
%             colormap(flipud(jet)); colorbar;
%             subplot 122;
%             imagesc(RR_EXACT(:,:,8),'AlphaData',st.maskimage);
%             colormap(flipud(jet)); colorbar;
%             pause
            
             
%             if name{1}=='RDI'
%                 figure(2)
%                 errorbar(1:Nfr,mean(harpi_C.segments_RR),std(harpi_C.segments_RR),'-','MarkerSize',3,'LineWidth',2)
%                 hold on
%                 errorbar(1:Nfr,mean(harpi_D.segments_RR),std(harpi_D.segments_RR),'-','MarkerSize',3,'LineWidth',2)
%                 hold on
%                 errorbar(1:Nfr,mean(harpi_C.segments_CC),std(harpi_C.segments_CC),'-','MarkerSize',3,'LineWidth',2)
%                 hold on
%                 errorbar(1:Nfr,mean(harpi_D.segments_CC),std(harpi_D.segments_CC),'-','MarkerSize',3,'LineWidth',2)
%                 hold off
%                 l = legend('Radial C','Radial D','Circumferential C','Circumferential D');
%                 l.Interpreter = 'latex';
%                 pause
%             end
            
            % Write exact displacement and strain
            mask_exact = st.maskimage(:,:,1);
            save([outputs{1},filename],'dxe','dye','RR_EXACT','CC_EXACT','mask_exact');

        end

    end
end


%% HARP analysis
if RUN_HARP
    for view={'base','mid','apex'}

        % Debug
        fprintf('\n  Processing CSPAMM data using HARP (%s view)',view{1})

        % Load data
        filename = sprintf('CI_%s.mat',view{1});
        IPath = [input_folder,'3D_experiments/',filename];
        load(IPath);
        M = squeeze(I.M(:,:,1,1,:));
        I = (squeeze(I.I1 - I.I2));
        
        % Re-format the images
        Nfr = size(I,4);
        for i=1:Nfr
            M(:,:,i) = M(:,:,i)';
            I(:,:,1,i) = I(:,:,1,i)';
            I(:,:,2,i) = I(:,:,2,i)';
        end

%         figure(1)
%         imagesc(abs(itok(I(:,:,10))))
%         colormap gray
%         caxis([0 2500])
%         axis equal off
%         pause        
        
        % SPAMM encoding frequency
        ke = [ke_spamm ke_spamm];

        % Image size
        Isz = size(I);

        % HARP displacements
        try
            args = struct(...
                'Mask',             M,...
                'EncFreq',          ke,...
                'FOV',              Isz(1:2).*pxszC,...
                'PixelSize',        [1 1].*pxszC,...
                'Frames',           fr,...
                'tol',              1e-2,...
                'maxiter',          30,...
                'GradientEval',     5,...
                'SearchWindow',     [0,0],...
                'PhaseWindow',      [2,2],...
                'Show',             false,...
                'ShowConvergence',  false,...
                'Seed',             'auto',...
                'theta',            [0 pi/2],...
                'Connectivity',     8,...
                'KspaceFilter',     KspaceFilter,...
                'BTW_cutoff',       BTW_cutoff,...
                'BTW_order',        BTW_order,...
                'KspaceFollowing',  KspaceFollowing);
            [dxh, dyh] = HARPTrackingOsman(I, args);

            % Prepare displacements
            displacements = permute(cat(4,dxh,dyh),[1 2 4 3]);

            % Temporal fitting
            tfit_args.Mask = M(:,:,1);
            [dxh,dyh] = TemporalFitting(displacements,tfit_args);

            % HARP strain
            [X, Y] = meshgrid(1:size(dxh,2), 1:size(dyh,1));
            options = struct(...
                'X', X,...
                'Y', Y,...
                'mask',M(:,:,1),...
                'times',1:Nfr,...
                'dx', dxh,...
                'dy', dyh,...
                'Origin', [],...
                'checknans',  false,...
                'Orientation', []);
            st = mypixelstrain(options);
            RR_HARP = NaN([Isz(1) Isz(2) Nfr]);
            CC_HARP = NaN([Isz(1) Isz(2) Nfr]);
            RR_HARP(repmat(st.maskimage(:,:,1),[1 1 Nfr])) = st.RR(:);
            CC_HARP(repmat(st.maskimage(:,:,1),[1 1 Nfr])) = st.CC(:);

%             figure(1)
%             subplot 121
%             imagesc(CC_HARP(:,:,8),'AlphaData',st.maskimage); colorbar
%             subplot 122
%             imagesc(RR_HARP(:,:,8),'AlphaData',st.maskimage); colorbar
%             pause(0.1)

            % Write displacement and strain
            mask_harp = st.maskimage(:,:,1);
            save([outputs{2},filename],...
                  'dxh','dyh','RR_HARP','CC_HARP','mask_harp');

        catch exception
            fprintf('  Could not process HARP data: %s\n', exception.message);
        end

    end
end


%% SinMod analysis
if RUN_SinMod
    for view={'base','mid','apex'}

        % Debug
        fprintf('\n  Processing CSPAMM data using SinMod (%s view)',view{1})

        % Load data
        filename = sprintf('CI_%s.mat',view{1});
        IPath = [input_folder,'3D_experiments/',filename];
        load(IPath);
        M = squeeze(I.M(:,:,1,1,:));
        I = (squeeze(I.I1 - I.I2));
        
        % Re-format the images
        Nfr = size(I,4);
        for i=1:Nfr
            M(:,:,i) = M(:,:,i)';
            I(:,:,1,i) = I(:,:,1,i)';
            I(:,:,2,i) = I(:,:,2,i)';
        end

        % SPAMMM encoding frequency
        ke = [ke_spamm ke_spamm];

        % Image size
        Isz = size(I);

        % SinMod displacements
        try

            options = struct(...
                'Mask',              M,...
                'EncFreq',           ke.*pxszC,...
                'FOV',               Isz(1:2),...
                'PixelSize',         [1 1],...
                'SearchWindow',      [2,2],...
                'Frames',            1:Nfr,...
                'show',              true,...
                'theta',             deg2rad([0 90]),...
                'UnwrapPhase',       false,...
                'Seed',              'auto',...
                'Connectivity',      8,...
                'CheckQuality',      false,...
                'QualityPower',      8,...
                'QualityFilterSize', 15,...
                'Window',            false,...
                'Frame2Frame',       true,...
                'KspaceFilter',      KspaceFilter,...
                'BTW_cutoff',        BTW_cutoff,...
                'BTW_order',         BTW_order,...
                'KspaceFollowing',   KspaceFollowing);
            [us] = get_SinMod_motion(I, options);
            dxs = squeeze(us(:,:,1,:));
            dys = squeeze(us(:,:,2,:));

            % Prepare displacements
            displacements = permute(cat(4,dxs,dys),[1 2 4 3]);

            % Temporal fitting
            tfit_args.Mask = M(:,:,1);
            [dxs,dys] = TemporalFitting(displacements,tfit_args);

            % SinMod Strain
            [X, Y] = meshgrid(1:size(dxs,2), 1:size(dys,1));
            options = struct(...
                'X', X,...
                'Y', Y,...
                'mask',M(:,:,1),...
                'times',1:Nfr,...
                'dx', dxs,...
                'dy', dys,...
                'Origin', [],...
                'Orientation', []);
            st = mypixelstrain(options);
            RR_SinMod = NaN([Isz(1) Isz(2) Nfr]);
            CC_SinMod = NaN([Isz(1) Isz(2) Nfr]);
            RR_SinMod(repmat(st.maskimage(:,:,1),[1 1 Nfr])) = st.RR(:);
            CC_SinMod(repmat(st.maskimage(:,:,1),[1 1 Nfr])) = st.CC(:);

%             figure(1)
%             imagesc(CC_SinMod(:,:,8),'AlphaData',st.maskimage)
%             colorbar
%             colormap jet
%             pause

            % Write displacement and strain
            mask_sinmod = st.maskimage(:,:,1);
            save([outputs{3},filename],'dxs','dys','RR_SinMod','CC_SinMod','mask_sinmod');

        catch exception
            fprintf('  Could not process SinMod data: %s\n', exception.message);
        end                

    end
end

%% DENSE analysis
if RUN_DENSE
    for view={'base','mid','apex'}

        % Debug
        fprintf('\n  Processing DENSE data using DENSEanalysis (%s view)',view{1})

        % Load data
        filename = sprintf('DI_%s.mat',view{1});
        IPath = [input_folder,'3D_experiments/',filename];
        load(IPath);
        M = squeeze(I.M(:,:,1,1,:));
        I = squeeze(I.I1 - I.I2);

        % Re-format the images
        Nfr = size(I,4);
        for i=1:Nfr
            M(:,:,i) = M(:,:,i)';
            I(:,:,1,i) = I(:,:,1,i)';
            I(:,:,2,i) = I(:,:,2,i)';
        end

        % Image size
        Isz = size(I);

        % Image filtering
        h = ButterworthFilter(Isz(1:2),Isz(1:2)/2,40,5);
%         figure(1)
%         subplot 131
%         imagesc(abs(itok(I(:,:,1,10))))
%         subplot 132
%         imagesc(h.*abs(itok(I(:,:,1,10))))
%         subplot 133
%         imagesc(h)
%         pause
        I = ktoi(h.*itok(I));

        % Displacement
        u = angle(I);

        % Displacement
        try

            args = struct(...
                'PixelSpacing',         pxszD,...
                'EncFreq',              [ke_dense, ke_dense, nan],...
                'Mask',                 M,...
                'FramesForAnalysis',    [1 Nfr],...
                'ResampleMethod',       'gridfit',...
                'SpatialSmoothing',     0.8,...
                'SeedFrame',            1,...
                'TemporalOrder',        7,...
                'Seed',                 'auto',...
                'OptionsPanel',         false,...
                'UnwrapConnectivity',   8);
            try
                [dxd, dyd] = GetDenseDisplacement(u, args);
            catch
                args.OptionsPanel = true;
                [dxd, dyd] = GetDenseDisplacement(u, args);
            end

            % Prepare displacements
            displacements = permute(cat(4,dxd,dyd),[1 2 4 3]); 

            % Temporal fitting
            tfit_args.Mask = M(:,:,1);
            [dxd,dyd] = TemporalFitting(displacements,tfit_args);

            % Strain
            [X, Y] = meshgrid(1:size(dxd,2), 1:size(dyd,1));
            options = struct(...
                'X', X,...
                'Y', Y,...
                'mask',M(:,:,1),...
                'times',1:Nfr,...
                'dx', dxd,...
                'dy', dyd,...
                'Origin', [],...
                'Orientation', []);
            st = mypixelstrain(options);
            RR_DENSE = NaN([Isz(1) Isz(2) Nfr]);
            CC_DENSE = NaN([Isz(1) Isz(2) Nfr]);
            RR_DENSE(repmat(st.maskimage,[1 1 Nfr])) = st.RR(:);
            CC_DENSE(repmat(st.maskimage,[1 1 Nfr])) = st.CC(:);

%             filename_exact = sprintf('RDI_%s.mat',view{1});
%             load([outputs{1},filename_exact])
%             figure(1)
%             subplot 221
%             imagesc(st.maskimage,'AlphaData',st.maskimage); hold on;
%             q = quiver(X,Y,dxe(:,:,6),dye(:,:,6)); hold off
%             q.AutoScale = 'off';
%             set(gca,'YDir','normal')
%             subplot 222
%             imagesc(st.maskimage,'AlphaData',st.maskimage); hold on;
%             q = quiver(X,Y,dxd(:,:,6),dyd(:,:,6)); hold off
%             q.AutoScale = 'off';
%             set(gca,'YDir','normal')
%             subplot 223
%             imagesc(CC_EXACT(:,:,8)); colorbar
%             set(gca,'YDir','normal')
%             subplot 224
%             imagesc(CC_DENSE(:,:,8)); colorbar
%             set(gca,'YDir','normal')
%             pause

%             figure(1)
%             subplot 121
%             imagesc(CC_DENSE(:,:,8),'AlphaData',st.maskimage)
%             colorbar; colormap(flipud(jet))
%             subplot 122
%             imagesc(RR_DENSE(:,:,8),'AlphaData',st.maskimage)
%             colorbar; colormap(flipud(jet))
%             pause

            % Write displacement and strain
            mask_dense = st.maskimage(:,:,1);
            save([outputs{4},filename],'dxd','dyd','RR_DENSE','CC_DENSE','mask_dense');

        catch exception
            fprintf('  Could not process DENSE data: %s\n', exception.message);
        end                    

    end
end

%% STRAIN CURVES
% Strain data
tmp = struct('base',[],'mid',[],'apex',[]);
strain = struct('ExactD',tmp,'ExactC',tmp,'HARP',tmp,'SinMod',tmp,'DENSE',tmp);

% Arguments
Nseg = 6;
api = struct(...
    'CC',               [],...
    'RR',               [],...
    'Mask',             [],...
    'Nseg',             Nseg,...
    'ClockWise',        false,...
    'Frames',           1:Nfr);

if STRAIN_CURVES
    for view={'base','mid','apex'}

        % Load data
        filename_1 = sprintf('DI_%s.mat',view{1});
        filename_2 = sprintf('RDI_%s.mat',view{1});
        filename_3 = sprintf('CI_%s.mat',view{1});
        filename_4 = sprintf('RCI_%s.mat',view{1});
        load([outputs{2},filename_3],'RR_HARP','CC_HARP','mask_harp');
        load([outputs{3},filename_3],'RR_SinMod','CC_SinMod','mask_sinmod');
        load([outputs{4},filename_1],'RR_DENSE','CC_DENSE','mask_dense'); 
        
        % Change number of segments
        if strcmp(view{1},'apex')
            api.Nseg = 4;
        end
        api.Nseg
        
        % Exact strains
        load([outputs{1},filename_2],'RR_EXACT','CC_EXACT','mask_exact');
        api.CC = CC_EXACT;
        api.RR = RR_EXACT;
        api.Mask = mask_exact;
        strain.ExactD.(view{1}) = getStrainBySegments(api);

        load([outputs{1},filename_4],'RR_EXACT','CC_EXACT','mask_exact');
        api.CC = CC_EXACT;
        api.RR = RR_EXACT;
        api.Mask = mask_exact;
        strain.ExactC.(view{1}) = getStrainBySegments(api);
        
        % HARP strain
        api.CC = CC_HARP;
        api.RR = RR_HARP;
        api.Mask = mask_harp;
        strain.HARP.(view{1}) = getStrainBySegments(api);        

        % SinMod strain
        api.CC = CC_SinMod;
        api.RR = RR_SinMod;
        api.Mask = mask_sinmod;
        strain.SinMod.(view{1}) = getStrainBySegments(api);

        % DENSE strain
        api.CC = CC_DENSE;
        api.RR = RR_DENSE;
        api.Mask = mask_dense;
        strain.DENSE.(view{1}) = getStrainBySegments(api);        
        
    end

end

%% PLOTS

% Plot settings
api = struct(...
    'AxesFontSize',  15,...
    'AxesLineWidth', 2,...
    'LegendFontSize', 8,...
    'Axis', [0 21 -30 30],...
    'XLabel', false,...
    'YLabel', true,...
    'XLabelStr', 'Cardiac phase',...
    'YLabelStr', 'Strain (\%)',...
    'YAxisTickValues', -30:10:30);
plot_line_width = 1.5;
plot_marker_size = 3.0;

% Y labels
labels = struct('base','Basal','mid','Mid','apex','Apical');

% Colors
co = [[0 0 0];
      [1 1 1]*0.6];

% Plot strain curves
figure(1)
t = tiledlayout(3,3);
t.TileSpacing = 'compact';
t.Padding = 'compact';

for view={'base','mid','apex'}

%     sgtitle(view{1})

    nexttile;
    plot(1:Nfr,100*mean(strain.HARP.(view{1}).segments_RR),'-^','Color',co(1,:),'MarkerSize',plot_marker_size,'LineWidth',plot_line_width)
    hold on
    plot(1:Nfr,100*mean(strain.ExactC.(view{1}).segments_RR),'-^','Color',co(2,:),'MarkerSize',plot_marker_size,'LineWidth',plot_line_width)
    hold on
    plot(1:Nfr,100*mean(strain.HARP.(view{1}).segments_CC),'-s','Color',co(1,:),'MarkerSize',plot_marker_size,'LineWidth',plot_line_width)
    hold on
    plot(1:Nfr,100*mean(strain.ExactC.(view{1}).segments_CC),'-s','Color',co(2,:),'MarkerSize',plot_marker_size,'LineWidth',plot_line_width)
    hold off
    if strcmp(view{1},'apex')
        api.XLabel = true;
    end
    api.YLabelStr = sprintf('%s strain (\\%%)',labels.(view{1}));
    api.YLabel = true;
    nice_plot_3(api);
    l = legend('RR HARP', 'RR Exact', 'CC HARP', 'CC Exact');
    l.Location = 'southeast';
    l.LineWidth = 1.5;
    l.FontSize = api.LegendFontSize;
    l.Interpreter = 'latex';
    
    nexttile;
    plot(1:Nfr,mean(100*strain.SinMod.(view{1}).segments_RR),'-^','Color',co(1,:),'MarkerSize',plot_marker_size,'LineWidth',plot_line_width)
    hold on
    plot(1:Nfr,mean(100*strain.ExactC.(view{1}).segments_RR),'-^','Color',co(2,:),'MarkerSize',plot_marker_size,'LineWidth',plot_line_width)
    hold on
    plot(1:Nfr,mean(100*strain.SinMod.(view{1}).segments_CC),'-s','Color',co(1,:),'MarkerSize',plot_marker_size,'LineWidth',plot_line_width)
    hold on
    plot(1:Nfr,mean(100*strain.ExactC.(view{1}).segments_CC),'-s','Color',co(2,:),'MarkerSize',plot_marker_size,'LineWidth',plot_line_width)
    hold off
    if strcmp(view{1},'apex')
        api.XLabel = true;
    end
    api.YLabel = false;
    nice_plot_3(api);
    l = legend('RR SinMod', 'RR Exact', 'CC SinMod', 'CC Exact');
    l.Location = 'southeast';
    l.LineWidth = 1.5;
    l.FontSize = api.LegendFontSize;
    l.Interpreter = 'latex';

    nexttile;
    plot(1:Nfr,100*mean(strain.DENSE.(view{1}).segments_RR),'-^','Color',co(1,:),'MarkerSize',plot_marker_size,'LineWidth',plot_line_width)
    hold on
    plot(1:Nfr,100*mean(strain.ExactD.(view{1}).segments_RR),'-^','Color',co(2,:),'MarkerSize',plot_marker_size,'LineWidth',plot_line_width)
    hold on
    plot(1:Nfr,100*mean(strain.DENSE.(view{1}).segments_CC),'-s','Color',co(1,:),'MarkerSize',plot_marker_size,'LineWidth',plot_line_width)
    hold on
    plot(1:Nfr,100*mean(strain.ExactD.(view{1}).segments_CC),'-s','Color',co(2,:),'MarkerSize',plot_marker_size,'LineWidth',plot_line_width)
    hold off
    if strcmp(view{1},'apex')
        api.XLabel = true;
    end
    api.YLabel = false;
    nice_plot_3(api);
    l = legend('RR DENSE', 'RR Exact', 'CC DENSE', 'CC Exact');
    l.Location = 'southeast';
    l.LineWidth = 1.5;
    l.FontSize = api.LegendFontSize;
    l.Interpreter = 'latex';
    
    % Errors
    err_CC = 100*sum(abs(mean(strain.ExactD.(view{1}).segments_CC) - ...
                     mean(strain.DENSE.(view{1}).segments_CC))/...
                     (Nfr*max(abs(mean(strain.ExactD.(view{1}).segments_CC)))));
    err_RR = 100*sum(abs(mean(strain.ExactD.(view{1}).segments_RR) - ...
                     mean(strain.DENSE.(view{1}).segments_RR))/...
                     (Nfr*max(abs(mean(strain.ExactD.(view{1}).segments_RR)))));
    fprintf('\n  [%s] DENSE errors in the E_CC and E_RR components: (%.1f,%.1f)',view{1},err_CC,err_RR)
    err_CC = 100*sum(abs(mean(strain.ExactC.(view{1}).segments_CC) - ...
                     mean(strain.HARP.(view{1}).segments_CC))/...
                     (Nfr*max(abs(mean(strain.ExactC.(view{1}).segments_CC)))));
    err_RR = 100*sum(abs(mean(strain.ExactC.(view{1}).segments_RR) - ...
                     mean(strain.HARP.(view{1}).segments_RR))/...
                     (Nfr*max(abs(mean(strain.ExactC.(view{1}).segments_RR)))));
    fprintf('\n  [%s] HARP errors in the E_CC and E_RR components: (%.1f,%.1f)',view{1},err_CC,err_RR)
    err_CC = 100*sum(abs(mean(strain.ExactC.(view{1}).segments_CC) - ...
                     mean(strain.SinMod.(view{1}).segments_CC))/...
                     (Nfr*max(abs(mean(strain.ExactC.(view{1}).segments_CC)))));
    err_RR = 100*sum(abs(mean(strain.ExactC.(view{1}).segments_RR) - ...
                     mean(strain.SinMod.(view{1}).segments_RR))/...
                     (Nfr*max(abs(mean(strain.ExactC.(view{1}).segments_RR)))));
    fprintf('\n  [%s] SinMod errors in the E_CC and E_RR components: (%.1f,%.1f)\n',view{1},err_CC,err_RR)
    
end

pause
print('-depsc','-r600','regional_strain')