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
    for i=1:Nfr
        h = figure('Visible','off');
        t = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
        nexttile
        imagesc(abs(itok(IC(:,:,1,i)))); colormap gray
        axis equal off; caxis([0 2.3091e+03])
        nexttile
        imagesc(abs((IC(:,:,1,i)))); colormap gray
        axis equal off;
        nexttile
        imagesc(abs(itok(ID(:,:,1,i)))); colormap gray
        axis equal off; caxis([0 2.0098e+03])
        nexttile
        imagesc(angle((ID(:,:,1,i)))); colormap gray
        axis equal off
        print('-depsc','-r1200',sprintf('images/kspaces_%s_%2d',view{1},i))
    end

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

        % SPAMM encoding frequency
        ke = [ke_spamm ke_spamm];        
        
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

        % Filtered image
        ke_norm = ke.*pxszC;
        H = HARPFilter(struct('Image',I,'CentralFreq',ke_norm,'Direction',deg2rad([0 90]),...
                            'FilterType','Butterworth','Butterworth_cuttoff',ke(1)/30,'Butterworth_order',5));
        If = H.filter(I);        
        
%         figure(1),
%         subplot 121
%         imagesc(abs(itok(I(:,:,1,10))))
%         subplot 122
%         imagesc(abs(itok(If(:,:,1,10))))
%         pause        

        % Image size
        Isz = size(I);

        % ROI
        [X,Y] = meshgrid(1:Isz(2),1:Isz(1));
        m = M(:,:,1);
        ROI = [min(Y(m))-1, max(Y(m))+1,...
               min(X(m))-1, max(X(m))+1];                        

        % HARP displacements
        try
            args = struct(...
                    'Mask',             M,...
                    'EncFreq',          pxszC.*ke,...
                    'FOV',              Isz(1:2),...
                    'PixelSize',        [1 1],....
                    'Frames',           1:Nfr,...
                    'Show',             false,...
                    'TagSpacing',       2*pi/ke(1),...
                    'SeedPoint',        [],...
                    'ROI',              ROI,...
                    'TemporalFitting',  false,...
                    'TemporalFittingOrder', 10);
            harp = HARP_SPHR(If, args);
            dxh = squeeze(harp.RawMotion(:,:,1,:));
            dyh = squeeze(harp.RawMotion(:,:,2,:));

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
%             imagesc(CC_HARP(:,:,8),'AlphaData',st.maskimage); colorbar; caxis([-0.5 0.2])
%             subplot 122
%             imagesc(RR_HARP(:,:,8),'AlphaData',st.maskimage); colorbar
%             colormap jet
%             pause

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

        % SPAMMM encoding frequency
        ke = [ke_spamm ke_spamm];       

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
        
        % Filtered image
        ke_norm = ke.*pxszC;
        H = HARPFilter(struct('Image',I,'CentralFreq',ke_norm,'Direction',deg2rad([0 90]),...
                            'FilterType','Butterworth','Butterworth_cuttoff',ke(1)/30,'Butterworth_order',5));
        If = H.filter(I);          

        % Image size
        Isz = size(I);

        % SinMod displacements
        try

            options = struct(...
                'Mask',              M,...
                'EncFreq',           ke.*pxszC,...
                'FOV',               Isz(1:2),...
                'PixelSize',         [1 1],...
                'SearchWindow',      [0,0],...
                'Frames',            1:Nfr,...
                'Show',              false,...
                'CheckQuality',      true,...
                'QualityPower',      8,...
                'QualityFilterSize', 15,...
                'Filter',            H,...
                'FrameToFrame',      true,...
                'TemporalFitting',   false,...
                'TemporalFittingOrder', 10);
            sinmod = SinMod(If, options);
            dxs = squeeze(sinmod.RawMotion(:,:,1,:));
            dys = squeeze(sinmod.RawMotion(:,:,2,:));

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

        % Image filtering
        ke_norm = ke_dense.*pxszD;
        H = HARPFilter(struct('Image',I,'CentralFreq',[0 0],'Direction',deg2rad([90 0]),...
                            'FilterType','Butterworth','Butterworth_cuttoff',10,...
                            'Butterworth_order',10));
        If = H.filter(I);        
        
        % Re-format the images
        Nfr = size(I,4);
        for i=1:Nfr
            M(:,:,i) = M(:,:,i)';
            I(:,:,1,i) = I(:,:,1,i)';
            I(:,:,2,i) = I(:,:,2,i)';
        end

        % Image size
        Isz = size(I);

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
                'TemporalOrder',        -1,...
                'Seed',                 'auto',...
                'OptionsPanel',         false,...
                'UnwrapConnectivity',   4);
            try
                [dxd, dyd] = GetDenseDisplacement(u, args);
            catch
                args.OptionsPanel = true;
                [dxd, dyd] = GetDenseDisplacement(u, args);
            end

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
plot_marker_size = 5.0;

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
    plot(1:Nfr,100*mean(strain.HARP.(view{1}).segments_RR),'-^','Color',co(1,:),'MarkerSize',plot_marker_size,'LineWidth',plot_line_width,'MarkerFaceColor',co(1,:))
    hold on
    plot(1:Nfr,100*mean(strain.ExactC.(view{1}).segments_RR),'-^','Color',co(2,:),'MarkerSize',plot_marker_size,'LineWidth',plot_line_width,'MarkerFaceColor',co(2,:))
    hold on
    plot(1:Nfr,100*mean(strain.HARP.(view{1}).segments_CC),'-s','Color',co(1,:),'MarkerSize',plot_marker_size,'LineWidth',plot_line_width,'MarkerFaceColor',co(1,:))
    hold on
    plot(1:Nfr,100*mean(strain.ExactC.(view{1}).segments_CC),'-s','Color',co(2,:),'MarkerSize',plot_marker_size,'LineWidth',plot_line_width,'MarkerFaceColor',co(2,:))
    hold off
    if strcmp(view{1},'apex')
        api.XLabel = true;
    end
    api.YLabelStr = sprintf('%s strain (\\%%)',labels.(view{1}));
    api.YLabel = true;
    nice_plot_3(api);
%     l = legend('RR SP-HR', 'RR Ref.', 'CC SP-HR', 'CC Exact');
%     l.Location = 'southeast';
%     l.LineWidth = 1.5;
%     l.FontSize = api.LegendFontSize;
%     l.Interpreter = 'latex';
    
    nexttile;
    plot(1:Nfr,mean(100*strain.SinMod.(view{1}).segments_RR),'-^','Color',co(1,:),'MarkerSize',plot_marker_size,'LineWidth',plot_line_width,'MarkerFaceColor',co(1,:))
    hold on
    plot(1:Nfr,mean(100*strain.ExactC.(view{1}).segments_RR),'-^','Color',co(2,:),'MarkerSize',plot_marker_size,'LineWidth',plot_line_width,'MarkerFaceColor',co(2,:))
    hold on
    plot(1:Nfr,mean(100*strain.SinMod.(view{1}).segments_CC),'-s','Color',co(1,:),'MarkerSize',plot_marker_size,'LineWidth',plot_line_width,'MarkerFaceColor',co(1,:))
    hold on
    plot(1:Nfr,mean(100*strain.ExactC.(view{1}).segments_CC),'-s','Color',co(2,:),'MarkerSize',plot_marker_size,'LineWidth',plot_line_width,'MarkerFaceColor',co(2,:))
    hold off
    if strcmp(view{1},'apex')
        api.XLabel = true;
    end
    api.YLabel = false;
    nice_plot_3(api);
%     l = legend('RR SinMod', 'RR Ref.', 'CC SinMod', 'CC Exact');
%     l.Location = 'southeast';
%     l.LineWidth = 1.5;
%     l.FontSize = api.LegendFontSize;
%     l.Interpreter = 'latex';

    nexttile;
    plot(1:Nfr,100*mean(strain.DENSE.(view{1}).segments_RR),'-^','Color',co(1,:),'MarkerSize',plot_marker_size,'LineWidth',plot_line_width,'MarkerFaceColor',co(1,:))
    hold on
    plot(1:Nfr,100*mean(strain.ExactD.(view{1}).segments_RR),'-^','Color',co(2,:),'MarkerSize',plot_marker_size,'LineWidth',plot_line_width,'MarkerFaceColor',co(2,:))
    hold on
    plot(1:Nfr,100*mean(strain.DENSE.(view{1}).segments_CC),'-s','Color',co(1,:),'MarkerSize',plot_marker_size,'LineWidth',plot_line_width,'MarkerFaceColor',co(1,:))
    hold on
    plot(1:Nfr,100*mean(strain.ExactD.(view{1}).segments_CC),'-s','Color',co(2,:),'MarkerSize',plot_marker_size,'LineWidth',plot_line_width,'MarkerFaceColor',co(2,:))
    hold off
    if strcmp(view{1},'apex')
        api.XLabel = true;
    end
    api.YLabel = false;
    nice_plot_3(api);
%     l = legend('RR DENSE', 'RR Ref.', 'CC DENSE', 'CC Exact');
%     l.Location = 'southeast';
%     l.LineWidth = 1.5;
%     l.FontSize = api.LegendFontSize;
%     l.Interpreter = 'latex';
    
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
set(gca,'Position',[0.680000000000000,0.068383404889673,0.272500000000000,0.274705531703442])
set(gcf,'Position',[1,67,1134,932])
drawnow
print('-depsc','-r600','regional_strain')

% Plot legends
view = 'base';
figure(2)
plot(1:Nfr,100*mean(strain.HARP.(view).segments_RR),'-^','Color',co(1,:),'MarkerSize',3*plot_marker_size,'LineWidth',3*plot_line_width)
hold on
plot(1:Nfr,100*mean(strain.ExactC.(view).segments_RR),'-^','Color',co(2,:),'MarkerSize',3*plot_marker_size,'LineWidth',3*plot_line_width)
hold on
plot(1:Nfr,100*mean(strain.HARP.(view).segments_CC),'-s','Color',co(1,:),'MarkerSize',3*plot_marker_size,'LineWidth',3*plot_line_width)
hold on
plot(1:Nfr,100*mean(strain.ExactC.(view).segments_CC),'-s','Color',co(2,:),'MarkerSize',3*plot_marker_size,'LineWidth',3*plot_line_width)
hold off

l = legend('RR estimated', 'RR reference', 'CC estimated', 'CC reference');
l.Location = 'northoutside';
l.LineWidth = 1.5;
l.FontSize = 3*api.LegendFontSize;
l.Interpreter = 'latex';
l.Orientation = 'horizontal';
drawnow