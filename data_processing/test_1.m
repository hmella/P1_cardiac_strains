%% PATHS
% Inputs path
input_folder = '../data_generation/inputs/';
addpath(input_folder)

% functions paths
addpath(genpath('utils/'))
addpath(genpath('/home/hernan/cardiac_motion'))
addpath(genpath('/home/hernan/denseanalysis'))


% Output folder
outputs = {'outputs/noise_free/Exact/',...
           'outputs/noise_free/HARP/',...
           'outputs/noise_free/SinMod/',...
           'outputs/noise_free/DENSE/',};
for i=1:numel(outputs)
   mkdir(outputs{i});
end
   
%% INPUT DATA
% Analysis to be performed
RUN_EXACT   = false;
RUN_HARP    = false;
RUN_SinMod  = false;
RUN_DENSE   = false;
RUN_ERROR   = true;

% Errors to be estimated
HARP_ERROR = true;
SinMod_ERROR = true;
DENSE_ERROR = true;

% Number of cardiac phases
fr  = 1:20;
Nfr = numel(fr);


%% FILTERS SPECS (for image processing)
% Filter specs
KspaceFilter = 'Transmission';
BTW_cutoff = 1;
BTW_order  = [];
KspaceFollowing = false;


%% IMAGING PARAMETERS
% Resolutions
resolutions = [33 39 51 67 101];
FOV = [0.1 0.1];
pxsz = 0.1./resolutions;

% Encoding frequencies
tag_spac = [0.0080, 0.0100, 0.0120, 0.0140, 0.0160]; % [m]
ke_spamm = 2*pi./tag_spac;                           % [rad/m]
ke_dense = 0.12*1000;                                %


%% FIELD INHOMOGENEITY
add_inhomogeneity = true;
phi = @(X,Y) (X+Y)/FOV(1);

% Modified output folders
if add_inhomogeneity
    outputs = {'outputs/noise_free/Exact/',...
               'outputs/noise_free/HARP_in/',...
               'outputs/noise_free/SinMod_in/',...
               'outputs/noise_free/DENSE/',};
    for i=1:numel(outputs)
       mkdir(outputs{i});
    end
end


%% MAIN CODE
% Data to be analyzed
data = 0:99;
bias_EXACT  = [];                                      % Corrupted exact data
bias_SinMod = [];                                      % Corrupted C-SPAMM data
bias_HARP   = [];                                      % Corrupted C-SPAMM data
bias = [bias_SinMod, bias_HARP, bias_EXACT];
data(bias) = [];

% Sizes
nos = numel(tag_spac);      % number of tag spacings
nod = numel(data);          % number of data
nor = numel(resolutions);   % number of resolutions

%% ERRORS
tmp = NaN([nod nos nor Nfr]);
nRMSE = struct(...
    'DENSE',    tmp,...
    'SinMod',   tmp,...
    'HARP',     tmp);
MDE = struct(...
    'DENSE',    tmp,...
    'SinMod',   tmp,...
    'HARP',     tmp);
nRMSE_CC = struct(...
    'DENSE',    tmp,...
    'SinMod',   tmp,...
    'HARP',     tmp);
nRMSE_RR = struct(...
    'DENSE',    tmp,...
    'SinMod',   tmp,...
    'HARP',     tmp);

%% EXACT ANALYSIS
if RUN_EXACT
    for d=ini:fin%1:nod
        for f=1
            for r=1:nor

                % Load data
                filename = sprintf('EI_%03d_%02d_%02d.mat',d-1,f-1,r-1);
                IPath = [input_folder,'noise_free_images/',filename];
                MPath = [input_folder,'masks/',filename];
                [I,M] = P1_read_EXACT(IPath,MPath,0);
                
                % Debug
                fprintf('\n [EXACT] Processing data %d, resolution %d',d,r)

                % Image size
                Isz = size(M);

                % Get displacements in pixels
                ue = 0.01*angle(I)/pxsz(r);
                
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
                st = pixelstrain(options);
                RR_EXACT = NaN([Isz(1) Isz(2) Nfr]);
                CC_EXACT = NaN([Isz(1) Isz(2) Nfr]);
                RR_EXACT(repmat(st.maskimage(:,:,1),[1 1 Nfr])) = st.RR(:);
                CC_EXACT(repmat(st.maskimage(:,:,1),[1 1 Nfr])) = st.CC(:);

                % figure(1)
                % subplot 121;
                % imagesc(CC_EXACT(:,:,8),'AlphaData',st.maskimage);
                % colormap(flipud(jet)); colorbar;
                % subplot 122;
                % imagesc(RR_EXACT(:,:,8),'AlphaData',st.maskimage);
                % colormap(flipud(jet)); colorbar;
                % pause(0.05)

                % Write exact displacement and strain
                mask_exact = st.maskimage(:,:,1);
                save([outputs{1},filename],'dxe','dye','RR_EXACT','CC_EXACT','mask_exact');

            end
        end
    end
end


%% HARP analysis
if RUN_HARP
    for d=ini:fin%1:nod
        for f=1:nos

            % SPAMM encoding frequency
            ke = [ke_spamm(f) ke_spamm(f)];        

            for r=1:nor

                % Load data
                filename = sprintf('CI_%03d_%02d_%02d.mat',d-1,f-1,r-1);
                IPath = [input_folder,'noise_free_images/',filename];
                MPath = [input_folder,'masks/',filename];
                [I,M] = P1_read_CSPAMM(IPath,MPath,0);
                I = permute(I,[2 1 3 4]);
                M = permute(M,[2 1 3 4]);
                
                % Image size
                Isz = size(I);                
                
                % Add field inhomogeneity
                if add_inhomogeneity
                    [X,Y] = meshgrid(pxsz(r)*(1:Isz(2)),pxsz(r)*(1:Isz(1)));
                    X = X - pxsz(r)*(Isz(2)+1)/2;
                    Y = Y - pxsz(r)*(Isz(1)+1)/2;
                    phi_in = phi(X,Y);
                    
                    I = I.*exp(1j*phi_in);
                    
%                     figure(1)
%                     subplot 131
%                     imagesc(real(I(:,:,1,1)))
%                     subplot 132
%                     imagesc(real(I_in(:,:,1,1)))                    
%                     subplot 133
%                     imagesc(wrap(phi_in))
%                     colormap gray
%                     pause
                    
                end

                % Filtered image
                ke_norm = ke.*[pxsz(r) pxsz(r)];
                H = HARPFilter(struct('Image',I,'CentralFreq',ke_norm,'Direction',deg2rad([0 90]),...
                                    'FilterType','Butterworth','Butterworth_cuttoff',ke(1)/100,'Butterworth_order',5));
                If = H.filter(I);
                
%                 figure(1),
%                 subplot 121
%                 imagesc(abs(itok(I(:,:,1,1))))
%                 subplot 122
%                 imagesc(abs(itok(If(:,:,1,1))))
%                 pause
                
                % Debug
                fprintf('\n [HARP] Processing data %d, tag spacing %d, resolution %d',d,f,r)
                
                % ROI
                [X,Y] = meshgrid(1:Isz(2),1:Isz(1));
                m = M(:,:,1);
                ROI = [min(Y(m))-1, max(Y(m))+1,...
                       min(X(m))-1, max(X(m))+1];                
                
                % HARP displacements
                args = struct(...
                        'Mask',             M,...
                        'EncFreq',          [pxsz(r) pxsz(r)].*ke,...
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
                st = pixelstrain(options);
                RR_HARP = NaN([Isz(1) Isz(2) Nfr]);
                CC_HARP = NaN([Isz(1) Isz(2) Nfr]);
                RR_HARP(repmat(st.maskimage(:,:,1),[1 1 Nfr])) = st.RR(:);
                CC_HARP(repmat(st.maskimage(:,:,1),[1 1 Nfr])) = st.CC(:);

%                 filename_exact = sprintf('EI_%03d_%02d_%02d.mat',d-1,0,r-1);
%                 load([outputs{1},filename_exact])
%                 figure(1)
%                 subplot 221
%                 imagesc(st.maskimage,'AlphaData',st.maskimage); hold on;
%                 q = quiver(X,Y,dxe(:,:,8),dye(:,:,8)); hold off
%                 q.AutoScale = 'off';
%                 subplot 222
%                 imagesc(st.maskimage,'AlphaData',st.maskimage); hold on;
%                 q = quiver(X,Y,dxh(:,:,8),dyh(:,:,8)); hold off
%                 q.AutoScale = 'off';
%                 subplot 223
%                 imagesc(CC_EXACT(:,:,8)); colorbar; caxis([-0.35 -0.1])
%                 subplot 224
%                 imagesc(CC_HARP(:,:,8)); colorbar; caxis([-0.35 -0.1])
                
%                 figure(1)
%                 subplot 121
%                 imagesc(CC_HARP(:,:,8)); colorbar
%                 subplot 122
%                 imagesc(RR_HARP(:,:,8)); colorbar
%                 pause(0.1)

                % Write displacement and strain
                mask_harp = st.maskimage(:,:,1);
                save([outputs{2},filename],...
                      'dxh','dyh','RR_HARP','CC_HARP','mask_harp');

            end
        end
    end
end


%% SinMod analysis
if RUN_SinMod
    for d=ini:fin%1:nod
        for f=1:nos

            % SPAMMM encoding frequency
            ke = [ke_spamm(f) ke_spamm(f)];        

            for r=1:nor

                % Load data
                filename = sprintf('CI_%03d_%02d_%02d.mat',d-1,f-1,r-1);
                mask_filename = sprintf('EI_%03d_%02d_%02d.mat',d-1,0,r-1);
                IPath = [input_folder,'noise_free_images/',filename];
                MPath = [input_folder,'masks/',mask_filename];
                [I,M] = P1_read_CSPAMM(IPath,MPath,0);
                I = permute(I,[2 1 3 4]);
                M = permute(M,[2 1 3 4]);

                % Image size
                Isz = size(I);                
                
                % Add field inhomogeneity
                if add_inhomogeneity
                    [X,Y] = meshgrid(pxsz(r)*(1:Isz(2)),pxsz(r)*(1:Isz(1)));
                    X = X - pxsz(r)*(Isz(2)+1)/2;
                    Y = Y - pxsz(r)*(Isz(1)+1)/2;
                    phi_in = phi(X,Y);

                    I = I.*exp(1j*phi_in);
                    
%                     figure(1)
%                     subplot 131
%                     imagesc(real(I(:,:,1,1)))
%                     subplot 132
%                     imagesc(real(I_in(:,:,1,1)))                    
%                     subplot 133
%                     imagesc(wrap(phi_in))
%                     colormap gray
%                     pause
                    
                end

                % Filtered image
                ke_norm = ke.*[pxsz(r) pxsz(r)];
                H = HARPFilter(struct('Image',I,'CentralFreq',ke_norm,'Direction',deg2rad([0 90]),...
                                    'FilterType','Butterworth','Butterworth_cuttoff',ke(1)/100,'Butterworth_order',5));
                If = H.filter(I);
                
                % Debug
                fprintf('\n [SinMod] Processing data %d, tag spacing %d, resolution %d',d,f,r)

                % Image size
                Isz = size(If);

                % SinMod displacements
                options = struct(...
                    'Mask',              M,...
                    'EncFreq',           ke.*[pxsz(r) pxsz(r)],...
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
                st = pixelstrain(options);
                RR_SinMod = NaN([Isz(1) Isz(2) Nfr]);
                CC_SinMod = NaN([Isz(1) Isz(2) Nfr]);
                RR_SinMod(repmat(st.maskimage(:,:,1),[1 1 Nfr])) = st.RR(:);
                CC_SinMod(repmat(st.maskimage(:,:,1),[1 1 Nfr])) = st.CC(:);

%                 filename_exact = sprintf('EI_%03d_%02d_%02d.mat',d-1,0,r-1);
%                 load([outputs{1},filename_exact])
%                 figure(1)
%                 subplot 221
%                 imagesc(st.maskimage,'AlphaData',st.maskimage); hold on;
%                 q = quiver(X,Y,dxe(:,:,8),dye(:,:,8)); hold off
%                 q.AutoScale = 'off';
%                 subplot 222
%                 imagesc(st.maskimage,'AlphaData',st.maskimage); hold on;
%                 q = quiver(X,Y,dxs(:,:,8),dys(:,:,8)); hold off
%                 q.AutoScale = 'off';
%                 subplot 223
%                 imagesc(CC_EXACT(:,:,8)); colorbar; caxis([-0.35 -0.1])
%                 subplot 224
%                 imagesc(CC_SinMod(:,:,8)); colorbar; caxis([-0.35 -0.1])                
                
%                 figure(1)
%                 imagesc(CC_SinMod(:,:,8),'AlphaData',st.maskimage)
%                 colormap jet

                % Write displacement and strain
                mask_sinmod = st.maskimage(:,:,1);
                save([outputs{3},filename],'dxs','dys','RR_SinMod','CC_SinMod','mask_sinmod');

            end
        end
    end
end

%% DENSE analysis
if RUN_DENSE
    for d=ini:fin%1:nod
        for f=1
            for r=1:nor

                % Load data
                filename = sprintf('DI_%03d_%02d_%02d.mat',d-1,f-1,r-1);
                IPath = [input_folder,'noise_free_images/',filename];
                MPath = [input_folder,'masks/',filename];
                [I,M] = P1_read_CSPAMM(IPath,MPath,0);
                I = permute(I,[2 1 3 4]);
                M = permute(M,[2 1 3 4]);                
                
                % Image size
                Isz = size(I,[1 2]);

                % Image filtering
                ke_norm = ke.*[pxsz(r) pxsz(r)];
                H = HARPFilter(struct('Image',I,'CentralFreq',[0 0],'Direction',deg2rad([90 0]),...
                                    'FilterType','Butterworth','Butterworth_cuttoff',10,...
                                    'Butterworth_order',10));
                If = H.filter(I);
                              
                % Displacement
                u = -angle(If);
                
                % Debug
                fprintf('\n [DENSE] Estimating strain from DENSE data %d, resolution %d',d,r)

                % Displacement
                args = struct(...
                    'PixelSpacing',         [1 1],...
                    'EncFreq',              [pxsz(r) pxsz(r) 1].*[ke_dense, ke_dense, nan],...
                    'Mask',                 M,...
                    'FramesForAnalysis',    [1 Nfr],...
                    'ResampleMethod',       'gridfit',...
                    'SpatialSmoothing',     0.5,...
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
                [X, Y] = meshgrid(1:Isz(2), 1:Isz(1));
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

%                 filename_exact = sprintf('EI_%03d_%02d_%02d.mat',d-1,f-1,r-1);
%                 load([outputs{1},filename_exact])
%                 figure(1)
%                 subplot 221
%                 imagesc(st.maskimage,'AlphaData',st.maskimage); hold on;
%                 q = quiver(X,Y,dxe(:,:,8),dye(:,:,8)); hold off
%                 q.AutoScale = 'off';
%                 subplot 222
%                 imagesc(st.maskimage,'AlphaData',st.maskimage); hold on;
%                 q = quiver(X,Y,dxd(:,:,8),dyd(:,:,8)); hold off
%                 q.AutoScale = 'off';
%                 subplot 223
%                 imagesc(CC_EXACT(:,:,8)); colorbar
%                 subplot 224
%                 imagesc(CC_DENSE(:,:,8)); colorbar
%                 pause
                
%                 figure(1)
%                 subplot 121
%                 imagesc(CC_DENSE(:,:,8),'AlphaData',st.maskimage); colorbar
%                 subplot 122
%                 imagesc(RR_DENSE(:,:,8),'AlphaData',st.maskimage); colorbar
%                 pause
                
                % Write displacement and strain
                mask_dense = st.maskimage(:,:,1);
                save([outputs{4},filename],'dxd','dyd','RR_DENSE','CC_DENSE','mask_dense');

            end
        end
    end
end


%% ERROR ANALYSIS
if RUN_ERROR
    for d=ini:fin%1:nod

        % Debug
        fprintf('\n Estimating error metrics in data %d',d-1)
      
        for f=1:nos
            for r=1:nor

                % Load global data
                filename_exact  = sprintf('EI_%03d_%02d_%02d.mat',d-1,0,r-1);
                filename_cspamm = sprintf('CI_%03d_%02d_%02d.mat',d-1,f-1,r-1);
                filename_dense  = sprintf('DI_%03d_%02d_%02d.mat',d-1,0,r-1);
                load([outputs{1},filename_exact]) % EXACT
                
                %% Errors estimation
                % Load data
                if HARP_ERROR
                    load([outputs{2},filename_cspamm])
                end
                if SinMod_ERROR
                    load([outputs{3},filename_cspamm])
                end
                if DENSE_ERROR
                    load([outputs{4},filename_dense])
                end

                
                % Reference mask
                m = and(and(and(mask_exact,mask_dense),mask_harp),mask_sinmod);

                % Error loop
                N = sum(m(:));           
                for l=1:Nfr

                    % References
                    dx_exact = dxe(:,:,l);  % x-displacement
                    dy_exact = dye(:,:,l);  % y-displacement
                    CC_EXACT_ = CC_EXACT(:,:,l);  % CC strain component
                    RR_EXACT_ = RR_EXACT(:,:,l);  % RR strain component
                    me = sqrt(dx_exact(m).^2 + dy_exact(m).^2);   % displacement magnitude
                    
                    % HARP errors
                    if HARP_ERROR
                        dx_harp = dxh(:,:,l);  % x-displacement
                        dy_harp = dyh(:,:,l);  % y-displacement
                        CC_HARP_ = CC_HARP(:,:,l);  % CC strain component
                        RR_HARP_ = RR_HARP(:,:,l);  % RR strain component
                        mh = sqrt(dx_harp(m).^2 + dy_harp(m).^2);   % displacement magnitude
                        angle_HARP = rad2deg(acos(abs(dx_harp(m).*dx_exact(m) + ...
                                    dy_harp(m).*dy_exact(m))./(mh.*me)));  % directional error
                        MDE.HARP(d,f,r,l)  = (1/N)*sum(angle_HARP(:));    % mean directional error
                        nRMSE.HARP(d,f,r,l)  = 1/(max(me)*sqrt(N))*sqrt(sum((dx_harp(m)-dx_exact(m)).^2 + (dy_harp(m)-dy_exact(m)).^2));
                        nRMSE_CC.HARP(d,f,r,l)  = 1/(max(abs(CC_EXACT_(m)))*sqrt(N))*sqrt(sum((CC_HARP_(m) - CC_EXACT_(m)).^2));
                        nRMSE_RR.HARP(d,f,r,l)  = 1/(max(abs(RR_EXACT_(m)))*sqrt(N))*sqrt(sum((RR_HARP_(m) - RR_EXACT_(m)).^2));

                    end

                    % SinMod errors
                    if SinMod_ERROR
                        dx_sinmod = dxs(:,:,l);  % x-displacement
                        dy_sinmod = dys(:,:,l);  % y-displacement
                        CC_SinMod_ = CC_SinMod(:,:,l);  % CC strain component
                        RR_SinMod_ = RR_SinMod(:,:,l);  % RR strain component
                        ms = sqrt(dx_sinmod(m).^2 + dy_sinmod(m).^2);   % displacement magnitude
                        angle_SinMod  = rad2deg(acos(abs(dx_sinmod(m).*dx_exact(m) + ...
                                        dy_sinmod(m).*dy_exact(m))./(ms.*me)));  % directional error
                        MDE.SinMod(d,f,r,l)  = (1/N)*sum(angle_SinMod(:));    % mean directional error
                        nRMSE.SinMod(d,f,r,l)  = 1/(max(me)*sqrt(N))*sqrt(sum((dx_sinmod(m)-dx_exact(m)).^2 + (dy_sinmod(m)-dy_exact(m)).^2));
                        nRMSE_CC.SinMod(d,f,r,l)  = 1/(max(abs(CC_EXACT_(m)))*sqrt(N))*sqrt(sum((CC_SinMod_(m) - CC_EXACT_(m)).^2));
                        nRMSE_RR.SinMod(d,f,r,l)  = 1/(max(abs(RR_EXACT_(m)))*sqrt(N))*sqrt(sum((RR_SinMod_(m) - RR_EXACT_(m)).^2));

                    end

                    % DENSE errors
                    if DENSE_ERROR
                        dx_DENSE = dxd(:,:,l);  % x-displacement
                        dy_DENSE = dyd(:,:,l);  % y-displacement
                        CC_DENSE_ = CC_DENSE(:,:,l);  % CC strain component
                        RR_DENSE_ = RR_DENSE(:,:,l);  % RR strain component
                        mi = sqrt(dx_DENSE(m).^2 + dy_DENSE(m).^2);   % displacement magnitude
                        angle_DENSE  = rad2deg(acos(abs(dx_DENSE(m).*dx_exact(m) + ...
                                        dy_DENSE(m).*dy_exact(m))./(mi.*me)));  % directional error
                        MDE.DENSE(d,f,r,l)  = (1/N)*sum(angle_DENSE(:));    % mean directional error
                        nRMSE.DENSE(d,f,r,l)  = 1/(max(me)*sqrt(N))*sqrt(sum((dx_DENSE(m)-dx_exact(m)).^2 + (dy_DENSE(m)-dy_exact(m)).^2));
                        nRMSE_CC.DENSE(d,f,r,l)  = 1/(max(abs(CC_EXACT_(m)))*sqrt(N))*sqrt(sum((CC_DENSE_(m) - CC_EXACT_(m)).^2));
                        nRMSE_RR.DENSE(d,f,r,l)  = 1/(max(abs(RR_EXACT_(m)))*sqrt(N))*sqrt(sum((RR_DENSE_(m) - RR_EXACT_(m)).^2));
                    end
                    
                end

            end % for loop
        end % for loop
    end % for loop

    %% Errors
    [mean_HARP_mag, std_HARP_mag]     = meanstd(100*nRMSE.HARP,1);
    [mean_SinMod_mag, std_SinMod_mag] = meanstd(100*nRMSE.SinMod,1);
    [mean_DENSE_mag, std_DENSE_mag]   = meanstd(100*nRMSE.DENSE,1);

    [mean_HARP_ang, std_HARP_ang]     = meanstd(MDE.HARP,1);    
    [mean_SinMod_ang, std_SinMod_ang] = meanstd(MDE.SinMod,1);
    [mean_DENSE_ang, std_DENSE_ang]   = meanstd(MDE.DENSE,1);

    [mean_HARP_CC, std_HARP_CC]     = meanstd(100*nRMSE_CC.HARP,1);
    [mean_SinMod_CC, std_SinMod_CC] = meanstd(100*nRMSE_CC.SinMod,1);
    [mean_DENSE_CC, std_DENSE_CC]   = meanstd(100*nRMSE_CC.DENSE,1);

    [mean_HARP_RR, std_HARP_RR]     = meanstd(100*nRMSE_RR.HARP,1);    
    [mean_SinMod_RR, std_SinMod_RR] = meanstd(100*nRMSE_RR.SinMod,1);
    [mean_DENSE_RR, std_DENSE_RR]   = meanstd(100*nRMSE_RR.DENSE,1);

    %% Save workspace
    save('outputs/noise_free/workspace.mat');

end

   
% %% plots
% option=1;
% if option==1
%     load('outputs/noise_free/workspace.mat')
% else
%     load('outputs.bak/noise_free/workspace.mat')
% end
% 
% fr = 8;
% spa = 4;
% res = 1:5;
% 
% figure,
% subplot 221
% errorbar(squeeze(mean_HARP_mag(spa,res,fr)),squeeze(std_HARP_mag(spa,res,fr)),'LineWidth',2); hold on
% errorbar(squeeze(mean_DENSE_mag(spa,res,fr)),squeeze(std_DENSE_mag(spa,res,fr)),'LineWidth',2); hold on
% errorbar(squeeze(mean_SinMod_mag(spa,res,fr)),squeeze(std_SinMod_mag(spa,res,fr)),'LineWidth',2); hold off
% legend('HARP','DENSE','SinMod')
% axis([0 6 0 20])
% xlabel('displacement (in wavelengths)', 'interpreter', 'LaTeX');
% ylabel('nRMSE (\%)', 'interpreter', 'LaTeX')
% ax = gca;
% ax.XAxis.TickLabels = [0.1 0.2 0.3 0.4 0.5];
% ax.XAxis.TickValues = [1 2 3 4 5];
% ax.YAxis.TickValues = [0 5 10 15 20 25];
% 
% subplot 222
% errorbar(squeeze(mean_HARP_ang(spa,res,fr)),squeeze(std_HARP_ang(spa,res,fr)),'LineWidth',2); hold on
% errorbar(squeeze(mean_DENSE_ang(spa,res,fr)),squeeze(std_DENSE_ang(spa,res,fr)),'LineWidth',2); hold on
% errorbar(squeeze(mean_SinMod_ang(spa,res,fr)),squeeze(std_SinMod_ang(spa,res,fr)),'LineWidth',2); hold off
% legend('HARP','DENSE','SinMod')
% axis([0 6 0 6])
% xlabel('displacement (in wavelengths)', 'interpreter', 'LaTeX');
% ylabel('DE ($^o$)', 'interpreter', 'LaTeX')
% ax = gca;
% ax.XAxis.TickLabels = [0.1 0.2 0.3 0.4 0.5];
% ax.XAxis.TickValues = [1 2 3 4 5];
% ax.YAxis.TickValues = [0 2 4 6 8 10];
% 
% subplot 223
% errorbar(squeeze(mean_HARP_CC(spa,res,fr)),squeeze(std_HARP_CC(spa,res,fr)),'LineWidth',2); hold on
% errorbar(squeeze(mean_DENSE_CC(spa,res,fr)),squeeze(std_DENSE_CC(spa,res,fr)),'LineWidth',2); hold on
% errorbar(squeeze(mean_SinMod_CC(spa,res,fr)),squeeze(std_SinMod_CC(spa,res,fr)),'LineWidth',2); hold off
% legend('HARP','DENSE','SinMod')
% axis([0 6 0 20])
% xlabel('displacement (in wavelengths)', 'interpreter', 'LaTeX');
% ylabel('CC nRMSE (\%)', 'interpreter', 'LaTeX')
% ax = gca;
% ax.XAxis.TickLabels = [0.1 0.2 0.3 0.4 0.5];
% ax.XAxis.TickValues = [1 2 3 4 5];
% ax.YAxis.TickValues = [0 5 10 15 20 25];
% 
% subplot 224
% errorbar(squeeze(mean_HARP_RR(spa,res,fr)),squeeze(std_HARP_RR(spa,res,fr)),'LineWidth',2); hold on
% errorbar(squeeze(mean_DENSE_RR(spa,res,fr)),squeeze(std_DENSE_RR(spa,res,fr)),'LineWidth',2); hold on
% errorbar(squeeze(mean_SinMod_RR(spa,res,fr)),squeeze(std_SinMod_RR(spa,res,fr)),'LineWidth',2); hold off
% legend('HARP','DENSE','SinMod')
% axis([0 6 0 100])
% xlabel('displacement (in wavelengths)', 'interpreter', 'LaTeX');
% ylabel('RR nRMSE (\%)', 'interpreter', 'LaTeX')
% ax = gca;
% ax.XAxis.TickLabels = [0.1 0.2 0.3 0.4 0.5];
% ax.XAxis.TickValues = [1 2 3 4 5];
% ax.YAxis.TickValues = [0 20 40 60 80 100];