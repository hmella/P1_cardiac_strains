close all; clear;

%% PATHS
% Inputs path
input_folder = '../data_generation/inputs/';
addpath(input_folder)

% functions paths
addpath(genpath('utils/'))
addpath(genpath('/home/hernan/cardiac-motion'))

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
RUN_DENSE   = false;
RUN_HARP    = false;
RUN_SinMod  = true;
RUN_ERROR   = true;

% Errors to be estimated
HARP_ERROR = true;
SinMod_ERROR = true;
DENSE_ERROR = true;

% Number of cardiac phases
fr  = 1:18;
Nfr = numel(fr);

%% FILTERS SPECS (for image processing)
% Filter specs
KSpaceFilter = 'Transmission';
BTW_cutoff = 1;
BTW_order  = [];
KSpaceFollowing = false;

%% IMAGING PARAMETERS
% Resolutions
resolutions = [33 40 50 67 100];
FOV = [0.1 0.1];
pxsz = 0.1./resolutions;

% Encoding frequencies
tag_spac = [0.0080, 0.0100, 0.0120, 0.0140, 0.0160]; % [m]
ke_spamm = 2*pi./tag_spac;                           % [rad/m]

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
    for d=1:nod
        for f=1:nos
            for r=1:nor


                % Load data
                filename = sprintf('I_d%02d_f%01d_r%01d.mat',d-1,f-1,r-1);
                IPath = [input_folder,'reference_images/',filename];
                MPath = [input_folder,'masks/',filename];
                [I,M] = P1_read_EXACT(IPath,MPath,1);
                
                % Debug
                fprintf('\n Processing EXACT data (%d/%d)',d-1,nod)
                
                % Image size
                Isz = size(M);            
                
                % Get displacements in pixels
                ue = -0.01*angle(I)/0.001;
                
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

                % figure(1)
                % subplot 121;
                % imagesc(CC_EXACT(:,:,end)); colorbar; %caxis([-0.2 0.2])
                % subplot 122;
                % imagesc(RR_EXACT(:,:,end)); colorbar; %caxis([-0.2 0.2])
                % pause(0.1)
                
                % Write exact displacement and strain
                mask_exact = st.maskimage(:,:,1);
                save([outputs{2},filename],'dxe','dye','RR_EXACT','CC_EXACT','mask_exact');

            end
        end
    end
end


%% HARP analysis
if RUN_HARP
    for d=1:nod
        for f=1:nos

            % SPAMM encoding frequency
            ke = [ke_spamm(f) ke_spamm(f)];        

            for r=1:nor

                % Load data
                filename = sprintf('CI_%03d_%02d_%02d.mat',d-1,f-1,r-1);
                IPath = [input_folder,'noise_free_images/',filename];
                MPath = [input_folder,'masks/',filename];
                [I,M] = P1_read_CSPAMM(IPath,MPath,1);
                
                % Debug
                fprintf('\n Processing data %d, tag spacing %d, resolution %d',d,f,r)

                % Image size
                Isz = size(I);

                % HARP displacements
                args = struct(...
                    'Mask',             M,...
                    'EncFreq',          ke*pxsz(r),...
                    'FOV',              Isz(1:2),...
                    'PixelSize',        [1 1],...
                    'Frames',           fr,...
                    'tol',              1e-2,...
                    'maxiter',          30,...
                    'GradientEval',     5,...
                    'SearchWindow',     [3,3],...
                    'PhaseWindow',      [2,2],...
                    'show',             false,...
                    'ShowConvergence',  false,...
                    'Seed',             'auto',...
                    'theta',            [0 pi/2],...
                    'Connectivity',     8,...
                    'KSpaceFilter',     KSpaceFilter,...
                    'BTW_cutoff',       BTW_cutoff,...
                    'BTW_order',        BTW_order,...
                    'KSpaceFollowing',  KSpaceFollowing);
                [ux_HARP, uy_HARP] = HARPTrackingOsman(I, args);
                dxh = ux_HARP;    % pixels
                dyh = uy_HARP;    % pixels

                % HARP strain
                % TODO: ELIMINAR OPCION ADICIONAL AÑADIDA A mypixelstrain
                [X, Y] = meshgrid(1:size(ux_HARP,2), 1:size(ux_HARP,1));
                options = struct(...
                    'X', X,...
                    'Y', Y,...
                    'mask',M(:,:,1),...
                    'times',1:Nfr,...
                    'dx', ux_HARP,...
                    'dy', uy_HARP,...
                    'Origin', [],...
                    'checknans',  false,...
                    'Orientation', []);
                st = mypixelstrain(options);
                RR_HARP = NaN([Isz(1) Isz(2) Nfr]);
                CC_HARP = NaN([Isz(1) Isz(2) Nfr]);
                RR_HARP(repmat(st.maskimage(:,:,1),[1 1 Nfr])) = st.RR(:);
                CC_HARP(repmat(st.maskimage(:,:,1),[1 1 Nfr])) = st.CC(:);

%                 figure(1)
%                 subplot 121
%                 imagesc(CC_HARP(:,:,8)); colorbar
%                 subplot 122
%                 imagesc(RR_HARP(:,:,8)); colorbar
%                 pause(0.1)

                % Write displacement and strain
                mask_harp = st.maskimage(:,:,1);
                save([outputs{3},filename],...
                      'dxh','dyh','RR_HARP','CC_HARP','mask_harp');

            end
        end
    end
end


%% SinMod analysis
if RUN_SinMod
    for d=1:nod
        for f=1:nos

            % SPAMM encoding frequency
            ke = [ke_spamm(f) ke_spamm(f)];        

            for r=1:nor

                % Load data
                filename = sprintf('CI_%03d_%02d_%02d.mat',d-1,f-1,r-1);
                IPath = [input_folder,'noise_free_images/',filename];
                MPath = [input_folder,'masks/',filename];
                [I,M] = P1_read_CSPAMM(IPath,MPath,1);

                % Debug
                fprintf('\n Processing data %d, tag spacing %d, resolution %d',d,f,r)

                % Image size
                Isz = size(I);

                % SinMod displacements
                options = struct(...
                    'Mask',              M,...
                    'EncFreq',           ke*pxsz(r),...
                    'FOV',               Isz(1:2),...
                    'PixelSize',         [1 1],...
                    'SearchWindow',      [2,2],...
                    'Frames',            1:Nfr,...
                    'show',              false,...
                    'theta',             deg2rad([0 90]),...
                    'UnwrapPhase',       false,...
                    'Seed',              'auto',...
                    'Connectivity',      8,...
                    'CheckQuality',      true,...
                    'QualityPower',      8,...
                    'QualityFilterSize', 15,...
                    'Window',            false,...
                    'Frame2Frame',       true,...
                    'KSpaceFilter',      KSpaceFilter,...
                    'BTW_cutoff',        BTW_cutoff,...
                    'BTW_order',         BTW_order,...
                    'KSpaceFollowing',   KSpaceFollowing);
                [us] = get_SinMod_motion(I, options);
                dxs = squeeze(us(:,:,1,:));
                dys = squeeze(us(:,:,2,:));

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

                figure(1)
                imagesc(CC_SinMod(:,:,8),'AlphaData',st.maskimage)
                colormap jet
                pause

                % Write displacement and strain
                mask_sinmod = st.maskimage(:,:,1);
                save([outputs{4},filename],'dxs','dys','RR_SinMod','CC_SinMod','mask_sinmod');

            end
        end
    end
end

%% DENSE analysis
if RUN_DENSE
    for d=1:nod
        for f=[0]
            for r=1:nor

                % Load data
                filename = sprintf('DI_%03d_%02d_%02d.mat',d-1,f-1,r-1);
                IPath = [input_folder,'noise_free_images/',filename];
                MPath = [input_folder,'masks/',filename];
                [I,M] = P1_read_CSPAMM(IPath,MPath,1);

                % Displacement
                u = -angle(I);
                
                % figure(2)
                % imagesc(mask(:,:,10).*u(:,:,1,10))
                % pause
                
                % Debug
                fprintf('\n Estimating strain from DENSE data %d, resolution %d',i,j)

                % Pixel size
                pxsz = [X(1,2)-X(1,1) Y(2,1)-Y(1,1)];

                % Displacement
                args = struct(...
                    'PixelSpacing',         pxsz,...
                    'EncFreq',              [ke_dense, ke_dense, nan],...
                    'Mask',                 M,...
                    'FramesForAnalysis',    [1 Nfr],...
                    'ResampleMethod',       'gridfit',...
                    'SpatialSmoothing',     0.8,...
                    'SeedFrame',            10,...
                    'TemporalOrder',        10,...
                    'Seed',                 'auto',...
                    'OptionsPanel',         false,...
                    'UnwrapConnectivity',   8);
                try
                    [ux_DENSE, uy_DENSE] = GetDenseDisplacement(u, args);
                catch
                    args.OptionsPanel = true;
                    [ux_DENSE, uy_DENSE] = GetDenseDisplacement(u, args);
                end
                dxd = ux_DENSE;    % millimeters
                dyd = uy_DENSE;    % millimeters

                % Strain
                [X, Y] = meshgrid(1:size(ux_DENSE,2), 1:size(ux_DENSE,1));
                options = struct(...
                    'X', X,...
                    'Y', Y,...
                    'mask',M(:,:,1),...
                    'times',1:Nfr,...
                    'dx', ux_DENSE,...
                    'dy', uy_DENSE,...
                    'Origin', [],...
                    'Orientation', []);
                st = mypixelstrain(options);
                RR_DENSE = NaN([Isz(1) Isz(2) Nfr]);
                CC_DENSE = NaN([Isz(1) Isz(2) Nfr]);
                RR_DENSE(repmat(st.maskimage,[1 1 Nfr])) = st.RR(:);
                CC_DENSE(repmat(st.maskimage,[1 1 Nfr])) = st.CC(:);

            end
        end
    end
end


%% ERROR ANALYSIS
if RUN_ERROR
    for f=[1 2 4]%1:nos         % spacing
        
        % Mean values
        mean_h = zeros([nod Nfr]);
        mean_s = zeros([nod Nfr]); 
        mean_i = zeros([nod Nfr]); 
        mean_e = zeros([nod Nfr]);            

        for d=1:nod

            % Load global data
            filename = sprintf('I_d%02d_f%01d_r%01d.mat',d-1,f-1,0);
            load(['outputs/noise_free/Exact/',filename]);
            load([input_folder,'masks/',filename]);

            % Squeeze mask
            M = squeeze(M(:,:,1,1,:));
            for i=1:size(M,3)
                M(:,:,i) = M(:,:,i)';
            end

            % Debug
            fprintf('\n Estimating error metrics in data %d, tag spacing %d',d-1,f-1)
            
            %% Errors estimation
            % Load data
            if SinMod_ERROR
                load(['outputs/noise_free/SinMod/',filename]);
            end
            if HARP_ERROR
                load(['outputs/noise_free/HARP/',filename]);
            end
            if HARPI_ERROR
                  load([harpi_output,filename]);    
            end

            
            % Reference mask
            m = mask_exact; 
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
                    MDE.HARP(d,f,l)  = (1/N)*sum(angle_HARP(:));    % mean directional error
                    nRMSE.HARP(d,f,l)  = 1/(max(me)*sqrt(N))*sqrt(sum((dx_harp(m)-dx_exact(m)).^2 + (dy_harp(m)-dy_exact(m)).^2));
                    nRMSE_CC.HARP(d,f,l)  = 1/(max(abs(CC_EXACT_(m)))*sqrt(N))*sqrt(sum((CC_HARP_(m) - CC_EXACT_(m)).^2));
                    nRMSE_RR.HARP(d,f,l)  = 1/(max(abs(RR_EXACT_(m)))*sqrt(N))*sqrt(sum((RR_HARP_(m) - RR_EXACT_(m)).^2));

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
                  MDE.SinMod(d,f,l)  = (1/N)*sum(angle_SinMod(:));    % mean directional error
                  nRMSE.SinMod(d,f,l)  = 1/(max(me)*sqrt(N))*sqrt(sum((dx_sinmod(m)-dx_exact(m)).^2 + (dy_sinmod(m)-dy_exact(m)).^2));
                  nRMSE_CC.SinMod(d,f,l)  = 1/(max(abs(CC_EXACT_(m)))*sqrt(N))*sqrt(sum((CC_SinMod_(m) - CC_EXACT_(m)).^2));
                  nRMSE_RR.SinMod(d,f,l)  = 1/(max(abs(RR_EXACT_(m)))*sqrt(N))*sqrt(sum((RR_SinMod_(m) - RR_EXACT_(m)).^2));

              end

              % HARPI errors
              if HARPI_ERROR
                  dx_harpi = dxi(:,:,l);  % x-displacement
                  dy_harpi = dyi(:,:,l);  % y-displacement
                  CC_HARPI_ = CC_HARPI(:,:,l);  % CC strain component
                  RR_HARPI_ = RR_HARPI(:,:,l);  % RR strain component
                  mi = sqrt(dx_harpi(m).^2 + dy_harpi(m).^2);   % displacement magnitude
                  angle_HARPI  = rad2deg(acos(abs(dx_harpi(m).*dx_exact(m) + ...
                                dy_harpi(m).*dy_exact(m))./(mi.*me)));  % directional error
                  MDE.HARPI(d,f,l)  = (1/N)*sum(angle_HARPI(:));    % mean directional error
                  nRMSE.HARPI(d,f,l)  = 1/(max(me)*sqrt(N))*sqrt(sum((dx_harpi(m)-dx_exact(m)).^2 + (dy_harpi(m)-dy_exact(m)).^2));
                  nRMSE_CC.HARPI(d,f,l)  = 1/(max(abs(CC_EXACT_(m)))*sqrt(N))*sqrt(sum((CC_HARPI_(m) - CC_EXACT_(m)).^2));
                  nRMSE_RR.HARPI(d,f,l)  = 1/(max(abs(RR_EXACT_(m)))*sqrt(N))*sqrt(sum((RR_HARPI_(m) - RR_EXACT_(m)).^2));
                  
                  if l==6 && f==2 && d==4
                      r = 61:140;
                      CA = [min(CC_EXACT_(:)) max(CC_EXACT_(:))];

                      figure('visible','off')
                      imagesc(CC_EXACT_(r,r),'alphadata',m(r,r));
                      cb = colorbar;
                      cb.Label.Interpreter = 'latex';
                      caxis(CA)
                      colormap(flipud(jet))
                      axis off square
                      print('-depsc','-r600', sprintf('%01d_CC_EXACT',d));
%                       print('-dpng','-r600', sprintf('%01d_CC_EXACT',d));

                      figure('visible','off')
                      imagesc(CC_HARP_(r,r),'alphadata',m(r,r));
                      cb = colorbar;
                      cb.Label.Interpreter = 'latex';
                      caxis(CA)
                      colormap(flipud(jet))
                      axis off square  
                      print('-depsc','-r600', sprintf('%01d_CC_HARP',d));
%                       print('-dpng','-r600', sprintf('%01d_CC_HARP',d));
                      
                      figure('visible','off')
                      imagesc(CC_SinMod_(r,r),'alphadata',m(r,r));
                      cb = colorbar;
                      cb.Label.Interpreter = 'latex';
                      caxis(CA)
                      colormap(flipud(jet))
                      axis off square
                      print('-depsc','-r600', sprintf('%01d_CC_SinMod',d));
%                       print('-dpng','-r600', sprintf('%01d_CC_SinMod',d));

                      figure('visible','off')
                      imagesc(CC_HARPI_(r,r),'alphadata',m(r,r)); 
                      cb = colorbar;
                      cb.Label.Interpreter = 'latex';
                      caxis(CA)
                      colormap(flipud(jet))
                      axis off square                      
                      print('-depsc','-r600', sprintf('%01d_CC_HARPI',d));
%                       print('-dpng','-r600', sprintf('%01d_CC_HARPI',d));
                  end
                  
              end
                
            end
        end
    end

    %% Errors
    [mean_HARPI_mag, std_HARPI_mag]   = meanstd(100*nRMSE.HARPI,1);
    [mean_SinMod_mag, std_SinMod_mag] = meanstd(100*nRMSE.SinMod,1);
    [mean_HARP_mag, std_HARP_mag]     = meanstd(100*nRMSE.HARP,1);

    [mean_HARPI_ang, std_HARPI_ang]   = meanstd(MDE.HARPI,1);
    [mean_SinMod_ang, std_SinMod_ang] = meanstd(MDE.SinMod,1);
    [mean_HARP_ang, std_HARP_ang]     = meanstd(MDE.HARP,1);    

    [mean_HARPI_CC, std_HARPI_CC]   = meanstd(100*nRMSE_CC.HARPI,1);
    [mean_SinMod_CC, std_SinMod_CC] = meanstd(100*nRMSE_CC.SinMod,1);
    [mean_HARP_CC, std_HARP_CC]     = meanstd(100*nRMSE_CC.HARP,1);

    [mean_HARPI_RR, std_HARPI_RR]   = meanstd(100*nRMSE_RR.HARPI,1);
    [mean_SinMod_RR, std_SinMod_RR] = meanstd(100*nRMSE_RR.SinMod,1);
    [mean_HARP_RR, std_HARP_RR]     = meanstd(100*nRMSE_RR.HARP,1);    

    %% Save workspace
    save([harpi_output,'workspace.mat']);

end


%% plots
spa = 4;
figure,
subplot 221
errorbar(mean_HARP_mag(spa,2:end),std_HARP_mag(spa,2:end),'LineWidth',2); hold on
errorbar(mean_HARPI_mag(spa,2:end),std_HARPI_mag(spa,2:end),'LineWidth',2); hold on
errorbar(mean_SinMod_mag(spa,2:end),std_SinMod_mag(spa,2:end),'LineWidth',2); hold off
legend('HARP','HARPI','SinMod')
axis([0 6 0 20])
xlabel('displacement (in wavelengths)', 'interpreter', 'LaTeX');
ylabel('nRMSE (\%)', 'interpreter', 'LaTeX')
ax = gca;
ax.XAxis.TickLabels = [0.1 0.2 0.3 0.4 0.5];
ax.XAxis.TickValues = [1 2 3 4 5];
ax.YAxis.TickValues = [0 5 10 15 20 25];

subplot 222
errorbar(mean_HARP_ang(spa,2:end),std_HARP_ang(spa,2:end),'LineWidth',2); hold on
errorbar(mean_HARPI_ang(spa,2:end),std_HARPI_ang(spa,2:end),'LineWidth',2); hold on
errorbar(mean_SinMod_ang(spa,2:end),std_SinMod_ang(spa,2:end),'LineWidth',2); hold off
legend('HARP','HARPI','SinMod')
axis([0 6 0 6])
xlabel('displacement (in wavelengths)', 'interpreter', 'LaTeX');
ylabel('DE ($^o$)', 'interpreter', 'LaTeX')
ax = gca;
ax.XAxis.TickLabels = [0.1 0.2 0.3 0.4 0.5];
ax.XAxis.TickValues = [1 2 3 4 5];
ax.YAxis.TickValues = [0 2 4 6 8 10];

subplot 223
errorbar(mean_HARP_CC(spa,2:end),std_HARP_CC(spa,2:end),'LineWidth',2); hold on
errorbar(mean_HARPI_CC(spa,2:end),std_HARPI_CC(spa,2:end),'LineWidth',2); hold on
errorbar(mean_SinMod_CC(spa,2:end),std_SinMod_CC(spa,2:end),'LineWidth',2); hold off
legend('HARP','HARPI','SinMod')
axis([0 6 0 20])
xlabel('displacement (in wavelengths)', 'interpreter', 'LaTeX');
ylabel('CC nRMSE (\%)', 'interpreter', 'LaTeX')
ax = gca;
ax.XAxis.TickLabels = [0.1 0.2 0.3 0.4 0.5];
ax.XAxis.TickValues = [1 2 3 4 5];
ax.YAxis.TickValues = [0 5 10 15 20 25];

subplot 224
errorbar(mean_HARP_RR(spa,2:end),std_HARP_RR(spa,2:end),'LineWidth',2); hold on
errorbar(mean_HARPI_RR(spa,2:end),std_HARPI_RR(spa,2:end),'LineWidth',2); hold on
errorbar(mean_SinMod_RR(spa,2:end),std_SinMod_RR(spa,2:end),'LineWidth',2); hold off
legend('HARP','HARPI','SinMod')
axis([0 6 0 100])
xlabel('displacement (in wavelengths)', 'interpreter', 'LaTeX');
ylabel('RR nRMSE (\%)', 'interpreter', 'LaTeX')
ax = gca;
ax.XAxis.TickLabels = [0.1 0.2 0.3 0.4 0.5];
ax.XAxis.TickValues = [1 2 3 4 5];
ax.YAxis.TickValues = [0 20 40 60 80 100];