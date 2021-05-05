%% PATHS
% Inputs path
input_folder = '../data_generation/inputs/';
addpath(input_folder)

% functions paths
addpath(genpath('utils/'))
addpath(genpath('/home/hernan/cardiac_motion'))
addpath(genpath('/home/hernan/denseanalysis'))


% Output folder
outputs = {'outputs/noisy/Exact/',...
           'outputs/noisy/HARP/',...
           'outputs/noisy/SinMod/',...
           'outputs/noisy/DENSE/',};
for i=1:numel(outputs)
   mkdir(outputs{i});
end
   
%% INPUT DATA
% Analysis to be performed
RUN_EXACT   = false;
RUN_HARP    = true;
RUN_DENSE   = true;

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
% Resolutions
resolutions = [33 39 51 67 101];
FOV = [0.1 0.1];
pxsz = 0.1./resolutions;

% Encoding frequencies
tag_spac = [0.0080, 0.0100, 0.0120, 0.0140, 0.0160]; % [m]
ke_spamm = 2*pi./tag_spac;                           % [rad/m]
ke_dense = 0.12*1000;                                %

%% MAIN CODE
% Data to be analyzed
data = 0:99;
bias_EXACT  = [];                                      % Corrupted exact data
bias_SinMod = [];                                      % Corrupted C-SPAMM data
bias_HARP   = [];                                      % Corrupted C-SPAMM data
bias = [bias_SinMod, bias_HARP, bias_EXACT];
data(bias) = [];

% Sizes
non = 5;                    % number of noise levels
nod = numel(data);          % number of data
nor = numel(resolutions);   % number of resolutions

%% ERRORS
tmp = NaN([nod non nor Nfr]);
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

% Success indicator
tmp = true([nod, non, nor]);
success = struct(...
    'HARP',    tmp,...
    'SinMod',   tmp,...
    'DENSE',     tmp);


%% EXACT ANALYSIS
if RUN_EXACT
    for d=ini:fin%1:nod
        for f=1
            for r=1:nor

                % Load data
                filename = sprintf('EI_%03d_%02d_%02d.mat',d-1,f-1,r-1);
                IPath = [input_folder,'noisy_images/',filename];
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
                st = mypixelstrain(options);
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


%% CSPAMM image
if RUN_HARP
    for d=1
        for f=2

            % SPAMM encoding frequency
            ke = [ke_spamm(f) ke_spamm(f)];        

            for r=4

                % Load data
                filename = sprintf('HI_%03d_%02d_%02d.mat',d-1,f-1,r-1);
                mask_filename = sprintf('CI_%03d_%02d_%02d.mat',d-1,0,r-1);
                IPath = [input_folder,'noisy_images/',filename];
                MPath = [input_folder,'masks/',mask_filename];
                [I,M] = P1_read_CSPAMM(IPath,MPath,0);
                
                % Debug
                fprintf('\n [HARP] Processing data %d, noise level %d, resolution %d',d,f,r)

                % Image size
                Isz = size(I);
                
                % Plot
                figure('Visible','off')
                subplot 121
                imagesc(abs(itok(I(:,:,1,10)))); colormap gray; axis equal off; box on
                subplot 122
                imagesc(abs(I(:,:,1,10))); colormap gray; axis equal off; box on
                print('-depsc','-r600','im_CSPAMM')
                
            end
        end
    end
end


%% DENSE image
if RUN_DENSE
    for d=1
        for f=2
            for r=3

                % Load data
                filename = sprintf('DI_%03d_%02d_%02d.mat',d-1,f-1,r-1);
                mask_filename = sprintf('DI_%03d_%02d_%02d.mat',d-1,0,r-1);
                IPath = [input_folder,'noisy_images/',filename];
                MPath = [input_folder,'masks/',mask_filename];
                [I,M] = P1_read_CSPAMM(IPath,MPath,0);

                % Image size
                Isz = size(I);

                % Image filtering
%                 h = ButterworthFilter(Isz(1:2),Isz(1:2)/2,10,10);                
%                 I = ktoi(h.*itok(I));

                % Displacement
                u = angle(I);
                
                % Plot
                figure('Visible','off')
                subplot 121
                imagesc(abs(itok(I(:,:,1,10)))); colormap gray; axis equal; box on
                set(gca,'YDir','normal');
                subplot 122
                imagesc(angle(I(:,:,1,10))); colormap gray; axis equal off; box on
                set(gca,'YDir','normal');
                print('-depsc','-r600','im_DENSE')
                
            end
        end
    end
end