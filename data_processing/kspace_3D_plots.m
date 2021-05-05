close all


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
    ID = squeeze(I.I1);

    % Re-format the images
    Nfr = size(ID,4);
    for i=1:Nfr
        M(:,:,i) = M(:,:,i)';
        ID(:,:,1,i) = ID(:,:,1,i)';
        ID(:,:,2,i) = ID(:,:,2,i)';
    end

    % Ranges for plotting
    c_max_k = max(max(abs(itok(IC(:,:,1,1)))));
    c_max_i = max(max(abs(IC(:,:,1,1))));
    d_max_k = max(max(abs(itok(ID(:,:,1,1)))));
    d_max_i = max(max(abs(ID(:,:,1,1))));
    
    % Concatenated images
    multi_KC = zeros(resC);
    multi_IC = zeros([108, 108]);
    multi_KD = zeros(resD);
    multi_ID = zeros([54 54]);
    c = 1;
    for i=1:3:Nfr
        multi_IC(:,((c-1)*108 + c):(c*108 + (c-1))) = IC(75:182,75:182,1,i);
        multi_IC(:,end+1) = NaN;
        multi_KC(:,((c-1)*resC(1) + c):(c*resC(1) + (c-1))) = itok(IC(:,:,1,i));
        multi_KC(:,end+1) = NaN;
        multi_ID(:,((c-1)*54 + c):(c*54 + (c-1))) = ID(38:91,38:91,1,i);
        multi_ID(:,end+1) = NaN;
        multi_KD(:,((c-1)*resD(1) + c):(c*resD(1) + (c-1))) = itok(ID(:,:,1,i));
        multi_KD(:,end+1) = NaN;        
        c = c + 1;
    end
    
    figure,
    subplot(2,6,1:6)
    imagesc(abs(multi_KC))
    axis off equal; colormap gray
    caxis([0 0.75*c_max_k])   
    subplot(2,6,7:12)
    imagesc(abs(multi_IC))
    axis off equal; colormap gray
    set(gca,'Position',[0.130000000000000,0.110000000000000,0.775000000000000,0.341162790697674])
    set(gcf,'Position',[1,41,1920,958])
    print('-depsc','-r1600',sprintf('images/SPAMM_kspaces_%s',view{1}))

    figure,
    subplot(2,6,1:6)
    imagesc(abs(multi_KD))
    axis off equal; colormap gray
    caxis([0 0.75*d_max_k])
    subplot(2,6,7:12)
    imagesc(angle(multi_ID))
    axis off equal; colormap gray
    set(gca,'Position',[0.130000000000000,0.110000000000000,0.775000000000000,0.341162790697674])
    set(gcf,'Position',[1,41,1920,958])
    print('-depsc','-r1600',sprintf('images/DENSE_kspaces_%s',view{1}))    
    
end


% Plot complmentary kspaces
for view={'base','mid','apex'}
    % Load SPAMM image
    filename = sprintf('CI_%s.mat',view{1});
    IPath = [input_folder,'3D_experiments/',filename];
    load(IPath);
    M = squeeze(I.M(:,:,1,1,:));
    IC = (squeeze(I.I1 - I.I2));

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

    % Re-format the images
    Nfr = size(ID,4);
    for i=1:Nfr
        M(:,:,i) = M(:,:,i)';
        ID(:,:,1,i) = ID(:,:,1,i)';
        ID(:,:,2,i) = ID(:,:,2,i)';
    end

    % Ranges for plotting
    c_max_k = max(max(abs(itok(IC(:,:,1,1)))));
    c_max_i = max(max(abs(IC(:,:,1,1))));
    d_max_k = max(max(abs(itok(ID(:,:,1,1)))));
    d_max_i = max(max(abs(ID(:,:,1,1))));
    
    % Concatenated images
    multi_KC = zeros(resC);
    multi_IC = zeros([108, 108]);
    multi_KD = zeros(resD);
    multi_ID = zeros([54 54]);
    c = 1;
    for i=1:3:Nfr
        multi_IC(:,((c-1)*108 + c):(c*108 + (c-1))) = IC(75:182,75:182,1,i);
        multi_IC(:,end+1) = NaN;
        multi_KC(:,((c-1)*resC(1) + c):(c*resC(1) + (c-1))) = itok(IC(:,:,1,i));
        multi_KC(:,end+1) = NaN;
        multi_ID(:,((c-1)*54 + c):(c*54 + (c-1))) = ID(38:91,38:91,1,i);
        multi_ID(:,end+1) = NaN;
        multi_KD(:,((c-1)*resD(1) + c):(c*resD(1) + (c-1))) = itok(ID(:,:,1,i));
        multi_KD(:,end+1) = NaN;        
        c = c + 1;
    end
    
    figure,
    subplot(2,6,1:6)
    imagesc(abs(multi_KC))
    axis off equal; colormap gray
    caxis([0 0.75*c_max_k])
    subplot(2,6,7:12)
    imagesc(abs(multi_IC))
    axis off equal; colormap gray
    set(gca,'Position',[0.130000000000000,0.110000000000000,0.775000000000000,0.341162790697674])
    set(gcf,'Position',[1,41,1920,958])
    print('-depsc','-r1600',sprintf('images/CSPAMM_kspaces_%s',view{1}))

    figure,
    subplot(2,6,1:6)
    imagesc(abs(multi_KD))
    axis off equal; colormap gray
    caxis([0 0.75*d_max_k])
    subplot(2,6,7:12)
    imagesc(angle(multi_ID))
    axis off equal; colormap gray
    set(gca,'Position',[0.130000000000000,0.110000000000000,0.775000000000000,0.341162790697674])
    set(gcf,'Position',[1,41,1920,958])
    print('-depsc','-r1600',sprintf('images/CDENSE_kspaces_%s',view{1}))    
    
end


