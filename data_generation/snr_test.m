clear; close all; clc

% Add paths
addpath(genpath('inputs/'))

% Load DENSE image
d = 0;
n = 1;

% SNR
for d=1:10
    for r=0:4

        % Load images
        D_filename = sprintf('DI_%03d_%02d_%02d.mat',d,n,r);
        M_filename = sprintf('DI_%03d_%02d_%02d.mat',d,0,r);
        N_filename = sprintf('SI_%03d_%02d_%02d.mat',d,n,r);
        IPath = ['inputs/noisy_images/',D_filename];
        MPath = ['inputs/masks/',M_filename];
        NPath = ['inputs/noisy_images/',N_filename];
        [I,M1] = P1_read_CSPAMM(IPath,MPath,1);
        N = P1_read_CSPAMM(NPath,MPath,1);

        % Image size
        Isz = size(I);

        % Image filtering
        h = ButterworthFilter(Isz(1:2),Isz(1:2)/2,10,10);                
        I = ktoi(h.*itok(I));    

        % SNR estimation    
        noise_std = zeros([1 size(I,4)]);
        signal_mean = zeros([1 size(I,4)]);
        for i=1:size(I,4)
            IM_mag = abs(I(:,:,1,i));
            noise = N(:,:,1,i);
            noise_std(i) = std(noise,0,[1 2]);
            signal_mean(i) = mean(IM_mag(M1(:,:,i)));
        end

        % Plots
        figure(2)
        plot(signal_mean./noise_std,'LineWidth',2); hold on

    end
end
xlabel('Frame')
ylabel('SNR')
hold off
