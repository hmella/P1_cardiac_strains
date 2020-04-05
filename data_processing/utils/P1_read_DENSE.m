function [I,M] = P1_read_DENSE(IPath,MPath,NbPixels)

    if nargin < 3
        NbPixels = 0;
    end

    % Read data
    load(IPath);
    load(MPath);
    
    % Rescale data
    m = I.complex.RescaleSlope;
    b = I.complex.RescaleIntercept;
    Ir = m*double(I.real.Image) + double(b);
    Ii = 1j*(m*double(I.complex.Image) + double(b));
    I = Ir + Ii;
    
    % Squeeze data
    I = squeeze(I(:,:,1,:,:));
    M = squeeze(M(:,:,1,1,:));

    %
    for i=1:Nfr
        M(:,:,i) = M(:,:,i)';
        I(:,:,1,i) = I(:,:,1,i)';
        I(:,:,2,i) = I(:,:,2,i)';
    end

    % Remove outer pixels
    M = removeOuterPixels(M, NbPixels);

    % Image filtering
    h = ButterworthFilter(Isz(1:2),Isz(1:2)/2,10,10);

    for e=1:size(I,4)
        for f=size(I,3)
            I(:,:,f,e) = ktoi(h.*itok(I(:,:,f,e)));
        end
    end

end