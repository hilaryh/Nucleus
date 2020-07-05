% This function computes the power spectrum for the given time-series data
% and applies an appropriate adaptive Butterworth filter.
% The output is the filtered signal.
% 
% Written by Greg Bass on 18/06/15 

function [filt] = filterNoisyLines(mylines,fs)

    % Pre-allocate the output matrix
    filt = zeros(size(mylines));

    % Loop through each noisy line
    parfor ii = 1:size(mylines,1)
        myline = mylines(ii,:);
        
        % For the power spectrum analysis, I only want to use the negative
        % control (approximately the first 25 seconds of the protocol) 
        % 25sec * 1000 ms/sec * 1 increment/0.25 ms = 10^5 increments
        mylinePowSpec = myline;
        
        % Adapted from the tutorial at http://sigpromu.org/brett/elec2400/matlab3.pdf (method 2)
        % Fourier transform
        X = fft(mylinePowSpec);
        
        % Sampling and window parameters
        N = length(mylinePowSpec);
        ws = 2*pi/N;
        wnorm = -pi:ws:pi;
        wnorm = wnorm(1:length(mylinePowSpec));
        w = wnorm*fs + ws*fs/2;
        
        % Plot the power spectrum
        powSpec = abs(fftshift(X));
        %figure, plot(w,powSpec);
        %xlabel('Frequency (Hz)'); ylabel('Magnitude')
        %title('Power spectral density');
        
        % Find the DC component and the next highest (pair of) peaks
        % (The DC component is the tallest peak, which is centered around zero)
        [pks,locs] = findpeaks(powSpec);
        DCloc = locs(pks==max(pks));
        %DCfreq = w(DCloc);
        DCidx = find(pks==max(pks));
        lowCutoff = abs(w(DCloc+1)-w(DCloc-1));  % take the difference of the next (lower) frequencies on either side of the DC component
        highCutoff = abs(w(locs(DCidx+1))-w(locs(DCidx-1)));  % take the difference of the positive and negative high frequency peaks


        % Set the cutoff frequencies for the filter
        % If the cutoffs are below a threshold, override them
        if highCutoff < 8
            f_high = 10; f_low = 2; iter = 10;
        else
            f_high = highCutoff; f_low = lowCutoff; iter = 10;
        end

        % Apply the adaptive Butterworth filter
        m = [repmat(mean(myline(1:500)),[1 1000]), myline, repmat(mean(myline((end-500):end)), [1 1000])];  % padding
        [tempfilt,~,~]=adaptiveB(m,fs,f_low,f_high,iter);  % input data must be a row vector
        filt(ii,:) = tempfilt((1+1000):(end-1000));  % remove padding
    
    end
             
end


