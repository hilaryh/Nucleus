function newmap = bluewhitered(m)
%BLUEWHITERED   Blue, white, and red color map.
%   BLUEWHITERED(M) returns an M-by-3 matrix containing a blue to white
%   to red colormap, with white corresponding to the CAXIS value closest
%   to zero.  This colormap is most useful for images and surface plots
%   with positive and negative values.  BLUEWHITERED, by itself, is the
%   same length as the current colormap.
%
%   Examples:
%   ------------------------------
%   figure
%   imagesc(peaks(250));
%   colormap(bluewhitered(256)), colorbar
%
%   figure
%   imagesc(peaks(250), [0 8])
%   colormap(bluewhitered), colorbar
%
%   figure
%   imagesc(peaks(250), [-6 0])
%   colormap(bluewhitered), colorbar
%
%   figure
%   surf(peaks)
%   colormap(bluewhitered)
%   axis tight
%
%   See also HSV, HOT, COOL, BONE, COPPER, PINK, FLAG, 
%   COLORMAP, RGBPLOT.


if nargin < 1
   m = size(get(gcf,'colormap'),1);
end

%blue white red
% bottom = [68 114 196]/255;%[0 0 0.5];
% botmiddle = [160 184 226]/255;%[0 0.5 1];
% middle = [1 1 1];
% topmiddle = [250,132,138]/255;%[1 0 0];
% top = [249 0 1]/255;%[0.5 0 0];

% pink yellow green
% bottom = [157 191 0]/255;%[0 0 0.5];
% botmiddle = [205 202 35]/255;%[0 0.5 1];
% middle = [253 213 59]/255;
% topmiddle = [240,153,86]/255;%[1 0 0];
% top = [226 94 106]/255;%[0.5 0 0];
% 
% 
% % tangerine offwhite pear
% bottom = [118 158 52]/255;%[0 0 0.5];
% botmiddle = [186 206 156]/255;%[0 0.5 1];
% middle = [252 253 254]/255;
% topmiddle = [248,200,168]/255;%[1 0 0];
% top = [242 147 77]/255;%[0.5 0 0];


% % cool
% bottom = [64 255 255]/255;%[0 0 0.5];
% botmiddle = [97 196 255]/255;%[0 0.5 1];
% middle = [122 137 255]/255;
% topmiddle = [181 86 246]/255;%[1 0 0];
% top = [249 28 255]/255;%[0.5 0 0];

% blue white pink
bottom = [91 129 173]/255;%[0 0 0.5];
botmiddle = [171 191 213]/255;%[0 0.5 1];
middle = [252 253 254]/255;
topmiddle = [245 184 177]/255;%[1 0 0];
top = [240 113 100]/255;%[0.5 0 0];

% Find middle
lims = get(gca, 'CLim');

% Find ratio of negative to positive
if (lims(1) < 0) & (lims(2) > 0)
    % It has both negative and positive
    % Find ratio of negative to positive
    ratio = abs(lims(1)) / (abs(lims(1)) + lims(2));
    neglen = round(m*ratio);
    poslen = m - neglen;
    
    % Just negative
    new = [bottom; botmiddle; middle];
    len = length(new);
    oldsteps = linspace(0, 1, len);
    newsteps = linspace(0, 1, neglen);
    newmap1 = zeros(neglen, 3);
    
    for i=1:3
        % Interpolate over RGB spaces of colormap
        newmap1(:,i) = min(max(interp1(oldsteps, new(:,i), newsteps)', 0), 1);
    end
    
    % Just positive
    new = [middle; topmiddle; top];
    len = length(new);
    oldsteps = linspace(0, 1, len);
    newsteps = linspace(0, 1, poslen);
    newmap = zeros(poslen, 3);
    
    for i=1:3
        % Interpolate over RGB spaces of colormap
        newmap(:,i) = min(max(interp1(oldsteps, new(:,i), newsteps)', 0), 1);
    end
    
    % And put 'em together
    newmap = [newmap1; newmap];
    
elseif lims(1) >= 0
    % Just positive
    new = [middle; topmiddle; top];
    len = length(new);
    oldsteps = linspace(0, 1, len);
    newsteps = linspace(0, 1, m);
    newmap = zeros(m, 3);
    
    for i=1:3
        % Interpolate over RGB spaces of colormap
        newmap(:,i) = min(max(interp1(oldsteps, new(:,i), newsteps)', 0), 1);
    end
    
else
    % Just negative
    new = [bottom; botmiddle; middle];
    len = length(new);
    oldsteps = linspace(0, 1, len);
    newsteps = linspace(0, 1, m);
    newmap = zeros(m, 3);
    
    for i=1:3
        % Interpolate over RGB spaces of colormap
        newmap(:,i) = min(max(interp1(oldsteps, new(:,i), newsteps)', 0), 1);
    end
    
end
% 
% m = 64;
% new = [bottom; botmiddle; middle; topmiddle; top];
% % x = 1:m;
% 
% oldsteps = linspace(0, 1, 5);
% newsteps = linspace(0, 1, m);
% newmap = zeros(m, 3);
% 
% for i=1:3
%     % Interpolate over RGB spaces of colormap
%     newmap(:,i) = min(max(interp1(oldsteps, new(:,i), newsteps)', 0), 1);
% end
% 
% % set(gcf, 'colormap', newmap), colorbar