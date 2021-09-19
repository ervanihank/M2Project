% My first try for Drifting Gabor
%Erva July/21....

close all;
clearvars; 
sca;  
Screen('Preference', 'SkipSyncTests', 1);

% Setting up Psychtoolbox
PsychDefaultSetup(2);

% Seed the random number generator.
rand('seed', sum(100 * clock));

% Screen Number
screenid = max(Screen('Screens'));   
 
% Colour codes
white = [1 1 1];
black = [0 0 0];  
grey = white / 2;  

% Open the screen 
[window, rect] = PsychImaging('OpenWindow', screenid, grey);

% Refresh rate of the monitor
ifi = Screen('GetFlipInterval', window);

%size of the window
[w, h] = RectSize(rect);

% Center of the screen
[centerx, centery] = RectCenter(rect);
 
%Define Gabor Parameters
gabordimension = 60;
% Sigma of Gaussian//sigma is the standard deviation of the Gaussian
% function, decide the roundness of the gabor. Large numbers = more rounded.
sigma = gabordimension / 6 ;
contrast = 0.5 ;
aspectRatio = 1; %specifies the ellipticity of the support of the Gabor function. For γ = 1, the support is circular. For γ < 1 the support is elongated in gaborlinesorientation of the parallel stripes of the function. 

% Spatial Frequency
numCycles = 5;
freq = numCycles / gabordimension;

backgroundgabor = [0.5 0.5 0.5 0.5];
disableNorm = 1; %‘disableNorm’ Optional, defaults to 0. If set to a value of 1, the special multiplicative normalization term normf = 1/(sqrt(2*pi) * sc) will not be applied to the computed gabor.
preContrastMultiplier = 0.5; %‘contrastPreMultiplicator’ Optional, defaults to 1. This value is multiplied as a scaling factor to the requested contrast value

% Build a procedural gabor texture
gabor = CreateProceduralGabor(window, gabordimension, gabordimension,[],backgroundgabor , disableNorm, preContrastMultiplier);

%Position of Gabors
dim =8; 
[x, y] = meshgrid(-dim:dim, -dim:dim); %create a matrix!!!!!

%distance of each gabor from the center of the array
dist = sqrt(x.^2 + y.^2);

% Cut out an inner annulus
innerDist = 2;
x(dist <= innerDist) = nan;
y(dist <= innerDist) = nan;

%Cut out an outer annulus
outerDist = 8; 
x(dist >= outerDist) = nan;
y(dist >= outerDist) = nan;

%Select only the finite values
x = x(isfinite(x));
y = y(isfinite(y));

% Center the annulus coordinates in the centre of the screen
xPos = x .* gabordimension + centerx;
yPos = y .* gabordimension + centery;

% Count how many Gabors there are
nGabors = numel(xPos);

% Make the destination rectangles for all the Gabors in the array
baseRect = [0 0 gabordimension gabordimension];
allRects = nan(4, nGabors);
for i = 1:nGabors
    allRects(:, i) = CenterRectOnPointd(baseRect, xPos(i), yPos(i));
end
 

%Global Motion Gabors 
percentglobalgabor=80;
nGlobalGabors=round((nGabors *percentglobalgabor)/100);
numgenerator=1:nGlobalGabors;
posglobalmotiongabor=Shuffle(numgenerator);
% posglobalmotiongabor=randi(nGabors,nGlobalGabors,1); 
globalgabors=allRects(1:4,posglobalmotiongabor); 

% Drift speed for global motion
degPerSec = 360 * 4;  
globaldriftspeed =  degPerSec * ifi;

globalmotiondirection= 0 ; %this gives the direction of all gabors. up/down direction. 0 is leftward, 90 would be downward
gaborlinesorientation=rand(1, nGlobalGabors) .* 180 - 90;
speedofeachgabors = cosd(gaborlinesorientation-globalmotiondirection) .* globaldriftspeed; %speed of all gabors 

%Properties matrix with the same phase for each gabor
phase=speedofeachgabors;
propertiesMat = repmat([NaN, freq, sigma, contrast, aspectRatio, 0, 0, 0],nGlobalGabors, 1);
propertiesMat(:, 1) = phase'; 

%Noise Gabors
nNoiseGabor=nGabors-nGlobalGabors;
posnoisegabor=[1:nGabors];
posnoisegabor(posglobalmotiongabor)=[];
noisegabors=allRects(1:4,posnoisegabor); 
noiselinesorientation=rand(1,nNoiseGabor) .* 360;
noisegaborspeed=globaldriftspeed;
phasenoise=noisegaborspeed;
noisepropertiesMat = repmat([NaN, freq, sigma, contrast, aspectRatio, 0, 0, 0],nNoiseGabor, 1);
noisepropertiesMat(:, 1) = phasenoise';

% Perform initial flip/start of animation
vbl = Screen('Flip', window);
 
% Numer of frames to wait before re-drawing
waitframes = 1 ; 
   
while ~KbCheck  
    Screen('DrawTextures', window, gabor, [], globalgabors, gaborlinesorientation' ,[], [], [], [], kPsychDontDoRotation, propertiesMat');
    Screen('DrawDots', window, [centerx; centery], 5, black, [], 2); % Draw the fixation point
    vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi); 
    phase = phase + speedofeachgabors; %The phase of each gabor is changing by the same amount on each frame(they each have the same speed in different directions).
    propertiesMat(:, 1) = phase'; 
    Screen('DrawTextures', window, gabor, [], noisegabors, noiselinesorientation' ,[], [], [], [], kPsychDontDoRotation, noisepropertiesMat');
% %     vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
    phasenoise=phasenoise+noisegaborspeed;
    noisepropertiesMat(:, 1) = phasenoise';
end
   
sca; 
close all; 
 

  