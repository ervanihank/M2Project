% M1 Project Stimulus
%Erva September/21

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
[center_x, center_y] = RectCenter(rect);

%Define gabor_text Parameters
gabor_text_dim = 30;
% Sigma of Gaussian//sigma is the standard deviation of the Gaussian
% function, decide the roundness of the gabor_text. Large numbers = more rounded.
sigma = gabor_text_dim / 6 ;
contrast = 0.5 ;
aspectRatio = 1; %specifies the ellipticity of the support of the gabor_text function. For γ = 1, the support is circular. 
%For γ < 1 the support is elongated in orientation_gabors of the parallel stripes of the function. 

% Spatial Frequency
num_cycles = 5;
freq = num_cycles / gabor_text_dim;

backgroundgabor_text = [0.5 0.5 0.5 0.5];
disable_norm = 1; %‘disable_norm’ Optional, defaults to 0. If set to a value of 1, the special multiplicative normalization term normf = 1/(sqrt(2*pi) * sc) will not be applied to the computed gabor_text.
pre_contrast_multiplier = 0.5; %‘contrastPreMultiplicator’ Optional, defaults to 1. This value is multiplied as a scaling factor to the requested contrast value

% Build a procedural gabor_text texture
gabor_text = CreateProceduralGabor(window, gabor_text_dim, gabor_text_dim,[],backgroundgabor_text , disable_norm, pre_contrast_multiplier);

%Position of gabor_texts
dim =8; 
[x, y] = meshgrid(-dim:dim, -dim:dim); %create a matrix!!!!!

%distance of each gabor_text from the center of the array
dist = sqrt(x.^2 + y.^2);

%Inner annulus
inner_dist = 1;
x(dist <= inner_dist) = nan;
y(dist <= inner_dist) = nan;

%Outer annulus
outer_dist = 4; 
x(dist >= outer_dist) = nan;
y(dist >= outer_dist) = nan;

%Select only the finite values
x = x(isfinite(x));
y = y(isfinite(y));

% Center the annulus coordinates in the centre of the screen
x_pos = x .* gabor_text_dim + center_x;
y_pos = y .* gabor_text_dim + center_y;

% Count how many gabor_texts there are
n_gabors = numel(x_pos);

% Make the destination rectangles for all the gabor_texts in the array
base_rect = [0 0 gabor_text_dim gabor_text_dim];
all_rects = nan(4, n_gabors);
for i = 1:n_gabors
    all_rects(:, i) = CenterRectOnPointd(base_rect, x_pos(i), y_pos(i));
end

%percentages and directions for n_trials
n_trials= 10;
percentages= randi([60,90],1,n_trials);
directions= [0,180];
directions= repmat(directions,1,n_trials/2);
directions= Shuffle(directions);

%Positions
numgenerator=1:n_gabors;
pos_gabors=Shuffle(numgenerator); 
pos_each_gabor=all_rects(1:4,pos_gabors); 

% Drift speed 
deg_per_sec = 360 * 4;  
drift_speed =  deg_per_sec * ifi;

respMat = nan(4, n_trials);

for trial=1:n_trials
%Percentages 
perc_signal_gabor=percentages(trial);
direction=directions(trial);
n_signal_gabor=round((n_gabors *perc_signal_gabor)/100);
n_noise_gabor=n_gabors-n_signal_gabor;

%Directions of Gabors
glob_direc_signal= repmat(direction,n_signal_gabor,1); %global direction=0 (leftward)
glob_direc_noise= randi(360,1,n_noise_gabor);

all_gabors_direc=[glob_direc_signal',glob_direc_noise];
all_gabors_direc=Shuffle(all_gabors_direc);

%Orientations of Gabors
orientation_gabors=rand(1, n_gabors) .* 180;

%Speed of each Gabor
speed_each_gabor = cosd(orientation_gabors-all_gabors_direc) .* drift_speed;

%Properties matrix with the same phase for each gabor_text
phase=speed_each_gabor;
propertiesMat = repmat([NaN, freq, sigma, contrast, aspectRatio, 0, 0, 0],n_gabors, 1);
propertiesMat(:, 1) = phase'; 

% Perform initial flip/start of animation
vbl = Screen('Flip', window);

% Interstimulus interval time in seconds and frames
isiTimeSecs = 0.3;
isiTimeFrames = round(isiTimeSecs / ifi);


% Numer of frames to wait before re-drawing
waitframes = 1 ; 

%Keyboard Info
escapeKey = KbName('ESCAPE');
leftKey = KbName('LeftArrow');
rightKey = KbName('RightArrow');

for frame = 1:isiTimeFrames - 1
   
Screen('DrawDots', window, [center_x; center_y], 5, black, [], 2); %fixation point
vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
end

tStart = GetSecs; %to count the RT
respToBeMade = true;
    while respToBeMade 
    Screen('DrawTextures', window, gabor_text, [], pos_each_gabor, orientation_gabors' ,[], [], [], [], kPsychDontDoRotation, propertiesMat');
    Screen('DrawDots', window, [center_x; center_y], 5, black, [], 2); %fixation point
    vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi); 
    phase = phase + speed_each_gabor; %The phase of each gabor_text is changing by the same amount on each frame(they each have the same speed in different all_gabors_direc).
    propertiesMat(:, 1) = phase'; 
        [keyIsDown,secs, keyCode] = KbCheck;
        if keyIsDown
             vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi); 
        if keyCode(escapeKey)
            ShowCursor;
            sca;
            return
        elseif keyCode(leftKey)
            if directions(trial)==0
            response = 1;
            elseif directions(trial)==180
                response=0;
            end
            respToBeMade = false;
        elseif keyCode(rightKey)
            if directions(trial)==180
            response = 1;
            elseif directions(trial)==0
            response=0;
            end
            respToBeMade = false;
        end
        end 
        while keyIsDown
            [keyIsDown] = KbCheck;
        end  
        vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
    end
    rt = secs - tStart;   
       
    respMat(1, trial) = perc_signal_gabor;
    respMat(2, trial) = direction;
    respMat(3, trial) = response;
    respMat(4, trial) = rt;
end

sca; 
% close all;