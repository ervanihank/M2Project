% M1 Project Stimulus
%Erva September/21

close all;
clearvars; 
sca;  
Screen('Preference', 'SkipSyncTests', 1);

commandwindow;

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

monitorPos = get(0,'MonitorPositions');

pixelperdeg=25;
%pixelperdeg=36;
% for the computer downstairs
% ScreenHypotenuse = 50.5;
% pixHypotenuse = sqrt((monitorPos(3)^2)+(monitorPos(4)^2));
% pixelperdeg = round(pixHypotenuse/ScreenHypotenuse);

% Open the screen 
[window, rect] = PsychImaging('OpenWindow', screenid, grey,[0,800;0,800]);

%Screen('BlendFunction',window, 'GL_SRC_ALPHA','GL_ONE_MINUS_SRC_ALPHA')
Screen('BlendFunction',window, 'GL_ONE','GL_ZERO');

% Refresh rate of the monitor
ifi = Screen('GetFlipInterval', window);

%size of the window
[w, h] = RectSize(rect);

% Center of the screen
[center_x, center_y] = RectCenter(rect);

%fixation cross
fixCrossDimPix = 4;
xCoords = [-fixCrossDimPix fixCrossDimPix 0 0];
yCoords = [0 0 -fixCrossDimPix fixCrossDimPix];
allCoords = [xCoords; yCoords];
lineWidthPix = 4;

%% Define gabor_text Parameters

gabor_text_dim = round(pixelperdeg*0.75);

% Sigma of Gaussian//sigma is the standard deviation of the Gaussian
% function, decide the roundness of the gabor_text. Large numbers = more rounded.
sigma = gabor_text_dim / 6;

contrast = 0.25 ;

aspectRatio = 1; %specifies the ellipticity of the support of the gabor_text function. For γ = 1, the support is circular. 
%For γ < 1 the support is elongated in orientation_gabors of the parallel stripes of the function. 

% Spatial Frequency
num_cycles = 3;
freq = num_cycles / gabor_text_dim;

backgroundgabor_text = [0.5 0.5 0.5 0.5];

disable_norm = 1; %‘disable_norm’ Optional, defaults to 0. If set to a value of 1, the special multiplicative normalization term normf = 1/(sqrt(2*pi) * sc) will not be applied to the computed gabor_text.

pre_contrast_multiplier = 0.5; %‘contrastPreMultiplicator’ Optional, defaults to 1. This value is multiplied as a scaling factor to the requested contrast value

%% Build a procedural gabor_text texture
gabor_text = CreateProceduralGabor(window, gabor_text_dim, gabor_text_dim,[],backgroundgabor_text , disable_norm, pre_contrast_multiplier);

%Position of gabor_texts
dim =(16/2)*pixelperdeg; 
[x, y] = meshgrid(-dim:gabor_text_dim:dim, -dim:gabor_text_dim:dim); %create a matrix!!!!!

%distance of each gabor_text from the center of the array
dist = sqrt(x.^2 + y.^2);

%Inner annulus
inner_dist = 1*pixelperdeg;
x(dist <= inner_dist) = nan;
y(dist <= inner_dist) = nan;

%Outer annulus
outer_dist = 8*pixelperdeg; 
x(dist >= outer_dist) = nan;
y(dist >= outer_dist) = nan;

%Select only the finite values
x = x(isfinite(x));
y = y(isfinite(y));

% Center the annulus coordinates in the centre of the screen
% x_pos = x .* gabor_text_dim + center_x;
% y_pos = y .* gabor_text_dim + center_y;
x_pos = x + center_x;
y_pos = y + center_y;

% Count how many gabor_texts there are
n_gabors = numel(x_pos);

% Make the destination rectangles for all the gabor_texts in the array
base_rect = [0 0 gabor_text_dim gabor_text_dim];
all_rects = nan(4, n_gabors);
for i = 1:n_gabors
    all_rects(:, i) = CenterRectOnPointd(base_rect, x_pos(i), y_pos(i));
end


%%
%run either contra gabor or signal gabor code depending on which one's
%percentage varry 

% %contra gabor
% trialpercontracoh=50;
% listcohcontra= linspace(0,40,9);
% alllistcoh= repmat(listcohcontra,1,trialpercontracoh);
% n_trials=length(alllistcoh);

%signal gabor
trialpersignalcoh=50;
listcohsignal= linspace(5,80,16);
alllistcoh= repmat(listcohsignal,1,trialpersignalcoh);
n_trials=length(alllistcoh);

directions= [0,180];
directions= repmat(directions,1,round(n_trials/2));
%directions= Shuffle(directions);


% edit by TB: make equal number of trials for each direction and percentage
alllistcoh = sort(alllistcoh); 
% now there should be 50 5s, followed by 50 10s, etc...
% the directions list goes [0,180,0,180... so these should be aligned with
% the sorted percentages to make equal number

trialOrder = randperm(n_trials); % this is a list from 1 to ntrials in random order

directions = directions(trialOrder);
alllistcoh = alllistcoh(trialOrder);
% now the pairings should stay aligned, though they have been sorted into a
% random order

%Positions
numgenerator=1:n_gabors;
pos_gabors=Shuffle(numgenerator); 
pos_each_gabor=all_rects(1:4,pos_gabors); 

% Drift speed 
deg_per_sec = 3;  
cyclesPerDeg= num_cycles/gabor_text_dim*pixelperdeg;
phasesPerCycle = 360;

% drift_speed =  deg_per_sec*pixelperdeg * ifi;
drift_speed =  phasesPerCycle*deg_per_sec * ifi;

%%

respMat = nan(n_trials,6); % tb edit: save response and correct response

%%
%kbDev = -1; %for comp downstairs and for my own computer =2 or 4 for the
%computer in lab
kbDev = -1; 
KbName('UnifyKeyNames');
%Keyboard Info
escapeKey = KbName('ESCAPE');
leftKey = KbName('LeftArrow');
rightKey = KbName('RightArrow');

queueList = zeros(1,256);
queueList([escapeKey,leftKey,rightKey ]) = 1;

KbQueueCreate(kbDev,queueList);

%RestrictKeysForKbCheck([escapeKey,leftKey,rightKey]);

%%
% Interstimulus interval time in seconds and frames??????
isiTimeSecs = 0.5;
isiTimeFrames = round(isiTimeSecs / ifi);

% Numer of frames to wait before re-drawing
waitframes = 1 ;
nextTrialStart =0;
vbl = Screen('Flip', window); 

%%
stimDur = 2;
stimFlips = round(stimDur/ifi);
stimParams = struct([]);

% perc_noise=10;

%%
for trial=1:n_trials
    %Percentages
    
    %for contra gabor percentage variation
    %     perc_signal_gabor=50;
    %     perc_cont_gabor= alllistcoh(trial);
    
    %for signal gabor percentage variation
    perc_signal_gabor = alllistcoh(trial);
    perc_cont_gabor = 20;
    
    %if both the percentage of signal and contra predefined
    %     perc_signal_gabor= 90;
    %     perc_cont_gabor= 0;
    
    
    direction=directions(trial);
    
    %Properties matrix and gabors orientation
    [propertiesMat,orientation_gabors] = genPropertiesMat([NaN, freq, sigma, contrast, aspectRatio, 0, 0, 0], n_gabors, 4, 32, stimFlips, direction, drift_speed, perc_signal_gabor/100, perc_cont_gabor/100);

    
    KbReleaseWait(kbDev);
    KbQueueFlush(kbDev);
    KbQueueStart(kbDev);
    
    while GetSecs <= nextTrialStart
    end
    
    tStart = GetSecs; %to count the RT
    
    for fi = 1:stimFlips   
        Screen('DrawTextures', window, gabor_text, [], pos_each_gabor, orientation_gabors' ,[], [], [], [], ...
        kPsychDontDoRotation, propertiesMat(:,:,fi)');
%       Screen('DrawLines', window, allCoords,lineWidthPix, white, [center_x; center_y], 2);
       Screen('DrawDots', window, [center_x; center_y], 2, black, [], 2); %fixation point
        Screen('DrawingFinished',window);
        vbl = Screen('Flip', window);        
    end
    
    Screen('DrawDots', window, [center_x; center_y], 2, black, [], 2); %fixation point
    Screen('DrawingFinished',window);
    vbl = Screen('Flip', window);
    
    pressed = 0;
    doContinue = 0;
    while ~doContinue
        while ~pressed
            [pressed, keys] = KbQueueCheck(kbDev);
        end
        
        thisKey = min(find(keys,1,'first'));
        keyTime = min(keys(thisKey));
        
        if thisKey == escapeKey
            fprintf('User Quit')
            ShowCursor;
            sca;
            return
        elseif thisKey == leftKey
            response = 1;
            if directions(trial)==0
                correctResponse = 1;
            elseif directions(trial)==180
                correctResponse=0;
            end
            doContinue = 1;
        elseif thisKey == rightKey
            response = 2;
            if directions(trial)==180
                correctResponse = 1;
            elseif directions(trial)==0
                correctResponse=0;
            end
            doContinue = 1;
        else
            pressed = 0;
        end
        
    end
    
    KbQueueStop(kbDev);
    
    
    rt = keyTime - tStart;
    
    
    respMat(trial,1) = perc_cont_gabor;
    respMat(trial,2) = direction;
    respMat(trial,3) = response;
    respMat(trial,4) = rt;
    respMat(trial,5) = correctResponse; % TB edit: save response and correct response
    respMat(trial,6) = perc_signal_gabor;
    
    
    % TB edit: add space for a break for us
    if rem(trial,100) == 0 % every 100 trials
        % if we don't wait for the key to be lifted above then we need to
        % do so here
        KbReleaseWait(kbDev);
        KbQueueFlush(kbDev);
        KbQueueStart(kbDev);
        
        breaktext = ['You have reached trial ' num2str(trial) ' of ' num2str(n_trials) '\nPress left to continue'];
        DrawFormattedText(window, breaktext,'center','center', black);
        Screen('Flip',window);
        
        continueNextTrial = 0;
        while ~continueNextTrial
            [continueNextTrial, keys] = KbQueueCheck(kbDev);
        end
        KbQueueStop(kbDev);
        
    end
       
    nextTrialStart = GetSecs+isiTimeSecs;
    
    stimParams(trial).positions = pos_each_gabor;
    stimParams(trial).orientation = orientation_gabors;
    stimParams(trial).propertiesmat = propertiesMat;  
    
end

respMatTable=array2table(respMat);
respMatTable.Properties.VariableNames(1:5) ={'PercContra','Direction','Response','ReactionTime','Correct','PercSignal'};

RestrictKeysForKbCheck([])
sca; 
% close all;