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

% for the computer downstairs
% ScreenHypotenuse = 50.5;
% pixHypotenuse = sqrt((monitorPos(3)^2)+(monitorPos(4)^2));
% pixelperdeg = round(pixHypotenuse/ScreenHypotenuse);


%pixelperdeg=36;
pixelperdeg=25;


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

%Define gabor_text Parameters
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

% Build a procedural gabor_text texture
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

%percentages and directions for n_trials
trialpercond=50;
% listpercentages=linspace(5,80,16);
% allpercentages=repmat(listpercentages,1,trialpercond);
% n_trials= length(allpercentages);
%percentages= Shuffle(allpercentages);


%contra gabor
trialpercontracoh=50;
listcoherence= linspace(0,40,9);
alllistcoh= repmat(listcoherence,1,trialpercontracoh);
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

respMat = nan(n_trials,5); % tb edit: save response and correct response

%kbDev = -1; %for comp downstairs
kbDev = 4;
KbName('UnifyKeyNames');
%Keyboard Info
escapeKey = KbName('ESCAPE');
leftKey = KbName('LeftArrow');
rightKey = KbName('RightArrow');

queueList = zeros(1,256);
queueList([escapeKey,leftKey,rightKey ]) = 1;

KbQueueCreate(kbDev,queueList);

%RestrictKeysForKbCheck([escapeKey,leftKey,rightKey]);
% Interstimulus interval time in seconds and frames
isiTimeSecs = 0.5;
isiTimeFrames = round(isiTimeSecs / ifi);

% Numer of frames to wait before re-drawing
waitframes = 1 ;
nextTrialStart =0;
vbl = Screen('Flip', window); 

%maxwait=5;

stimDur = 5;
stimFlips = stimDur/ifi;
stimParams = struct([]);

% perc_noise=10;


for trial=1:n_trials
    %Percentages
    % edit TB: change to allpercentages
    %perc_signal_gabor=percentages(trial);
% %     perc_signal_gabor=allpercentages(trial);
    perc_signal_gabor=50;
    direction=directions(trial);
    perc_cont_gabor= alllistcoh(trial);
    n_signal_gabor=round((n_gabors *perc_signal_gabor)/100);
    n_contra_gabor=round((n_gabors *perc_cont_gabor)/100);
    n_noise_gabor=n_gabors-n_signal_gabor-n_contra_gabor;
    
    %Directions of Gabors
    glob_direc_signal= repmat(direction,n_signal_gabor,1); %global direction=0 (leftward)
    glob_direc_contra= repmat(180-direction,n_contra_gabor,1); 
    glob_direc_noise= randi(360,1,n_noise_gabor);
    
    all_gabors_direc=[glob_direc_signal',glob_direc_contra',glob_direc_noise];
    all_gabors_direc=Shuffle(all_gabors_direc);
    
    %Orientations of Gabors
    orientation_gabors=rand(1, n_gabors) .* 180;
    
    %Speed of each Gabor
    speed_each_gabor = cosd(orientation_gabors-all_gabors_direc) .* drift_speed;
    
    %Properties matrix with the same phase for each gabor_text
    % phase=speed_each_gabor;
    phase=randi(360,1,n_gabors);
    propertiesMat = repmat([NaN, freq, sigma, contrast, aspectRatio, 0, 0, 0],n_gabors, 1);
    propertiesMat(:, 1) = phase';
    
    % Perform initial flip/start of animation
    % vbl = Screen('Flip', window);
    % reset queue
    KbReleaseWait(kbDev);
    KbQueueFlush(kbDev);
    KbQueueStart(kbDev);
    
    while GetSecs <= nextTrialStart
    end
    
    
    tStart = GetSecs; %to count the RT
    %stoptime=GetSecs+maxwait;
    % [keyIsDown,secs, keyCode] = KbCheck;
    % while keyIsDown
    %     Screen('DrawTextures', window, gabor_text, [], pos_each_gabor, orientation_gabors' ,[], [], [], [], kPsychDontDoRotation, propertiesMat');
    %     Screen('DrawDots', window, [center_x; center_y], 5, black, [], 2); %fixation point
    %     vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
    %     phase = phase + speed_each_gabor;
    %     propertiesMat(:, 1) = phase';
    %     respToBeMade = false;
    %     [keyIsDown]=KbCheck;
    % end
    
    for fi = 1:stimFlips
        
        Screen('DrawTextures', window, gabor_text, [], pos_each_gabor, orientation_gabors' ,[], [], [], [], kPsychDontDoRotation, propertiesMat');
        Screen('DrawDots', window, [center_x; center_y], 5, black, [], 2); %fixation point
        Screen('DrawingFinished',window);
        vbl = Screen('Flip', window);
        phase = phase + speed_each_gabor; %The phase of each gabor_text is changing by the same amount on each frame(they each have the same speed in different all_gabors_direc).
        propertiesMat(:, 1) = phase';
        
    end
    
    Screen('DrawDots', window, [center_x; center_y], 5, black, [], 2); %fixation point
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
    
    % TB edit: this doesn't need to be done each trial - can wait till
    % after
    %     respMatTable=array2table(respMat);
    %     respMatTable.Properties.VariableNames(1:4) ={'Percentage','Direction','Response','Reaction Time'};
    
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
    % for frame = 1:isiTimeFrames - 1
    %
    % Screen('DrawDots', window, [center_x; center_y], 5, black, [], 2); %fixation point
    % vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
    % end
    
    stimParams(trial).positions = pos_each_gabor;
    stimParams(trial).orientation = orientation_gabors;
    stimParams(trial).speed = speed_each_gabor;
    stimParams(trial).direction = all_gabors_direc;
    
    
end

respMatTable=array2table(respMat);
respMatTable.Properties.VariableNames(1:5) ={'Percentage','Direction','Response','ReactionTime','Correct'};

RestrictKeysForKbCheck([])
sca; 
% close all;