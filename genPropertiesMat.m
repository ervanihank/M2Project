% function to generate the properties matrix for limited lifetime gabor
% stimuli

% ----------- INPUT --------------

% properties:   a 1x8 array of the basic properties for procedural gabors
%               should contain phase, frequency, sigma, contrast, aspect
%               ratio, followed by 3 zeros. phase and contrast will be
%               reset with this code.
% ngabors:      the number of gabors 
% ngroups:      the number of groups (the members of one group have their
%               lifetime reset at the same time)
% lifetime:     the number of screen flips in the total lifetime of the
%               groups. With 4 groups and a lifetime of 8, the gabors will 
%               display for 6 frames with 2 frames off.
% nflips:       the number of screen flips of stimulus display
% direction:    the direction of the signal gabors
% speed:        the speed at which each gabor is moving
% pSignal:      the proportion of gabors moving in the signal direction.
%               Can be a single value (always the same pSignal each 
%               lifetime, or an mxn matrix, where m is sufficient to cover 
%               the number of lifetimes, and n is the number of groups.
% pContra:      the proportion of gabors moving in the opposite direction
%               to the signal gabors (pSignal+pContra<=1). Same size as
%               pSignal.


function [propertiesMat, orientations] = genPropertiesMat(properties,ngabors,ngroups,lifetime,nflips,direction,speed,pSignal,pContra)

% check and define the number of flips on in the lifetime

onFlips = (lifetime/ngroups)*(ngroups-1);

if rem(onFlips,1)>0
    fprintf('\n\nLifetime and number of groups poorly specified, rounding number of flips');
    onFlips = round(onFlips);
end

offFlips = lifetime-onFlips;

if offFlips<1
    fprintf('\n\nLifetime does not allow for off flips, quitting')
    return
end

totalLifetimes = ceil(nflips/lifetime);

% check pSignal and pContra
if any(size(pSignal)~=size(pContra))
    fprintf('\n\npSignal and pContra must be the same size')
    return
end
if any((pSignal+pContra)>1)
    fprintf('\n\npSignal+pContra>1, must be less than or equal to 1')
    return
end
if length(pSignal)>1
    if size(pSignal,2)~=ngroups
        fprintf('\n\nWhen pSignal changes, it must be specified for each group')
        return
    end
end

% asign groups

if rem(ngabors/ngroups,1)
    
    % one group will have more/fewer gabors than the rest
    nGaborsPerGroup = zeros(1,ngroups);
    nGaborsPerGroup(1:ngroups-1) = round(ngabors/ngroups);
    nGaborsPerGroup(ngroups) = ngabors-sum(nGaborsPerGroup);
       
else
    
    % same number of members in each group
    nGaborsPerGroup = ones(1,ngroups)*ngabors/ngroups;

end

groupCount = [0,cumsum(nGaborsPerGroup)];

groupIndex = zeros(ngabors,1);
for gi = 1:ngroups
    groupIndex(groupCount(gi)+1:groupCount(gi+1)) = gi;
end

% asign groups to random positions
groupIndex = Shuffle(groupIndex);

% to update the gabor phases we first need to asign the directions
if length(pSignal)==1
    totalSignal = round(ngabors*pSignal);
    totalContra = round(ngabors*pContra);
    
    nSignalPerGroup = round(ones(1,ngroups).*nGaborsPerGroup.*pSignal);
    nContraPerGroup = round(ones(1,ngroups).*nGaborsPerGroup.*pContra);
    
    % make sure we end up with the right amount of signal and contra in total
    if sum(nSignalPerGroup)~=totalSignal
        nSignalPerGroup(ngroups) = nSignalPerGroup(ngroups) + (sum(nSignalPerGroup)-totalSignal);
        nContraPerGroup(ngroups) = nContraPerGroup(ngroups) + (sum(nContraPerGroup)-totalContra);
    end
    nNoisePerGroup = nGaborsPerGroup - (nSignalPerGroup+nContraPerGroup);
    
    nSignalPerGroup = repmat(nSignalPerGroup,(totalLifetimes+ngroups),1);
    nContraPerGroup = repmat(nContraPerGroup,(totalLifetimes+ngroups),1);
    nNoisePerGroup = repmat(nNoisePerGroup,(totalLifetimes+ngroups),1);
    
else % new signal ratio for each lifetime
    
    if size(pSignal,1)<(totalLifetimes+ngroups)
        pSignal = repmat(pSignal,ceil((totalLifetimes+ngroups)/size(pSignal,1)),1);
        pContra = repmat(pContra,ceil((totalLifetimes+ngroups)/size(pContra,1)),1);
    end
    pSignal = pSignal(1:(totalLifetimes+ngroups),:);
    pContra = pContra(1:(totalLifetimes+ngroups),:);
    
    nSignalPerGroup = zeros((totalLifetimes+ngroups),ngroups);
    nContraPerGroup = zeros((totalLifetimes+ngroups),ngroups);
    
    for li = 1:(totalLifetimes+ngroups)
        
        nSignalPerGroup(li,:) = round(ones(1,ngroups).*nGaborsPerGroup.*pSignal(li,:));
        nContraPerGroup(li,:) = round(ones(1,ngroups).*nGaborsPerGroup.*pContra(li,:));
        
    end
    
    nNoisePerGroup = nGaborsPerGroup - (nSignalPerGroup+nContraPerGroup);
        
end
contraDir = direction+180;

% make a big matrix, with more than enough lifetimes to cover the number of
% flips and the offset in lifetime time

phases = zeros(ngabors,1,(totalLifetimes+ngroups)*lifetime);
% start with a random phase
phases(:,1,1) = randi([1,360],ngabors,1,1);

% set up the Gabor orientations
orientations = randi([1,360],ngabors,1);

% go through the groups to asign the phases each lifetime
lifeUpdate = 1:lifetime;

for gi = 1:ngroups
    
    % the orientations of the gabors in this group
    thisgrouporis = orientations(groupIndex==gi);
    
    for li = 1:(totalLifetimes+ngroups)
        
        % the list of directions in this group
        thisgroupdirs = [ones(1,nSignalPerGroup(li,gi))*direction,...
            ones(1,nContraPerGroup(li,gi))*contraDir,...
            randi([1, 360],1,nNoisePerGroup(li,gi))];
        
        % shuffle the directions
        thisgroupdirs = Shuffle(thisgroupdirs);
        
        % calculate the phase update each frame
        phaseUpdate = cosd(thisgrouporis-thisgroupdirs') .* speed;
        
        phaseUpdate = phaseUpdate.*lifeUpdate;
        
        % allocate to the appropriate identity, at the appropriate time
        startphase = phases(groupIndex==gi,1,((li-1)*lifetime)+1);
        phases(groupIndex==gi,1,((li-1)*lifetime)+1:(li*lifetime)) = startphase+phaseUpdate;
        % for the next round
        phases(groupIndex==gi,1,(li*lifetime)+1) = phases(groupIndex==gi,1,(li*lifetime));
    end
end

    
% set up the propertiesMat
propertiesMat = repmat(properties,ngabors,1,round(nflips));

% allocate the phases and the contrast for each group
contrastMultiplier = repmat([ones(1,onFlips),zeros(1,offFlips)],1,(totalLifetimes+ngroups));

for gi = 1:ngroups
    
    % each groups lifetime starts at a different time
    startFlip= lifetime-((gi-1)*offFlips)+1;
    thisGroupInds = startFlip:(nflips+startFlip-1);
    
    % phase
    propertiesMat(groupIndex==gi,1,:) = phases(groupIndex==gi,1,thisGroupInds);
    
    % contrast
    thiscontrast = contrastMultiplier(thisGroupInds);
    thiscontrast = repmat(thiscontrast,nGaborsPerGroup(gi),1);
    thiscontrast = reshape(thiscontrast,nGaborsPerGroup(gi),1,round(nflips));
    
    propertiesMat(groupIndex==gi,4,:) = propertiesMat(groupIndex==gi,4,:).*thiscontrast;
    
end
    


end

