function [correctedInitialEpoches, correctedInitialStates, exitflag] = MultipleShooting(dynamicFcn, initialEpoches, initialStates, positionTolerance, velocityTolerance, odeOptions)
%multiple shooting method in the given dynamic system
%
% This function will automatically use 'parfor' (after R2017a) if a parpool has been open outside. 
%
% References:
%   Marchand, Belinda G., Kathleen C. Howell, and Roby S. Wilson. 2007. “Improved Corrections Process for Constrained Trajectory Design in the N-Body Problem.” Journal of Spacecraft and Rockets 44 (4):884–97. https://doi.org/10.2514/1.27205.
%   Parker, Jeffrey S., and Rodney L. Anderson. 2014. Low-Energy Lunar Trajectory Design. 1st ed. JPL Deep-Space Communications and Navigation Series, July. Wiley. http://descanso.jpl.nasa.gov/Monograph/series12_chapter.cfm?force_external=0.
%
% see also: PositionShooting

if nargin < 4
    positionTolerance = 1e-9;
    velocityTolerance = 1e-7;
    odeOptions = odeset('AbsTol',1e-12,'RelTol',1e-12);
    disp('Default position, velocity, and propagation tolerance is used.');
end
flagAutoUseParallel = 1; % 1: use 'parfor' if find one parallel pool
initialPositionFixed = false; % 为 true 时较容易收敛
finalPositionFixed = false; % 为 true 时较难收敛

iterationNumberLevelTwoMax = 20;

segmentNumber = length(initialEpoches)-1; % the last element is the target epoch, similar to initialStates

iterationNumberLevelTwo = 1;
correctedInitialEpoches = initialEpoches; % will be corrected by multiple shooting method
correctedInitialStates = initialStates;   % will be corrected by multiple shooting method
while 1
    
    
    %---- level-1 shooting ----
    
    
    % simple shooting to connect position of every segment
    stateTransitionMatrixes = zeros(segmentNumber, 6, 6);
    exitflag = zeros(segmentNumber,1);
    correctedFinalStates = zeros(segmentNumber,6); % used to calculate Delta V
    p = gcp('nocreate');
    if flagAutoUseParallel && ~isempty( p ) % this will test if a parpool is open, and suppress gcp to creat one
        % use parallel
        %   usually faster, but your dynamicFcn must support parallel computing
        disp(['#    level-1 shooting begin... (using ' num2str(p.NumWorkers) ' workers)']);
        parforCorrectedInitialEpoches = correctedInitialEpoches;  % auxiliary variables for parfor, will be discarded after parfor
        parforCorrectedInitialStates = correctedInitialStates;    % auxiliary variables for parfor, will be discarded after parfor
        parfor iiSegment = 1:segmentNumber
            [correctedInitialStates(iiSegment,:), correctedFinalStates(iiSegment,:), stateTransitionMatrixes(iiSegment,:,:), exitflagLevelOne(iiSegment)]...
                = PositionShooting(...
                dynamicFcn,...
                correctedInitialEpoches(iiSegment),         correctedInitialStates(iiSegment,:),...
                parforCorrectedInitialEpoches(iiSegment+1), parforCorrectedInitialStates(iiSegment+1,:),...
                positionTolerance, odeOptions);
            if exitflagLevelOne(iiSegment) ~= 1
                disp(['# !!! fail: segment ' num2str(iiSegment) ' fails at level-1 shooting.']);
            else
                disp(['#      segment ' num2str(iiSegment) ' done.']);
            end
        end
    else
        % use serial
        %   slower, but easier for testing, debugging, and profiling
        disp(['#    level-1 position shooting begin...']);
        for iiSegment = 1:segmentNumber
            [correctedInitialStates(iiSegment,:), correctedFinalStates(iiSegment,:), stateTransitionMatrixes(iiSegment,:,:), exitflagLevelOne(iiSegment)]...
                = PositionShooting(...
                dynamicFcn,...
                correctedInitialEpoches(iiSegment),   correctedInitialStates(iiSegment,:),...
                correctedInitialEpoches(iiSegment+1), correctedInitialStates(iiSegment+1,:),...
                positionTolerance, odeOptions);
            if exitflagLevelOne(iiSegment) ~= 1
                disp(['# !!! fail: segment ' num2str(iiSegment) ' fails at the level-1 shooting.']);
                break;
            else
                disp(['#      segment ' num2str(iiSegment) ' done.']);
            end
        end
    end
    disp(['#    level-1 done.']);
    % draw level-1 shooting results
    figure(99); clf; PlotInitialState(dynamicFcn, correctedInitialEpoches, correctedInitialStates);
    
    % test failure
    if any(exitflagLevelOne ~= 1)
        exitflag = -2;
        break; % exitflag already set to -2 in above code
    end
    
    
    
    %---- level-2 shooting ----
    
    
    
    % collcect the target error
    deltaVelocity = correctedFinalStates(1:segmentNumber-1,4:6) - correctedInitialStates(2:segmentNumber,4:6);
    deltaVelocity = reshape( deltaVelocity.', [], 1 );
    disp(['#  level-2 iter ' num2str(iterationNumberLevelTwo) ' norm: ' num2str(norm(deltaVelocity))]);
    
    % test early stop
    if norm(deltaVelocity) < velocityTolerance
        exitflag = 1;
        disp(['# multiple-shooting success. norm(dV) = ' num2str(norm(deltaVelocity))]);
        break;
    end
    
    % modify epoch and velocity of all segments at once
    % after one modification, shooting position again
    stateRelationshipMatrix = zeros( (segmentNumber-2)*3+3, (segmentNumber-2)*4+12 );
    for ii = 2:segmentNumber
        % generate state relationship matrix
        stm21 = squeeze( stateTransitionMatrixes(ii-1, :, :) );
        stm12 = inv(stm21);
        stm32 = squeeze( stateTransitionMatrixes(ii, :, :) );
        %stm23 = inv(stm32); % seems do not need this
        %
        v1plus  = correctedInitialStates(ii-1, 4:6).';
        v2minus = correctedFinalStates(ii-1, 4:6).';
        v2plus  = correctedInitialStates(ii, 4:6).';
        v3minus = correctedFinalStates(ii, 4:6).';
        %
        a2minus = dynamicFcn( correctedInitialEpoches(ii), correctedFinalStates(ii-1, :).' ); a2minus = a2minus(4:6);
        a2plus  = dynamicFcn( correctedInitialEpoches(ii), correctedInitialStates(ii, :).' ); a2plus = a2plus(4:6);
        %
        stateRelationshipMatrix( (ii-2)*3+1:(ii-2)*3+3, (ii-2)*4+1:(ii-2)*4+12 ) = [...
            -inv(stm12(1:3,4:6)),...
            stm12(1:3,4:6)\v1plus,...
            -stm32(1:3,4:6)\stm32(1:3,1:3) + stm12(1:3,4:6)\stm12(1:3,1:3),...
            (a2plus-a2minus) + (stm32(1:3,4:6)\stm32(1:3,1:3)*v2plus-stm12(1:3,4:6)\stm12(1:3,1:3)*v2minus),...
            inv(stm32(1:3,4:6)),...
            -stm32(1:3,4:6)\v3minus];
    end
    %disp(stateRelationshipMatrix);
    
    % solve for the correction
    if initialPositionFixed
        stateRelationshipMatrix = stateRelationshipMatrix(:, 5:end);
    end
    if finalPositionFixed
        stateRelationshipMatrix = stateRelationshipMatrix(:, 1:end-2);
    end
    correctionAtPositionAndEpoch = pinv(stateRelationshipMatrix) * deltaVelocity;
    
    % updata initial epoches and states for next iteration
    updateSegments = 1:segmentNumber+1;
    indexOffset = 0;
    if initialPositionFixed
        updateSegments = setdiff(updateSegments, 1); % do not update first state if it should be fixed
        indexOffset = -1; % to compensate that first element of correctionAtPositionAndEpoch is for the 2nd point
    end
    if finalPositionFixed
        updateSegments = setdiff(updateSegments, segmentNumber+1); % do not update last state if it should be fixed
    end
    sigma = 1;
    for ii = updateSegments
        correctedInitialEpoches(ii)    = correctedInitialEpoches(ii)    + sigma * correctionAtPositionAndEpoch((ii+indexOffset-1)*4+4);
        correctedInitialStates(ii,1:3) = correctedInitialStates(ii,1:3) + sigma * correctionAtPositionAndEpoch((ii+indexOffset-1)*4+1:(ii+indexOffset-1)*4+3).';
    end
    iterationNumberLevelTwo = iterationNumberLevelTwo + 1;
    
    % stop after too many iterations
    if iterationNumberLevelTwo > iterationNumberLevelTwoMax
        exitflag = -1;
        disp(['#  !!! fail: level-2 shooting exceeds maximum iteration number ' num2str(iterationNumberLevelTwoMax) '.']);
        break;
    end
    
end

% return corrected epoches and states
% correctedInitialEpoches = correctedInitialEpoches;
% correctedInitialStates = correctedInitialStates;

end








function [correctedInitialState, correctedFinalState, correctedStateTransitionMatrix, exitflag] = PositionShooting(dynamicFcn, initialEpoch, initialState, targetEpoch, targetState, positionTolerance, odeOptions)
%[correctedInitialState, correctedFinalState, correctedStateTransitionMatrix, exitflag] = PositionShooting(dynamicFcn, initialEpoch, initialState, targetEpoch, targetState, positionTolerance, odeOptions)
%correctedInitialState = PositionShooting(dynamicFcn, initialEpoch, initialState, targetEpoch, targetState, positionTolerance)
%simple shooting method for position in the given dynamic system
% This is called by MultipleShooting.m for the level-1 shooting
%
%   fixed variables: initial epoch, initial position, target epoch, target position
%    free variables: initlal velocity
%
% inputs:
%   dynamicFcn
%   initialEpoch
%   initialState
%   targetEpoch
%   targetState
%   positionTolerance
%   [odeOptions]: has a default value
% outputs:
%   correctedInitialState: shooting result
%   correctedFinalState: final state is used by MultipleShooting, but not for most other situations.
%   correctedStateTransitionMatrix: required by MultipleShooting to avoid recalculation.
%   exitflag: 1 for suceess; -1 for too many iterations. (need to check more cases)
%
%   see also: MultipleShooting
%
% References:
%   Marchand, Belinda G., Kathleen C. Howell, and Roby S. Wilson. 2007. “Improved Corrections Process for Constrained Trajectory Design in the N-Body Problem.” Journal of Spacecraft and Rockets 44 (4):884–97. https://doi.org/10.2514/1.27205.
%   Parker, Jeffrey S., and Rodney L. Anderson. 2014. Low-Energy Lunar Trajectory Design. 1st ed. JPL Deep-Space Communications and Navigation Series, July. Wiley. http://descanso.jpl.nasa.gov/Monograph/series12_chapter.cfm?force_external=0.

% check inputs
if nargin < 7
    odeOptions = odeset('AbsTol',1e-9,'RelTol',1e-9);
    disp('Default propagation accuracy of 1e-9 is used.');
end
iterationNumberMax = 50;

fixedPropagationInterval = [initialEpoch, (initialEpoch+targetEpoch)/2, targetEpoch]; % insert a middle point to eliminate too much redundant results
iterationNumber = 1;
while 1
    % calculate state transition matrix
    [~,allState] = ode113( dynamicFcn, fixedPropagationInterval, [initialState,reshape(eye(6),1,36)], odeOptions );
    finalState = allState(end, 1:6);
    stateTransitionMatrix = reshape( allState(end,7:42).', 6, 6 );
    
    % check if target is reached
    errorFinalState = targetState(1:3) - finalState(1:3);
    %disp(['debug: position shooting: iter ' num2str(iterationNumber) ': error is ' num2str(norm(errorFinalState(1:3)))]);

    % test early stop
    if norm(errorFinalState(1:3)) < positionTolerance
        exitflag = 1; % success
        break;
    end
    
    % solve for the correction at initial state
    L = stateTransitionMatrix( 1:3, 4:6 );
    correctionAtInitialState = ( L \ errorFinalState.' ).';

    % update state for next iteration
    sigma = 0.618;
    initialState(4:6) = initialState(4:6) + sigma * correctionAtInitialState(1:3);
    iterationNumber = iterationNumber + 1;
    
    % stop after too many iterations
    if iterationNumber > iterationNumberMax
        exitflag = -1; % too much iterations
        disp(['position shooting: max iteration reached.']);
        break;
    end
end

% generate output
correctedInitialState = initialState;
correctedFinalState = finalState;
correctedStateTransitionMatrix = stateTransitionMatrix;

end