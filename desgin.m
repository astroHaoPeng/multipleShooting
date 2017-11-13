% 参考文献：
% books:
%   Parker, Jeffrey S., and Rodney L. Anderson. 2014. Low-Energy Lunar Trajectory Design. 1st ed. JPL Deep-Space Communications and Navigation Series, July. Wiley. http://descanso.jpl.nasa.gov/Monograph/series12_chapter.cfm?force_external=0.
% thesis:
%   Pernicka, Henry John. 1990. “The Numerical Determination of Nominal Libration Point Trajectories and Development of a Station-Keeping Strategy.” Ph.D., Purdue University. http://adsabs.harvard.edu/abs/1990PhDT........26P.
%   Parker, Jeffrey S. 2007. “Low-Energy Ballistic Lunar Transfers.” Doctor of Philosophy, University of Colorado.





function dotState = DynamicRTBP(time, state, massParameter, eccentricity, varargin)
% dynamic function of general Restricted Three-Body Problem (RTBP)
%
% nd - number of states dimension, should be 3 (default) or 2
% ns - number of segment
%
% inputs:
%   time: [normalized] 
%   state: [normalized] 
%   massParameter: mu = m2 / (m1+m2), with m1 > m2.
%   eccentricity: 0 for CRTBP
% outputs:
%   dotState: ns 
end



function dotState = DynamicRTBPCeterAtLibrationPoint(time, state, massParameter, eccentricity, varargin)
% dynamic fucntions tranformed to libration point
% HOW TO WRTIE THIS?
end



function dotState = DynamicEphemerisMice(time, state, frame, celestialBodyIDs, varargin)
% 
% inputs:
%   time: [?]
%   state: [?]
%   frame: spice kernal frame in which state is defined 
%   celestialBodyIDs: ID of celestial bodies that considered in this function
%   varargin: tbd
% outputs:
%   dotState:
end




function [finalEpoch, finalState, exitflag] = PositionShootingVaryFinalEpoch(dynamicFcn, initialEpoch, initialState, targetEpoch, targetState, positionTolerance, odeOptions)
% [status: design, done]
%simple shooting method for position in the given dynamic system
% This is called by MultipleShooting.m for the level-1 shooting
%
%   fixed variables: initial epoch, initial position, target position
%   free  variables: initlal velocity, final epoch (or propagation duration)
%
%   see also: MultipleShooting
%
% References:
%   Howell, Kathleen C., and Henry John Pernicka. 1987. “Numerical Determination of Lissajous Trajectories in the Restricted Three-Body Problem.” Celestial Mechanics 41 (1–4):107–24. https://doi.org/10.1007/BF01238756.

% check inputs
if nargin < 7
    odeOptions = odeset('AbsTol',1e-9,'RelTol',1e-9);
    warning('Default propagation accuracy of 1e-9 is used.');
end

iterationNumber = 1;
propagationDuration = targetEpoch - initialEpoch; 
while 1
    % calculate state transition matrix
    tspan = initialEpoch + propagationDuration * [0, 1/2, 1]; % insert a middle point to eliminate too much redundant results
    [allEpoch,allState] = ode113( dynamicFcn, tspan, initialState, odeOptions );
    finalEpoch = allEpoch(end);
    finalState = allState(end, 1:6);
    stateTransitionMatrix = reshape( allState(end,7:42).', 6, 6 );
    
    % check if target is reached
    errorFinalState = targetState(1:3) - finalState(1:3);
    if norm(errorFinalState(1:3)) < positionTolerance
        exitflag = 1; % success
        break;
    end
    if iterationNumber > 20
        exitflag = -1; % too much iterations
        break;
    end

    % solve for the correction at initial state
    finalStateDot = dynamicFcn( finalState, finalState );
    L = [ stateTransitionMatrix( 1:3, 4:6 ), finalStateDot(1:3) ]; % ？是否也可以是 finalState(4:6)，可以节省一次调用
    correctionAtInitialStateAndDuration = pinv(L) * errorFinalState;

    % update state and duration for next iteration
    propagationDuration = propagationDuration + correctionAtInitialStateAndDuration(4);
    initialState(4:6) = initialState(4:6) + correctionAtInitialStateAndDuration(1:3);
end

% generate output
finalEpoch = initialEpoch;
finalState = initialState;

end



function [newTime, newState] = Convertion(directionOldToNew, oldDynamicFcn, oldTIme, oldState, newDynamic, newTime, newState)
%convertion function that will be required by shooting method, to covert between old and new dynamic models
% direction: 1. old -->> new
%           -1. old <<-- new
end


% function [newTime, newState, exitflag] = RefineOrbit()
% nd - number of states dimension, should be 3 (default) or 2
% ns - number of segment
%
% variables:
%   numberLevel: 1. simple shooting; 2. connection shooting then smoothness optimization
%   numberSegment
%   
%   oldDynamicFcn: 
%   oldTime: nd x 1
%   oldState: nd x ns
%   oldFixedIndex: index of varaibles that should not be changed
%   oldInitialConstraintFcn
%   oldFinalConstraintFcn
%
%   newDynamicFcn: 
%   newTime: nd x 1 (output)
%   newState: nd x ns (output)
%   newFixedIndex: index of varaibles that should not be changed
%   newInitialConstraint
%   newFinalConstraintFcn
%
%   convertionFcn: convert between old and new state
%   
%   tolerance = [tolPosition, tolVelocity, tolInitial, tolFinal] [new units]
%
%   exitflag: -1. fail
%              0. success
%              other 
%   

%% check parameters
% end



% frames which state should be freely converted from and to 
%   barycenter-centered rotating frame (ERTBP) (pulsating)
%   barycenter-centered rotating frame (CRTBP) (none-pulsating)
%   m1-centered rotating (pulsating or none-pulsating)
%   m2-centered rotating (pulsating or none-pulsating)
%   m1-centered inertial
%   m2-centered inertial
%   reference ephemeris frame (abbr. REF)
%       this frame connects nominal orbits and refined orbits
%       it's easy to be converted to any ephemeris frame after then
%           e.g. X0 should be converted to the REF, then refined in ephemeris model using shooting method
%       ? Earth- or m2-centered? Synodic one seems reasonable.
%       ? Barycenter-centered? Should be synodic. 
%
% Need to design the convertion sequence, choose a root one.
% Prefer the m2- or m1-centered inertial frame
% because for Sun-Earth/Moon system m2-centered inertial frame is at Earth
% and for Earth-Moon system m1-center frame is at Earth.
% Besides, it is easy to convert to ephemeris model.
%
% Need to check what I already have, and clean up the codes.













