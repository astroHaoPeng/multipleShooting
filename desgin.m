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



function [newTime, newState, exitflag] = ShootingCorrection()
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

%% numberLevel = 1, simple shooting
while ~done
    % 算 Jacobian 矩阵
    % ode113 ?? 貌似比 ode45 效果更好
    
    % 求伪逆
    % The linear system given in Eq. (2.20) is underdetermined; it is common practice to use the smallest Euclidean norm to produce a good solution [134] (Parker2014, p58)
    
    % 判断 done，达到则退出
    
    % 更新状态
    
end

%% numberLevel = 2, connection shooting then smoothness optimization 
%
%（还不太确定如何编程，需要查文献）
% 暂时搁置
%
while ~done
end

%% if fails
exitflag = 0;

end




function [newTime, newState] = Convertion(directionOldToNew, oldDynamicFcn, oldTIme, oldState, newDynamic, newTime, newState)
%convertion function that will be required by shooting method, to covert between old and new dynamic models
% direction: 1. old -->> new
%           -1. old <<-- new
end



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













