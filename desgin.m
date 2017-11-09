


function dotState = dynamicRTBP(time, state, massParameter, eccentricity, varargin)
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



function dotState = dynamicRTBPAtLibrationPoint(time, state, massParameter, eccentricity, varargin)
% dynamic fucntions tranformed to libration point
% HOW TO WRTIE THIS?
end



function dotState = dynamicEphemerisMice(time, state, frame, celestialBodyIDs, varargin)
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



function [newTime, newState, exitflag] = shootingMethod()
% nd - number of states dimension, should be 3 (default) or 2
% ns - number of segment
%
% variables:
%   numberLayer: 1. simple shooting; 2. connection shooting then smoothness optimization
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

%% numberLayer = 1, simple shooting
while ~done
    % 算 Jacobian 矩阵
    % ode113 ?? 貌似比 ode45 效果更好
    
    % 求伪逆
    
    % 判断 done，达到则退出
    
    % 更新状态
    
end

%% numberLayer = 2, connection shooting then smoothness optimization 
%
%（还不太确定如何编程，需要查文献）
% 暂时搁置
%
while ~done
end

%% if fails
exitflag = 0;

end








function [newTime, newState] = convertion(directionOldToNew, oldDynamicFcn, oldTIme, oldState, newDynamic, newTime, newState)
% direction: 1. old -->> new
%           -1. old <<-- new
end
















