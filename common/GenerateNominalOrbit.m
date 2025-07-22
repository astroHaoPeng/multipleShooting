function [outT, outX] = GenerateNominalOrbit(dynamicFcn,t,X,odeOptions)
%% Generate segmented trajectory with initial condition for each segment saved in `t` and `X`.
%   The difference to `PlotInitialState.m` is this function doesn't plot.

for ii = 1:size(X,1)-1
    if nargin==4
        [plotT, plotX] = ode113(dynamicFcn,[t(ii),t(ii+1)],X(ii,:),odeOptions);
    else
        [plotT, plotX] = ode113(dynamicFcn,[t(ii),t(ii+1)],X(ii,:));
    end
    if ii == 1
        outT = plotT;
        outX = plotX;
    else
        outT = [outT; plotT];
        outX = [outX; plotX];
    end
end

end
