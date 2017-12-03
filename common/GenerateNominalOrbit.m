function [outT, outX] = GenerateNominalOrbit(dynamicFcn,t,X,odeOptions)

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
