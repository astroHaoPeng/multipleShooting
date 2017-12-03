function [outT, outX] = PlotInitialState(dynamicFcn,t,X,odeOptions)

for ii = 1:size(X,1)-1
    if nargin==4
        [plotT,plotX] = ode113(dynamicFcn,[t(ii),t(ii+1)],X(ii,:),odeOptions);
    else
        [plotT,plotX] = ode113(dynamicFcn,[t(ii),t(ii+1)],X(ii,:));
    end
    
    if nargout == 0
        plot3(plotX(:,1),plotX(:,2),plotX(:,3),'-'); hold on;
        plot3(plotX(1,1),plotX(1,2),plotX(1,3),'b*'); hold on;
        plot3(plotX(end,1),plotX(end,2),plotX(end,3),'ro'); hold on;
    end
    if ii == 1
        outT = plotT;
        outX = plotX;
    else
        outT = [outT; plotT(2:end,:)];
        outX = [outX; plotX(2:end,:)];
    end
end

if nargout == 0
    axis equal;
    xlabel('x');
    ylabel('y');
    drawnow update;
end

end
