function PlotInitialState(dynamicFcn,t,X,odeOptions)
for ii = 1:size(X,1)-1
    if nargin==4
        [~,plotX] = ode113(dynamicFcn,[t(ii),t(ii+1)],X(ii,:),odeOptions);
    else
        [~,plotX] = ode113(dynamicFcn,[t(ii),t(ii+1)],X(ii,:));
    end
    plot3(plotX(:,1),plotX(:,2),plotX(:,3),'-'); hold on;
    if ii == 1
        plot3(plotX(1,1),plotX(1,2),plotX(1,3),'b*'); hold on;
        %         disp([t(ii),plotX(1,1:6)]);
    elseif ii == size(X,1)-1
        plot3(plotX(end,1),plotX(end,2),plotX(end,3),'ro'); hold on;
        %         disp([t(end), plotX(end,1:6)]);
    end
end
axis equal;
xlabel('x');
ylabel('y');
drawnow update;
end
