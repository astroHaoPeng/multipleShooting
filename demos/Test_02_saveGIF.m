movietimer = tic;

% MyFigure = figure(7);
% scrsz = get(0,'ScreenSize');
% set(MyFigure,'OuterPosition',[200, 100, scrsz(4)/1.2 ,scrsz(4)/1.2], 'MenuBar', 'none');

M = 121;

a = linspace( 0, 360, M );
b = [linspace( 5, 30, round(M/2) ), linspace( 30, 0, M-round(M/2) )];

set(gca,'FontSize',15);

kk = 1;
clear figX;
for ii = 1:M % 画三维
    
    %
    % 插入画图代码
    %
    set(gca,'view',[a(ii),b(ii)]);
    
    GifData = getframe(gcf);
    if ii == 1
        [figX(:,:,1,kk),Map] = rgb2ind(GifData.cdata, 256);
    else
        figX(:,:,1,kk) = rgb2ind(GifData.cdata, Map);
    end
    kk = kk + 1;
    
end
imwrite(figX,Map,'Test_02_VL_Ephemeris.gif','LoopCount',Inf,'DelayTime',0);

toc(movietimer);