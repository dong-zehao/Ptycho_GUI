function imshow(axis, img, clb)
    if nargin < 3
        clb = 'parula';
    end
    imagesc(axis, img); 
    xlim(axis, 'tight');
    ylim(axis, 'tight')
    axis.DataAspectRatio = [1 1 1];
    axis.XTick = [];
    axis.YTick = [];
    axis.XColor = 'none';  
    axis.YColor = 'none';  
    axis.Box = 'off';
    colormap(axis, clb);
end