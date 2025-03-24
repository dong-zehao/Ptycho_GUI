function plot_object(app, object, cache)
    object_rot = zeros(size(object{1},1),size(object{1},2),size(object,2));
    for ilayer = 1:size(object,2)
        object_rot(:,:,ilayer) = object{ilayer};
    end
    object_rot = imrotate(object_rot, cache.init_rot, 'bilinear', 'crop');
    object_plot = angle(object_rot(cache.object_ROI{:},app.LayertoplotSpinner.Value));
    utils.imshow(app.UIAxes_2, object_plot, 'gray');
end