function [cbed, parameters] = read_h5_arina(filename, ind_x, ind_y, totx, toty)
    % read a smaller region from large dataset on arina, with a faster speed
    % 
    % written by Zehao Dong

    % arguments: 
    % filename: full directory towards the master .h5 file
    % (optional) ind_x: selected data index x, e.g. [1:100]
    % (optional) ind_y: selected data index y, e.g. [1:100]

    % outputs:
    % cbed: the 4D-STEM data
    % patameters: all detector patameters read from h5 file

% 根据master HDF5文件结构，可以看到它包含以下主要组：
% 
% % /entry 组：
% % 
% % % 这是数据的入口点，具有NX_class属性为'NXentry'。
% % % 包含一个名为'definition'的数据集，该数据集存储一个字符串，表示版本信息。
% % % 包含一个名为'data'的组，其中存储了数据集链接，指向另外的HDF5文件。
% % /entry/data 组：
% % 
% % % 这是数据的主要存储位置，具有NX_class属性为'NXdata'。
% % % 包含多个外部链接（external link），这些链接指向其他HDF5文件中的数据集，例如'data_000001'、'data_000002'等等。
% % % 这些链接指向的数据集包含实际的数据。
% % /entry/instrument 组：
% % 
% % % 这是与仪器相关的信息的存储位置，具有NX_class属性为'NXinstrument'。
% % % 包含子组/entry/instrument/beam和/entry/instrument/detector，分别存储了与光束和探测器相关的信息。
% % /entry/instrument/detector 组：
% % 
% % % 包含有关探测器的详细信息。
% % % 包含了探测器的几何信息、数据处理参数（如像素大小、计数时间等）以及探测器相关的属性。
% % % 包含一个子组/entry/instrument/detector/detectorSpecific，其中存储了特定于探测器的信息，如自动求和、压缩、校正等。
% % /entry/instrument/detector/geometry 组：
% % 
% % % 存储了探测器的几何信息，包括位置和方向。
% % /entry/sample 组：
% % 
% % % 存储了与样品相关的信息，具有NX_class属性为'NXsample'。
% % % 包含了样品的旋转信息等。

    % filename = 'G:\Data_Arina\15\Test_STO_good_dataset_015_master.h5';
    % ind_x = [1:100];
    % ind_y = [1:100];
    
    % get the chunk size for each slave dataset
    [dp_size, chunk_size] = get_chunk_size(filename);
    mod_chunks = max(chunk_size); %假设每个分组大小都是相同的，可能会出bug？

    if nargin < 5
        att_filename = strrep(filename, '_master', '');
        parameters = get_h5_att_all(att_filename, '/STEM Metadata');
        initial_size = [parameters.NumPixelsX parameters.NumPixelsY];
        detector_info = get_h5_data_all(filename, '/entry/instrument/detector');
        parameters.detector_info = detector_info;
        parameters.ScanDwellTime = detector_info.count_time;
        parameters.dp_size = dp_size;
        parameters.chunk_size = chunk_size;
    end

    if nargin == 5
        initial_size = [totx toty];
    end

    if nargin < 2
        ind_x = 1:initial_size(1);
        ind_y = 1:initial_size(2);
    end

    % get the index for the selected region
    [ind_xx, ind_yy] = meshgrid(ind_x, ind_y);
    image_index = sub2ind(initial_size, ind_xx, ind_yy);
    

    data_group_info = h5info(filename, '/entry/data');
    parent_folder = fileparts(filename);
    cbed = zeros(dp_size(1), dp_size(2),size(ind_y,2),size(ind_x,2), 'single');

    for ix = 1:size(image_index,1)
        if mod(ix, size(image_index,1)/20) == 0
            disp(ix/size(image_index,1));
        end

        %fprintf('%d/%d\n', ix, size(image_index,1))
        %由于image_index中每一行必然是连续的，所以只要对每一行的数据分别读取即可
        ind_file = ceil(image_index(ix, :)/mod_chunks); %利用坐标取模来得知这一行需要读取哪一个h5文件
        ind_read = mod(image_index(ix, :),mod_chunks);
        ind_uniq = unique(ind_file);

        if length(ind_uniq) == 1
            link_path = data_group_info.Links(ind_uniq).Value{1};
            slave_data = h5read(fullfile(parent_folder, link_path), '/entry/data/data', [1 1 ind_read(1)], [dp_size length(ind_read)]);
            cbed(:,:,ix,:) = slave_data;
        elseif length(ind_uniq) == 2
            ind_file1 = (ind_file == ind_uniq(1));
            ind_file2 = (ind_file == ind_uniq(2));
            i_to_read1 = ind_read(ind_file1);
            i_to_read2 = ind_read(ind_file2);
            link_path1 = data_group_info.Links(ind_uniq(1)).Value{1};
            link_path2 = data_group_info.Links(ind_uniq(2)).Value{1};
            slave_data1 = h5read(fullfile(parent_folder, link_path1), '/entry/data/data', [1 1 i_to_read1(1)], [dp_size length(i_to_read1)]);
            slave_data2 = h5read(fullfile(parent_folder, link_path2), '/entry/data/data', [1 1 i_to_read2(1)], [dp_size length(i_to_read2)]);
            cbed(:,:,ix,:) = cat(3,slave_data1,slave_data2);
        else
            error('undefined behaviour: some row corresponds to more than 2 files!')
        end
    end
end
    
    


function [dp_size, chunk_size] = get_chunk_size(filename)
    dp_size = h5read(filename, '/entry/instrument/detector/module/data_size');
    dp_size = double(dp_size');

    chunk_size = [];
    parent_folder = fileparts(filename);
    data_group_info = h5info(filename, '/entry/data');
    
    % 获取链接文件的数量
    num_links = numel(data_group_info.Links);
    for i = 1:num_links
        % 获取链接文件路径
        link_path = data_group_info.Links(i).Value{1};
        slave_data_info = h5info(fullfile(parent_folder, link_path), '/entry/data/data');

        % 提取Size属性值
        chunk_size = [chunk_size slave_data_info.Dataspace.Size(3)];
    end
end


function res = get_h5_data_all(filename, path)
    res = struct();

    % 获取指定组中的所有数据集信息
    info = h5info(filename, path);
    
    % 获取数据集的数量
    num_datasets = numel(info.Datasets);
    
    % 逐个读取每个数据集的数据
    for i = 1:num_datasets
        dataset_name = info.Datasets(i).Name;
        dataset_path = [path, '/', dataset_name];
        
        % 读取数据集的数据
        data = h5read(filename, dataset_path);
        
        % 处理数据，例如打印数据集的大小和名称
        res.(dataset_name) = data;
    end
end

function res = get_h5_att_all(filename, path)
    res = struct();

    % 获取指定组中的所有数据集信息
    info = h5info(filename, path);
    
    % 获取数据集的数量
    num_datasets = numel(info.Attributes);
    
    % 逐个读取每个数据集的数据
    for i = 1:num_datasets
        dataset_name = info.Attributes(i).Name;

        % 读取数据集的数据
        data = h5readatt(filename, path, dataset_name);
        dataset_name = strrep(dataset_name, ' ', '_');
        % 处理数据，例如打印数据集的大小和名称
        res.(dataset_name) = data;
    end
end