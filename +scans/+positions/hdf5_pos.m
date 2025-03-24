%hdf5_pos loads ptycho scan positions from hdf5 files 
%Written by YJ

function [ p ] = hdf5_pos( p )

for ii = 1:p.numscans
    switch p.scan.type
        case 'default'
            pos_file = strcat(p.base_path,sprintf(p.scan.format, p.scan_number(ii)),'/data_roi',p.scan.roi_label,'_para.hdf5');   
            if exist(pos_file,'file')
                ppX = h5read(pos_file,'/ppX');
                ppY = h5read(pos_file,'/ppY');
                ppX = ppX(:);
                ppY = ppY(:);
                positions_real = zeros(length(ppX),2); 

                positions_real(:,1) = -ppY;
                positions_real(:,2) = -ppX;                
            else
                disp(strcat(p.base_path,sprintf(p.scan.format, p.scan_number(ii)),'/data_roi',p.scan.roi_label,'_para.hdf5'))
            	error('Could not find function or data file %s', pos_file);
            end
        case 'custom'
            % modified by dzh to match the data position file
            if isempty(p.scan.custom_positions_source) %guess the position file name from base path
                pos_file = strcat(p.base_path,sprintf(p.scan.format, p.scan_number(ii)),'/data_position.hdf5');
            end
            if exist(pos_file,'file')
                
            	ps = h5read(pos_file,'/probe_positions_0');
                ppX = ps(:,1); % Angstrom to meter
                ppY = ps(:,2);
                positions_real = zeros(length(ppX),2); 

                positions_real(:,1) = -ppY(:);
                positions_real(:,2) = -ppX(:);                
            else
            	error('Could not find function or data file %s', pos_file);
            end         
        otherwise
            error('Unknown scan type %s.', p.scan.type);
    end
    utils.verbose(2, strcat('Loaded scan positions from:', pos_file))
    p.numpts(ii) = size(positions_real,1);
    p.positions_real = [p.positions_real ; positions_real]; %append position
end
    
end

