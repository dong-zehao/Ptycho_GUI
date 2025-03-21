function OnPrepare(app)
    % This script prepares experimental electron ptycho. data for PtychoShelves
    app.UITable_3.Data = [];
    cd(app.PtychopathEditField.Value);
    init_dir = fullfile(app.DatapathEditField.Value, app.DfolderpatternEditField.Value);
    filename = app.FilenameEditField.Value;
    fileList = dir(fullfile(init_dir, filename));
    
    for i_run = 1:size(app.UITable_2.Data,1)
        index = app.UITable_2.Data(i_run,1);
        % read data
        data_dir = fullfile(fileparts(init_dir),'All_Data');
        if ~exist(data_dir)
            mkdir(data_dir)
        end 
        fname_here = fullfile(fileList(index).folder, fileList(index).name);
        LogMessage(app, fname_here);
        [dp, rot_gatan] = ReadDM4File_4D(fname_here);

        % update UITable for preparation
        app.UITable_3.Data = [app.UITable_3.Data; index size(dp, [3 4])];

        % how much to bin?
        Np_binto = app.Np_bintoSpinner.Value;
        k3_rotation = 50.5 + 180 + rot_gatan;
        binsize = size(dp,1)/Np_binto;
        cbed = bin_cbed(dp, binsize);
        switch app.RotateCBEDSwitch.Value
            case 'Yes'
                cbed = imrotate(cbed,k3_rotation,'bilinear');
            case 'No'

            otherwise
                LogMessage(app, 'RotateCBED cannot determine')
        end
        clear dp
        
        % how much to pad?
        Np_padto = app.Np_padtoSpinner.Value+mod(size(cbed,1),2); 
        % Np_padto = size(cbed,1);
    
        % generate mask
        [~, maskname] = mygen_mask(app.PtychopathEditField.Value, Np_padto, Np_binto, k3_rotation, 'crop');
    
    % Step 2: prepare data
        scan_number = index;
        crop_idx0 = [2 size(cbed,3) 2 size(cbed,4)];
        ADU = 1;
        alpha = app.alphaEditField.Value;
        rbf = app.rbfEditField.Value/binsize;
        voltage= app.AccVoltEditField.Value;
        
        df=app.UITable_2.Data(i_run, 2);
    
        Ninter=1;           % skip some scan points
        icenter = 1;        % recenter disk
    
        % dk = 0.086327 ./ 10;
    
        % Step 3: go back to .../fold_slice/ptycho and pre-process data
        Np_p=[Np_padto,Np_padto];
        cbed=cbed(:,:,crop_idx0(1):crop_idx0(2),crop_idx0(3):crop_idx0(4));
        
        if exist('Ninter','var') && Ninter > 1
            cbed=cbed(:,:,1:Ninter:end,1:Ninter:end);
        end
        
        % recenter
        if icenter==1
            pacbed=sum(cbed, [3 4]);%sum(cbed(:,:,ind_for_center(1):ind_for_center(2),ind_for_center(3):ind_for_center(4)),[3 4]);
            [sftx0,sfty0]=shiftmeasure(pacbed,0.5);
            cbed=shift_pix(cbed,[sfty0,sftx0]);
        else
            sftx0=0;
            sfty0=0;
        end
        [ndpy,ndpx,npy,npx]=size(cbed);
        
        if ndpy < Np_p(1) % pad zeros
            cbed=padarray(cbed,[(Np_p(1)-ndpy)/2,(Np_p(2)-ndpx)/2,0,0],0,'both');
        else
            cbed=crop_pad(cbed,Np_p);
        end
        
        %dp_max = unique(max(max(cbed(:,:,1,1))));
    
        %cbed = cbed / dp_max / 2;
        cbed = cbed / ADU; % convert to electron count
        %cbed(cbed<0)=0; % set background,
        cbed=reshape(cbed,Np_p(1),Np_p(2),[]);
        Itot=mean(squeeze(sum(sum(cbed,1),2))); %need this for normalizting initial probe
        
        % calculate pxiel size (1/A) in diffraction plane
        [~,lambda]=electronwavelength(voltage);
        dk = alpha/1e3./rbf/lambda; %%% PtychoShelves script needs this %%%
        % dk = 0.086327 ./ 10;
        
        % Step 4: save CBED in a .hdf5 file (needed by Ptychoshelves)
        % scan_number = 3; %Ptychoshelves needs
       
        
        save_dir = fullfile(data_dir,num2str(scan_number,'%02d'), '/');
        if exist('save_dir','dir')
            error('folder exists, change folder name')
        else
            mkdir(save_dir)
        end
        imagesc(cbed(:,:,1)); colorbar; axis equal;
        roi_label = strcat('0_Ndp',num2str(Np_p(1)));
        
        switch app.SavedpSwitch.Value
            case 'Yes'
                saveName = strcat('data_roi',roi_label,'_dp.hdf5');
                h5create(fullfile(save_dir,saveName), '/dp', size(cbed),'ChunkSize',[size(cbed,1), size(cbed,2), 1],'Deflate',4,'Datatype','single')
                h5write(fullfile(save_dir,saveName), '/dp', cbed)
            case 'No'
            otherwise
                LogMessage(app, 'Save dp cannot determine')
        end
    
        % Step 5: prepare initial probe
        dx=1./Np_p./dk; %% pixel size in real space (angstrom)
        cs = 0;
        probe=generateProbeFunction2(dx,Np_p,0,0,df,cs,voltage,alpha,0);
        probe=probe/sqrt(sum(sum(abs(probe.^2))))*sqrt(Itot)/sqrt(Np_p(1)*Np_p(2));
        probe=single(probe);
        % add parameters for PtychoShelves_electron
        p = {};
        p.binning = false;
        p.detector.binning = false;
        
        % Step 6: save initial probe
        save(fullfile(save_dir,'/init_probe.mat'),'probe','p')
        clear cbed;
        save(fullfile(save_dir,'/params_backup.mat'),'ADU','alpha','rbf','dx','voltage',...
            'df','Np_p','dk','cs','sft*', 'crop_idx0', 'maskname', 'binsize');
        copyfile(strcat(mfilename('fullpath'),'.m'),save_dir);
    end
end

function [mask, maskname] = mygen_mask(path, in, out, rot, methods)
    if nargin < 5
        methods = 'loose';
    end
    mask = imrotate(gen_mask(in, out), rot, "nearest", methods);
    maskname = fullfile(path, ['det_mask' num2str(in) '_outside' num2str(out) '_rot' num2str(rot) '.mat']);
    save(maskname, 'mask');
end