% FIND_GEOM_CORRECTION use current probe positions estimates to update geometry model and
% improve the  new probe positions 
%
% [self] = find_geom_correction(self,cache, par, iter,best_mode_id)
% 
%
% ** self      structure containing inputs: e.g. current reconstruction results, data, mask, positions, pixel size, ..
% ** cache     structure with precalculated values to avoid unnecessary overhead
% ** par       structure containing parameters for the engines 
% ** iter      current iteration 
% ** best_mode_id   strongest mode id
%
% returns:
% ++ self        self-like structure with final reconstruction
%
%
function [self] = find_geom_correction(self,cache, par, iter,best_mode_id, update_position_weight)

    import engines.GPU_MS.GPU_wrapper.*
    import engines.GPU_MS.shared.*
    import utils.*
    import math.*
    
    mode = self.modes{best_mode_id};
   
    %% constrain the detector rotation 
     
    % store only the single update per scan 
    if iter > par.detector_rotation_search
        for ii = 1:length(self.reconstruct_ind)
            ind = self.reconstruct_ind{ii};
            mrot(ii) = mean(mode.probe_rotation_all(ind)); 
        end
        if par.mirror_objects
            % enforce zero average rotation if two mirror scans are provided
            mrot = mrot - mean(mrot);
        end
        for ii = 1:length(self.reconstruct_ind)
            mode.probe_rotation_all(ind) = mrot(ii) ;
        end
        mode.probe_rotation(end+1,:) = mrot;
    end
    
    
    if iter <= par.probe_position_search
        self.modes{best_mode_id} = mode;
        return
    end

    
    pos = mode.probe_positions;
    pos_0 = mode.probe_positions_0;
       
    
    if all(isnan(mode.probe_positions_weight(:))) || all(mode.probe_positions_weight(:)==0) || update_position_weight % by ZC
%     if update_position_weight %modified by YJ. add this option so users can update poisition weights more than once
        %% EMPIRICAL ESTIMATION OF POSITION RELIABILITY
        verbose(1,'EMPIRICAL ESTIMATION OF POSITION RELIABILITY\n')
        illum = utils.crop_pad(abs(self.probe{1}(:,:,1)).^2, self.Np_p/2); 
        total_variation = zeros(self.Npos,2, 'single'); 
        
        for ii = 1:par.Nscans
            best_layer = fix(par.Nlayers/2)+1; % change to middle by ZC
            o_tmp =  self.object{min(end,ii), best_layer}; 
            o_tmp = o_tmp ./ max2(abs(o_tmp(cache.object_ROI{:})));
            % keep it more memory effecient (important for GPU !! )
            Npos = length(self.reconstruct_ind{ii});
            for jj = 1:ceil(Npos/par.grouping)
                ind = 1+(jj-1)*par.grouping:min(Npos, jj*par.grouping);
                obj_proj = get_views(o_tmp,[],1,1,self.reconstruct_ind{ii}(ind),cache);
                obj_proj = utils.crop_pad(obj_proj, self.Np_p/2); 

                [nx, ny,~]  = size(obj_proj); 
                [X,Y] = meshgrid(-ny/2:ny/2-1, -nx/2:nx/2-1);
                % suppress edge effects of the FFT derivatives 
                spatial_filter = exp(-(X.^16+Y.^16)/(min(nx,ny)/2.2)^16);
                obj_proj = obj_proj.* spatial_filter;
                [dX, dY] = get_img_grad(obj_proj);
                clear obj_proj 
                illum_proj = get_views(utils.imgaussfilt2_fft(cache.illum_sum_0{min(ii,end)},self.Np_p/10),[],1,1,self.reconstruct_ind{ii}(ind),cache);
                illum_proj = utils.crop_pad(illum_proj, self.Np_p/2); 

                dX = abs(dX) .* illum_proj.* illum; 
                dY = abs(dY) .* illum_proj.* illum;
                clear illum_proj
                total_variation(self.reconstruct_ind{ii}(ind),:) = Ggather(sqrt(squeeze([mean2(dX),mean2(dY)]))');
                clear dX dY
            end
        end
        mode.probe_positions_weight = total_variation.^4./mean(total_variation.^4); 
    end
    % adjust weight close to 1 by ZC
    if par.probe_positions_weight_rescale ~= 1
        mode.probe_positions_weight=(mode.probe_positions_weight-1)/par.probe_positions_weight_rescale+1;
    end

    if par.position_weight_no  % no weighting but keep geom, by ZC
        mode.probe_positions_weight(:,:) = 1;
    end
    probe_positions_weight = double(mode.probe_positions_weight);
    jj = size(mode.scales,1)+1; 
        
    % find geometry for each scan separatelly: by dzh
    for ii = 1:par.Nscans
        ind = self.reconstruct_ind{ii};
        C0 = mode.affine_matrix(:,:,ii) - eye(2);
        C0 = C0(:);

        if par.Nscans > 1 && par.share_object 
            % it the case of multiple scans allow also freedom of coordinates shifts
            pos_fun = @(C)(( [1+C(1), C(2); C(3), 1+C(4)]*pos_0(ind,:)')' + C([5,6])' );
            if isfield(mode, 'shift_scans' ) && size(mode.shift_scans,2)>=ii
                C0(5:6) = mode.shift_scans(:,ii);
            else
                C0(5:6) = 0;
            end
        else
            pos_fun = @(C)(( [1+C(1), C(2); C(3), 1+C(4)]*pos_0(ind,:)')' );
        end
            
        err_fun = @(C)( probe_positions_weight(ind,:) .* (pos(ind,:) - pos_fun(C))); 

        options = optimoptions('lsqnonlin','Display','off');
        %C(:,ii) = lsqnonlin( err_fun, C0,[],[],options) ; 
        %% modified by YJ to avoid fitting error
        try 
            C(:,ii) = lsqnonlin( err_fun, C0,[],[],options);
        catch
            disp('Fitting error during geom correction...')
            C(:,ii) = [0,0,0,0,0,0];
        end
        %% restrict the geometry model only to the allowed degreed of freedom 
        % ===================================================================
        M{ii} = reshape(C(1:4,ii),2,2)+eye(2); 
        M_fit = M; %added by YJ to keep track of estimated affine matrix
        [scale, asymmetry, rotation, shear] = decompose_affine_matrix(M{ii}); 
        
        if ~ismember('scale', par.probe_geometry_model)
            scale = 1;
        end
        if ~ismember('asymmetry', par.probe_geometry_model)
            asymmetry = 0;
        end
        if ~ismember('rotation', par.probe_geometry_model)
            rotation = 0;
        end
        if ~ismember('shear', par.probe_geometry_model)
            shear = 0;
        end
        if ~ismember('shiftx', par.probe_geometry_model)
            shiftx = 0;
        end
        if ~ismember('shifty', par.probe_geometry_model)
            shifty = 0;
        end
        M{ii} = compose_affine_matrix(scale, asymmetry, rotation, shear);
        
        % ===================================================================
        mode.scales(jj,ii) =  scale;
        mode.asymmetry(jj,ii) = asymmetry;
        mode.rotation(jj,ii) =  rotation;
        mode.shear(jj,ii)    =  shear;

        
        if par.Nscans > 1 && par.share_object 
            mode.shift_scans(:,ii) =   C(5:6,ii);
            mode.shiftx(jj,ii)    =  C(5,ii);
            mode.shifty(jj,ii)    =  C(6,ii);
            verbose(0,sprintf('scan %d: global shift corrected to [%.4f, %.4f]',ii, C(5,ii), C(6,ii)));
        else
            mode.shift_scans(:,ii) =   [0,0];
            mode.shiftx(jj,ii)    =  0;
            mode.shifty(jj,ii)    =  0;
        end
        
        % store initial guess 
        mode.affine_matrix(:,:,ii) = M{ii}; 

        % calculate ideal model positions 
        pos_model(ind,:) = pos_fun([reshape(M{ii} - eye(2), [],1); mode.shift_scans(:,ii)]); 

    end

    self.affine_matrix = M;
    self.affine_matrix_fit = M_fit; %added by YJ to keep track of estimated affine matrix

    verbose(2,['-----  Geom. correction  ', repmat('%3.3g ', 1,length(C))], C)

    % use average 
    resid_pos= pos - pos_model; 
    
    % ignore errors in the global shift of the positions 
    for ii = 1:par.Nscans
        ind = self.reconstruct_ind{ii};
        resid_pos(ind,:)  = resid_pos(ind,:) - mean(resid_pos(ind,:)); 
    end
    
    err = abs(resid_pos); 

    max_err =  par.probe_position_error_max ./ self.pixel_size .* self.relative_pixel_scale; 

    verbose(1, '==== AVG position error %3.2g px MAX error %3.2g LIMIT %3.2g px ', mean(err(:)), max(err(:)),  max(max_err))
    
    %modified by YJ: add par.apply_relaxed_position_constraint to allow position update without constraints from geometry model.
    %Useful if there are big jumps in positions
    
%     if par.apply_relaxed_position_constraint 
    if par.apply_relaxed_position_constraint && iter >= par.geometry_model_Niter(1) && mod(iter-par.geometry_model_Niter(1), par.geometry_model_Niter(2)) == 0 % By ZC
        %% apply only relaxed constrain on the probe positions  !!! 
%         relax = 0.1;
        if isfield(par,'apply_relaxed_position_constraint_relax' ) % added relax factor by ZC
            relax= par.apply_relaxed_position_constraint_relax; 
        else 
            relax=0.1;
        end
        % constrain more the probes in flat regions 
        W = relax*(1-  (probe_positions_weight./ (1+probe_positions_weight)));  
        % penalize positions that are further than max_err from origin 
        W = min(10*relax, W+max(0,err - max_err).^2 ./ max_err.^2  );  % avoid travel larger than max error
    else
        W=0; %no geom model imposed to regularize positions
    end
   
   % allow free movement in depenence on realibility and max allowed error 
   % test only
%     t=reshape(W(:,1),64,[]);
%     figure;imagesc(t);axis image;
    %%%
    
    self.modes{best_mode_id} = mode;
    
    % added by dzh: correct tilt angle in TCMEP reconstructions
    if par.precession_mode && par.correct_angle && ismember(iter, par.iter_correct_angle)
        layer4pos = par.layer4pos;
        tilt_angle_temp = zeros(1, par.Nscans*2);
        for iscan = 1:par.Nscans
            % ipos saves the average position update for each layer in layer4pos
            ipos = zeros(length(layer4pos),2);
            for ilayer = 1:length(layer4pos)
                pos_upd = self.modes{layer4pos(ilayer)}.probe_positions_update{1}(self.reconstruct_ind{iscan},:);
                ipos(ilayer,:) = mean(pos_upd);
            end

            % linear fit the ipos to find the tilt angle update
            xcoord = 0:(length(layer4pos)-1);
            fitx = polyfit(xcoord,ipos(:,1),1);
            fity = polyfit(xcoord,ipos(:,2),1);
            intershift_x = fitx(1);
            intershift_y = fity(1);
            layermiddle = round(mean(layer4pos));

            % update tilt angle: dzh
            ipos_angle = zeros(par.Nlayers,2);
            for ilayer = 1:par.Nlayers
                self.modes{ilayer}.probe_positions(self.reconstruct_ind{iscan},:) = self.modes{ilayer}.probe_positions(self.reconstruct_ind{iscan},:) ...
                + self.modes{layermiddle}.probe_positions_update{1}(self.reconstruct_ind{iscan},:) + [intershift_x intershift_y]*(ilayer-layermiddle);
                % [intershift_x intershift_y]*(ilayer-layermiddle + imiddle - 1) - ipos(imiddle,:) + [shift_x shift_y];
                
                % save the average position in this layer for fitting
                ipos_angle(ilayer,:) = mean(self.modes{ilayer}.probe_positions(self.reconstruct_ind{iscan},:));
            end

            % verbose tilt angle in this iteration: dzh
            xcoord_angle = 1:par.Nlayers;
            fitx = polyfit(xcoord_angle,ipos_angle(:,1),1);
            fity = polyfit(xcoord_angle,ipos_angle(:,2),1);
            tilt_x_pixel = fitx(1);
            tilt_y_pixel = fity(1);
            angle_temp = [tilt_x_pixel tilt_y_pixel]/self.z_distance(1)/pi*180.*self.pixel_size;
            tilt_angle_temp((iscan-1)*2+1:iscan*2) = angle_temp;
            verbose(0,sprintf('scan %d: tilt angle corrected to [%.4f, %.4f]',iscan,angle_temp(1),angle_temp(2)))
        end
        
        % save tilt angle correction result
        if ~isfield(self, 'tilt_angle')
            self.tilt_angle = [];
        end
        self.tilt_angle = [self.tilt_angle; tilt_angle_temp];

        if any(isnan(mode.probe_positions(:)))
            keyboard
        end
        return
    elseif par.precession_mode && (~par.correct_angle || ~ismember(iter, par.iter_correct_angle))
        for ilayer = 1:par.Nlayers
            self.modes{ilayer}.probe_positions=self.modes{ilayer}.probe_positions+self.modes{ilayer}.probe_positions_update{1}; % - mean(self.modes{ilayer}.probe_positions_update{1});
        end
        upd_global_shift = zeros(par.Nscans, 2);
        for iscan = 1:par.Nscans
            % dzh: save the temporary shift to shift back the probes
            upd_global_shift(iscan,:) = mean(self.modes{best_mode_id}.probe_positions_update{1}(self.reconstruct_ind{iscan},:));
        end
        self.upd_global_shift = upd_global_shift;
        if any(isnan(mode.probe_positions(:)))
            keyboard
        end
        return
    elseif ~par.precession_mode
    
        layer4pos = par.layer4pos;
        pos = self.modes{layer4pos}.probe_positions;
        pos_new =  pos .*(1-W)+W.*pos_model;
        mode.probe_positions = pos_new;
        mode.probe_positions_model = pos_model;
        if any(isnan(mode.probe_positions(:)))
                keyboard
        end
        for ilayer = 1:size(self.modes,1)
            self.modes{ilayer}.probe_positions = pos_new;
        end
        global_shift = mean(self.modes{1}.probe_positions)-mean(self.modes{1}.probe_positions_0);
        verbose(0,sprintf('global shift corrected to [%.4f, %.4f]',global_shift(1),global_shift(2)))
        
    end
    self.modes{best_mode_id} = mode;
end

