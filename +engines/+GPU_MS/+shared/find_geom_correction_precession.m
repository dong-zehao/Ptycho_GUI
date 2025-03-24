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
function [self] = find_geom_correction_precession(self,cache, par, iter,best_mode_id, update_position_weight)

    import engines.GPU_MS.GPU_wrapper.*
    import engines.GPU_MS.shared.*
    import utils.*
    import math.*
    
    struc_fit_res = struct();
    struc_fit_res.scales = zeros(length(par.layer4pos),par.Nscans) ;
    struc_fit_res.asymmetry = zeros(length(par.layer4pos),par.Nscans);
    struc_fit_res.rotation =  zeros(length(par.layer4pos),par.Nscans);
    struc_fit_res.shear    =  zeros(length(par.layer4pos),par.Nscans);
    struc_fit_res.shift_scansx =   zeros(length(par.layer4pos),par.Nscans);
    struc_fit_res.shift_scansy =   zeros(length(par.layer4pos),par.Nscans);

    for ilayer = par.layer4pos
        mode = self.modes{ilayer}; 
%{
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
            self.modes{ilayer} = mode;
            return
        end
    
%}      
        pos = mode.probe_positions;
        pos_0 = mode.probe_positions_0;
           
%{
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

%}
        probe_positions_weight = double(mode.probe_positions_weight);
        jj = size(mode.scales,1)+1; 

        % find geometry for each scan separatelly 
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
            M{ii} = compose_affine_matrix(scale, asymmetry, rotation, shear);
            
            % ===================================================================
            struc_fit_res.scales(ilayer,ii) = scale ;
            struc_fit_res.asymmetry(ilayer,ii) = asymmetry;
            struc_fit_res.rotation(ilayer,ii) =  rotation;
            struc_fit_res.shear(ilayer,ii)    =  shear;

            if par.Nscans > 1 && par.share_object 
                struc_fit_res.shift_scansx(ilayer,ii) =   C(5,ii);
                struc_fit_res.shift_scansy(ilayer,ii) =   C(6,ii);
            else
                struc_fit_res.shift_scansx(ilayer,ii) =   0;
                struc_fit_res.shift_scansy(ilayer,ii) =   0;
            end
        end
    end
            
        scale = mean(struc_fit_res.scales);
        asymmetry = mean(struc_fit_res.asymmetry);
        rotation = mean(struc_fit_res.rotation);
        shear = mean(struc_fit_res.shear);
        inter_shift_scansx = mean(diff(struc_fit_res.shift_scansx));
        inter_shift_scansy = mean(diff(struc_fit_res.shift_scansy));
        start_shift_scansx = (struc_fit_res.shift_scansx());
        start_shift_scansy = (struc_fit_res.shift_scansy);

        for ilayer = 1:par.Nlayers
            mode = self.modes{ilayer};
            mode.scales(jj) =  scale;
            mode.asymmetry(jj) = asymmetry;
            mode.rotation(jj) =  rotation;
            mode.shear(jj)    =  shear;
            
            if par.Nscans > 1 && par.share_object 
                shift_scanx = middle_shift_scansx - 
                shift_scany
                mode.shift_scans(:) = [shift_scanx; shift_scany];
            else
                mode.shift_scans(:,ii) =   [0,0];
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
        pos_new =  pos .*(1-W)+W.*pos_model;
        
        mode.probe_positions = pos_new;
        mode.probe_positions_model = pos_model;
    
        if any(isnan(mode.probe_positions(:)))
            keyboard
        end
        
        self.modes{ilayer} = mode;
    % par.precession_vacuum_thickness
end

