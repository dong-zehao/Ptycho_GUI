%LSQML least squared maximum likelihood solver 
%
%[self, cache, fourier_error] =  LSQML(self,par,cache,fourier_error,iter)
%
% ** self      structure containing inputs: e.g. current reconstruction results, data, mask, positions, pixel size, ..
% ** par       structure containing parameters for the engines 
% ** cache     structure with precalculated values to avoid unnecessary overhead
% ** fourier_error  array [Npos,1] containing evolution of reconstruction error 
% ** iter      current iteration number 
% returns:
% ++ self        self-like structure with final reconstruction
% ++ cache     structure with precalculated values to avoid unnecessary overhead
% ++ fourier_error  array [Npos,1] containing evolution of reconstruction error 
%

function [self, cache, fourier_error] =  LSQML(self,par,cache,fourier_error,iter)
    import engines.GPU_MS.shared.*
    import math.*
    import utils.*
    import plotting.*
    import engines.GPU_MS.GPU_wrapper.*
    import engines.GPU_MS.LSQML.*

    assert( ~(par.Nscans > 1 &&  par.object_modes > 1), 'Multiobject + multiscan not supported')

    % define some useful variables 
    object_modes = size(self.object,1);
    object_upd_sum = cell(object_modes,par.Nlayers);
    obj_illum_sum = cell(object_modes,par.Nlayers);
    obj_proj = cell(par.object_modes,1);
    apply_subpx_shift(); % reset persistent values 

    for ll = 1:object_modes
        for layer = 1:par.Nlayers
            obj_illum_sum{ll,layer} = Gzeros(self.Np_o);
            object_upd_sum{ll,layer} = Gzeros(self.Np_o, true);
        end
    end

    probe_update_sum = Gzeros(size(self.probe{1})); 
    probe_amp_corr = [0,0]; 
    %%%%%%%%%%%%%%%%%% LSQ-ML algorithm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
    beta_probe = ones(self.Npos,par.Nlayers); 
    beta_object = ones(self.Npos,par.Nlayers); 

   %for ML is more useful to get close/overlapping positions 
   % use already precalculated indices 
    if is_method(par, 'MLs')
        % load sparse indices
        rand_ind = randi(length( cache.preloaded_indices_sparse));
        indices = cache.preloaded_indices_sparse{rand_ind}.indices;
        scan_ids = cache.preloaded_indices_sparse{rand_ind}.scan_ids;
    else % load compact indices 
        indices = cache.preloaded_indices_compact{1}.indices;
        scan_ids = cache.preloaded_indices_compact{1}.scan_ids;
    end
    Nind = length(indices);

    for jj = 1:Nind
         layer_ids{jj} = 1:par.Nlayers; 
    end
    
    
    % shuffle order but keep same over iterations 
    if is_method(par, 'MLs')
        % shuffle order
        ind_range = randperm(Nind);
    else
        % MLc, call groups in given order, they are sorted by size to make
        % execution on GPU more effecient + stable convergence 
        ind_range = 1:Nind;
    end
    
   N_amp_out=0; % count the No of times (group*layers) whose amplitude outside the range defined.
   min_amp=inf; % minimum ampltide 
   max_amp=-inf; % maximum ampltide 
   %% apply updated in parallel over sets indices{jj}
   for  jj = ind_range
       % list of positions solved in current subiteration 
        g_ind = indices{jj}; 
        
        for ll = 1:max(par.probe_modes, par.object_modes) 
            % generate indices of the used probes 
             % single probe only
             if par.share_probe  % share incoherent modes 
                p_ind{ll} = 1;
             else
                 if all(scan_ids{jj} == scan_ids{jj}(1))
                    p_ind{ll} =  scan_ids{jj}(1);
                 else
                    % scan positions from multiple scans are processed in a single bunch 
                    p_ind{ll} =  scan_ids{jj};
                 end
             end
        end
        
        %% suppress large or small amplitude of object, quick try by ZC, inner loop
        if iter >= par.object_change_start && par.Nlayers > 1  && (par.amplitude_threshold_object < inf || par.amplitude_threshold_object_lower > 0 )
            for ll=1:par.Nlayers
                temp = abs(self.object{ll});
                if any(temp(:) > par.amplitude_threshold_object) || any(temp(:) < par.amplitude_threshold_object_lower)
                    N_amp_out=N_amp_out+1; % counting 
                    min_amp=min([min_amp,min(temp(:))]);
                    max_amp=max([max_amp,max(temp(:))]);
%                     if iter > 10 && jj == ind_range(end) % beginning iterations not monitor
%                         warning(strcat('amplitude range=[',num2str(min(temp(:)),'%5.2f'),', ',num2str(max(temp(:)),'%5.2f'),...
%                     '], force outliers to 1.0 (inner loop)!'));
%                     end
                    temp (temp > par.amplitude_threshold_object) =1;
                    temp (temp < par.amplitude_threshold_object_lower) = 1;
                    self.object{ll} = temp.* exp(1i* angle(self.object{ll}));                    
                end
            end
        end
        N_phs_out=0; % count the No of times (group*layers) whose amplitude outside the range defined.
        min_phs=inf; % minimum phase 
        max_phs=-inf; % maximum phase 
        %% suppress large or small phase of object, quick try by ZC, inner loop, be careful to use
        if iter >= par.object_change_start && par.Nlayers > 1  && (par.phase_threshold_object < inf || par.phase_threshold_object_lower > -inf )
            for ll=1:par.Nlayers
                temp = angle(self.object{ll});
                temp(isnan(temp))=0;
                temp(isinf(temp))=0;
                temp = temp - mean(temp(:),'omitnan'); % zero mean                
                if any(temp(:) > par.phase_threshold_object) || any(temp(:) < par.phase_threshold_object_lower)
                    N_phs_out=N_phs_out+1; % counting
                    min_phs=min([min_phs,min(temp(:))]);
                    max_phs=max([max_phs,max(temp(:))]);
                    temp (temp > par.phase_threshold_object) =0;
                    temp (temp < par.phase_threshold_object_lower) = 0;
                    self.object{ll} = abs(self.object{ll}).* exp(1i*temp);                    
                end
            end
        end
        
        % estimate forward model, ie wavefront behind the sample 
        clear probe
        [self, probe, obj_proj, psi] = get_forward_model(self, obj_proj, par,cache, g_ind, p_ind, scan_ids{jj}, layer_ids{jj});
        
        %% load data to GPU 
        modF = get_modulus(self, cache, g_ind,true,jj);
        mask = get_mask(self, cache.mask_indices, g_ind, par.damped_mask);
        noise = get_noise(self, par, g_ind); 

        drawnow;

        % get intensity (modulus) on detector including different corrections
        [aPsi, ~, cache, self] = get_reciprocal_model(self, psi, modF, mask,iter, g_ind, par, cache);

        %%%%%%%%%%%%%%%%%%%%%%% LINEAR MODEL CORRECTIONS END %%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%            
        if iter > 0 && (verbose >= -1 || par.number_iterations > par.plot_results_every) && ... 
                ((mod(iter,min(20, 2^(floor(2+iter/50)))) == 0 || iter < 10) ...
              || verbose()> -1) || any(iter == max(1,[1,par.number_iterations])) %  || par.accelerated_gradients_start-1 >= iter
            % calculate only sometimes to make it faster
            [fourier_error(iter,g_ind)] = get_fourier_error(modF, aPsi, noise,mask, par.likelihood);
            if any(~isfinite(fourier_error(iter,g_ind)))
                %disp(any(~isfinite(fourier_error(iter,g_ind))))

                if par.accelerated_gradients_start < par.number_iterations || par.momentum > 0
                    keyboard
                    error('Convergence failed, error contains NaNs, quitting ... \n%s',  'If repeated, try to set eng.accelerated_gradients_start = inf or eng.momentum = 0;  ')
                else 
                     keyboard
                    error('Convergence failed, error contains NaNs, quitting ... \n%s',  'If repeated, try to set eng.beta_LSQ = 0.5 or less ')
                end
            end
        end
        
        if strcmp(par.likelihood, 'poisson')   || iter == 0
            [chi,R] = modulus_constraint(modF,aPsi,psi, mask, noise, par, 1);
            if iter == 0
                % in the first iteration only find optimal scale for the probe
                probe_amp_corr(1) = probe_amp_corr(1) + Ggather(sum(modF(:).^2));
                probe_amp_corr(2) = probe_amp_corr(2) + Ggather(sum(aPsi(:).^2));
                continue
            end
        else
            chi = modulus_constraint(modF,aPsi,psi, mask, noise, par, 1);
        end
        
        if ~strcmp(par.likelihood, 'poisson')    % soft memory cleanup 
            mask = []; aPsi = []; noise = []; modF = [];  R=[];
        end
        
        if iter > par.estimate_NF_distance
            % update estimation of the nearfield propagation distance 
            [self, cache] = gradient_NF_propagation_solver(self,psi(:,end),chi, cache, g_ind);
        end
        if ~strcmp(par.likelihood, 'poisson')    % soft memory cleanup 
           psi = [];
        end
        
        if strcmp(par.likelihood, 'poisson')    % calculate only for the first mode
            %ll = 1;
            %% automatically find  optimal step-size, note that for Gauss it is 1 !! 
            aPsi2 = sumsq_cell(Psi);
            beta_xi  =  gradient_descent_xi_solver(self,modF, aPsi2, R,mask, g_ind, mean(cache.beta_xi(g_ind)), cache);
            R = [];
            for ll = 1:max(par.probe_modes, par.object_modes)
                chi{ll} = chi{ll} .*  beta_xi ;
            end
        end
        
        % At this point, chi is at the far-field (detector) plane
        % Now propagate chi back to the last object layer 
        % This is same as using back_fourier_proj w. distance = inf and no camera angle refinement 
        for ll= 1:max(par.object_modes, par.probe_modes)
        	chi{ll} = ifft2_safe(chi{ll});  % fully farfield inverse fft
        end
        
        %% %%%%%%%%%% LSQ optimization code (probe & object updates)%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for layer = par.Nlayers:-1:1
            for ll = 1:max(par.probe_modes, par.object_modes)
                object_reconstruct = iter >= par.object_change_start && (par.apply_multimodal_update || is_used(par, 'fly_scan') || ll <= par.object_modes );
                probe_reconstruct = (iter >= par.probe_change_start && iter < par.probe_change_end);

                llo = min(par.object_modes, ll);
                llp = min(par.probe_modes, ll);
                
                if layer < par.Nlayers
                    % modified by dzh: take in interlayer subpixel shift into account (only for the first mode, why?)
                    if (ll == 1 && par.apply_subpix_shift && isinf(self.z_distance(end)))  || is_used(par,'fly_scan')
                        chi{ll} = apply_subpx_shift(chi{ll} , -self.modes{layer+1}.sub_px_shift(g_ind,:)+self.modes{layer}.sub_px_shift(g_ind,:));
                    end
                    chi{ll} = back_fourier_proj(chi{ll}, self.modes{layer},g_ind);
                end
                
                if layer ~= par.Nlayers 
                    % if only single layer is used, reuse obj_proj already
                    % loaded, but avoid storing obj_proj for each layer, rather load it again 
                    obj_proj{llo} = get_views(self.object, obj_proj{llo},layer_ids{jj}(layer),llo, g_ind, cache, scan_ids{jj},[]);
                end

                % get update directions for each scan positions 
                if ( probe_reconstruct || layer > 1) && object_reconstruct
                    [probe_update, object_update_proj] = Gfun(@get_update_both, chi{ll}, obj_proj{llo}, probe{llp,layer});
                elseif ( probe_reconstruct || layer > 1)
                    probe_update       = Gfun(@get_update, chi{ll}, obj_proj{llo});
                    object_update_proj = 0;
                else
                    probe_update       = 0; m_probe_update = 0;
                    object_update_proj = Gfun(@get_update, chi{ll}, probe{llp,layer});
                end
                
                % debug by ZC
%                  if any(isnan(object_update_proj(:))) || any(isinf(object_update_proj(:))) 
%                      keyboard
%                      error('Convergence failed')
%                  end
                 % end debug

                % refine single optimal probe update direction (use overlap constraint)
                if probe_reconstruct || layer > 1
                   [self,m_probe_update, probe_update, cache]  = refine_probe_update(self, obj_proj{llo}, probe_update, chi{ll},layer,ll,p_ind{ll},g_ind, par, cache);
                end
                if layer == 1
                   probe_update = [] ; % soft memory clean 
                end

                % refine single optimal object update direction (use overlap constraint)
                if object_reconstruct                    
                    [object_upd_sum,object_update_proj, cache] = refine_object_update(self,  ...
                        object_update_proj,object_upd_sum,layer_ids{jj}(layer),scan_ids{jj},g_ind, par, cache);
                end

                % debug by ZC
%                  if any(isnan(object_update_proj(:))) || any(isinf(object_update_proj(:))) 
%                      keyboard
%                      error('Convergence failed')
%                  end
                 % end debug

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%  calculate the optimal step %%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if  ll == 1   % && layer == par.Nlayers
                    if par.beta_LSQ > 0 && object_reconstruct && probe_reconstruct && par.Nlayers == 1  % for Nlayers > 1 seems to be better using gradient_projection_solver in order to keep behaviour the same for each layer !!
                        % calculate the optimal step using LSQ method 
                        [beta_probe(g_ind,layer), beta_object(g_ind,layer)] = ...
                            get_optimal_LSQ_step(self,chi{ll},object_update_proj,m_probe_update,obj_proj{llo},probe{llp,layer},p_ind{ll} , par); 

                    elseif par.beta_LSQ > 0 && (object_reconstruct || probe_reconstruct )
                        % computationally cheaper method that assumes only
                        % diagonal terms of the AA matrix 
                        % debug by ZC
%                          if any(isnan(object_update_proj(:))) || any(isinf(object_update_proj(:))) 
%                              keyboard
%                              error('Convergence failed')
%                          end
                         % end debug
                        [beta_probe(g_ind,layer),beta_object(g_ind,layer)] = ...
                            gradient_projection_solver(self,chi{ll},obj_proj{llo},probe{llp,layer},...
                                object_update_proj, m_probe_update,p_ind{ll}, par, cache); 
                    elseif par.beta_LSQ == 0
                        % use constant step 
                        beta_probe(g_ind,layer)  = par.beta_probe;
                        beta_object(g_ind,layer)  = par.beta_object;
                    end

                    beta_object(g_ind,layer) =  beta_object(g_ind,layer) / par.Nlayers;

                    if (is_used(par, 'fly_scan') || par.apply_multimodal_update)
                         beta_object(g_ind,layer) =  beta_object(g_ind,layer) / par.probe_modes;
                    end
                end
                object_update_proj = [];   % soft memory clean 

                
                beta_object = max(0.2, beta_object); 
                beta_probe  = max(0, beta_probe); 

                if any(jj == ind_range(end)) && ll == 1   
                   %  show the optimal steps calculation 
    %                if verbose() > 2
    %                     if mod(iter, 5) == 0
    %                         plotting.smart_figure(121122)
    %                         plot(squeeze(real(beta_object)))
    %                         hold all 
    %                         plot(squeeze(real(beta_probe)))
    %                         hold off 
    %                         legend({'O', 'P'})
    %                         title('Update step for each scan position')
    %                         xlabel('Scan positions #')
    %                         ylabel('Step')
    %                         drawnow
    %                     end
    %                end
                   verbose(3,'Average step p%i: %3.3g  o%i: %3.3g layer %i', llp,mean(beta_probe(g_ind,layer)),llo,mean(beta_object(g_ind,layer)), layer);
                end

                %%%%%%%%%%%%%%% apply update with the optimal LSQ step %%%%%%%%%%%%%%%%%
                if probe_reconstruct && layer == 1 && ll <= max(par.probe_modes)    % multilayer extension -> update probe only from the first layer
                    self.probe{ll} = update_probe(self.probe{ll}, m_probe_update, par, p_ind{ll}, g_ind, beta_probe, Nind); % finally update also the probe                   
                end

                if object_reconstruct && is_method(par, 'MLs') 
                    self.object = update_object(self, self.object, object_upd_sum, layer_ids{jj}(layer), llo, {g_ind}, scan_ids(jj), par, cache, beta_object);
                    
%                      % suppress large amplitude of object, Added by ZC
%                     if par.amplitude_threshold_object < inf 
%                         for ilay=1:par.Nlayers
%                             temp = abs(self.object{ilay});
%                             if any(temp(:) > par.amplitude_threshold_object)
%                                 warning('amplitude too large (inner loop), force to 1');
%                                 temp (temp > par.amplitude_threshold_object) = 1;
%                                 self.object{ilay} = temp.* exp(1i* angle(self.object{ilay}));
%                             end                            
%                         end
%                     end
                    % debug by ZC
%                     if any(abs(self.object{30}(:))>1.5 | abs(self.object{30}(:))<0.5)
%                         keyboard;
%                     end
                end

%                 if par.probe_position_multiprobe  &&  layer == par.layer4pos                          
%                     [pos_update, pos_rotate_upd,probe_scale_upd, cache] = gradient_position_solver(self, chi{1}, obj_proj{1},probe{ll,layer}, g_ind, iter, cache, par);
%                     self.modes{1}.probe_positions(g_ind,:)=self.modes{1}.probe_positions(g_ind,:)+pos_update;                       
                if  ll == 1  &&  ismember(layer, par.layer4pos) %layer == ceil(par.Nlayers/2)  % assume that sample in center is best constrained . Changed to variable layer by Zhen Chen
                    %%%%%%%%%%%%% update other parameters of the ptychography model
                    if  iter >= par.probe_position_search || iter >= par.detector_rotation_search || iter >= par.detector_scale_search
                        % find optimal position shift that minimize chi{1} in current iteration 
                        if par.precession_mode && par.correct_angle && ismember(iter, par.iter_correct_angle - 1)
                            %modified by dzh to take in tilt angle as a free parameter
                            [pos_update, ~,~, cache] = gradient_position_solver(self, chi{1}, obj_proj{1},probe{1,layer}, g_ind, iter, cache, par);
                            self.modes{layer}.probe_positions_update{1}(g_ind,:) = pos_update;
                            % self.modes{layer}.probe_scale_upd(end+1)=self.modes{layer}.probe_scale_upd(end)+mean(probe_scale_upd); 
                            % self.modes{layer}.probe_positions(g_ind,:)=self.modes{layer}.probe_positions(g_ind,:)+pos_update;
                            % self.modes{layer}.probe_rotation_all(g_ind)=self.modes{layer}.probe_rotation_all(g_ind)+squeeze(pos_rotate_upd);
                        elseif par.precession_mode && (~par.correct_angle || ~ismember(iter, par.iter_correct_angle - 1))
                            % precession mode but not correcting angle
                            if layer==round(mean(par.layer4pos))
                                [pos_update, ~,~, cache] = gradient_position_solver(self, chi{1}, obj_proj{1},probe{1,layer}, g_ind, iter, cache, par);
                                % dzh: save the desired update for each layer; apply the position update after this iteration
                                for ilayer = 1:par.Nlayers
                                    self.modes{ilayer}.probe_positions_update{1}(g_ind,:) = pos_update;
                                end
                            end
                        elseif ~par.precession_mode %normal mode
                            [pos_update, pos_rotate_upd,probe_scale_upd, cache] = gradient_position_solver(self, chi{1}, obj_proj{1},probe{1,layer}, g_ind, iter, cache, par);
                            self.modes{layer}.probe_scale_upd(end+1)=self.modes{layer}.probe_scale_upd(end)+mean(probe_scale_upd); 
                            self.modes{layer}.probe_positions(g_ind,:)=self.modes{layer}.probe_positions(g_ind,:)+pos_update;
                            self.modes{layer}.probe_rotation_all(g_ind)=self.modes{layer}.probe_rotation_all(g_ind)+squeeze(pos_rotate_upd);
                        end
                    end

                    if iter > par.probe_fourier_shift_search 
                        % search position corrections in the Fourier space, use
                        % only informatiom from the first mode, has to be after
                        % the probes updated , SEARCH ONLY FOR THE FIRST MODE 
                        self.modes{1} = gradient_fourier_position_solver(chi{1}, obj_proj{1},probe{1,layer},self.modes{1}, g_ind);
                    end
                end

                if par.Nlayers > 1
                    % get update direction for the next layer 
                     chi{ll} = probe_update; % .* median(beta_probe(g_ind,layer));  
                else
                     chi{ll} = []; % soft memory clean 
                end
            end 
        end

        if iter > par.estimate_NF_distance
            %% correct propagation distance if updated
            for i = 1:par.Nlayers %par.Nmodes % layers by Zhen Chen   
                ASM = exp( self.modes{i}.distances(end)* cache.ASM_difference);
                self.modes{i}.ASM_factor = ASM;
                self.modes{i}.cASM_factor = conj(ASM);
            end
        end
        % to be used for momentum calculation, use only the last layer 
        if par.momentum
            uniq_p_ind = unique(p_ind{ll});
            probe_update_sum(:,:,uniq_p_ind,1) = probe_update_sum(:,:,uniq_p_ind,1) + m_probe_update / Nind; 
        end
        
        %% constraint periodic along propagation,Added by ZC
        % regularize_layers works better (set regularize_layers > 0), then this is unnecessary.
%         if par.Nlayers > 1 && isfield(par.p,'samelayer') && par.p.samelayer && (iter >= par.p.Nst_samelayer && iter <= par.p.Nend_samelayer)
%             object_avg=0;
%             for layer=1:par.Nlayers-1 % last layer is fixed to 1s.
%                 object_avg = object_avg + self.object{layer};
%             end
%             object_avg = object_avg ./(par.Nlayers-1);
% %             object_avg = (object_avg .* exp(1i*angle(self.object{par.Nlayers})))./(par.Nlayers-1);
%             for layer=1:par.Nlayers-1
%                 self.object{layer} = object_avg ;
%             end  
%         end
%         % self.object{par.Nlayers}(:)=1;
%         

   end

   % output warning message for amplitude too large or small.
    if N_amp_out > 0 % beginning iterations not monitor
          warning(strcat('amplitude range=[',num2str(min([min_amp,par.amplitude_threshold_object_lower]),'%5.2f'),', ',num2str(max([max_amp,par.amplitude_threshold_object]),'%5.2f'),...
                    '], ',num2str(N_amp_out,'%6d'),' times,  force outliers to 1 (inner loop)!'));
    end
    if N_phs_out > 0 % beginning iterations not monitor
          warning(strcat('phase range=[',num2str(min([min_phs,par.phase_threshold_object_lower]),'%5.2f'),', ',num2str(max([max_phs,par.phase_threshold_object]),'%5.2f'),...
                    '], ',num2str(N_phs_out,'%6d'),' times,  force outliers to 1 (inner loop)!'));
    end
        
   if iter == 0
       % apply initial correction for the probe intensity and return
       probe_amp_corr = sqrt(probe_amp_corr(1) / probe_amp_corr(2)); %% calculate ratio between modF^2 and aPsi^2
       for ii = 1:par.probe_modes
           self.probe{ii} = self.probe{ii}*probe_amp_corr;
       end
       verbose(2,'Probe amplitude corrected by %.3g',probe_amp_corr)
       return
   end
  
     
   % applying single update emulates behaviour of the original ML method ->
   % provides better noise robustness 
   % advantage is that less memory is needed and no linesearch is required
   
    object_reconstruct = iter >= par.object_change_start; % && (par.apply_multimodal_update || is_used(par, 'fly_scan') || ll <= par.object_modes );
    probe_reconstruct = (iter >= par.probe_change_start && iter < par.probe_change_end);
    % if true, caclulate momentum and use is for acceleration 
    momentum_acceleration = isfield(par, 'momentum')  && par.momentum && par.number_iterations < par.accelerated_gradients_start; 
    if object_reconstruct && is_method(par, 'MLc')   
        for ll = 1:par.object_modes
            for layer = 1:par.Nlayers
                [self.object, object_upd_sum] = update_object(self, self.object, object_upd_sum, layer, ll, indices, scan_ids, par, cache, beta_object(:,layer));
            end
        end
        %% apply momentum update on the object
        if momentum_acceleration
            [self, cache] = add_momentum_object(self, cache, par, object_upd_sum, iter, fourier_error, beta_object);
        end
    end
   
    if probe_reconstruct && momentum_acceleration
        %% apply momentum update on the probe
        [self, cache] = add_momentum_probe(self, cache, par, {probe_update_sum}, iter, fourier_error, beta_probe);
    end
    %if iter>0
    %   disp(any(~isfinite(fourier_error(iter,:))))
    %end
    %% FLY-SCAN: join all subprobes
    if iter >= par.probe_change_start && iter < par.probe_change_end
       if is_used(par,'fly_scan')
           probe_new = 0;
           for ll = 1:par.Nmodes
                probe_new = probe_new + self.probe{ll}/par.Nmodes;
           end
           for ll = 1:par.Nmodes
               % allow variation of the modes intensity 
                proj(ll,1,:) = real(sum2( self.probe{ll} .* conj( probe_new)) ./ sum2(abs(probe_new).^2)) ;
                self.probe{ll} = proj(ll,1,:) .* probe_new;
                % assume constant intensity 
%                  self.probe{ll} =  probe_new;
           end
       end
    end
%     if iter >= par.probe_change_start 
%         for ll = 1:par.probe_modes
%             self.probe{ll} = apply_probe_contraints(self.probe{ll}, self.modes{min(ll,end)});
%         end
%     end
    
   
end


%% merged CUDA kernels for faster calculations 
function update = get_update(chi, proj)
    update = chi .* conj(proj); 
end

function [update_1, update_2] = get_update_both(chi, proj_1, proj_2)
    update_1 = chi .* conj(proj_1); 
    update_2 = chi .* conj(proj_2); 
end




