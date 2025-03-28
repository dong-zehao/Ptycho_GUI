% PTYCHO_SOLVER  the main loop of ptychography. Calls the selected engine 
% apply additional constraints, and tries to remove ambiguities
% 
% [outputs, fourier_error, fsc_score] = ptycho_solver(self, par, cache)
%
% ** self      structure containing inputs: e.g. current reconstruction results, data, mask, positions, pixel size, ..
% ** par       structure containing parameters for the engines 
% ** cache     structure with precalculated values to avoid unnecessary overhead
%
% returns:
% ++ outputs        self-like structure with final reconstruction
% ++ fourier_error  array [Npos,1] containing evolution of reconstruction error 
% ++ fsc_score      [] or a structure with outputs from online estimation of FSC curve
%

function [outputs, fourier_error, fsc_score] = ptycho_solver(self, par, cache)

import engines.GPU_MS.analysis.*
import engines.GPU_MS.shared.*
import engines.GPU_MS.initialize.*
import engines.GPU_MS.GPU_wrapper.*
import math.*
import utils.*



% precalculate the parallel block sizes and sets 
[cache, par] = get_parallel_blocks(self, par, cache); 
 
global gpu use_gpu
verbose( par.verbose_level )    

if  (nargout == 1 && verbose()  == 0 && isinf(par.plot_results_every))
    par.get_error = false;
elseif ~isfield(par, 'get_error')
    par.get_error = true;
end

if use_gpu 
    %modified by YJ to print out more info
    if isfield(par,'extraPrintInfo')
        verbose(struct('prefix',['GPU-',num2str(gpu.Index),'_', par.method,'_',par.extraPrintInfo]))
    else
        verbose(struct('prefix',['GPU-',num2str(gpu.Index),'_', par.method]))
    end
    %verbose(struct('prefix',['GPU-', par.method]))
    verbose(0,'Started solver using %s method on GPU %i', par.method, gpu.Index )
    
else
    verbose(struct('prefix',['CPU-', par.method]))
    verbose(0,'Started solver using %s method on CPU', par.method)
end


lastwarn('') 
fsc_score = cell(1,0);
par.Nscans = length(self.reconstruct_ind);

%% move everything on GPU if needed; 
if par.use_gpu
    % if not sparse solvers as ePIE, hPIE, MLs are used presplit data into
    % bunches (allow larger data to be processed )
    %%%%split_data = is_method(par, {'MLc', 'DM'});
    split_data = false; % modified by YJ, seems to avoid some errors from GPU
    [self, cache] =  move_to_gpu(self,cache, par.keep_on_gpu, split_data);
end 
if par.share_object && par.object_modes == 1
    % enforce only a single object 
    self.object = self.object(1,:);
    cache.illum_sum_0 = cache.illum_sum_0(1);
end

%% allocate memory 
fourier_error = Garray( nan(par.number_iterations, self.Npos));
if is_method(par, {'DM'})
    psi_dash = cell(max(par.probe_modes,par.object_modes), length(cache.preloaded_indices_simple{1}.indices));
end
if is_method(par, {'PIE', 'ML'})
    if par.beta_LSQ
        cache.beta_object = ones(self.Npos,par.Nlayers,'single')*par.beta_object;
        cache.beta_probe  = ones(self.Npos,par.Nlayers,'single')*par.beta_probe;
    else
        cache.beta_object = single(par.beta_object);
        cache.beta_probe  = single(par.beta_probe);
    end
    switch lower(par.likelihood)
        case 'l1', cache.beta_xi = 1; % optimal step for gauss
        case 'poisson', cache.beta_xi = 0.5*ones(1,1,self.Npos,'single'); % for poisson it will be further refined 
    end
end

%% in case of multilayer extension assume that the provided probe is positioned in middle of the sample -> shift it at the beginning
if par.Nlayers > 1 && par.preshift_ML_probe
   probe_offset = -sum(self.z_distance(1:end-1))/2; 
   for ii = 1:par.probe_modes
        self.probe{ii} = utils.prop_free_nf(self.probe{ii}, self.lambda , probe_offset ,self.pixel_size) ;
   end
end

%% in case of tilted plane ptychography, tilt the provided probe 
if any(par.p.sample_rotation_angles(1:2)) && check_option(par.p, 'apply_tilted_plane_correction', 'propagation') 
    % apply propagators to the tilted plane 
    for ii = 1:par.probe_modes % bugs! because modes are layers
        self.probe{ii} = self.modes{ii}.tilted_plane_propagate_fwd(self.probe{ii});
    end
end

global pprev;
pprev = -1;

mode_id = round(mean(par.layer4pos));  % main mode (assume single most important mode for approchimations) 

t0 = tic;
t_start = tic;
iter_start = 1;

% object averaging for DM code 
for ll = 1:length(self.object)
    object_avg{ll} = 0; 
end
N_object_avg = 0;
%par.initial_probe_rescaling = true or false
%% added by YJ
update_position_weight = true; %calculate position weight the first time find_geom_correction() is called
fourier_error_mean_previous = 1e10; %for convergence check
if par.fourier_error_threshold<inf
    verbose(0, 'Convergence check is enabled. Threshold = %3.3g.', par.fourier_error_threshold)
end
%% added by YJ: update figures after initialization
if par.p.use_display
    if verbose() <= 0
    	ptycho_plot_wrapper(self, par, 0)
    end
end
%% added by YJ: save initial probe (after init_solver.m's pre-processing)
probe_temp = Ggather(self.probe{1});
probe_init = zeros(size(probe_temp,1),size(probe_temp,2),par.probe_modes,par.variable_probe_modes+1);
for ll = 1:par.probe_modes
   probe_temp = Ggather(self.probe{ll});
   probe_init(:,:,ll,:) = probe_temp(:,:,1,:);
end
%%
for iter =  (1-par.initial_probe_rescaling):par.number_iterations
    if iter == par.probe_position_search || iter==par.detector_scale_search+1 || iter == par.detector_rotation_search+1
        t0 = tic; %reset time for estimating avgTimePerIter. added by YJ
        iter_start = iter;
    end
    if iter > 0
        %{
        if verbose() == 0
            progressbar(iter, par.number_iterations, max(20,round(sqrt(par.number_iterations))))
        else
            verbose(1,'Iteration %s: %i / %i  (time %3.3g  avg:%3.3g)', par.method, iter, par.number_iterations, toc(t_start), toc(t0)/(iter-1))
        end
        %}
        %modified by YJ to print out more details

        %verbose(0,'Iteration %s: %i / %i  (time %3.3g  avg:%3.3g)', par.method, iter, par.number_iterations, toc(t_start), toc(t0)/(iter-1))
        avgTimePerIter = toc(t0)/(iter-iter_start);
        timeLeft = (par.number_iterations-iter+1)*avgTimePerIter;
        %verbose(0, 'Method: %s, GPU id: %i',par.method, gpu.Index)
        if timeLeft>3600
            verbose(0, 'Iteration: %i / %i  (Time left:%3.3g hour. avg:%3.3g sec)', iter, par.number_iterations, timeLeft/3600, avgTimePerIter)
            utils.LogMessage(par.p.apphandle, sprintf('GPU-%i%s: Iteration: %i / %i  (Time left:%3.3g hour. avg:%3.3g sec)',par.gpu_id, par.extraPrintInfo, iter, par.number_iterations, timeLeft/3600, avgTimePerIter))
        elseif timeLeft>60
            verbose(0,'Iteration: %i / %i  (Time left:%3.3g min. avg:%3.3g sec)', iter, par.number_iterations, timeLeft/60, avgTimePerIter)
            utils.LogMessage(par.p.apphandle, sprintf('GPU-%i%s: Iteration: %i / %i  (Time left:%3.3g min. avg:%3.3g sec)',par.gpu_id, par.extraPrintInfo, iter, par.number_iterations, timeLeft/60, avgTimePerIter))
        else
            verbose(0,'Iteration: %i / %i  (Time left:%3.3g sec. avg:%3.3g sec)', iter, par.number_iterations, timeLeft, avgTimePerIter)
            utils.LogMessage(par.p.apphandle, sprintf('GPU-%i%s: Iteration: %i / %i  (Time left:%3.3g sec. avg:%3.3g sec)',par.gpu_id, par.extraPrintInfo, iter, par.number_iterations, timeLeft, avgTimePerIter))
        end

    end

    t_start = tic;
    if  iter > 0.9*par.number_iterations && is_method(par, 'DM')
        for ll = 1:length(self.object)
            object_avg{ll} = object_avg{ll}  + self.object{ll}; 
        end
        N_object_avg = N_object_avg +1 ; 
        verbose(1,'==== Averaging DM result ======')
    end
     
    %% GEOMETRICAL CORRECTIONS        
    if (iter > par.probe_position_search || iter > par.detector_rotation_search) && is_method(par, {'PIE', 'ML'})
        self = find_geom_correction(self,cache,par,iter,mode_id, update_position_weight);
        % Added by YJ
        if mod(iter-par.probe_position_search+1, par.update_pos_weight_every) == 0
            update_position_weight = true;
        else
            update_position_weight = false;
        end
    end  

    %% in TCMEP mode, shift back the probe according to the global shift: dzh test
    % if isfield(self, 'upd_global_shift')
    %     for imode = 1:par.Nmodes
    %         self.probe{imode} = apply_subpx_shift(self.probe{imode}, self.upd_global_shift);
    %     end
    % end

    %% remove extra degree of freedom for OPRP and other optimizations
    if iter > par.probe_fourier_shift_search
        for kk = 1:par.Nscans
            ind = self.reconstruct_ind{kk};
            self.modes{1}.probe_fourier_shift(ind,:) = self.modes{1}.probe_fourier_shift(ind,:) - mean(self.modes{1}.probe_fourier_shift(ind,:));
        end    
    end
        
    %% remove ambiguity related to the variable probe 
    if par.variable_probe && iter > par.probe_change_start && is_method(par, 'ML')
         self = remove_variable_probe_ambiguities(self,par); 
    end

   %% remove the ambiguity in the probe / object reconstruction => keep average object transmission around 1
   if mod(iter,10)==1 &&  par.remove_object_ambiguity  && ~is_used(par, {'fly_scan'}) &&  ~is_method(par, {'DM', 'PIE'})  % too slow for variable probe 
        self = remove_object_ambiguity(self, cache, par) ; 
   end
   
    if (mod(iter, 10) == 1 || iter  < 5) && check_option(par, 'get_fsc_score')   && ...
       (((par.Nscans > 1 ) && size(self.object,1) == par.Nscans) || ... 
       ( check_option(self, 'object_orig') ))
            
        %% Fourier ring correlation between two scans with independend objects 
        aux = online_FSC_estimate(self, par, cache, fsc_score(end,:), iter); 
        fsc_score(end+1,1:size(aux,1), 1:size(aux,2)) = aux; 
    end
    
    %% ADVANCED FLY SCAN. disabled by YJ. it's already calculated init_solver.m . probably useful with position correction?
    if  is_used(par, 'fly_scan')
        if iter == 1
           disp(['== AVG fly scan step ', num2str( median(sqrt(sum(diff(self.probe_positions_0).^2,1)))  )]) 
        end
        self = prepare_flyscan_positions(self, par); 
    end
    
    %% update current probe positions (views)
    if iter <= 1 || iter >= par.probe_position_search
        %%%%%%%%  crop only ROI of the full image for subsequent calculations  %%%%%%%%%% 
        for ll = 1:par.Nlayers
            if ~is_used(par, 'fly_scan') && par.precession_mode
                %% modified by dzh to add oROI for every layer; different scans are included in each layer
                [cache.oROI_s{ll},cache.oROI{ll},sub_px_shift] = find_reconstruction_ROI( self.modes{ll}.probe_positions,self.Np_o, self.Np_p); 
                self.modes{ll}.sub_px_shift = sub_px_shift; 
            elseif ~is_used(par, 'fly_scan') && ~par.precession_mode
                %% conventional farfield/nearfield ptycho -> update view coordinates and keep subpixels shift < 1 px 
                [cache.oROI_s{ll},cache.oROI{ll},sub_px_shift] = find_reconstruction_ROI( self.modes{1}.probe_positions,self.Np_o, self.Np_p); 
                self.modes{ll}.sub_px_shift = sub_px_shift; 
            else 
                % if is_used(par, 'fly_scan')   
                %% flyscan farfield ptycho -> update view coordinates and keep subpixels shift < 1 px 
                [cache.oROI_s{ll},cache.oROI{ll},sub_px_shift] = find_reconstruction_ROI( self.modes{1}.probe_positions,self.Np_o, self.Np_p); 
                %% use the fftshift  for much larger corrections in the case of the fly scan 
                self.modes{ll}.sub_px_shift = self.modes{ll}.probe_positions -  self.modes{1}.probe_positions_0 + sub_px_shift;                
            end
        end
    end
    
    %% update probe fft support window  if mode.probe_scale_upd(end) ~= 0, important to avoid issues during subpixel probe rescaling (is pixel scale search)
    if self.modes{1}.probe_scale_upd(end) > 0
        self.modes{1}.probe_scale_window = get_window(self.Np_p, 1+self.modes{1}.probe_scale_upd(end), 1) .* get_window(self.Np_p, 1+self.modes{1}.probe_scale_upd(end), 2); 
    elseif self.modes{1}.probe_scale_upd(end) < 0
        self.modes{1}.probe_scale_window = fftshift(get_window(self.Np_p, 1-self.modes{1}.probe_scale_upd(end), 1) .* get_window(self.Np_p, 1-self.modes{1}.probe_scale_upd(end), 2)) ; 
    else
        self.modes{1}.probe_scale_window = [];
    end
   
    %% updated illumination; modified by dzh to take TCMEP into account
    if iter <= 1 || ( iter > par.probe_change_start && iter < par.probe_change_end && (mod(iter, 10) == 1 || iter < par.probe_change_start+10 ))
        if par.precession_mode
            layerid = ceil(par.precession_vacuum_thickness/self.z_distance(1));
        else
            layerid = 1;
        end
        for ll = 1:size(self.object,1)
            for iscan = 1:par.Nscans
                aprobe2 = abs(self.probe{1}(:,:,iscan)).^2; 
                illum_sum_temp = Gzeros(self.Np_o);
                ind = [self.reconstruct_ind{iscan}];
                % divide the ind into parts of 60000 to keep within uint16: 
                % debug by dzh
                ind_parts = ceil(ind/60000); % which part does all inds belong to
                for ipart = 1:max(ind_parts)
                    ind_temp = ind(ind_parts == ipart);
                    illum_sum_temp = set_views(illum_sum_temp, Garray(aprobe2), layerid,1, ind_temp, cache);
                end
                cache.illum_sum_0{ll} = cache.illum_sum_0{ll} + illum_sum_temp/par.Nscans;
                % avoid oscilations by adding momentum term 
                cache.illum_norm(ll) = norm2(cache.illum_sum_0{ll});
                cache.MAX_ILLUM(ll) = max2(cache.illum_sum_0{ll});
            end
        end
    end
    
    %% improve convergence speed by gradient acceleration 
    if is_method(par, 'MLc') && iter >= par.accelerated_gradients_start
        [self, cache] = accelerate_gradients(self, par, cache, iter); 
    end
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%  PERFORM ONE ITERATION OF THE SELECTED METHOD %%%%%%%%%%%%%%%%%%%%%
    switch  lower(par.method)
        case {'epie', 'hpie'}
            [self, cache, fourier_error] = engines.GPU_MS.PIE(self,par,cache,fourier_error,iter);
        case { 'mls','mlc'}
            [self, cache, fourier_error] = engines.GPU_MS.LSQML(self,par,cache,fourier_error,iter);
        case 'dm'
            [self, cache,psi_dash,fourier_error] =  engines.GPU_MS.DM(self,par,cache,psi_dash,fourier_error,iter);
        otherwise
            error('Not implemented method')
    end
    
    % check stop
    if par.p.apphandle.StopRequested
        utils.LogMessage(app, 'Stopping received ...');
        drawnow;
        return;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if iter == 0; continue; end  % interation 0 is used only to calibrate iinitial probe intensity
    
    if any(~isnan(fourier_error(iter,:)))       
        switch lower(par.likelihood)
            case  'l1'
                verbose(0,'=====  Fourier error = %3.4g ', nanmedian(fourier_error(iter,:)) ); 
                utils.LogMessage(par.p.apphandle, sprintf('=====  Fourier error = %3.4g ', nanmedian(fourier_error(iter,:))))
            case 'poisson'
                err = fourier_error(iter,:) - fourier_error(1,:);
                verbose(0,'=====  Log likelihood = %3.5g ', nanmedian(err)); 
        end
    end
    %% %%%%%%%%%%%%%%%%%%%%%%%%%% convergence check added by YJ %%%%%%%%%%%%%%%%%%%%%%%%%%
    if any(~isnan(fourier_error(iter,:)))
        switch lower(par.likelihood)
            case  'l1'
                fourier_error_mean = mean(fourier_error(iter,:)); %maybe median is better?
                fourier_error_diff = (fourier_error_mean - fourier_error_mean_previous)/fourier_error_mean_previous;

                if fourier_error_diff > par.fourier_error_threshold
                    %verbose(0, 'Convergence check is enabled.Threshold = %3.3g.', par.fourier_error_threshold)
                	error('Current Fourier error exceeds the previous one. Terminate reconstruction.')
                	%error('Current Fourier error exceeds the previous one.')
                end
                fourier_error_mean_previous = fourier_error_mean;
            case 'poisson'
            %    err = fourier_error(iter,:) - fourier_error(1,:);
            %    verbose(1,'=====  Log likelihood = %3.5g ', nanmedian(err)); 
        end 
    end
    
       %% suppress large or small amplitude of object, Added by ZC
    if iter >= par.object_change_start && (par.amplitude_threshold_object < inf || par.amplitude_threshold_object_lower > 0 )
        for ll=1:par.Nlayers
            temp = abs(self.object{ll});
            if any(temp(:) > par.amplitude_threshold_object) || any(temp(:) < par.amplitude_threshold_object_lower)
                warning(strcat('amplitude range=[',num2str(min(temp(:)),'%5.2f'),', ',num2str(max(temp(:)),'%5.2f'),...
                    '], force outliers to 1.0!'));
                temp (temp > par.amplitude_threshold_object) = 1;
                temp (temp < par.amplitude_threshold_object_lower) = 1;
                self.object{ll} = temp.* exp(1i* angle(self.object{ll}));
            end 
        end
    end
    %% %%%%%%%%%%%%%%%%%%%%%%%%%% CORRECTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % apply probe constraints 
    % if iter >= par.probe_change_start
    %     if (check_option(par,'probe_support_tem') && check_option(par,'probe_support_tem_Nend') && iter < par.probe_support_tem_Nend)  % quick implement for TEM aperture constrain by Zhen Chen, bugs for probe rescale.
    % 
    %         % propagate probe on the detector 
    %         probe_temp = fft2_safe(self.probe{mode_id});  
    % 
    %         if ~isempty(self.modes{mode_id}.probe_support_fft)
    %             probe_temp = probe_temp .* self.modes{mode_id}.probe_support_fft;
    %         end
    % 
    %         % propagate probe back to the sample plane 
    %         self.probe{mode_id} = ifft2_safe(probe_temp);  
    %     elseif ~ check_option(par,'probe_support_tem') % force not use near field propagation if probe_support_tem exists, by Zhen Chen
    %         self.probe{mode_id} = apply_probe_contraints(self.probe{mode_id}, self.modes{mode_id});
    %     end
    % end
    
  if iter >= par.probe_change_start && iter < par.probe_change_end % do correction for all modes by ZC
         for mode_ii=1:length(self.probe) % modes fixed to 1, correct bug for more modes than Nlayers
                if (check_option(par,'probe_support_tem') && check_option(par,'probe_support_tem_Nend') && iter < par.probe_support_tem_Nend)  % quick implement for TEM aperture constrain by Zhen Chen, bugs for probe rescale.
                    
                    % propagate probe on the detector 

                    probe_temp = fft2_safe(self.probe{mode_ii});  

                    if ~isempty(self.modes{1}.probe_support_fft)
                        probe_temp = probe_temp .* self.modes{1}.probe_support_fft;
                    end

                    % propagate probe back to the sample plane 
                    self.probe{mode_ii} = ifft2_safe(probe_temp);  
                elseif ~ check_option(par,'probe_support_tem') % force not use near field propagation if probe_support_tem exists, by Zhen Chen
                    self.probe{mode_ii} = apply_probe_contraints(self.probe{mode_ii}, self.modes{1});
                end
         end
  end
    
    % push low illum regions of object to zero  
    if par.delta > 0 
        if  iter > par.object_change_start && iter < par.probe_change_end
            %for ll = 1:max(par.object_modes, par.Nscans) %old code
            for ll = 1:par.Nlayers % For multilayer. Added by ZC
                % push everywhere to zero, even out of the object region
                self.object{ll} = Gfun(@regular_object_out_regions, self.object{ll}, cache.illum_sum_0{1}, cache.MAX_ILLUM(1),par.delta); 
            end
        end
        if  iter > par.probe_change_start && iter < par.probe_change_end && ~par.variable_probe
            for ll = 1:par.probe_modes
                % push everywhere to zero, even out of the object region
                self.probe{ll} = self.probe{ll}  .* (1-par.delta);
            end
        end
    end

    
    %% suppress uncontrained values !! 
    if iter > par.object_change_start && par.object_regular(1) > 0
        %for ll = 1:max(par.object_modes, par.Nscans) %old          
        for ll = 1:par.Nlayers % Added by ZC. quick fix for multilayer, but does not work for multiple scans or modes,  
            self.object{ll} = apply_smoothness_constraint(self.object{ll},par.object_regular(1)); % blur only intensity, not phase 
        end
    end
    
    %% constrain amplitude to phase *0.1 by Zhen Chen
        reg_amp=par.reg_amp;
        if reg_amp ~= 0 && iter < par.reg_amp_iter
            Wphase = min(1, 10*(cache.illum_sum_0{1}/cache.MAX_ILLUM(1))); 
            for ii = 1:par.Nlayers   
                temp_phs= math.unwrap2D_fft2(self.object{ii},[],0,Wphase,-1);
                temp_phs=temp_phs-min(min(temp_phs.*Wphase,[],1),[],2);
                temp_amp1=abs(self.object{ii});
                temp_amp=(1-temp_phs*0.1)*reg_amp+temp_amp1*(1-reg_amp);
                offset=mean(temp_amp1(:))-mean(temp_amp(:)); % make the mean not change
                self.object{ii} = exp(1i*angle(self.object{ii})).*(temp_amp+offset); 
            end
        end
     %% Added by YJ. TV regularization on object
    if iter > par.object_change_start  && isfield(par,'TV_lambda') && par.TV_lambda > 0
        N_tv_iter = 10;
        for ll = 1:par.Nlayers
            self.object{ll} = local_TV2D_chambolle(self.object{ll}, par.TV_lambda, N_tv_iter);
         %self.object{ll}(cache.object_ROI{:}) = ...
         %    Gfun(@local_TV2D_chambolle,self.object{ll}(cache.object_ROI{:}), par.TV_lambda, N_tv_iter);
        end
    end
        
    %% weak positivity object 
    if iter > par.object_change_start  && any(par.positivity_constraint_object)
        for ll = 1:par.Nlayers 
        %for ll = 1:par.object_modes
             self.object{ll}(cache.object_ROI{:}) = ...
                 Gfun(@positivity_constraint_object,self.object{ll}(cache.object_ROI{:}), par.positivity_constraint_object);
         end
    end
    
    %% probe orthogonalization 
    if iter>par.probe_change_start && iter < par.probe_change_end && par.probe_modes > par.Nrec && (~is_method(par, 'DM') || iter == par.number_iterations)
        %  orthogonalization of incoherent probe modes 
        if is_used(par, 'fly_scan')
            probes = self.probe;
            % orthogonalize the modes with all the other shifted modes 
            for i = 1:par.Nrec
                dx = mean(self.modes{i}.sub_px_shift - self.modes{1}.sub_px_shift); 
                % apply average shift 
                probes{i} = imshift_fft(probes{i}, dx); 
            end
            probes = ortho_modes(probes);  % perform othogonalization 
            % update only the incoherent orthogonal modes 
            self.probe(1+par.Nrec:par.probe_modes) = probes(1+par.Nrec:par.probe_modes);
        else
            %% orthogonalize the incoherent probe modes ; modified by dzh to take in multiple scans
            for iscan = 1:par.Nscans
                for ii = 1:par.probe_modes
                    P(:,:,ii) = self.probe{ii}(:,:,iscan); 
                end
                P = core.probe_modes_ortho(P);
                for ii = 1:par.probe_modes
                    self.probe{ii}(:,:,iscan) = P(:,:,ii); 
                end
            end
        end
    end

  %% regularize multilayer reconstruction 
  if par.regularize_layers > 0 && par.Nlayers > 1 %  && mod(iter, 2) == 1, ## change from Nlayers > 1 to Nlayers > 2 by Zhen Chen
        self = regulation_multilayers(self, par, cache);
  end
  
  % constraint periodic along propagation, by Zhen Chen.
        % regularize_layers works better (set regularize_layers > 0), then this is unnecessary.
%     if par.Nlayers > 1 && isfield(par.p,'forcelayer') && par.p.forcelayer && (iter >= par.p.Nst_forcelayer && iter <= par.p.Nend_forcelayer)
%         object_avg=0;
%         for layer=1:par.Nlayers 
%             object_avg = object_avg + self.object{layer};
%         end
%         object_avg = object_avg ./(par.Nlayers);
% %             object_avg = (object_avg .* exp(1i*angle(self.object{par.Nlayers})))./(par.Nlayers-1);
%         for layer=1:par.Nlayers
%             self.object{layer} = object_avg ;
%         end  
%     end
  
    if mod(iter, par.plot_results_every)==0 || iter == 1
        par.p.apphandle.object_temp = self.object;
        par.p.apphandle.cache_temp.init_rot = cache.init_rot;
        par.p.apphandle.cache_temp.object_ROI = cache.object_ROI;
        utils.plot_object(par.p.apphandle, self.object, cache)
        utils.plot_probe(par.p.apphandle, self.probe)
        utils.plot_error(par.p.apphandle, fourier_error)
    end
  %% PLOTTING 
    %%%% plot  results %%%%%%%%%%%%
    if par.p.use_display && mod(iter, par.plot_results_every)==0 &&  par.plot_results_every~=0     
        try
            if verbose() <= 0
                % use cSAXS plorring rutines 
                ptycho_plot_wrapper(self, par, fourier_error)
            else
                %   use more detailed plotting rutines 
                if (par.probe_modes > 1)
                    %% probe incoherent modes 
                    plot_probe_modes(self,par);                        
                end
                if (par.Nlayers > 1) || (par.Nscans > 1 && ~par.share_object)
                    %% object incoherent modes 
                    plot_object_modes(self, cache)
                end

                plot_results(self,cache, par, Ggather(fourier_error), ...
                self.modes{mode_id}.probe_positions)
            end

            % show variable modes 
            if (par.variable_probe && par.variable_probe_modes > 0)
                plot_variable_probe(self, par)
            end

            % show position correction in the fourier plane 
            if iter > par.probe_fourier_shift_search
                plotting.smart_figure(24654)
                clf
                hold all
                for ll = 1:par.Nmodes    
                    plot(self.modes{ll}.probe_fourier_shift)
                end
                hold off
                grid on 
                axis tight
                xlabel('Position #')
                ylabel('Corrected probe shift in Fourier plane [px]')
                title('Fourier space probe shift')
            end                

            % show position correction 
            if iter > min([par.probe_position_search, par.estimate_NF_distance, par.detector_rotation_search, par.detector_scale_search]) ...
                    && is_method(par, {'PIE', 'ML'}) 
                     plot_geom_corrections(self, self.modes{mode_id}, Ggather(self.object{1}),iter, par, cache)
                if iter > min(par.detector_rotation_search, par.detector_scale_search)
                    for i = 1:max(1,par.Nrec)
                        verbose(1,sprintf(' Reconstruction id: %i =============  Detector pixel scale: %0.5g  Detector rotation: %0.5g deg', ...
                                            i,  1-self.modes{i}.probe_scale_upd(end) , self.modes{i}.probe_rotation(end,1)))
                    end
                end
            end

            if par.get_fsc_score && ~isempty(fsc_score)
                % plot score estimated by the fourier ring correlation 
                plot_frc_analysis(fsc_score, par)
            end

            drawnow
   
        catch err
            warning(err.message)
            if verbose()  > 1
            	keyboard
            end
        end
    end
    %% save intermediate images, added by YJ
    if isfield(par,'save_results_every') && (mod(iter, par.save_results_every ) == 0 &&  par.save_results_every ~=0) || iter == par.number_iterations
         if ~exist(par.fout, 'dir')
            mkdir(par.fout);
         end
        % save temporary outputs
        if iter == par.number_iterations
            outputs = self;
            outputs.diffraction = [];
        else
            outputs = struct();
        end
        %% ZC's format
        %{
        for ll = 1:par.probe_modes
            outputs.probe{ll} = Ggather(self.probe{ll});
        end
        for ll = 1:max(par.object_modes,par.Nlayers)
            outputs.object{ll} = Ggather(self.object{ll});
            object_roi = Ggather(self.object{ll});
            outputs.object_roi{ll} = object_roi(cache.object_ROI{:});
        end
        %}
        %% YJ's format
        probe_temp = Ggather(self.probe{1});
        %dzh
        probe = zeros(size(probe_temp,1),size(probe_temp,2),par.probe_modes, par.Nscans,'single');
        for ll = 1:par.probe_modes
           for kk = 1:par.Nscans
               probe_temp = Ggather(self.probe{ll});
               probe(:,:,ll,kk) = probe_temp(:,:,kk);
           end
        end
        
        object_temp = Ggather(self.object{1});
        object = zeros(size(object_temp,1),size(object_temp,2),par.Nlayers,'single');
        
        for ll = 1:par.Nlayers %for multislice recon
        	object_temp = Ggather(self.object{ll});
        	object(:,:,ll) = object_temp(:,:,1,1);
        end
        
        % data error
        if strcmp(par.likelihood, 'poisson')
            fourier_error_out = fourier_error- fourier_error(1,:);
            fourier_error_out = Ggather(fourier_error_out);
        else
            fourier_error_out = Ggather(fourier_error);
        end
        outputs.fourier_error_out = mean(fourier_error_out,2,'omitnan'); % omit nan by ZC
        %save the lastest error of all dp
        try
            iter_error = find(~isnan(mean(fourier_error_out,2)));
            iter_error = iter_error(end);
            outputs.fourier_error_dp = fourier_error_out(iter_error,:);
            clear fourier_error_out probe_temp object_temp iter_error %free up memory
        catch
        end
        
        if par.variable_probe && isfield(self,'probe_evolution')
            %self.probe_evolution=(S*V')'. dimension:[Nposition,eng.variable_probe_modes]
            %to obtain probes at each scan position: compute probe*self.probe_evolution'
            outputs.probe_evolution = self.probe_evolution; %
        end
        
        if isfield(self,'sub_objects') %save sub objects corresponding to each dp
            %To be implemented later
            %{
            outputs.sub_objects = self.sub_objects; %
            
            %crop to smaller size 
            if isfield(par,'save_sub_objects_N') && par.save_sub_objects_N < size(outputs.sub_objects,1)
                N_temp = par.save_sub_objects_N;
                cen_temp = floor(size(outputs.sub_objects,1)/2)+1;
                outputs.sub_objects = self.sub_objects(cen_temp-N_temp/2:cen_temp+N_temp/2-1,cen_temp-N_temp/2:cen_temp+N_temp/2-1,:); %
                %disp(size(outputs.sub_objects))
            end
            %}
        end
        
        if par.get_fsc_score && ~isempty(fsc_score)
            outputs.fsc_score = fsc_score;
        end
        
        % scan positions
        outputs.probe_positions = self.modes{mode_id}.probe_positions;
        try %probe_positions_model may not exist if geom refinement is not enabled
            outputs.probe_positions_model = self.modes{mode_id}.probe_positions_model;
        catch
        end
        outputs.affine_matrix = self.affine_matrix{1};
        if isfield(self,'affine_matrix_fit') %save input affine matrix
            outputs.affine_matrix_fit = self.affine_matrix_fit{1};
        end
        if isfield(par,'affine_matrix_init') %save input affine matrix
            outputs.affine_matrix_init = par.affine_matrix_init;
        end
        outputs.probe_positions_0 = self.probe_positions_0;
        outputs.probe_positions_weight = double(self.modes{mode_id}.probe_positions_weight);
        if isfield(self,'tilt_angle')
            outputs.tilt_angle = self.tilt_angle;
        end
        %global errors refinement results
        outputs.detector_rotation = self.modes{mode_id}.probe_rotation(end,:);
        outputs.detector_scale = 1+self.modes{mode_id}.probe_scale_upd(end);
        outputs.relative_pixel_scale = self.modes{mode_id}.scales;
        outputs.rotation =  self.modes{mode_id}.rotation;
        outputs.shear =   self.modes{mode_id}.shear;
        outputs.asymmetry =   self.modes{mode_id}.asymmetry; % added by ZC
        if isfield(self.modes{mode_id}, 'shiftx')
            outputs.shiftx =   self.modes{mode_id}.shiftx; % added by ZC
            outputs.shifty =   self.modes{mode_id}.shifty; % added by ZC
        end
        outputs.pixel_size=self.pixel_size; % added by ZC
        outputs.z_distance = self.modes{mode_id}.distances;

        % 
        outputs.avgTimePerIter = avgTimePerIter;
        outputs.fourier_error_threshold = par.fourier_error_threshold;
        
        % additional parameter for PSI's IO code
        p = {};
        p.object_ROI{1} = cache.object_ROI{1};
        p.object_ROI{2} = cache.object_ROI{2};
        p.binning = false;
        p.detector.binning = false;
        p.dx_spec = self.pixel_size;
        p.lambda = self.lambda;
        
        % store regularization parameters
        if par.amplitude_threshold_object<inf
            p.reg.amplitude_threshold_object = par.amplitude_threshold_object;
        end

        % store parameters related to multi-slice
        p.multi_slice_param.z_distance = self.z_distance;
        p.multi_slice_param.regularize_layers = par.regularize_layers;
        p.multi_slice_param.preshift_ML_probe = par.preshift_ML_probe;
        p.multi_slice_param.layer4pos = par.layer4pos;
        
        % store parameters for object initialization
        p.obj_init_param.init_layer_select = par.init_layer_select;
        p.obj_init_param.init_layer_preprocess = par.init_layer_preprocess;
        p.obj_init_param.init_layer_interp = par.init_layer_interp;
        p.obj_init_param.init_layer_append_mode = par.init_layer_append_mode;
        p.obj_init_param.init_layer_scaling_factor = par.init_layer_scaling_factor;
        if isfield(par.p,'initial_iterate_object_file') && ~isempty(par.p.initial_iterate_object_file{1}) % if initial object is from a given file
        	p.init_object_file = par.p.initial_iterate_object_file;
            if isfield(par.p,'multiple_layers_obj')
                p.obj_init_param.multiple_layers_obj = par.p.multiple_layers_obj;
            end
        end
        % store initial position file
        if isfield(par.p.scan,'custom_positions_source') && ~isempty(par.p.scan.custom_positions_source)
            p.init_position_file = par.p.scan.custom_positions_source;
        end
        % store initial probe
        if isfield(par.p,'initial_probe_file') && ~isempty(par.p.initial_probe_file) % if initial probe is given by a file
        	p.init_probe_file = par.p.initial_probe_file;
        end
        if isfield(par.p,'normalize_init_probe')
            p.normalize_init_probe = par.p.normalize_init_probe;
        else
            p.normalize_init_probe = true;
        end
        if par.save_init_probe
            p.init_probe = probe_init; %store initial probe (after init_solver.m's pre-processing)
        end
        [scale,~,rot,~] = decompose_affine_matrix(par.affine_matrix_init);
        if scale < 1
            rot = rot + 180;
        end
        outputs.init_rot = rot;
        %object_rot = imrotate(object, rot,"bilinear", 'crop' );
        save(strcat(par.fout,'Niter',num2str(iter),'.mat'),'outputs','p','probe','object');
        %% save object phase
        if isfield(par,'save_phase_image') && par.save_phase_image
            %%
            %{
            for ll=1:par.Nlayers
                O_phase_roi = phase_unwrap(angle(outputs.object_roi{ll}));
                %O_phase_roi = angle(outputs.object_roi);
                saveName = strcat('O_phase_roi_','Niter',num2str(iter),'_Layer',num2str(ll),'.tiff');
                saveDir = strcat(par.fout,'/O_phase_roi/');
                if ~exist(saveDir, 'dir')
                    mkdir(saveDir)
                end
                %imwrite(mat2gray(O_phase_roi),strcat(saveDir,fileName),'tiff');

                % save as single tiff
                tagstruct.ImageLength     = size(O_phase_roi,1);
                tagstruct.ImageWidth      = size(O_phase_roi,2);
                tagstruct.Photometric     = Tiff.Photometric.MinIsBlack;
                tagstruct.BitsPerSample   = 16;
                tagstruct.SamplesPerPixel = 1;
                tagstruct.RowsPerStrip    = 16;
                tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
                tagstruct.Software        = 'MATLAB';
                t = Tiff(strcat(saveDir,saveName),'w');
                t.setTag(tagstruct)
                t.write(uint16(mat2gray(O_phase_roi)*2^16));
                t.close();
                %}
            [~,~,rot,~] = decompose_affine_matrix(par.affine_matrix_init);
            object_temp = Ggather(self.object{1});
            object_roi_temp = object_temp(cache.object_ROI{:});        
            O_phase_roi = zeros(size(object_roi_temp,1), size(object_roi_temp,2)*par.Nlayers);
            
            for ll=1:par.Nlayers
                object_temp = Ggather(self.object{ll});
                object_temp = imrotate(object_temp, rot,"bilinear", 'crop' );
                object_roi = object_temp(cache.object_ROI{:});
                x_lb = (ll-1)*size(object_roi,2)+1;
                x_ub = ll* size(object_roi,2);
                O_phase_roi(:,x_lb:x_ub) = phase_unwrap(angle(object_roi));
            end
            
            saveName = strcat('O_phase_roi_Niter',num2str(iter),'.tiff');
            saveDir = strcat(par.fout,'/O_phase_roi/');
            if ~exist(saveDir, 'dir')
                mkdir(saveDir)
            end
            %imwrite(mat2gray(O_phase_roi),strcat(saveDir,fileName),'tiff');

            % save as single tiff
            tagstruct.ImageLength     = size(O_phase_roi,1);
            tagstruct.ImageWidth      = size(O_phase_roi,2);
            tagstruct.Photometric     = Tiff.Photometric.MinIsBlack;
            tagstruct.BitsPerSample   = 16;
            tagstruct.SamplesPerPixel = 1;
            tagstruct.RowsPerStrip    = 16;
            tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
            tagstruct.Software        = 'MATLAB';
            t = Tiff(strcat(saveDir,saveName),'w');
            t.setTag(tagstruct)
            t.write(uint16(mat2gray(O_phase_roi)*2^16));
            t.close();
        end
        
        %% save probe mage
        if isfield(par,'save_probe_mag') && par.save_probe_mag
            %%
            %{
            for ii=1:size(outputs.probe{1},3)
                probe_mag = zeros(size(outputs.probe{1}(:,:,ii,1),1),size(outputs.probe{1}(:,:,ii,1),2)*length(outputs.probe));
                for jj=1:length(outputs.probe)
                    x_lb = (jj-1)*size(outputs.probe{1}(:,:,ii,1),2)+1;
                    x_ub = jj* size(outputs.probe{1}(:,:,ii,1),2);
                    probe_mag(:,x_lb:x_ub) = abs(outputs.probe{jj}(:,:,ii,1));
                end

                saveName = strcat('probe_mag_Niter',num2str(iter),'_',num2str(ii),'th_probe','.tiff');
                saveDir = strcat(par.fout,'/probe_mag/');
                if ~exist(saveDir, 'dir')
                    mkdir(saveDir)
                end
                imwrite(mat2gray(probe_mag)*64,parula,strcat(saveDir,saveName),'tiff')
            end
            %}

            %probe_mag = zeros(size(probe,1),size(probe,2)*par.probe_modes);
            probe_mag = zeros(size(probe,1)*par.Nscans,size(probe,2)*par.probe_modes);

            % modified by dzh: save probe image for every scans
            for i=1:par.Nscans
                y_lb = (i-1)*size(probe,1)+1;
                y_ub = i* size(probe,1);
                for jj=1:par.probe_modes
                    x_lb = (jj-1)*size(probe,2)+1;
                    x_ub = jj* size(probe,2);
                    probe_mag(y_lb:y_ub,x_lb:x_ub) = abs(probe(:,:,jj,i));
                end
            end
            
            saveName = strcat('probe_mag_Niter',num2str(iter),'.tiff');
            saveDir = strcat(par.fout,'/probe_mag/');
            if ~exist(saveDir, 'dir')
                mkdir(saveDir)
            end
            imwrite(mat2gray(probe_mag)*64,parula,strcat(saveDir,saveName),'tiff')
        end
    end
    
end

    %% return results 
    if is_method(par, 'DM')
        % return average from last 10% iterations 
        for ll = 1:length(self.object)
             self.object{ll} = object_avg{ll} /N_object_avg; 
        end
    end


    iter_time = toc(t0) / par.number_iterations; 
    verbose(1,' ====  Time per one iteration %3.3fs', iter_time)    
    verbose(1,' ====  Total time %3.2fs', toc(t0))    

    
    % clip outliers from the low illum regions 
    MAX_OBJ = 0;
    for ii = 1:length(self.object)
        MAX_OBJ = max(MAX_OBJ, max2(abs(self.object{ii}(cache.object_ROI{:}))));
    end
    for ii = 1:length(self.object)
        aobj = abs(self.object{ii});
        self.object{ii} = min(MAX_OBJ, aobj) .* self.object{ii} ./ max(aobj, 1e-3);
    end
    
    %% in case of tilted plane ptychography, back-tilt to the detector probe 
    if any(par.p.sample_rotation_angles(1:2)) && check_option(par.p, 'apply_tilted_plane_correction', 'propagation') 
        % apply propagators to the tilted plane 
        for ii = 1:par.probe_modes
            self.probe{ii} = self.modes{ii}.tilted_plane_propagate_back(self.probe{ii});
        end
    end

    %% in case of multilayer extension assume return the reconstructed probe in middle of the sample !! 
    if par.Nlayers > 1 && par.preshift_ML_probe
       probe_offset = sum(self.z_distance(1:end-1))/2; 
       for ii = 1:par.probe_modes
            self.probe{ii} = utils.prop_free_nf(self.probe{ii}, self.lambda , probe_offset ,self.pixel_size) ;
       end
    end
    
    outputs = self;
    
    % avoid duplication in memory 
    self.noise = [];
    self.diffraction= [];
    self.mask = [];
    
    % store useful parameters back to the main structure only for the 1st mode 
    outputs.relative_pixel_scale = self.modes{mode_id}.scales(end,:);
    outputs.rotation =  self.modes{mode_id}.rotation(end,:);
    outputs.shear =   self.modes{mode_id}.shear(end,:);
    outputs.z_distance = self.modes{mode_id}.distances(end,:);
    outputs.shift_scans = self.modes{mode_id}.shift_scans;
    outputs.probe_fourier_shift = self.modes{mode_id}.probe_fourier_shift; 
    outputs.probe_positions = self.modes{mode_id}.probe_positions;
    outputs.detector_rotation = self.modes{mode_id}.probe_rotation(end,:);
    outputs.detector_scale = 1+self.modes{mode_id}.probe_scale_upd(end);

    if strcmp(par.likelihood, 'poisson')
        fourier_error = fourier_error- fourier_error(1,:);
    end
    fourier_error = Ggather(fourier_error);


    outputs.illum_sum = cache.illum_sum_0; 

    outputs.diffraction = [];
    outputs.noise = [];
    outputs.mask = [];

    %% move everything back to RAM from GPU 
    if par.use_gpu
        outputs =  move_from_gpu(outputs);
    end
    
    %% report results
    try
        verbose(0,'==== REPORT ==== \n SNR %4.3g  RES %3.2g (%3.3gnm) AuC: %3.3g\n\n', fsc_score{end,1}.SNR_avg,...
                    fsc_score{end,1}.resolution,mean(self.pixel_size)*1e9/fsc_score{end,1}.resolution, fsc_score{end,1}.AUC)
    end

end

function x = positivity_constraint_object(x, relax)
   x = relax.*abs(x) + (1-relax).*x;
end

function object = regular_object_out_regions(object, illum, max_illum, delta)
    % push everywhere to zero, even out of the object region
    W = illum / max_illum;
    W = W ./ (0.1+W);
    object = object .* (W + (1-W).*(1-delta));
end

function win = get_window(Np_p, scale, ax)
    % aux function 
    % apodize window for img to prevent periodic boundary errors 
    import engines.GPU_MS.GPU_wrapper.Garray
    win = ones(floor(Np_p(ax)/scale/2-2)*2); 
    win = utils.crop_pad(win, [Np_p(ax),1]);
    win = shiftdim(win, 1-ax);
    win = Garray(win);
end


 
