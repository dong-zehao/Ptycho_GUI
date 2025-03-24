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

function [outputs, fourier_error, fsc_score] = ptycho_solver_affine(self, par, cache)

import engines.GPU_MS.analysis.*
import engines.GPU_MS.shared.*
import engines.GPU_MS.initialize.*
import engines.GPU_MS.GPU_wrapper.*
import math.*
import utils.*

%% in case of precession mode ptychography, shift each layer by its shifting vector: dzh
if par.precession_mode
    for ilayer = 1:par.Nlayers
        for iscan = 1:par.Nscans
            ind_scan = self.reconstruct_ind{iscan};
            self.modes{ilayer}.probe_positions(ind_scan,:) = self.modes{ilayer}.probe_positions(ind_scan,:) + self.modes{ilayer}.interlayer_shift_position{iscan}./self.pixel_size - par.refine_shift{iscan}./self.pixel_size;
            self.modes{ilayer}.probe_positions_0(ind_scan,:) = self.modes{ilayer}.probe_positions_0(ind_scan,:) + self.modes{ilayer}.interlayer_shift_position{iscan}./self.pixel_size - par.refine_shift{iscan}./self.pixel_size;
        end
    end
end

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

mode_id = 1;  % main mode (assume single most important mode for approchimations) 

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
for iter =  1:par.number_iterations
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
        elseif timeLeft>60
            verbose(0,'Iteration: %i / %i  (Time left:%3.3g min. avg:%3.3g sec)', iter, par.number_iterations, timeLeft/60, avgTimePerIter)
        else
            verbose(0,'Iteration: %i / %i  (Time left:%3.3g sec. avg:%3.3g sec)', iter, par.number_iterations, timeLeft, avgTimePerIter)
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
   
    %% updated illumination
    if iter <= 1 || ( iter > par.probe_change_start && (mod(iter, 10) == 1 || iter < par.probe_change_start+10 ))
        aprobe2 = abs(self.probe{1}(:,:,1)).^2; 
        for ll = 1:size(self.object,1)
            if par.share_object
                ind = [self.reconstruct_ind{:}];
            else
                ind = self.reconstruct_ind{ll};
            end
            % avoid oscilations by adding momentum term 
            cache.illum_sum_0{ll} = set_views(cache.illum_sum_0{ll}, Garray(aprobe2), 1,1,ind, cache)/2;
            cache.illum_norm(ll) = norm2(cache.illum_sum_0{ll});
            cache.MAX_ILLUM(ll) = max2(cache.illum_sum_0{ll});
        end
    end
    
    %% improve convergence speed by gradient acceleration 
    if is_method(par, 'MLc') && iter >= par.accelerated_gradients_start
        [self, cache] = accelerate_gradients(self, par, cache, iter); 
    end
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%  PERFORM ONE ITERATION OF THE SELECTED METHOD %%%%%%%%%%%%%%%%%%%%%
    
    affine_matrix = compose_affine_matrix(scale, asymmetry, rot_ang, shear);
    positions_here = affine_matrix*self.probe_positions_0.';
    p.positions_real = p.positions_real.';
    [~, ~, fourier_error] = engines.GPU_MS.LSQML(self,par,cache,fourier_error,iter);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if iter == 0; continue; end  % interation 0 is used only to calibrate iinitial probe intensity
    
    if verbose() > 0  && any(~isnan(fourier_error(iter,:)))       
        switch lower(par.likelihood)
            case  'l1', verbose(1,'=====  Fourier error = %3.4g ', nanmedian(fourier_error(iter,:)) ); 
            case 'poisson'
                err = fourier_error(iter,:) - fourier_error(1,:);
                verbose(1,'=====  Log likelihood = %3.5g ', nanmedian(err)); 
        end
    end
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


 
