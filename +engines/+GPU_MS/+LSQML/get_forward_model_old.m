% GET_FORWARD_MODEL from the provided object and probe calculate the exit wave 
%
% [self, probe, obj_proj, psi] = get_forward_model(self, obj_proj, par, cache, g_ind, p_ind, scan_ids, layer_ids)
% 
% ** self      structure containing inputs: e.g. current reconstruction results, data, mask, positions, pixel size, ..
% ** obj_proj  [Nx,Ny,N] array, just a preallocated array on GPU, can be empty 
% ** par       structure containing parameters for the engines 
% ** cache     structure with precalculated values to avoid unnecessary overhead
% ** g_ind      indices corresponding to the current group that is solved in parallel 
% ** p_ind      indices containg corresponding probe id for each processed position
% ** scan_ids   determines to which scan correponds each of the position 
% ** layer_ids  id of the solved layer for multilayer ptycho 
%
%
% returns:
% ++ self      structure containing inputs: e.g. current reconstruction results, data, mask, positions, pixel size, ..
% ++ probe     either [Nx,Nx,1] or [Nx,Nx,N] aarray of shared probe or variable probe that differs for each position 
% ++ obj_proj  [Nx,Ny,N] array, views of the object for each scan position 
% ++ psi       [Nx,Ny,N] array, complex valued exit-wave (psi = P*O)
%
% see also: engines.GPU.LSQML 


% Academic License Agreement
%
% Source Code
%
% Introduction 
% •	This license agreement sets forth the terms and conditions under which the PAUL SCHERRER INSTITUT (PSI), CH-5232 Villigen-PSI, Switzerland (hereafter "LICENSOR") 
%   will grant you (hereafter "LICENSEE") a royalty-free, non-exclusive license for academic, non-commercial purposes only (hereafter "LICENSE") to use the cSAXS 
%   ptychography MATLAB package computer software program and associated documentation furnished hereunder (hereafter "PROGRAM").
%
% Terms and Conditions of the LICENSE
% 1.	LICENSOR grants to LICENSEE a royalty-free, non-exclusive license to use the PROGRAM for academic, non-commercial purposes, upon the terms and conditions 
%       hereinafter set out and until termination of this license as set forth below.
% 2.	LICENSEE acknowledges that the PROGRAM is a research tool still in the development stage. The PROGRAM is provided without any related services, improvements 
%       or warranties from LICENSOR and that the LICENSE is entered into in order to enable others to utilize the PROGRAM in their academic activities. It is the 
%       LICENSEE’s responsibility to ensure its proper use and the correctness of the results.”
% 3.	THE PROGRAM IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR 
%       A PARTICULAR PURPOSE AND NONINFRINGEMENT OF ANY PATENTS, COPYRIGHTS, TRADEMARKS OR OTHER RIGHTS. IN NO EVENT SHALL THE LICENSOR, THE AUTHORS OR THE COPYRIGHT 
%       HOLDERS BE LIABLE FOR ANY CLAIM, DIRECT, INDIRECT OR CONSEQUENTIAL DAMAGES OR OTHER LIABILITY ARISING FROM, OUT OF OR IN CONNECTION WITH THE PROGRAM OR THE USE 
%       OF THE PROGRAM OR OTHER DEALINGS IN THE PROGRAM.
% 4.	LICENSEE agrees that it will use the PROGRAM and any modifications, improvements, or derivatives of PROGRAM that LICENSEE may create (collectively, 
%       "IMPROVEMENTS") solely for academic, non-commercial purposes and that any copy of PROGRAM or derivatives thereof shall be distributed only under the same 
%       license as PROGRAM. The terms "academic, non-commercial", as used in this Agreement, mean academic or other scholarly research which (a) is not undertaken for 
%       profit, or (b) is not intended to produce works, services, or data for commercial use, or (c) is neither conducted, nor funded, by a person or an entity engaged 
%       in the commercial use, application or exploitation of works similar to the PROGRAM.
% 5.	LICENSEE agrees that it shall make the following acknowledgement in any publication resulting from the use of the PROGRAM or any translation of the code into 
%       another computing language:
%       "Data processing was carried out using the cSAXS ptychography MATLAB package developed by the Science IT and the coherent X-ray scattering (CXS) groups, Paul 
%       Scherrer Institut, Switzerland."
%
% Additionally, any publication using the package, or any translation of the code into another computing language should cite for difference map:
% P. Thibault, M. Dierolf, A. Menzel, O. Bunk, C. David, F. Pfeiffer, High-resolution scanning X-ray diffraction microscopy, Science 321, 379–382 (2008). 
%   (doi: 10.1126/science.1158573),
% for maximum likelihood:
% P. Thibault and M. Guizar-Sicairos, Maximum-likelihood refinement for coherent diffractive imaging, New J. Phys. 14, 063004 (2012). 
%   (doi: 10.1088/1367-2630/14/6/063004),
% for mixed coherent modes:
% P. Thibault and A. Menzel, Reconstructing state mixtures from diffraction measurements, Nature 494, 68–71 (2013). (doi: 10.1038/nature11806),
% and/or for multislice:
% E. H. R. Tsai, I. Usov, A. Diaz, A. Menzel, and M. Guizar-Sicairos, X-ray ptychography with extended depth of field, Opt. Express 24, 29089–29108 (2016). 
%   (doi: 10.1364/OE.24.029089).
% 6.	Except for the above-mentioned acknowledgment, LICENSEE shall not use the PROGRAM title or the names or logos of LICENSOR, nor any adaptation thereof, nor the 
%       names of any of its employees or laboratories, in any advertising, promotional or sales material without prior written consent obtained from LICENSOR in each case.
% 7.	Ownership of all rights, including copyright in the PROGRAM and in any material associated therewith, shall at all times remain with LICENSOR, and LICENSEE 
%       agrees to preserve same. LICENSEE agrees not to use any portion of the PROGRAM or of any IMPROVEMENTS in any machine-readable form outside the PROGRAM, nor to 
%       make any copies except for its internal use, without prior written consent of LICENSOR. LICENSEE agrees to place the following copyright notice on any such copies: 
%       © All rights reserved. PAUL SCHERRER INSTITUT, Switzerland, Laboratory for Macromolecules and Bioimaging, 2017. 
% 8.	The LICENSE shall not be construed to confer any rights upon LICENSEE by implication or otherwise except as specifically set forth herein.
% 9.	DISCLAIMER: LICENSEE shall be aware that Phase Focus Limited of Sheffield, UK has an international portfolio of patents and pending applications which relate 
%       to ptychography and that the PROGRAM may be capable of being used in circumstances which may fall within the claims of one or more of the Phase Focus patents, 
%       in particular of patent with international application number PCT/GB2005/001464. The LICENSOR explicitly declares not to indemnify the users of the software 
%       in case Phase Focus or any other third party will open a legal action against the LICENSEE due to the use of the program.
% 10.	This Agreement shall be governed by the material laws of Switzerland and any dispute arising out of this Agreement or use of the PROGRAM shall be brought before 
%       the courts of Zürich, Switzerland. 


function [self, probe, obj_proj, psi] = get_forward_model(self, obj_proj, par, cache, g_ind, p_ind, scan_ids, layer_ids)
    import engines.GPU.shared.*
    import engines.GPU.GPU_wrapper.*
    import engines.GPU.LSQML.*
    import math.*
    import utils.*
    import plotting.*

    if isempty(obj_proj{1})
        for ll = 1:par.object_modes
            obj_proj{ll} = Gzeros([self.Np_p, 0], true);
        end
    end
    
% % allocate memory first by ZC
% probe = self.probe;
    if par.Nlayers > 1
        probe=cell(par.probe_modes,par.Nlayers+1);
    end
% % % % %     probe{1,1}=self.probe;

    
    % get illumination probe 
    for ll = 1:par.probe_modes
        if (ll == 1 && (par.variable_probe || par.variable_intensity))
            % add variable probe (OPRP) part into the constant illumination 
            probe{ll,1} =  get_variable_probe(self.probe{ll}, self.probe_evolution(g_ind,:),p_ind{ll});
        else
            % store the normal (constant) probe(s)
            probe{ll,1} = self.probe{min(ll,end)}(:,:,min(end,p_ind{ll}),1);
        end

        if (ll == 1 && par.apply_subpix_shift && isinf(self.z_distance(end)))  || is_used(par,'fly_scan')
            % only in farfield mode 
            probe{ll,1} = apply_subpx_shift(probe{ll,1}, self.modes{min(end,ll)}.sub_px_shift(g_ind,:) );
        end
        if (ll == 1)
            probe{ll,1} = apply_subpx_shift_fft(probe{ll,1}, self.modes{1}.probe_fourier_shift(g_ind,:)); 
        end
    end
    
%     % for debug
%     temp_p=probe;
% allocate by ZC
   psi=cell(max(par.object_modes, par.probe_modes),1);

   % get projection of the object and probe 
   for layer = 1:par.Nlayers
       for ll = 1:max(par.object_modes, par.probe_modes)
           llo = min(ll, par.object_modes); 
           llp = min(ll, par.probe_modes); 
            % get objects projections 
            obj_proj{llo} = get_views(self.object, obj_proj{llo},layer_ids(layer),llo, g_ind, cache, scan_ids,[]);
            if (ll == 1 && par.apply_subpix_shift && ~isinf(self.z_distance(end)))
                % only in nearfield mode , apply shift in the opposite direction 
                obj_proj{ll} = apply_subpx_shift(obj_proj{ll} .* cache.apodwin, -self.modes{min(end,ll)}.sub_px_shift(g_ind,:) ) ./ cache.apodwin;
            end
            
            % get exitwave after each layer
            psi{ll} =probe{llp,layer} .* obj_proj{llo};
            % fourier propagation  
            [psi{ll}] = fwd_fourier_proj(psi{ll} , self.modes{layer}, g_ind);  
            if par.Nlayers > 1
                 probe{llp,layer+1} = psi{llp};
            end
       end
   end
   % multilayer but not count final inf layer, then an additional farfield
   % fft needed, by Zhen Chen
   if par.Nlayers > 1 && isinf(self.z_distance(end)) && length(self.z_distance) > par.Nlayers %&& par.rmvac 
       for ll= 1:max(par.object_modes, par.probe_modes)
           psi{ll} = fft2_safe(psi{ll});  % fully farfield 
           probe{ll,par.Nlayers+1} = psi{ll};
       end
   end       
   
            % debug by Zhen Chen
%         temp=gather(psi{1});
%         if any(isnan(temp(:))) || any(temp(:)> 1e3)
%             keyboard;
%         end
        

       
end
