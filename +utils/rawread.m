function [cbed]=rawread(fname,irot180)
% read raw dataset from EMPAD by Zhen Chen @ Cornell University
% output: cbed, reshape to 4D dataset after crop the edge bad pixels
%         data format: 124 x 124 x N_scan_vertical x N_scan_horizontal.
% Inputs:
%       fname: file name of the *.raw dataset
%       irot180 = 1, rotate scan and diffraction 180 degrees to make the
%       final image the same orientation as those generated from the PAD acquisition.

if nargin < 2
    irot180 = 1;
end
if nargin < 1
    [fname,fdir]=uigetfile('*.raw','Open a *.raw file');
    fname=fullfile(fdir,fname);
end

x_cbed = 128;
y_cbed = 130;  % default size of PAD diffraction

fid=fopen(fname,'r');
fin = fread(fid,'*float32');
fclose(fid);

fnm_parts=strsplit(fname,'_');
npy=str2double(fnm_parts{end}(2:end-4)); %real space size
npx=str2double(fnm_parts{end-1}(2:end));

fun= reshape(fin,x_cbed,y_cbed,npx,npy); % reshape for data
clear fin;
cbed = single(fun(3:126,3:126,:,:)); % crop edges
clear fun;
    
if irot180 == 1
    cbedrot=rot90(cbed,2); % rotate CBED 180d
    clear cbed;
    cbedrot=permute(cbedrot,[3,4,1,2]); 
    cbedrot=rot90(cbedrot,2);   % rotate scanning 180d
      
    cbed=permute(cbedrot,[3,4,2,1]); % transpose scan x/y, as there is a storage difference. scanning along  2nd dimension of matrix, horizontal if plotting by imagesc
    clear cbedrot;
else
    cbed=permute(cbed,[1,2,4,3]); % only transpose scan x/y
end
    

end