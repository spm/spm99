function Vo = spm_write_filtered(Z,XYZ,DIM,M,descrip)
% Writes the filtered SPM as an image
% FORMAT spm_write_filtered(Z,XYZ,DIM,M,descrip)
%
% Z       - {1 x ?} vector point list of SPM values for MIP
% XYZ     - {3 x ?} matrix of coordinates of points (voxel coordinates)
% DIM     - image dimensions {voxels}
% M       - voxels - > mm matrix (used for header & written in *.mat file)
%           [default spm_matrix(-(DIM+1)/2) ]
% descrip - description string [default 'SPM-filtered']
%
%-----------------------------------------------------------------------
%
% spm_write_filtered takes a pointlist image (parallel matrixes of
% co-ordinates and voxel intensities), and writes it out into an image
% file.
%
% It is intended for writing out filtered SPM's from the results
% section of SPM, but can be used freestanding.
%_______________________________________________________________________
% @(#)spm_write_filtered.m	2.7 FIL 99/10/29

%-Parse arguments
%-----------------------------------------------------------------------
if nargin<3, error('Insufficient arguments'), end
if nargin<4, M = spm_matrix(-(DIM+1)/2); end
if nargin<5, descrip='SPM-filtered'; end

%-Get filename
%-----------------------------------------------------------------------
Q       = spm_str_manip(spm_input('Output filename',1,'s'),'sdv');
spm('Pointer','Watch')

%-Set up header information
%-----------------------------------------------------------------------
Vo      = struct(...
		'fname',	Q,...
		'dim',		[DIM', spm_type('uint8')],...
		'mat',		M,...
		'descrip', 	descrip);

%-Reconstruct (filtered) image from XYZ & Z pointlist
%-----------------------------------------------------------------------
Y      = zeros(DIM(1:3)');
OFF    = XYZ(1,:) + DIM(1)*(XYZ(2,:)-1 + DIM(2)*(XYZ(3,:)-1));
Y(OFF) = Z.*(Z > 0);

%-Write the reconstructed volume
%-----------------------------------------------------------------------
Vo = spm_write_vol(Vo,Y);
spm('alert"',{'Written:',['    ',spm_get('CPath',Q)]},mfilename,sqrt(-1));

%-End
%-----------------------------------------------------------------------
spm('Pointer','Arrow');
