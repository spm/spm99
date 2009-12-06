function varargout = spm_getxyz(varargin)
% extracts the finite non-zero coordinates from a volume - a compiled routine
% FORMAT XYZ = spm_getxyz(V)
% V      -  is a memory mapped image volume.
%           Finite and non-zero values are considered to be part of
%           the search volume.
% XYZ    - The co-ordinates (in millimeters) of the voxels.
%_______________________________________________________________________
% @(#)spm_getxyz.m	2.2 John Ashburner 99/04/19

%-This is merely the help file for the compiled routine
error('spm_getxyz.c not compiled - see spm_MAKE.sh')
