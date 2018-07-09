% Copyright (c) 2018 NCATS, NIH

% Permission is hereby granted, free of charge, to any person obtaining a 
% copy of this software and associated documentation files (the "Software"),
% to deal in the Software without restriction, including without limitation
% the rights to use, copy, modify, merge, publish, distribute, sublicense, 
% and/or sell copies of the Software, and to permit persons to whom the 
% Software is furnished to do so, subject to the following conditions:

% The above copyright notice and this permission notice shall be included 
% in all copies or substantial portions of the Software.

% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL 
% THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
% FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
% DEALINGS IN THE SOFTWARE.

clear

% this script should be run from matlab command line
% the parfor loop in this script allows parallel analysis of data from each microwell in the plate (data set for one experiment)
% the script calls function 'Fnctn_Boutin2018_3d_Sphere_Seg_v02.m'. this function should be placed in active matlab directory. 

% below line must load matrix with variable name: 'Well_Row_Col_List'
% Well_Row_Col_List contains two columns where each row contains the microwell column and row data for individual microwells to test 
load('Z:\Ty\170829 Parallel Segmentation Scripts and Testing\180330 Rename Variables\Well_Row_Col_List_short_test.mat')

PlanesNum = 60; %number of z-slices in the imaged z-stacks

%input path where the TIFF image files exported from Acapella imaging
%software are stored:
In_Path_Im_1 = 'Z:\Ty\170901_Molly_T47_U87exp1_3\170328-T47-exp3-5DIV__2017-03-28T10_13_56-Measurement 1\Images'
%output path where you would like the generated .mat files for each
%spheroid to be stored:
Out_Path_Data_1 = 'Z:\Ty\170901_Molly_T47_U87exp1_3\170328-T47-exp3-5DIV__2017-03-28T10_13_56-Measurement 1\Wells_Segmented'

parfor WellCounter = 1:size(Well_Row_Col_List,1)
    
     Output_1  = Fnctn_Boutin2018_3d_Sphere_Seg_v02( WellCounter, In_Path_Im_1, Well_Row_Col_List, PlanesNum, Out_Path_Data_1)
    
end
    