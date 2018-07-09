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

%%% Run this script by pasting entire code into Matlab command line.

%%% This script loads previous 3D region information regarding individual
%%% nuclei regions, then this script segments the 'whole spheroid region',
%%% and creates measurements/statistics that relate the behavior of the individual nuclei
%%% to the overall spheroid structure.

%%%%% This would be a simple segmentation of whole spheroid, but it is
%%%%% complicated by the requirements of this specicific experiment. For
%%%%% this experiment, the same algorithm must segement 'whole spheroid
%%%%% region' from uncleared (highly blurred, dispersed light, wide and weak edge features) spheroids 
%%%%% and also segment 'whole spheroid region' from spheroids following the optical clearing procedure (much better signal to noise, strong edge features). 

%Historical notes
%%**Adds in the ability to calculate ROI area/z slice area from 170824. 
%**Adds in ability to determine "clearing metric", as 50%smoothed#ROI/area
%over spheroid height (2x z with max area) 170825
%**Adds in variable ClearingQualityMetric_Matrix which has all clearing
%metric calculated to be >1 set to 1 170901
%**Modifies clearing metric to be 50% threshold of ROIarea/area over
%spheroid height (2x z with max area) 170905
%**Adjustment such that the max ROI area/area can only be
%found between slice 0 and x% down the ROI area curve. Gets rid of
%artifacts from when max ROI area/area is found to be after the spheroid
%ends. 
%*Compensate for condition where sphere is too large for detection of 
% x% down from the ROI area/area postion 170907, sets those ROI height=inf
%Outputs a modified s9 variable that is cut at the threshold. 
 
clear  
% Directory where previously saved Matlab data, 3D nuclei segmentation, one
% Matlab data file for each micro-well, presumably one spheroid per well
% directory also contains matlab data file with matrix of mirco-well
% coordinates 'Well_Row_Col_List_short_test.mat', this is variable
% 'Well_Row_Col_List'

%%% NOTICE!!! 3rd column of variable 'Well_Row_Col_List' must contain number for experimental condition ID of each micro-well,
%%% This ID number may need to be manually added, or could be added by another informatics script 

Dir_Temp_1 = 'Z:\Ty\170810_Molly_T47exp2_U87exp4\170324-T47-exp2-3DIV__2017-03-24T09_43_22-Measurement 1\test_Molly_180604'
cd(Dir_Temp_1)

load('Well_Row_Col_Conditions_postseg_col3only.mat') %loads in Well_Row_Col_List variable (column 1 = microwell row, column 2 = microwell column, column 3 = experimental condition ID)

Condition_List = unique(Well_Row_Col_List(:,3));
idx = find(Condition_List(:,1) > 0); %determines the number of conditions

Condition_List = Condition_List(idx,1);

% Replicate_Number, defines which well to call when multiple wells for each experimental condition in experiment 
% C1 defines the position in the 'Condition List'

TotalFieldCounter = 0;
clear All_Wells_WholeSphereROI_7_Volume
for C1 = 1:size(Condition_List,1)
% C1 = 2 %can modify this to run by hand
    Condition = Condition_List(C1);
    clear Temp_Concatenate_s9
    
    idx = find(Well_Row_Col_List(:,3)==C1);
    Replicate_Number_List = Well_Row_Col_List(idx,:);
    
    for C2 = 1:size(Replicate_Number_List,1)
    % C2 = 13 %can modify this to run by hand
    Replicate_Number = C2;

    Well_Row_Col_List_Single_Conditon = Replicate_Number_List;

    Well_Row_Col = Well_Row_Col_List_Single_Conditon(Replicate_Number,:);

        Well_Row = Well_Row_Col(1,1);
        Well_Column = Well_Row_Col(1,2);
        % Well_Row = 2, Well_Column =11 %can modify this to run by hand
  
        % run below (), searches drive for available files
        file_list_1 = ls;
        file_list_1 = cellstr(file_list_1);

        search_key1 = strcat('_r',num2str(Well_Row), '_c', num2str(Well_Column),'_'); %MB made this edit '_' on 170908
        % search_key1 = strcat('_r',num2str(Well_Row), '_c', num2str(Well_Column));
        % search_key1 = '_r6_c3' %can modify this to run by hand

        index = strfind(file_list_1,search_key1);
        index = find(~cellfun(@isempty,index));
        file_list_2 = file_list_1(index,1);
            
            if size(file_list_2,1) > 0
        
                load(char(file_list_2(1,1))); % need to run this to load the image analysis file in
%                 figure, imshow(In_Stack_2(:,:,35),[])
%                 figure, imshow(In_Stack_3(:,:,22),[])
%                 figure, imshow(WholeSphereROI_7(:,:,50),[])

SubVol_Seg_15_bw1 = false(size(SubVol_Seg_15));
SubVol_Seg_15_bw1(SubVol_Seg_15 > 1) = 1;
SubVol_Seg_15_bw2 = imresize(SubVol_Seg_15_bw1, 0.125); % factor 1/8 XY pixel size vs Z slice interval
      
        SubVol_Seg_15_bw2_3D_strel_radius = 9;

        SubVol_Seg_15_bw2_3D_strel = strel('sphere', SubVol_Seg_15_bw2_3D_strel_radius);

        SubVol_Seg_15_bw2b = padarray(SubVol_Seg_15_bw2, [SubVol_Seg_15_bw2_3D_strel_radius+1 SubVol_Seg_15_bw2_3D_strel_radius+1 SubVol_Seg_15_bw2_3D_strel_radius+1]);

        SubVol_Seg_15_bw3 = imclose(SubVol_Seg_15_bw2b,SubVol_Seg_15_bw2_3D_strel); 

        SubVol_Seg_15_bw3b = SubVol_Seg_15_bw3((SubVol_Seg_15_bw2_3D_strel_radius+2):end-(SubVol_Seg_15_bw2_3D_strel_radius+1),(SubVol_Seg_15_bw2_3D_strel_radius+2):end-(SubVol_Seg_15_bw2_3D_strel_radius+1),(SubVol_Seg_15_bw2_3D_strel_radius+2):end-(SubVol_Seg_15_bw2_3D_strel_radius+1));
        % figure, imshow(SubVol_Seg_15_bw3(:,:,15),[])


%%%%% 171003 -- redefine WholeSphereROI_7, remove small debris at bottom of well and
%%%%% increase whole spheroid segmentation robustness in upper z planes
        In_Stack_2_reserve = In_Stack_2;
        In_Stack_2_resizeXY = imresize(In_Stack_2, 0.125);
        
         
        In_Stack_2_MedFltr = uint16(zeros(size(In_Stack_2_resizeXY)));
        for J1 = 1:size(In_Stack_2_resizeXY,3)
        In_Stack_2_MedFltr(:,:,J1) = medfilt2(In_Stack_2_resizeXY(:,:,J1),[21 21], 'zeros');
        end
        
          In_Stack_2_MedFltr_Adj = uint16(zeros(size(In_Stack_2_resizeXY)));
      for J1 = 1:size(In_Stack_2_resizeXY,3)
          In_Stack_2_MedFltr_StretchLim = stretchlim(In_Stack_2_MedFltr(:,:,J1),0.0001);
           In_Stack_2_MedFltr_Adj(:,:,J1) = imadjust(In_Stack_2_MedFltr(:,:,J1),In_Stack_2_MedFltr_StretchLim);
           In_Stack_2_MedFltr_Adj_MaxList(J1,1) = max(max(In_Stack_2_MedFltr_Adj(:,:,J1)));
      end
      
        In_Stack_2_MedFltr_2 = uint16(zeros(size(In_Stack_2_resizeXY)));
        In_Stack_2_MedFltr_2_Avg = uint16(zeros(size(In_Stack_2_resizeXY)));
        for J1 = 1:size(In_Stack_2_resizeXY,3)
        In_Stack_2_MedFltr_2(:,:,J1) = medfilt2(In_Stack_2_resizeXY(:,:,J1),[21 21], 'symmetric');
        In_Stack_2_MedFltr_2_Avg(:,:,J1) = imfilter(In_Stack_2_MedFltr_2(:,:,J1), fspecial('disk',3),'symmetric');
        end
       
       In_Stack_2_MedFltr_2_Avg(:,:,J1) = imfilter(In_Stack_2_MedFltr_2(:,:,J1), fspecial('disk',3) );
        
       E1 = imfilter(In_Stack_2_MedFltr_2, fspecial('laplacian'), 'symmetric');
        %figure, imshow(E1(:,:,3),[])
        
       E1_Avg = imfilter(In_Stack_2_MedFltr_2_Avg, fspecial('laplacian'), 'symmetric');
        %figure, imshow(E1_Avg(:,:,3),[])
        
       E1_Med3D = medfilt3(E1_Avg, [3 3 3], 'zeros');
       %figure, imshow(E1_Med3D(:,:,3),[])
      
      for J1 = 1:size(In_Stack_2_resizeXY,3)
          if size(E1_Med3D(E1_Med3D(:,:,J1)>0), 1)
           E1_StretchLim = stretchlim(E1_Med3D(E1_Med3D(:,:,J1) > 0),0.001);
           E1(:,:,J1) = imadjust(E1_Med3D(:,:,J1),E1_StretchLim);
          end 
      end
      
      % figure, imshow(E1(:,:,3),[])
 
       clear  E1_GlobalThresh
       for J1 = 1:size(In_Stack_2_resizeXY,3)
            E1_PixList = E1(:,:,J1);
           E1_GlobalThresh(J1,1) = graythresh(E1_PixList(E1_PixList > 0));

       end
    
        E2 = false(size(In_Stack_2_resizeXY));
        for J1 = 1:size(In_Stack_2_resizeXY,3)
            E2(:,:,J1) = imbinarize(E1(:,:,J1),(E1_GlobalThresh(J1,1)*1));
         %   Slice = J1, figure, imshow(In_Stack_2_resizeXY(:,:,Slice),[]),title(num2str(Slice)),figure, imshow(E2(:,:,Slice),[]),title(num2str(Slice))
        end
        
        E3 = imdilate(E2, strel('sphere',3));
        % Slice = 45, figure, imshow(In_Stack_2_resizeXY(:,:,Slice),[]),title(num2str(Slice)),figure, imshow(E3(:,:,Slice),[]),title(num2str(Slice))
        
        clear E5 E9
        E6 = uint16(zeros(size(In_Stack_2_resizeXY)));
        E4_Adj = uint16(zeros(size(In_Stack_2_resizeXY)));
        
        for J1 = 1:size(In_Stack_2_resizeXY,3)
            % J1 = 10 %can modify this to run by hand
            E4 = In_Stack_2_MedFltr_Adj(:,:,J1);
            
            if max(max(max(E3(:,:,J1)))) == 1
                
                
                [E5(J1,1) E5(J1,2)] = graythresh(E4(E3(:,:,J1)==1));
                E5(J1,3) = E5(J1,1) * 2^16; 
                E6(:,:,J1) = imbinarize(E4,E5(J1,1));
                
                 E6_2D_strel_radius = 5;

                    E6_2D_strel = strel('disk', E6_2D_strel_radius);

                E6_b = padarray(E6, [E6_2D_strel_radius+1 E6_2D_strel_radius+1 0],0);
                E4_b = padarray(E4, [E6_2D_strel_radius+1 E6_2D_strel_radius+1 0],'replicate');
                
                
                E7 = imdilate(E6_b(:,:,J1), E6_2D_strel);
                E8 = E7-E6_b(:,:,J1);
                
                E9(J1,1) = mean(E4_b(E6_b(:,:,J1)==1))/mean(E4_b(E8==1));
                
                % set infinities or NaNs to zero, these are invalid, only due to debris in image
                if E9(J1,1) == Inf || isnan(E9(J1,1)) ==1
                    E9(J1,1) = 0;
                end
            else
            end
          % Slice = J1, figure, imshow(In_Stack_2_resizeXY(:,:,Slice),[]),title(num2str(Slice)),figure, imshow(E6(:,:,Slice),[]),title(num2str(Slice))    
        end
        %%%% Slice = 30, figure, imshow(In_Stack_2_resizeXY(:,:,Slice),[]),figure, imshow(E6(:,:,Slice),[])
        
        WholeSpheroid_RatioIntensityThreshold = 2.6; %find peaks bigger than 2.6, this is an intensity ratio
  
        
         if max(E6(:)) > 0 && max(E9(:)) > WholeSpheroid_RatioIntensityThreshold % conditional statement makes sure that a potential 'whole spheroid' region is found
    
             
                E9(:,2) = 1:size(E9,1);
                

                E9_MedFilt = medfilt2(E9(:,1),[5 1]); % currently set at 5 slice window, may need to be adjusted. bigger number = more smoothing
                
            if max(E9_MedFilt(:)) > WholeSpheroid_RatioIntensityThreshold % conditional statement makes sure that a potential 'whole spheroid' region is found
                
                E9_RegionMax = imregionalmax(E9_MedFilt(:,1));
                E9b = E9(E9_RegionMax==1,:);
                E9c = E9b(E9b(:,1) > WholeSpheroid_RatioIntensityThreshold,:); 
              
                %%%% LA variables (below section) are for line analysis of diameter of
                %%%% whole sphere at each Z-plane. Mid plane of sphere has
                %%%% the largest diameter of associated 2D region.
                %%%% Automated graphical analysis finds this feature.
                
                %%%% Finding the mid plane of the sphere allows rough
                %%%% estimation of total sphere z height.
                
                LA_1 = max(E9c(:,1)); 
                LA_1= E9c(find(E9c(:,1)==LA_1),2);
                
                % added 171020,makes sure that if multiple values for regional max (a vector with size > 1) then select first value in the vector
                clear LA_1b
                LA_1b = LA_1(1,:);  

                
                LA_1_rangeMid = 3; % range of window
                LA_2 = LA_1b - LA_1_rangeMid;
                if LA_2 < 1
                    LA_2 = 1;
                end

               LA_3 = LA_1b + LA_1_rangeMid;
                if LA_3 > size(E9,1)
                    LA_3 = size(E9,1);
                end

                clear LA_SliceArea_1
                for LA_C1 = 1:size(E9,1)
                LA_SliceArea_1(LA_C1,1) = size(find(E6(:,:,LA_C1)> 0),1);
                LA_SliceArea_1(LA_C1,2) = LA_C1;
                end

                LA_4 = LA_SliceArea_1(LA_2:LA_3,:);

                % mid point center z plane of sphere 
                LA_5 = LA_4(find(LA_4(:,1) == max(LA_4(:,1))),:);
                LA_5 = LA_5(1,2);

                % Right tail (high z planes) of sphere, local minimum in area per slice
                LA_6 = LA_SliceArea_1(LA_5:size(LA_SliceArea_1,1),:);
                LA_7 = LA_6(find(LA_6(:,1) == min(LA_6(:,1))),:);
                LA_7 = LA_7(1,2);

                LA_9 = E6;

                    LA_8 = graythresh(In_Stack_2_MedFltr_Adj(LA_9(:,:,LA_7) == 1)); 
                    LA_9(:,:,LA_7+1) = imbinarize(In_Stack_2_MedFltr_Adj(:,:,LA_7),LA_8);
                    LA_9_SliceArea_1_Holder = size(find(LA_9(:,:,LA_7)> 0),1);
                    % figure, imshow(In_Stack_2_MedFltr_Adj(:,:,41),[]);
                clear LA_9_SliceArea_1 LA_9_SliceArea_1_Holder_Zstop_TopSlice LA_8list
                
                for LA_C1 = LA_7+1:size(In_Stack_2_resizeXY,3)
                    LA_8 = graythresh(In_Stack_2_MedFltr_Adj(LA_9(:,:,LA_C1-1) == 1)); %In_Stack_2_MedFltr_Adj(LA_9(:,:,LA_C1-1) == 1)
                    LA_8list(LA_C1,1) = LA_8;
                    LA_9(:,:,LA_C1) = imbinarize(In_Stack_2_MedFltr_Adj(:,:,LA_C1),LA_8);
                    LA_9_SliceArea_1(LA_C1,1) = size(find(LA_9(:,:,LA_C1)> 0),1);

                    if LA_9_SliceArea_1(LA_C1,1) < LA_9_SliceArea_1_Holder
                        LA_9_SliceArea_1(LA_C1,2) = 1;
                        LA_9_SliceArea_1_Holder = LA_9_SliceArea_1(LA_C1,1);
                        if LA_9_SliceArea_1(LA_7+1:LA_C1,2) > 0
                            LA_9_SliceArea_1_Holder_Zstop_TopSlice = LA_C1;
                        else
                            LA_9(:,:,LA_C1) = 0;
                        end
                    else
                        LA_9_SliceArea_1(LA_C1,2) = 0;
                        LA_9(:,:,LA_C1) = 0;
                    end

                    LA_9b = imdilate(LA_9(:,:,LA_C1), strel('disk',5));
                    LA_9c = LA_9b-LA_9(:,:,LA_C1);
                    LA_9d(LA_C1,1) = mean(In_Stack_2_MedFltr_Adj(LA_9(:,:,LA_C1)==1))/mean(In_Stack_2_MedFltr_Adj(LA_9c==1));

                end

                % left tail (low z planes) of sphere, local minimum in area per slice
                LA_10 = LA_SliceArea_1(1:LA_5(1,1),:);
                LA_11 = LA_10(LA_10(:,1)>0,:);
                LA_12 = LA_11(find(LA_11(:,1) == min(LA_11(:,1))),:);
                LA_12 = LA_12(1,2);

                % detected minimum z slice often contains much blur in image... so, start 1 slice above detected minimum and work backward toward the beginning slice
                LA_12 = LA_12+1;

                    LA_8 = graythresh(In_Stack_2_MedFltr_Adj(LA_9(:,:,LA_12) == 1)); 
                    LA_9(:,:,LA_12-1) = imbinarize(In_Stack_2_MedFltr_Adj(:,:,LA_12),LA_8);
                    LA_9_SliceArea_1_Holder = size(find(LA_9(:,:,LA_12)> 0),1);

                clear LA_9_SliceArea_1_Holder_Zstop_BottomSlice LA_8list
                for LA_C1 = LA_12-1:-1:1
                    LA_9_SinglePlane = LA_9(:,:,LA_C1+1);
                    LA_9_SinglePlane_Dilate = imdilate(LA_9_SinglePlane,strel('disk',9));

                    [LA_8, LA_8b] = graythresh(In_Stack_2_MedFltr_Adj(LA_9_SinglePlane_Dilate == 1));
                    LA_8list(LA_C1,1) = LA_8;
                    LA_8list(LA_C1,2) = LA_8b;
                    LA_9(:,:,LA_C1) = imbinarize(In_Stack_2_MedFltr_Adj(:,:,LA_C1),LA_8);
                    LA_9_SliceArea_1(LA_C1,1) = size(find(LA_9(:,:,LA_C1)> 0),1);

                    if LA_9_SliceArea_1(LA_C1,1) < LA_9_SliceArea_1_Holder
                        LA_9_SliceArea_1(LA_C1,2) = 1;
                        LA_9_SliceArea_1_Holder = LA_9_SliceArea_1(LA_C1,1);
                        if LA_9_SliceArea_1(LA_C1:LA_12-1,2) > 0
                            LA_9_SliceArea_1_Holder_Zstop_BottomSlice = LA_C1;
                        else
                            LA_9(:,:,LA_C1) = 0;
                        end
                    else
                        LA_9_SliceArea_1(LA_C1,2) = 0;
                        LA_9(:,:,LA_C1) = 0;
                    end
                    LA_9b = imdilate(LA_9(:,:,LA_C1), strel('disk',5));
                    LA_9c = LA_9b-LA_9(:,:,LA_C1);
                    LA_9d(LA_C1,1) = mean(In_Stack_2_MedFltr_Adj(LA_9(:,:,LA_C1)==1))/mean(In_Stack_2_MedFltr_Adj(LA_9c==1));
                end

                lg_1 = LA_9;

                % 171018 replace bottom planes
                lg_1(:,:,1:10) = SubVol_Seg_15_bw3b(:,:,1:10);

        % figure, imshow(lg_1(:,:,1),[]),figure, imshow(lg_1(:,:,3),[]),figure, imshow(lg_1(:,:,5),[])
        % figure, imshow(In_Stack_2_resizeXY(:,:,1),[]),figure, imshow(In_Stack_2_resizeXY(:,:,3),[]),figure, imshow(In_Stack_2_resizeXY(:,:,5),[])

                 lg_2D_strel_radius = 11;
                 lg_3D_strel_radius = 3;

                 lg_2D_strel = strel('disk', lg_2D_strel_radius);
                  lg_3D_strel = strel('sphere', lg_3D_strel_radius);

                lg_1_b = padarray(lg_1, [lg_2D_strel_radius+1 lg_2D_strel_radius+1 0]);

                lg_2 = imclose(lg_1_b, lg_2D_strel);
                lg_3 = imclose(lg_2,lg_3D_strel);

                lg_3b = lg_3((lg_2D_strel_radius+2):end-(lg_2D_strel_radius+1),(lg_2D_strel_radius+2):end-(lg_2D_strel_radius+1),:);

                WholeSphereROI_7_reserveFromFileLoad = WholeSphereROI_7;
                WholeSphereROI_7 = imresize(lg_3b, 8);
                
                if max(WholeSphereROI_7(:)) > 0 %  make sure that WholeSphereROI_7 has a least one defined region
                
            %WSC: whole spheroid check
                WSC_1 = regionprops(WholeSphereROI_7);
                
                if WSC_1.Area > 100000 % check features of potential whole spheroid, for now it at least needs to be big, can add more criteria if needed 
              

        %%%%% below section calculates new 3D distance to spheroid edge
        %%%%% values from each nucleus (nuc segmentation done in previous
        %%%%% step). Nuc stats were loaded from .mat file into variable 's9'
        
                s9 = rmfield(s9, 'DistCentroidToSphereEdge'); % remove old 3D distance values that were loaded from saved segmentation output file 
                WholeSphereROI_7 = logical(WholeSphereROI_7);
                s6b = imcomplement(WholeSphereROI_7);

                % enlarge and interpolate 3d matrix
                M = double(s6b);
                k = 8;

                   d=size(M);
                z=1:d(3);
                zi=1:1/k:d(3);
                M=permute(interp1(z,permute(M,[3 1 2]),zi),[2 3 1]);

                s6b = logical(M);  
                clear M

                s6b(:,:,1) = 1; % clearing bottom slice on purpose, need bottom spheroid edge to be reliable

                s7 = bwdist(s6b);

                %%% enlarge and interpolate 3d matrix
                M = double(s7);
                k = .125;

                   d=size(M);
                z=1:d(3);
                zi=1:1/k:d(3);
                M=permute(interp1(z,permute(M,[3 1 2]),zi),[2 3 1]);

                s8 = uint16(M);  
                clear M s7
                
                for J1 = 1:size(s9,1)
                s9(J1).DistCentroidToSphereEdge = s8(uint16([s9(J1).Centroid(2)]), uint16([s9(J1).Centroid(1)]), uint16([s9(J1).Centroid(3)]));
                end
                
        %%%%%
        %%%%%
        %%%%%
                    % Determine cut-off for spheroid z height

                        Centroid_Matrix = reshape([s9.Centroid], 3,size(s9,1))';
                       % MeasurementToColor = [s9.DistCentroidToSphereEdge]'; figure, scatter3(Centroid_Matrix(:,1),Centroid_Matrix(:,2),Centroid_Matrix(:,3),3,MeasurementToColor)

                        for V1 = 1: size(In_Stack_2,3)
                            % V1 = 1
                            SpheroidVolume_VoxelCounts_PerSlice(V1,1) = V1;
                            %sphere roi area per Z plane
                            SpheroidVolume_VoxelCounts_PerSlice(V1,2) = size(find(WholeSphereROI_7(:,:,V1)> 0),1);
                            %nuclei roi count per Z plane
                            SpheroidVolume_VoxelCounts_PerSlice(V1,3) = size(find(uint16(Centroid_Matrix(:,3)) == V1),1);
                            %nuclei roi count normalized for spheroid area per Z plane
                            SpheroidVolume_VoxelCounts_PerSlice(V1,4) = SpheroidVolume_VoxelCounts_PerSlice(V1,3)/SpheroidVolume_VoxelCounts_PerSlice(V1,2);
                            % ROI area per Z plane
                            SpheroidVolume_VoxelCounts_PerSlice(V1,5) = size(find(SubVol_Seg_15(:,:,V1) > 2),1);
                            % ROI area per Z plane normalized for spheroid area per Z plane
                            SpheroidVolume_VoxelCounts_PerSlice(V1,6) = SpheroidVolume_VoxelCounts_PerSlice(V1,5)/SpheroidVolume_VoxelCounts_PerSlice(V1,2);

                        end
                        
        %                 figure, plot(SpheroidVolume_VoxelCounts_PerSlice(:,1),SpheroidVolume_VoxelCounts_PerSlice(:,2)), title('sphere roi area per Z plane')        
        %                 figure, plot(SpheroidVolume_VoxelCounts_PerSlice(:,1),SpheroidVolume_VoxelCounts_PerSlice(:,3)), title('nuclei roi count per Z plane')
        %                 figure, plot(SpheroidVolume_VoxelCounts_PerSlice(:,1),SpheroidVolume_VoxelCounts_PerSlice(:,4)), title('nuclei roi count normalized for spheroid area per Z plane')  
        %                 figure, plot(SpheroidVolume_VoxelCounts_PerSlice(:,1),SpheroidVolume_VoxelCounts_PerSlice(:,5)), title('ROI area per Z plane')
        %                 figure, plot(SpheroidVolume_VoxelCounts_PerSlice(:,1),SpheroidVolume_VoxelCounts_PerSlice(:,6)), title('ROI area per Z plane normalized for spheroid area per Z plane')  

        
                        %%%% finding the maximum of the spheroid area per slice over z curve
                        SphereRegionAreaPerSlice = SpheroidVolume_VoxelCounts_PerSlice(:,2); 

                        SphereRegionAreaPerSlice_Max = max(SphereRegionAreaPerSlice(:));
                        
                        SphereRegionAreaPerSlice_Max_Idx = find(SphereRegionAreaPerSlice == SphereRegionAreaPerSlice_Max); %finding the z slice corresponding to the max area

                        SphereRegionAreaPerSlice_PostMaxEnd = SphereRegionAreaPerSlice(SphereRegionAreaPerSlice_Max_Idx(1,1):end,1);%makes vector with area values from z max to z end

                        %output of below 2 lines is not currently being used in metric
                        SphereRegionAreaPerSlice_CutOff = SphereRegionAreaPerSlice_Max * 0.1; %%where the size cutoff is set 10
                        SphereRegionAreaPerSlice_CutOff_Idx = find(SphereRegionAreaPerSlice_PostMaxEnd < SphereRegionAreaPerSlice_CutOff)+(SphereRegionAreaPerSlice_Max_Idx(1,1)-1);
                        
                        %%%% calculate value that creates deepest (z slice) limit where 
                        %%%% max nuc per spheroid value can be detected, this is based on a defined point (z slice) deeper than where max spheroid area per slice is detected 
                        %%%%  (finding the z slice that corresponds to x % down from the max area)
                        SphereRegionAreaPerSlice_CutOff_ForLimit_QualMetric = SphereRegionAreaPerSlice_Max * 0.97; %%where the size cutoff is set 97%
                        SphereRegionAreaPerSlice_CutOff_Idx_ForLimit_QualMetric = find(SphereRegionAreaPerSlice_PostMaxEnd < SphereRegionAreaPerSlice_CutOff_ForLimit_QualMetric)+(SphereRegionAreaPerSlice_Max_Idx(1,1)-1);

                       % compensate for condition where sphere is too large for
                       % detection of max size z slice postion
                         if size(SphereRegionAreaPerSlice_CutOff_Idx_ForLimit_QualMetric,1) == 0 
                             %%%% calculate value for height of sphere
                             Sphere_ApproxHeight = double(inf);
                             % set quality detection limit for full size z stack
                             SphereRegionAreaPerSlice_CutOff_Idx_ForLimit_QualMetric = size(In_Stack_2,3);
                         else
                         %%%% calculate value for height of sphere
                         Sphere_ApproxHeight = double(2*SphereRegionAreaPerSlice_Max_Idx(1,1));    
                         end

                    %%%%%%%%%%%
                        % column 6 "total nuc area per slice area"
                        NucCountPerSliceArea = SpheroidVolume_VoxelCounts_PerSlice(:,6);
                        %below loop smooths
                        for V1 = 2:size(In_Stack_2,3)-1
                            NucCountPerSliceArea_Smooth(V1+1,1) = mean(NucCountPerSliceArea((V1-1):(V1+1),1));
                        end
                        NucCountPerSliceArea_Smooth(1,1)= NucCountPerSliceArea(1,1);
                        NucCountPerSliceArea_Smooth(size(In_Stack_2,3),1) = NucCountPerSliceArea(size(In_Stack_2,3),1);

                        NucCountPerSliceArea_Smooth_Reserve = NucCountPerSliceArea_Smooth;

                         NucCountPerSliceArea_Smooth = NucCountPerSliceArea_Smooth(1:SphereRegionAreaPerSlice_CutOff_Idx_ForLimit_QualMetric(1,1)); % modify range based on z slice where max sphere area per slice is detected
                        
                        % below gives plot of smoothed up until the threshold (97% below max area)
                        % figure, plot(NucCountPerSliceArea_Smooth), title('smoothed nuclei roi count normalized for spheroid area per Z plane')  
                        
                        NucCountPerSliceArea_Smooth_Max = max(NucCountPerSliceArea_Smooth(:));
                        NucCountPerSliceArea_Smooth_Max_Idx = find(NucCountPerSliceArea_Smooth == NucCountPerSliceArea_Smooth_Max);

                        NucCountPerSliceArea_Smooth = NucCountPerSliceArea_Smooth_Reserve;
                        
                        %below gives full smoothed plot that the 50% down is pulled from
                        % figure, plot(NucCountPerSliceArea_Smooth), title('smoothed nuclei roi count normalized for spheroid area per Z plane')  

                        NucCountPerSliceArea_Smooth_PostMaxEnd = NucCountPerSliceArea_Smooth(NucCountPerSliceArea_Smooth_Max_Idx(1,1):end,1);
                        NucCountPerSliceArea_CutOff = NucCountPerSliceArea_Smooth_Max * 0.50; %%where the quality cutoff is set 50

                        NucCountPerSliceArea_CutOff_Idx = find(NucCountPerSliceArea_Smooth_PostMaxEnd < NucCountPerSliceArea_CutOff)+(NucCountPerSliceArea_Smooth_Max_Idx(1,1)-1);

                        if size(NucCountPerSliceArea_CutOff_Idx,1) > 0

                            NucCountPerSliceArea_Smooth_CutOff_SliceNum = NucCountPerSliceArea_CutOff_Idx(1,1);
                        else
                            NucCountPerSliceArea_Smooth_CutOff_SliceNum = inf;
                        end

                        if size(SphereRegionAreaPerSlice_CutOff_Idx,1) > 0
                            SphereRegionAreaPerSlice_Final_CutOff_SliceNum = SphereRegionAreaPerSlice_CutOff_Idx(1,1);
                        else
                            SphereRegionAreaPerSlice_Final_CutOff_SliceNum = inf;
                        end

                        %%%%%% creating a cutoff to output a cutoff s9 variable
                        %%%%%% (not calculating metric here)
                        if NucCountPerSliceArea_Smooth_CutOff_SliceNum < SphereRegionAreaPerSlice_Final_CutOff_SliceNum
                        PostStats_Analysis_Cutoff_Zslice = NucCountPerSliceArea_Smooth_CutOff_SliceNum;
                        else
                        PostStats_Analysis_Cutoff_Zslice =  SphereRegionAreaPerSlice_Final_CutOff_SliceNum; 
                        end


                        if PostStats_Analysis_Cutoff_Zslice < inf


                            for V1 = 1:size(s9,1)
                                s9(V1,1).OriginalIdNum = V1; 
                            end

                           s9_OriginalReserve = s9;
                           PostStats_GoodZ_AnalysisRange_Idx = find(Centroid_Matrix(:,3) < PostStats_Analysis_Cutoff_Zslice);
                           s9 = s9(PostStats_GoodZ_AnalysisRange_Idx);
                      
                        end

                % defining sub-regions of spheroid by slice range (3x sub-regions are defined here)
                    Z_RangeLimit_Low_Maximum = 12;
                    % middle range is defined by the above line of code and below line of code
                    Z_RangeLimit_Top_Minimum = 24;
                    
                % defining positive cells by threshold
                    CellThresh_RedChan = 2500; 
                    
                % defining inner vs outer cells by threshold
                    CellThresh_SphereDistInner = 20; % number of voxels distance from spheroid edge

                    Centroid_Matrix = reshape([s9.Centroid], 3,size(s9,1))';

                    idx = find(uint16(Centroid_Matrix(:,3)) <= Z_RangeLimit_Low_Maximum);
                    Parameters_PerSliceRange_Low = s9(idx);

                    idx = find(uint16(Centroid_Matrix(:,3)) > Z_RangeLimit_Low_Maximum & uint16(Centroid_Matrix(:,3)) < Z_RangeLimit_Top_Minimum);
                    Parameters_PerSliceRange_Mid = s9(idx);

                    idx = find(uint16(Centroid_Matrix(:,3)) >= Z_RangeLimit_Top_Minimum);
                    Parameters_PerSliceRange_Top = s9(idx);


                    % Red Channel measurement - uncomment all red channel
                    % code to analyze second channel of immunostaining,
                    % etc.

%                     idx = find([s9.NucIntRedChan] > CellThresh_RedChan); 
%                     Parameters_RedChanPosCells = s9(idx);

%                     idx = find([Parameters_PerSliceRange_Low.NucIntRedChan] > CellThresh_RedChan); 
%                     Parameters_PerSliceRange_Low_RedChanPosCells = Parameters_PerSliceRange_Low(idx);
% 
%                     idx = find([Parameters_PerSliceRange_Mid.NucIntRedChan] > CellThresh_RedChan); 
%                     Parameters_PerSliceRange_Mid_RedChanPosCells = Parameters_PerSliceRange_Mid(idx);
% 
%                     idx = find([Parameters_PerSliceRange_Top.NucIntRedChan] > CellThresh_RedChan); 
%                     Parameters_PerSliceRange_Top_RedChanPosCells = Parameters_PerSliceRange_Top(idx);

                    idx = find([s9.DistCentroidToSphereEdge] > CellThresh_SphereDistInner); 
                    Parameters_SphereDistInner = s9(idx);

                    idx = find([s9.DistCentroidToSphereEdge] < CellThresh_SphereDistInner); 
                    Parameters_SphereDistOuter = s9(idx);

                    % Red channel measurement
                    
%                     idx = find([Parameters_SphereDistInner.NucIntRedChan] > CellThresh_RedChan); 
%                     Parameters_SphereDistInner_RedChanPosCells = Parameters_SphereDistInner(idx);
% 
%                     idx = find([Parameters_SphereDistOuter.NucIntRedChan] > CellThresh_RedChan); 
%                     Parameters_SphereDistOuter_RedChanPosCells = Parameters_SphereDistOuter(idx);

                    idx = find([Parameters_PerSliceRange_Low.DistCentroidToSphereEdge] > CellThresh_SphereDistInner);
                    Parameters_PerSliceRange_Low_SphereDistInner = Parameters_PerSliceRange_Low(idx);
                    
                    % Red Channel

%                     idx = find([Parameters_PerSliceRange_Low_SphereDistInner.NucIntRedChan] > CellThresh_RedChan);
%                     Parameters_PerSliceRange_Low_SphereDistInner_RedChanPosCells = Parameters_PerSliceRange_Low_SphereDistInner(idx); 

                    idx = find([Parameters_PerSliceRange_Low.DistCentroidToSphereEdge] < CellThresh_SphereDistInner);
                    Parameters_PerSliceRange_Low_SphereDistOuter = Parameters_PerSliceRange_Low(idx);
                    
                    % red channel
% 
%                     idx = find([Parameters_PerSliceRange_Low_SphereDistOuter.NucIntRedChan] > CellThresh_RedChan);
%                     Parameters_PerSliceRange_Low_SphereDistOuter_RedChanPosCells = Parameters_PerSliceRange_Low_SphereDistOuter(idx); 

                    idx = find([Parameters_PerSliceRange_Mid.DistCentroidToSphereEdge] > CellThresh_SphereDistInner);
                    Parameters_PerSliceRange_Mid_SphereDistInner = Parameters_PerSliceRange_Mid(idx);

                    % red channel
%                     idx = find([Parameters_PerSliceRange_Mid_SphereDistInner.NucIntRedChan] > CellThresh_RedChan);
%                     Parameters_PerSliceRange_Mid_SphereDistInner_RedChanPosCells = Parameters_PerSliceRange_Mid_SphereDistInner(idx); 

                    idx = find([Parameters_PerSliceRange_Mid.DistCentroidToSphereEdge] < CellThresh_SphereDistInner);
                    Parameters_PerSliceRange_Mid_SphereDistOuter = Parameters_PerSliceRange_Mid(idx);
                    
                    % red channel
%                     idx = find([Parameters_PerSliceRange_Mid_SphereDistOuter.NucIntRedChan] > CellThresh_RedChan);
%                     Parameters_PerSliceRange_Mid_SphereDistOuter_RedChanPosCells = Parameters_PerSliceRange_Mid_SphereDistOuter(idx); 

                    idx = find([Parameters_PerSliceRange_Top.DistCentroidToSphereEdge] > CellThresh_SphereDistInner);
                    Parameters_PerSliceRange_Top_SphereDistInner = Parameters_PerSliceRange_Top(idx);

                    
                    % red channel
%                     idx = find([Parameters_PerSliceRange_Top_SphereDistInner.NucIntRedChan] > CellThresh_RedChan);
%                     Parameters_PerSliceRange_Top_SphereDistInner_RedChanPosCells = Parameters_PerSliceRange_Top_SphereDistInner(idx); 

                    %%%

                    idx = find([Parameters_PerSliceRange_Top.DistCentroidToSphereEdge] < CellThresh_SphereDistInner);
                    Parameters_PerSliceRange_Top_SphereDistOuter = Parameters_PerSliceRange_Top(idx);
                    
                    % red channel
%                     idx = find([Parameters_PerSliceRange_Top_SphereDistOuter.NucIntRedChan] > CellThresh_RedChan);
%                     Parameters_PerSliceRange_Top_SphereDistOuter_RedChanPosCells = Parameters_PerSliceRange_Top_SphereDistOuter(idx); 

                    %%%
                    Summary.Counts_Total(C1,C2) = size(s9,1);
                    Summary.Counts_PerSliceRange_Low(C1,C2) = size(Parameters_PerSliceRange_Low,1);
                    Summary.Counts_PerSliceRange_Mid(C1,C2) = size(Parameters_PerSliceRange_Mid,1);
                    Summary.Counts_PerSliceRange_Top(C1,C2) = size(Parameters_PerSliceRange_Top,1);
                    Summary.Counts_SphereDistInner(C1,C2) = size(Parameters_SphereDistInner,1);
                    Summary.Counts_SphereDistOuter(C1,C2) = size(Parameters_SphereDistOuter,1);            
                    % Summary.Counts_RedChanPosCells(C1,C2) = size(Parameters_RedChanPosCells,1);

                    %%%
                    Summary.PercentTotal_PerSliceRange_Low(C1,C2) = 100*size(Parameters_PerSliceRange_Low,1)/size(s9,1);
                    Summary.PercentTotal_PerSliceRange_Mid(C1,C2) = 100*size(Parameters_PerSliceRange_Mid,1)/size(s9,1);
                    Summary.PercentTotal_PerSliceRange_Top(C1,C2) = 100*size(Parameters_PerSliceRange_Top,1)/size(s9,1);
                    Summary.PercentTotal_SphereDistInner(C1,C2) = 100*size(Parameters_SphereDistInner,1)/size(s9,1);
                    Summary.PercentTotal_SphereDistOuter(C1,C2) = 100*size(Parameters_SphereDistOuter,1)/size(s9,1);
                    % Summary.PercentTotal_RedChanPosCells(C1,C2) = 100*size(Parameters_RedChanPosCells,1)/size(s9,1);

                    %%%
                    Summary.PercentSliceRange_Low_SphereDistInner(C1, C2) = 100*size(Parameters_PerSliceRange_Low_SphereDistInner,1)/size(Parameters_PerSliceRange_Low,1);
                    Summary.PercentSliceRange_Mid_SphereDistInner(C1, C2) = 100*size(Parameters_PerSliceRange_Mid_SphereDistInner,1)/size(Parameters_PerSliceRange_Mid,1);
                    Summary.PercentSliceRange_Top_SphereDistInner(C1, C2) = 100*size(Parameters_PerSliceRange_Top_SphereDistInner,1)/size(Parameters_PerSliceRange_Top,1);

                    %%%
                    %Summary.PercentSliceRange_Low_RedChanPosCells(C1, C2) = 100*size(Parameters_PerSliceRange_Low_RedChanPosCells,1)/size(Parameters_PerSliceRange_Low,1);
                    %Summary.PercentSliceRange_Mid_RedChanPosCells(C1, C2) = 100*size(Parameters_PerSliceRange_Mid_RedChanPosCells,1)/size(Parameters_PerSliceRange_Mid,1);
                    %Summary.PercentSliceRange_Top_RedChanPosCells(C1, C2) = 100*size(Parameters_PerSliceRange_Top_RedChanPosCells,1)/size(Parameters_PerSliceRange_Top,1);

                    %%%
                    % Summary.PercentSphereDist_Inner_RedChanPosCells(C1, C2) = 100*size(Parameters_SphereDistInner_RedChanPosCells,1)/size(Parameters_SphereDistInner,1);
                    % Summary.PercentSphereDist_Outer_RedChanPosCells(C1, C2) = 100*size(Parameters_SphereDistOuter_RedChanPosCells,1)/size(Parameters_SphereDistOuter,1);

                    %%%
                    % Summary.PercentSliceRangeAndSphereDist_Low_Inner_Red(C1, C2) = 100*size(Parameters_PerSliceRange_Low_SphereDistInner_RedChanPosCells,1)/size(Parameters_PerSliceRange_Low_SphereDistInner,1);
                    % Summary.PercentSliceRangeAndSphereDist_Low_Outer_Red(C1, C2) = 100*size(Parameters_PerSliceRange_Low_SphereDistOuter_RedChanPosCells,1)/size(Parameters_PerSliceRange_Low_SphereDistOuter,1);

                    % Summary.PercentSliceRangeAndSphereDist_Mid_Inner_Red(C1, C2) = 100*size(Parameters_PerSliceRange_Mid_SphereDistInner_RedChanPosCells,1)/size(Parameters_PerSliceRange_Mid_SphereDistInner,1);
                    % Summary.PercentSliceRangeAndSphereDist_Mid_Outer_Red(C1, C2) = 100*size(Parameters_PerSliceRange_Mid_SphereDistOuter_RedChanPosCells,1)/size(Parameters_PerSliceRange_Mid_SphereDistOuter,1);

                    % Summary.PercentSliceRangeAndSphereDist_Top_Inner_Red(C1, C2) = 100*size(Parameters_PerSliceRange_Top_SphereDistInner_RedChanPosCells,1)/size(Parameters_PerSliceRange_Top_SphereDistInner,1);
                    % Summary.PercentSliceRangeAndSphereDist_Top_Outer_Red(C1, C2) = 100*size(Parameters_PerSliceRange_Top_SphereDistOuter_RedChanPosCells,1)/size(Parameters_PerSliceRange_Top_SphereDistOuter,1);

                    %%%
                    %%%
                    %%%
                    Summary.Well_Row(C1,C2) = Well_Row;
                    Summary.Well_Column(C1,C2) = Well_Column;

                    Summary.TopEndZ_SpheroidRegion(C1,C2) = SphereRegionAreaPerSlice_Final_CutOff_SliceNum;
                    Summary.TopEndZ_ReliableNucSegmentation(C1,C2) = NucCountPerSliceArea_Smooth_CutOff_SliceNum;

                    
                    Summary.ClearingQualityMetric(C1, C2) = double(NucCountPerSliceArea_Smooth_CutOff_SliceNum)/Sphere_ApproxHeight;
                    Summary.Sphere_ApproxHeight(C1, C2) = Sphere_ApproxHeight;

                    Summary.PostStats_Analysis_Cutoff_Zslice(C1, C2) = PostStats_Analysis_Cutoff_Zslice;
                    end % conditional, check region properties of whole spheroid   
                 end % conditional, make sure that WholeSphereROI_7 contains a regiion 
              end  %conditional, median filter at E9 produces acceptable values 
           end % conditional, make sure that a potential whole spheroid region is detected at variable E6  
         end % conditional loop to make sure new data file was loaded, ie: input data exists for selected micro well coordinates
         
    end % C2 experimental replicate loop
end % C1 experimental condition loop 

if exist('Summary') == 1
        Summary_ImportReserve = Summary;

           Summary_Field_Names = fieldnames(Summary);

           Summary_Matrix_Well_Row = [Summary.Well_Row];
           Summary_Matrix_Well_Column = [Summary.Well_Column];

           clear Summary_Analyzed_Wells_PerCondition
           for J1 = 1:size(Condition_List,1)
               % J1 = 1
               idx = find(Summary_Matrix_Well_Row(J1,:)>0);
               Summary_Analyzed_Wells_PerCondition(J1).Well_ReplicateCount = size(idx,2);
               Summary_Analyzed_Wells_PerCondition(J1).Well_InternalIdx_List = idx;
               Summary_Analyzed_Wells_PerCondition(J1).Well_Row_List = Summary_Matrix_Well_Row(J1,idx);
               Summary_Analyzed_Wells_PerCondition(J1).Well_Column_List = Summary_Matrix_Well_Column(J1,idx);
           end

           % does not exclude zeros that were generated by real wells
           for J2 = 1:size(Summary_Field_Names,1)
               % J2 = 1
               Field_Name_CharArray = char(Summary_Field_Names(J2));
               Temp_Summary_SingleField_Matrix = Summary.(Field_Name_CharArray);
               for J3 = 1:size(Condition_List,1)
                   %J3 = 1
                   Temp_Summary_SingleField_SingleCondition_Matrix = Temp_Summary_SingleField_Matrix(J3, Summary_Analyzed_Wells_PerCondition(J3).Well_InternalIdx_List);
                    idx = find(isnan(Temp_Summary_SingleField_SingleCondition_Matrix) ==1);

                    Temp_Summary_SingleField_SingleCondition_Matrix_2 = Temp_Summary_SingleField_SingleCondition_Matrix;

                    if size(idx,1)>0 
                    Temp_Summary_SingleField_SingleCondition_Matrix_2(idx) = [];            
                    end

                   Field_Name_Concat = char(strcat(string(Summary_Field_Names(J2)),'_RepCount'));
                   Summary_Stats.(Field_Name_Concat)(J3,1) = size(Temp_Summary_SingleField_SingleCondition_Matrix_2,2);

                   if size(Temp_Summary_SingleField_SingleCondition_Matrix_2,2) > 0
                       Field_Name_Concat = char(strcat(string(Summary_Field_Names(J2)),'_Mean'));
                       Summary_Stats.(Field_Name_Concat)(J3,1) = mean(Temp_Summary_SingleField_SingleCondition_Matrix_2,2);

                       Field_Name_Concat = char(strcat(string(Summary_Field_Names(J2)),'_StDev'));
                       Summary_Stats.(Field_Name_Concat)(J3,1) = std(single(Temp_Summary_SingleField_SingleCondition_Matrix_2),0,2);
                   else
                        Field_Name_Concat = char(strcat(string(Summary_Field_Names(J2)),'_Mean'));
                       Summary_Stats.(Field_Name_Concat)(J3,1) = NaN;

                       Field_Name_Concat = char(strcat(string(Summary_Field_Names(J2)),'_StDev'));
                       Summary_Stats.(Field_Name_Concat)(J3,1) = NaN;              
                   end

               end
           end

           % excludes all zeros
         ClearingQualityMetric_Matrix = [Summary.ClearingQualityMetric];
         ClearingQualityMetric_Matrix(ClearingQualityMetric_Matrix > 1) = 0;
         for C1 = 1:size(ClearingQualityMetric_Matrix)
         ClearingQualityMetric_Matrix_SingleCondition_nonZero_idx = find(ClearingQualityMetric_Matrix(C1,:) > 0);
         Summary_Stats_ClearingQualityMetric_Mean(C1,1) = mean(ClearingQualityMetric_Matrix(C1,ClearingQualityMetric_Matrix_SingleCondition_nonZero_idx));
         Summary_Stats_ClearingQualityMetric_StDev(C1,1) = std(ClearingQualityMetric_Matrix(C1,ClearingQualityMetric_Matrix_SingleCondition_nonZero_idx));
         Summary_Stats_ClearingQualityMetric_Count(C1,1) = size(ClearingQualityMetric_Matrix_SingleCondition_nonZero_idx,2);
         end
 
end % conditional, ensure 'Summary' variable exists
   

 Summary_Stats_FieldNames = fieldnames(Summary_Stats);
