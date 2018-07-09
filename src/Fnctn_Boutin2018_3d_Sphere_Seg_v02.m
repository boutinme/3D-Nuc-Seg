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

function [ Output_1 ] = Fnctn_Boutin2018_3d_Sphere_Seg_v01( WellCounter, In_Path_Im_1, Well_Row_Col_List, PlanesNum, Out_Path_Data_1 )
%
%   This is main function for Boutin et al 2018, 3D segmentation of DAPI
%   stained nuclei within spheroids following optical clearing treatment.
%   Allows segmentation optimized for high-throughput imaging. Relatively
%   sparse Z-sampling is employed for increased sample throughput.
%
% This function segments 3D volumes for each nucleus based on shape of volumes following image filtering with Gaussian blur and Laplacian edge filters.
% Gaussian filter size should be set to approximate nucleus size for data set.
% Original parameter settings are for confocal mode on PerkinElmer Opera
% Phenix High-Throughput Imaging Instrument fitted with 20x water immersion objective lens,
% 2x digital camera binning. This imaging set-up provides XY pixel size of ~640 nm. 
% Z-slice interval during acquisition is 5 microns. 
%
% This function is called by ParFor loop Matlab script:
% 'a170829_ParFor_FunctionCall_Fnctn_Boutin2018_3d_Sphere_Seg.m'
% This combination allows parallel analysis of each micro-well within the
% plate.
%

% Ratio of physical size of Z slice to XY pixel size in image stack
ScaleFactor_ZtoXY = 8; 

% starting value for well column list
% for WellCounter = 1:size(Well_Row_Col_List,1)
    % WellCounter = 1 % can modify to run by hand
    Well_Row = Well_Row_Col_List(WellCounter,1);
    Well_Column = Well_Row_Col_List(WellCounter,2);

% location of input tiff image files (into this matlab script) exported
% from Acapella PerkinElmer Imaging Software... different tiff importing
% code will be needed for different systematic tif file naming schemes.

            cd(In_Path_Im_1)
            
            file_list_4 = ls;
            file_list_4 = cellstr(file_list_4);
              
            if Well_Row > 9 && Well_Column > 9
            search_key3 = strcat('r',num2str(Well_Row), 'c', num2str(Well_Column));
            end 
            
            if Well_Row < 10 && Well_Column > 9
            search_key3 = strcat('r',num2str(0),num2str(Well_Row), 'c', num2str(Well_Column));
            end
            
             if Well_Row > 9 && Well_Column < 10
            search_key3 = strcat('r',num2str(Well_Row), 'c',num2str(0), num2str(Well_Column));
             end
             
            if Well_Row < 10 && Well_Column < 10
            search_key3 = strcat('r',num2str(0),num2str(Well_Row), 'c',num2str(0), num2str(Well_Column));
             end
            
            index = strfind(file_list_4,search_key3);
            index = find(~cellfun(@isempty,index));
            file_list_5 = file_list_4(index,1);

                %%%
                % channel 1 file names - nuc marker channel
                search_key1 = '-ch1';
                index = strfind(file_list_5,search_key1);
                index = find(~cellfun(@isempty,index));
                file_list_6 = file_list_5(index,1);
                file_list_6_numbers = uint16(zeros(size(file_list_6,1),1));
                
                
                for C1 = 1:size(file_list_6,1)
                   Str_Flag_Start =findstr(cell2mat(file_list_6(C1,1)),'p');
                   Str_Flag_End =findstr(cell2mat(file_list_6(C1,1)),'-ch');
                   Str_FileName = cell2mat(file_list_6(C1,1));
                   file_list_6_numbers(C1,1) = str2num(Str_FileName(Str_Flag_Start+1:Str_Flag_End-1)); 
                end
                
                [file_list_6_numbers_sort file_list_6_numbers_sort_idx] = sort(file_list_6_numbers);
                
                file_list_6b = file_list_6(file_list_6_numbers_sort_idx,1);
  
 % below IF statement goes all the way to very end of code            
 if size(file_list_6b,1) > 0
                        clear In_Stack_2
                        C2 = 0;
                        for C1 = 1:PlanesNum
                            C2 = C2+1;
                           In_Stack_2(:,:,C2) = imread(cell2mat(file_list_6b(C1))); 
                        end
                        
                        % Use for channel 2 or 3 file names - immunostain or marker channel
                        search_key1 = '-ch1'; % change this accordingly if want to additionally analyze a second channel (is currently ch1, nuclear channel)
                        index = strfind(file_list_5,search_key1);
                        index = find(~cellfun(@isempty,index));
                        file_list_6 = file_list_5(index,1);
                        file_list_6_numbers = uint16(zeros(size(file_list_6,1),1)); 


                        for C1 = 1:size(file_list_6,1)

                           Str_Flag_Start =findstr(cell2mat(file_list_6(C1,1)),'p');
                           Str_Flag_End =findstr(cell2mat(file_list_6(C1,1)),'-ch');
                           Str_FileName = cell2mat(file_list_6(C1,1));
                           file_list_6_numbers(C1,1) = str2num(Str_FileName(Str_Flag_Start+1:Str_Flag_End-1)); 
                        end

                        [file_list_6_numbers_sort file_list_6_numbers_sort_idx] = sort(file_list_6_numbers);

                        file_list_6b = file_list_6(file_list_6_numbers_sort_idx,1);

                        clear In_Stack_3

                         for C1 = 1:PlanesNum
                            In_Stack_3(:,:,C1) = imread(cell2mat(file_list_6b(C1))); 
                        end

                        In_Stack_2_reserve = In_Stack_2;
                        In_Stack_3_reserve = In_Stack_3;
        %               In_Stack_4_reserve = In_Stack_4;

        % Below section: 
        % Define ROI for entire spheroid by simple thresholding 
        % and binary morphological filtering
        % more sophisticated methods could be introduced later
        % This preliminary region allows cropping of image stack, increases
        % computation speed of subsequent steps
        
        
        % reduce overal stack size to increase computation speed 
        In_Stack_2_resizeXY = imresize(In_Stack_2, 1/ScaleFactor_ZtoXY);

        In_Stack_2_resizeXY_BG = mode(In_Stack_2_resizeXY(:));
        WholeSphereROI_FoldBG_ThresholdFactor = 15;
        
        WholeSphereROI_FoldBG_Threshold = WholeSphereROI_FoldBG_ThresholdFactor * In_Stack_2_resizeXY_BG;

        WholeSphereROI_1 = false(size(In_Stack_2_resizeXY));
        WholeSphereROI_1_ImageEdgeFrame = WholeSphereROI_1;
        WholeSphereROI_1(In_Stack_2_resizeXY > WholeSphereROI_FoldBG_Threshold) = 1;

        WholeSphereROI_1_ImageEdgeFrame(:,[1,end],:) = 1; 
        WholeSphereROI_1_ImageEdgeFrame([1,end],:,:) = 1;

        % structuring element size for 2D image filtering
         WholeSphereROI_1_2D_strel_radius = 11;
        % structuring element size for 3D image filtering
        WholeSphereROI_1_3D_strel_radius = 3;

        % create structuring elements that will be used in image morph filtering
         WholeSphereROI_1_2D_strel = strel('disk', WholeSphereROI_1_2D_strel_radius);
          WholeSphereROI_1_3D_strel = strel('sphere', WholeSphereROI_1_3D_strel_radius);

        WholeSphereROI_1_b = padarray(WholeSphereROI_1, [WholeSphereROI_1_2D_strel_radius+1 WholeSphereROI_1_2D_strel_radius+1 0]);

        WholeSphereROI_2 = imclose(WholeSphereROI_1_b, WholeSphereROI_1_2D_strel);
        WholeSphereROI_3 = imclose(WholeSphereROI_2,WholeSphereROI_1_3D_strel);

        % remove padding from stack after filtering
        WholeSphereROI_3_b = WholeSphereROI_3((WholeSphereROI_1_2D_strel_radius+2):end-(WholeSphereROI_1_2D_strel_radius+1),(WholeSphereROI_1_2D_strel_radius+2):end-(WholeSphereROI_1_2D_strel_radius+1),:);

        % find multiple large volumes that do not overlab with image edge,
        % label them -- potentially fragmented spheroids or debris
        WholeSphereROI_3_c = bwlabeln(WholeSphereROI_3_b);
        for J2 = 1:max(WholeSphereROI_3_c(:))
            WholeSphereROI_3_d = WholeSphereROI_1_ImageEdgeFrame(WholeSphereROI_3_c==J2);

            if max(WholeSphereROI_3_d(:)) == 1
                WholeSphereROI_3_b(WholeSphereROI_3_c==J2)=0;
            end
        end
        
        % restore original stack size
        WholeSphereROI_4 = imresize(WholeSphereROI_3_b,ScaleFactor_ZtoXY);

        

        WholeSphereROI_5 = bwlabeln(WholeSphereROI_4);

        PreSegment_LargeObjects = regionprops(WholeSphereROI_5,'BoundingBox','Area');

        PreSegment_LargeObjects_SelectIdx = find([PreSegment_LargeObjects.Area] > 20000);
        PreSegment_LargeObjects = PreSegment_LargeObjects(PreSegment_LargeObjects_SelectIdx);

        WholeSphereROI_5b = ismember(WholeSphereROI_5,PreSegment_LargeObjects_SelectIdx);
        WholeSphereROI_5c = bwlabeln(WholeSphereROI_5b);

        tic
        % loop for large objects, spheroid scale.
        if size(PreSegment_LargeObjects,1) > 0
                for J1 = 1:size(size(PreSegment_LargeObjects,1))
                    % J1 = 1

                    clear In_Stack_2 In_Stack_3 In_Stack_4

                In_Stack_2 = In_Stack_2_reserve;
                In_Stack_3 = In_Stack_3_reserve;
        %       In_Stack_4 = In_Stack_4_reserve;

                % Image stacks cropped and used as input for further nuc level analysis
                In_Stack_2 = In_Stack_2(uint16(round(PreSegment_LargeObjects(J1).BoundingBox(2))):uint16(floor(PreSegment_LargeObjects(J1).BoundingBox(2))+floor(PreSegment_LargeObjects(J1).BoundingBox(5))),uint16(round(PreSegment_LargeObjects(J1).BoundingBox(1))):uint16(floor(PreSegment_LargeObjects(J1).BoundingBox(1))+floor(PreSegment_LargeObjects(J1).BoundingBox(4))),:);
                In_Stack_3 = In_Stack_3(uint16(round(PreSegment_LargeObjects(J1).BoundingBox(2))):uint16(floor(PreSegment_LargeObjects(J1).BoundingBox(2))+floor(PreSegment_LargeObjects(J1).BoundingBox(5))),uint16(round(PreSegment_LargeObjects(J1).BoundingBox(1))):uint16(floor(PreSegment_LargeObjects(J1).BoundingBox(1))+floor(PreSegment_LargeObjects(J1).BoundingBox(4))),:);
        %        In_Stack_4 = In_Stack_4(uint16(round(PreSegment_LargeObjects(J1).BoundingBox(2))):uint16(floor(PreSegment_LargeObjects(J1).BoundingBox(2))+floor(PreSegment_LargeObjects(J1).BoundingBox(5))),uint16(round(PreSegment_LargeObjects(J1).BoundingBox(1))):uint16(floor(PreSegment_LargeObjects(J1).BoundingBox(1))+floor(PreSegment_LargeObjects(J1).BoundingBox(4))),:);

        WholeSphereROI_6 = false(size(WholeSphereROI_5c));
        WholeSphereROI_6(WholeSphereROI_5c==J1) = 1;

        % saved variable: simple stack image label for whole spheroid
        % region
        WholeSphereROI_7 = WholeSphereROI_6(uint16(round(PreSegment_LargeObjects(J1).BoundingBox(2))):uint16(floor(PreSegment_LargeObjects(J1).BoundingBox(2))+floor(PreSegment_LargeObjects(J1).BoundingBox(5))),uint16(round(PreSegment_LargeObjects(J1).BoundingBox(1))):uint16(floor(PreSegment_LargeObjects(J1).BoundingBox(1))+floor(PreSegment_LargeObjects(J1).BoundingBox(4))),:);
        
        % code (below) to segment individual nuclei from each large object
                % 2D filtering to reduce fine detail, interior volume of
                % each nucleus should be come uniform intensity
                % applied only in 2D becuase 3D is sparsely sampled in stack
                Nuc3D_Seg_3 = imfilter(In_Stack_2, fspecial('gaussian', 21, 2.5));
                % figure, imshow(Nuc3D_Seg_3(:,:,10),[])

                % 3D filtering to detect edge regions surrounding each nucleus 
                Nuc3D_Seg_4 = imfilter(Nuc3D_Seg_3, fspecial3('laplacian'));
                % figure, imshow(Nuc3D_Seg_4(:,:,10),[])

                % expanding edge region to fuse changing intensity to solid
                % volume that surrounds each nucleus, this function
                % operates on each XY plane
                Nuc3D_Seg_4a = imdilate(Nuc3D_Seg_4, strel('disk', 1));
                % figure, imshow(Nuc3D_Seg_4a(:,:,10),[])

                % setting changing surrounding each nuc to infinity
                Nuc3D_Seg_4b =Nuc3D_Seg_4a;
                Nuc3D_Seg_4b(Nuc3D_Seg_4>0) = inf; 

                % 2D filtering to detect edge regions surrounding each nucleus, more detail in 2D compared to 3D 
                Nuc3D_Seg_5 = imfilter(Nuc3D_Seg_3, fspecial('laplacian'));
                % figure, imshow(Nuc3D_Seg_5(:,:,10),[])

                % setting changing surrounding each nuc to infinity
                Nuc3D_Seg_5b = Nuc3D_Seg_5;
                Nuc3D_Seg_5b(Nuc3D_Seg_5>0)= inf;
                % figure, imshow(Nuc3D_Seg_5b(:,:,10),[])

                % using 2D(XY) information, remove small isolated regions
                % that were previously detected
                Stack_2dSegmentedRegions = In_Stack_2;
                for C1 = 1:size(In_Stack_2,3)
                t1 = bwlabel(Nuc3D_Seg_5b(:,:,C1));
                t1_stats = regionprops(t1, 'area');
                idx = find([t1_stats.Area] > 1000);
                Stack_2dSegmentedRegions(:,:,C1) = ismember(t1,idx);
                end

                clear Nuc3D_Seg_3 Nuc3D_Seg_4 Nuc3D_Seg_4a Nuc3D_Seg_5 Nuc3D_Seg_5b  t1 t1_stats idx

                Stack_2dSegmentedRegions(Stack_2dSegmentedRegions > 0) = inf;

                % combine edges detected by 2D and 3D filtering operations
                Stack_2dSegmentedRegions_wVerticalEdges = Stack_2dSegmentedRegions + Nuc3D_Seg_4b;
                %figure, imshow(Stack_2dSegmentedRegions_wVerticalEdges(:,:,10),[])

                clear Stack_2dSegmentedRegions

                % calculate threshold value
                k1 = graythresh(In_Stack_2)*65536;

                % erode with 2D(XY) filter
                Stack_2dSegmentedRegions_wVerticalEdges = imerode(Stack_2dSegmentedRegions_wVerticalEdges, strel('disk', 1));

                % apply threshold value
                Stack_2dSegmentedRegions_wVerticalEdges(In_Stack_2<k1)=inf;
                Stack_2dSegmentedRegions_wVerticalEdges(:,:,1)=inf;
                Stack_2dSegmentedRegions_wVerticalEdges(:,:,end)=inf;

                clear k1
                % figure, imshow(Stack_2dSegmentedRegions_wVerticalEdges(:,:,10),[])
                % figure, imshow(In_Stack_2(:,:,10),[])

                Nuc3D_Seg_6 = imcomplement(Stack_2dSegmentedRegions_wVerticalEdges);

                %%% enlarge and interpolate 3d matrix
                M = double(Stack_2dSegmentedRegions_wVerticalEdges);
                k = ScaleFactor_ZtoXY;
                % tic
                   d=size(M);
                z=1:d(3);
                zi=1:1/k:d(3);
                M=permute(interp1(z,permute(M,[3 1 2]),zi),[2 3 1]);
                % toc
                %%%
                Nuc3D_Seg_7 = uint16(M);  
                clear M clear Stack_2dSegmentedRegions_wVerticalEdges

                % convert to black and white image stack
                Nuc3D_Seg_8 = logical(Nuc3D_Seg_7);
                Nuc3D_Seg_8(:,:,1) = 1; % set first plane to null background.
               %figure, imshow(Nuc3D_Seg_8(:,:,20),[])

                clear Nuc3D_Seg_7

                %%%% 
                
                DistMap_3D_1 = bwdist(Nuc3D_Seg_8);
                % figure, imshow(DistMap_3D_1(:,:,10),[])
                clear Nuc3D_Seg_8

                DistMap_3D_2 = imcomplement(DistMap_3D_1);
                clear DistMap_3D_1
           
                % supress small local minimums in distance map
                % this can be adjusted to correct for over-splitting and under-splitting of true nuclei 
                
                % this output distance map (DistMap_3D_3) should define 
                %centers of individual nuclei (large local minima) and boundaries between touching nuclei
                DistMap_3D_3 = imhmin(DistMap_3D_2, 1.5);
                clear DistMap_3D_2

                %%% enlarge and interpolate 3d matrix
                M = double(Nuc3D_Seg_6);
                k = ScaleFactor_ZtoXY;
                %tic
                   d=size(M);
                z=1:d(3);
                zi=1:1/k:d(3);
                M=permute(interp1(z,permute(M,[3 1 2]),zi),[2 3 1]);
                %toc
                %%%
                Nuc3D_Seg_6_expanded = uint16(M);  
                clear M Nuc3D_Seg_6

                % set background regions (non-potential nuclei regions) in distance map to null (-inf) 
                DistMap_3D_3(Nuc3D_Seg_6_expanded == 0) = -inf;
                clear Nuc3D_Seg_6_expanded
                
                % blanking (-inf) top and bottom planes in distance map
                DistMap_3D_3(:,:,1) = -inf;
                DistMap_3D_3(:,:,end) = -inf;

                % watershed segmentation based on corrected distance map
                % Nuc3D_Seg_9 is 3D label matrix, partially separated nuclei
                Nuc3D_Seg_9 = watershed(DistMap_3D_3,26);
                clear DistMap_3D_3
                

                %%%%
                %%%%

                %%% reduce 3d matrix
                M = double(Nuc3D_Seg_9);
                clear test16
                k = double(1/ScaleFactor_ZtoXY);
                % tic
                   d=size(M);
                z=1:d(3);
                zi=1:1/k:d(3);
                M=permute(interp1(z,permute(M,[3 1 2]),zi),[2 3 1]);
                % toc
                %%%
                Nuc3D_Seg_10 = uint16(M);  
                clear M

                % statistics for individual 3D regions, 3D regions are
                % partially separated nuclei
                Nuc3D_Seg_10b = regionprops(Nuc3D_Seg_10, 'Area', 'Centroid', 'BoundingBox');

                % Below loop - based on stats evaluate partially sepatarted nuclei regions and
                % attempt further segmentation if needed.
                if size(Nuc3D_Seg_10b,1) > 1
                for C1 = 2:size(Nuc3D_Seg_10b,1) % this counter starts at '2', omits large non-spheroid region
                    % C1 = 1
                    % define whole number XYZ coordinates for bounding box of each 3D region
                    BoundingBox_Mat = uint16(round(Nuc3D_Seg_10b(C1).BoundingBox));

                    % Y coordinates for bounding box
                    BB_Y1 = BoundingBox_Mat(2);
                    BB_Y2 = BoundingBox_Mat(2)+BoundingBox_Mat(5);
                    if BB_Y2 > size(Nuc3D_Seg_10,1)
                        BB_Y2 = size(Nuc3D_Seg_10,1);
                    end

                    % X coordinates for bounding box
                    BB_X1 = BoundingBox_Mat(1);
                    BB_X2 = BoundingBox_Mat(1)+BoundingBox_Mat(4);
                    if BB_X2 > size(Nuc3D_Seg_10,2)
                        BB_X2 = size(Nuc3D_Seg_10,2);
                    end

                    % Z coordinates for bounding box
                    BB_Z1 = BoundingBox_Mat(3);
                    BB_Z2 = BoundingBox_Mat(3)+BoundingBox_Mat(6);
                    if BB_Z2 > size(Nuc3D_Seg_10,3)
                        BB_Z2 = size(Nuc3D_Seg_10,3);
                    end

                    % binary bounding box subvolume containing a single 3D region
                    BoundingImStack_1 = Nuc3D_Seg_10(BB_Y1:BB_Y2,BB_X1:BB_X2,BB_Z1:BB_Z2);
                    BoundingImStack_1(Nuc3D_Seg_10(BB_Y1:BB_Y2,BB_X1:BB_X2,BB_Z1:BB_Z2) ~= C1) = 0;
                    BoundingImStack_1(Nuc3D_Seg_10(BB_Y1:BB_Y2,BB_X1:BB_X2,BB_Z1:BB_Z2) == C1) = 1;

                    % for testing purposes, create 3D morphology variables
                    
                    % which variables describe a subvolume with 
                    % single nuclei vs non-nuclei vs multiple nuclei??
                    BoundingImStack_1_SinglePlane_Stats_Sum_Solidity = 0;
                    BoundingImStack_1_SinglePlane_Stats_Sum_Eccentricity = 0;
                    BoundingImStack_1_SinglePlane_Stats_Sum_Objects = 0;


                    % additional 3D variables are calculated by summing
                    % stack info for subvolume region
                    for C2 = 1:size(BoundingImStack_1,3) % loop through planes for each object
                        %C2 = 1
                        BoundingImStack_1_SinglePlane = BoundingImStack_1(:,:,C2);
                        BoundingImStack_1_SinglePlane_Stats = regionprops(BoundingImStack_1_SinglePlane, 'Area', 'ConvexArea','Eccentricity', 'Extent','Solidity');

                        if size(BoundingImStack_1_SinglePlane_Stats,1) > 0
                            BoundingImStack_1_SinglePlane_Stats_AreaFraction = BoundingImStack_1_SinglePlane_Stats(1).Area/Nuc3D_Seg_10b(C1).Area;

                            BoundingImStack_1_SinglePlane_Stats_Sum_Solidity = BoundingImStack_1_SinglePlane_Stats_Sum_Solidity + (BoundingImStack_1_SinglePlane_Stats(1).Solidity * BoundingImStack_1_SinglePlane_Stats_AreaFraction);
                            BoundingImStack_1_SinglePlane_Stats_Sum_Eccentricity = BoundingImStack_1_SinglePlane_Stats_Sum_Eccentricity+ (BoundingImStack_1_SinglePlane_Stats(1).Eccentricity * BoundingImStack_1_SinglePlane_Stats_AreaFraction);

                            BoundingImStack_1_SinglePlane_Stats_Sum_Objects = BoundingImStack_1_SinglePlane_Stats_Sum_Objects + max(max(bwlabel(BoundingImStack_1_SinglePlane)));

                        else            
                            BoundingImStack_1_SinglePlane_Stats_Sum_Objects = BoundingImStack_1_SinglePlane_Stats_Sum_Objects + 1;

                        end


                    end
                    Nuc3D_Seg_10b(C1).Z_Count = C2;
                    Nuc3D_Seg_10b(C1).WeightedObjectsCount = BoundingImStack_1_SinglePlane_Stats_Sum_Objects/C2;

                    Nuc3D_Seg_10b(C1).WeightedSolidity = BoundingImStack_1_SinglePlane_Stats_Sum_Solidity;
                    Nuc3D_Seg_10b(C1).WeightedEccentricity = BoundingImStack_1_SinglePlane_Stats_Sum_Eccentricity;

                    Nuc3D_Seg_10b(C1).OriginalLabelNum = C1;
                end
                % toc
                %%%%
                %%%%
                %%%%
                %%%%

                Nuc3D_Seg_10c = Nuc3D_Seg_10b;
                Nuc3D_Seg_10c(1) = []; % clear large surounding region... no spheroid in this region


                Nuc3D_Seg_10d = Nuc3D_Seg_10c(1:size([Nuc3D_Seg_10b.Z_Count]',1))';
                % figure, histogram([test51.Area], [1:100:2000]), title('all ROIs area')
                % figure, histogram([test51.WeightedSolidity]) , title('all ROIs WeightedSolidity')
                % figure, histogram([test51.Z_Count]) , title('all ROIs z count')
                % figure, histogram([test51.WeightedObjectsCount]) , title('all ROIs WeightedObjectsCount')
                
                % values from emperical testing used to select subvolumes
                % for further analysis
                
                Nuc3D_Seg_10e = find([Nuc3D_Seg_10d.Area] > 200 & [Nuc3D_Seg_10d.Area] < 8000 & [Nuc3D_Seg_10d.Z_Count] > 3 & [Nuc3D_Seg_10d.WeightedSolidity] < .84);
                Nuc3D_Seg_10f = Nuc3D_Seg_10d(Nuc3D_Seg_10e);

                    % create temporary list of object id numbers, those
                    % objects in label matrix Nuc3D_Seg_10... after
                    % selection in above code
                 temp_L = [Nuc3D_Seg_10d(Nuc3D_Seg_10e).OriginalLabelNum]';

                %tic
                % holder for removal of clumped ROIs that will be 're-split', 
                % the re-split ROIs will be added back here, then all labels will be re-numbered 
                SubVol_Seg_1 = Nuc3D_Seg_10;
                SubVol_Seg_2 = Nuc3D_Seg_10;
                SubVol_Seg_2(:) = 0;

                ROI_Counter_SplitROIs_1stPass = max(Nuc3D_Seg_10(:));

                % loop roles through all potential nuc ROIs, attempts to
                % separate touching nuclei
               if size(temp_L,1) > 2
                % loop through list of clumped ROIs
                for C2 = 1:size(temp_L,1)-1
                % pull bounding box from single roi, #1: undersegmented region -- C2 = 11

                % create blank volume that will soon contain a single ROI for re-splitting
                SubVol_Seg_3 = Nuc3D_Seg_10;
                SubVol_Seg_3(:) = 0;


                SubVol_Seg_3(Nuc3D_Seg_10==temp_L(C2)) = 1; % single ROI for re-splitting
                SubVol_Seg_1(Nuc3D_Seg_10==temp_L(C2)) = 0; % clearing the ROI in the original list

                 if Nuc3D_Seg_10d(temp_L(C2)).BoundingBox(6) > 1
                    BoundingBox_Mat = uint16(round(Nuc3D_Seg_10d(temp_L(C2)).BoundingBox));

                        BB_Y1 = BoundingBox_Mat(2);
                        BB_Y2 = BoundingBox_Mat(2)+BoundingBox_Mat(5);
                        if BB_Y2 > size(Nuc3D_Seg_10,1)
                            BB_Y2 = size(Nuc3D_Seg_10,1);
                        end

                        BB_X1 = BoundingBox_Mat(1);
                        BB_X2 = BoundingBox_Mat(1)+BoundingBox_Mat(4);
                        if BB_X2 > size(Nuc3D_Seg_10,2)
                            BB_X2 = size(Nuc3D_Seg_10,2);
                        end

                        BB_Z1 = BoundingBox_Mat(3);
                        BB_Z2 = BoundingBox_Mat(3)+BoundingBox_Mat(6);
                        if BB_Z2 > size(Nuc3D_Seg_10,3)
                            BB_Z2 = size(Nuc3D_Seg_10,3);
                        end



                    SubVol_Seg_4 = SubVol_Seg_3(BB_Y1:BB_Y2,BB_X1:BB_X2,BB_Z1:BB_Z2);

                    SubVol_Seg_5 = imcomplement(logical(SubVol_Seg_4));


                    %%% enlarge and interpolate 3d matrix
                    M = double(SubVol_Seg_5);
                    k = ScaleFactor_ZtoXY;
                    %tic
                       d=size(M);
                    z=1:d(3);
                    zi=1:1/k:d(3);
                    M=permute(interp1(z,permute(M,[3 1 2]),zi),[2 3 1]);
                    % toc
                    %%%
                    SubVol_Seg_6 = uint16(M);  
                    clear M

                    % distance map for subvolume regions
                    SubVol_Seg_7 = bwdist(logical(SubVol_Seg_6));
                    SubVol_Seg_7(SubVol_Seg_7 == 0) = inf;

                    SubVol_Seg_8 = imcomplement(SubVol_Seg_7);

 
                    SubVol_Seg_8b = SubVol_Seg_8;
                    SubVol_Seg_8b(SubVol_Seg_8 == inf) = inf;

                    
                    % create custom filter, can be adjusted to different
                    % levels of Z influence (Z influence can be eliminated)
                    Strel_ImErode1 = strel('sphere',2);
                    
                      xxx = [Strel_ImErode1(1).Neighborhood];
                      xxx(:,:,1:2) = 0;
                      xxx(:,:,4:5) = 0;
                    Strel_ImErode1 = strel('arbitrary',xxx);

                    % filter image stack with custom filter
                   SubVol_Seg_9 = imerode(SubVol_Seg_8b, Strel_ImErode1);
                    % figure, imshow(SubVol_Seg_9(:,:,16),[])

                    SubVol_Seg_10 = SubVol_Seg_9;
                    SubVol_Seg_10(SubVol_Seg_7 == inf) = -inf;

                    % blank (-inf) top and bottom z planes
                    SubVol_Seg_10(:,:,1) = -inf;
                    SubVol_Seg_10(:,:,end) = -inf;

                    %tic 
                    SubVol_Seg_11 = watershed(SubVol_Seg_10);
                    % toc
                    
                    SubVol_Seg_11b = SubVol_Seg_11;
                    SubVol_Seg_11b(SubVol_Seg_11==1)=0;

                    SubVol_Seg_11d = imdilate(SubVol_Seg_11b, Strel_ImErode1);

                    %%% Reduce and interpolate 3d matrix
                    M = double(SubVol_Seg_11d);
                    k = double(1/ScaleFactor_ZtoXY);%
                    %tic
                       d=size(M);
                    z=1:d(3);
                    zi=1:1/k:d(3);
                    M=permute(interp1(z,permute(M,[3 1 2]),zi),[2 3 1]);
                    %toc
                    %%%
                    SubVol_Seg_11_Zshort = uint16(M);  
                    clear M

                    SubVol_Seg_11_Zshort_Stats = regionprops(SubVol_Seg_11_Zshort, 'Area', 'Centroid', 'BoundingBox', 'PixelList');

                    if max(SubVol_Seg_11_Zshort(:)) > 1
                        for C3 = 2:max(SubVol_Seg_11_Zshort(:)) % test38
                            % C3 = 2
                            ROI_Counter_SplitROIs_1stPass = ROI_Counter_SplitROIs_1stPass+1;
                            
                            % create pixel list for single roi from subvolume, linear format
                            LinearROI = sub2ind(size(SubVol_Seg_2),(uint16(floor([SubVol_Seg_11_Zshort_Stats(C3).PixelList(:,2)])) + (BB_Y1 -1)),(uint16(floor([SubVol_Seg_11_Zshort_Stats(C3).PixelList(:,1)])) + (BB_X1-1)),uint16(floor([SubVol_Seg_11_Zshort_Stats(C3).PixelList(:,3)])) + (BB_Z1-1));
                            
                           % new single nuc ROIs from subvolumes are
                           % extracted and inserted into full size image stack, labeled with new sequential numbers
                           SubVol_Seg_1(LinearROI) = ROI_Counter_SplitROIs_1stPass;
                           SubVol_Seg_2(LinearROI) = ROI_Counter_SplitROIs_1stPass; 
                        end
                    end

                end % conditional loop, checking for 3D bounding box

                    end % C2 loop, count through selected clumped ROIs - re-splitting loop
               end % size of Temp_L conditional loop
                   
                SubVol_Seg_12 = regionprops(SubVol_Seg_1, 'Area', 'Centroid', 'BoundingBox');

        
                % label list of larger nuclei volumes
                SubVol_Seg_13 = uint16(find([SubVol_Seg_12.Area]' > 200));
                
                SubVol_Seg_14 = SubVol_Seg_12(SubVol_Seg_13);


                SubVol_Seg_15 = SubVol_Seg_1;
                SubVol_Seg_15(:) = 0;
                for C1 = 1:size(SubVol_Seg_13,1)
                    %for C1 = 1:100
                 SubVol_Seg_15(SubVol_Seg_1==SubVol_Seg_13(C1)) = C1;
                end

                %%%%%%%%%%
                %%%%%%%%%%
        
                SphereDistance_6b = imcomplement(WholeSphereROI_7);



                % enlarge and interpolate 3d matrix
                M = double(SphereDistance_6b);
                k = ScaleFactor_ZtoXY;
                %tic
                   d=size(M);
                z=1:d(3);
                zi=1:1/k:d(3);
                M=permute(interp1(z,permute(M,[3 1 2]),zi),[2 3 1]);
                %toc
                %%%
                SphereDistance_6b = logical(M);  
                clear M


                % blank first image in stack... 
                % important for accurate 3D distance calc at bottom of spheroid
                SphereDistance_6b(:,:,1) = 1;

                % create 3D distance map from exterior to interior for whole spheroid
                SphereDistance_7 = bwdist(SphereDistance_6b);

                %%% reduce and interpolate 3d matrix
                M = double(SphereDistance_7);
                k =double(1/ScaleFactor_ZtoXY);
                %tic
                   d=size(M);
                z=1:d(3);
                zi=1:1/k:d(3);
                M=permute(interp1(z,permute(M,[3 1 2]),zi),[2 3 1]);
                %toc
                %%%
                SphereDistance_8 = uint16(M);  
                clear M


                % s9 is datastructure for final output, fluor intensity values and
                % spheroid edge distances for each nuclear volume
                s9 = regionprops(SubVol_Seg_15,'Area','Centroid','PixelList');

                for C1 = 1:size(s9,1)
                    % C1 = 2 
                   s9(C1).DistCentroidToSphereEdge = SphereDistance_8(uint16([s9(C1).Centroid(2)]), uint16([s9(C1).Centroid(1)]), uint16([s9(C1).Centroid(3)]));
                   s9(C1).NucIntBlueChan = mean(In_Stack_2(SubVol_Seg_15 == C1));

                   % s9(C1).NucIntRedChan = mean(In_Stack_3(SubVol_Seg_15 == C1)); %need to uncomment to run a second channel                 


                end
                
                toc
                    % ,'In_Stack_4'
                 cd(Out_Path_Data_1)
                 save(strcat('Test_a170830_AutoSeg_Parallel_v01_r',num2str(Well_Row),'_c',num2str(Well_Column),'_ObjectID',num2str(J1),'.mat'),'SubVol_Seg_15','In_Stack_2', 'In_Stack_3','PreSegment_LargeObjects_SelectIdx', 's9', 'WholeSphereROI_7','In_Stack_2_resizeXY_BG');
                    SavedOutputFileName = strcat('a170830_AutoSeg_Parallel_v01_r',num2str(Well_Row),'_c',num2str(Well_Column),'_ObjectID',num2str(J1),'.mat')
                end % conditional loop test50b must be greater than 1 
                 end % J1 counter --- large object, presegment loop
        end % conditional, ensure large objects are found
        
        toc

    end % conditional statement to ensure a file is found
    Output_1 = strcat('WellCounter = ',num2str(WellCounter))
    
 end



