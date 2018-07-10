# A high-throughput imaging and nuclear segmentation analysis protocol for cleared 3D culture models

The following four Matlab script files are provided as open source software under the conditions set forth by 'The MIT License' (details in bottom section of this ReadMe file).

```
a170829_ParFor_FunctionCall_Fnctn_Boutin2018_3d_Sphere_Seg.m
Fnctn_Boutin2018_3d_Sphere_Seg_v02.m
a170817_Boutin2018_PostStats_v02_ch1andch2.m
a170817_Boutin2018_PostStats_v02_ch1only.m
```

## Brief description of file function

The `Fnctn_Boutin2018_3d_Sphere_Seg_v02.m` Matlab file is the primary image analysis function for identifying and measuring individual nuclei volumes within the larger multi-cellular spheroid volume. The primary image analysis function is run on parallel cpu processing cores by the `a170829_ParFor_FunctionCall_Fnctn_Boutin2018_3d_Sphere_Seg.m` Matlab file.

`a170817_Boutin2018_PostStats_v02_ch1andch2.m` and `a170817_Boutin2018_PostStats_v02_ch1only.m` Matlab files are for secondary statistical analysis and load data files that are created by the primary image analysis function.

`Well_Row_Col_List_short_test.mat` Matlab data file contains list of well coordinates and experimental condition groups for analysis.

The `fspecial3.m` Matlab file is a previously published function. Full credits are provided within comments of this file.

## Academic credits

The proof of concept for the use of these Matlab files in cellular analysis is published in a peer-reviewed journal article, which is cited below.

A high-throughput imaging and nuclear segmentation analysis protocol for cleared 3D culture models. Molly E. Boutin, Ty C. Voss, Steven A. Titus, Kennie Cruz-Gutierrez, Sam Michael, Marc Ferrer. Scientific Reports. 2018. Pre-press at this time: volume/issue/page numbers to be detemined. 

Any future projects using this source code or derivative code should cite the above journal article. Full citation information for this journal article will be available at https://www.ncbi.nlm.nih.gov/pubmed

## License Details

The MIT License
https://opensource.org/licenses/MIT

Copyright 2018 National Center for Advancing Translational Sciences (NCATS), National Institutes of Health (NIH)
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
