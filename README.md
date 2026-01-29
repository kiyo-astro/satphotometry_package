# satphotometry_package
Python program and notebook set with satphotometry library including satellite orbit calculator, satellite identification and photometry tool

### Usage
Before use, edit SPICE toolkit myfile (./config/myfile.txt) and set your space-track.org user id and password at ./config/space-track-org.config<br>
Also, following kernel files for SPICE toolkit are required. You can download additional kernel file set from https://www.kiyoaki.jp/wp-content/uploads/kernel.zip and put the files into `./cspice/data/` folder.
* geophysical.ker
* naif0012.tls
* de440.bsp
* pck00011.tpc
* earth_assoc_itrf93.tf
* earth_720101_070426.bpc
* earth_latest_high_prec.bpc
* earth_200101_990628_predict.bpc
* myframe.tf

### Developer
This package is developed by Kiyoaki Okudaira and Kyushu University Hanada Lab (SSDL)<br>
(c) 2026 Kiyoaki Okudaira
