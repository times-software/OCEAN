#!/usr/bin/perl
use strict;
require JSON::PP;
use JSON::PP;


# X-RAY DATA BOOKLET
# Center for X-ray Optics and Advanced Light Source
my %allEdges = ('H-K',13.6,'He-K',24.6,'Li-K',54.7,'Li-L1',5.3,'Be-K',111.5,'Be-L1',8,
    'Be-L2',3,'Be-L3',3,'B-K',188,'B-L1',12.6,'B-L2',4.7,'B-L3',4.7,
    'C-K',284.2,'C-L1',18,'C-L2',7.2,'C-L3',7.2,'N-K',409.9,'N-L1',37.3,
    'N-L2',17.5,'N-L3',17.5,'O-K',543.1,'O-L1',41.6,'O-L2',18.2,'O-L3',18.2,
    'F-K',696.7,'F-L1',45,'F-L2',19.9,'F-L3',19.9,'Ne-K',870.2,'Ne-L1',48.5,
    'Ne-L2',21.7,'Ne-L3',21.6,'Na-K',1070.8,'Na-L1',63.5,'Na-L2',30.4,'Na-L3',30.5,
    'Mg-K',1303,'Mg-L1',88.6,'Mg-L2',49.6,'Mg-L3',49.21,'Mg-M1',2,'Mg-M2',1,
    'Mg-M3',1,'Al-K',1559,'Al-L1',117.8,'Al-L2',72.9,'Al-L3',72.5,'Al-M1',4,
    'Al-M2',2,'Al-M3',2,'Si-K',1839,'Si-L1',149.7,'Si-L2',99.8,'Si-L3',99.2,
    'Si-M1',8,'Si-M2',2,'Si-M3',2,'P-K',2145.5,'P-L1',189,'P-L2',136,
    'P-L3',135,'P-M1',12,'P-M2',7,'P-M3',6,'S-K',2472,'S-L1',230.9,
    'S-L2',163.6,'S-L3',162.5,'S-M1',14,'S-M2',8,'S-M3',7,'Cl-K',2822,
    'Cl-L1',270,'Cl-L2',202,'Cl-L3',200,'Cl-M1',18,'Cl-M2',10,'Cl-M3',10,
    'Ar-K',3205.9,'Ar-L1',326.3,'Ar-L2',250.6,'Ar-L3',248.4,'Ar-M1',29.3,'Ar-M2',15.9,
    'Ar-M3',15.7,'K-K',3608.4,'K-L1',378.6,'K-L2',297.3,'K-L3',294.6,'K-M1',34.8,
    'K-M2',18.3,'K-M3',18.3,'Ca-K',4038.5,'Ca-L1',438.4,'Ca-L2',349.7,'Ca-L3',346.2,
    'Ca-M1',44.3,'Ca-M2',25.4,'Ca-M3',25.4,'Sc-K',4492,'Sc-L1',498,'Sc-L2',403.6,
    'Sc-L3',398.7,'Sc-M1',51.1,'Sc-M2',28.3,'Sc-M3',28.3,'Ti-K',4966,'Ti-L1',560.9,
    'Ti-L2',460.2,'Ti-L3',453.8,'Ti-M1',58.7,'Ti-M2',32.6,'Ti-M3',32.6,'Ti-M4',2,
    'Ti-M5',2,'V-K',5465,'V-L1',626.7,'V-L2',519.8,'V-L3',512.1,'V-M1',66.3,
    'V-M2',37.2,'V-M3',37.2,'V-M4',2,'V-M5',2,'Cr-K',5989,'Cr-L1',696,
    'Cr-L2',583.8,'Cr-L3',574.1,'Cr-M1',74.1,'Cr-M2',42.2,'Cr-M3',42.2,'Cr-M4',2,
    'Cr-M5',2,'Mn-K',6539,'Mn-L1',769.1,'Mn-L2',649.9,'Mn-L3',638.7,'Mn-M1',82.3,
    'Mn-M2',47.2,'Mn-M3',47.2,'Mn-M4',2,'Mn-M5',2,'Fe-K',7112,'Fe-L1',844.6,
    'Fe-L2',719.9,'Fe-L3',706.8,'Fe-M1',91.3,'Fe-M2',52.7,'Fe-M3',52.7,'Fe-M4',2,
    'Fe-M5',2,'Co-K',7709,'Co-L1',925.1,'Co-L2',793.2,'Co-L3',778.1,'Co-M1',101,
    'Co-M2',58.9,'Co-M3',59.9,'Co-M4',3,'Co-M5',3,'Ni-K',8333,'Ni-L1',1008.6,
    'Ni-L2',870,'Ni-L3',852.7,'Ni-M1',110.8,'Ni-M2',68,'Ni-M3',66.2,'Ni-M4',4,
    'Ni-M5',4,'Cu-K',8979,'Cu-L1',1096.7,'Cu-L2',952.3,'Cu-L3',932.7,'Cu-M1',122.5,
    'Cu-M2',77.3,'Cu-M3',75.1,'Cu-M4',5,'Cu-M5',5,'Zn-K',9659,'Zn-L1',1196.2,
    'Zn-L2',1044.9,'Zn-L3',1021.8,'Zn-M1',139.8,'Zn-M2',91.4,'Zn-M3',88.6,'Zn-M4',10.2,
    'Zn-M5',10.1,'Zn-N2',1,'Zn-N3',1,'Ga-K',10367,'Ga-L1',1299,'Ga-L2',1143.2,
    'Ga-L3',1116.4,'Ga-M1',159.51,'Ga-M2',103.5,'Ga-M3',100,'Ga-M4',18.7,'Ga-M5',18.7,
    'Ga-N1',1,'Ga-N2',2,'Ga-N3',2,'Ge-K',11103,'Ge-L1',1414.6,'Ge-L2',1248.1,
    'Ge-L3',1217,'Ge-M1',180.1,'Ge-M2',124.9,'Ge-M3',120.8,'Ge-M4',29.8,'Ge-M5',29.2,
    'Ge-N1',5,'Ge-N2',3,'Ge-N3',3,'As-K',11867,'As-L1',1527,'As-L2',1359.1,
    'As-L3',1323.6,'As-M1',204.7,'As-M2',146.2,'As-M3',141.2,'As-M4',41.7,'As-M5',41.7,
    'As-N1',8,'As-N2',3,'As-N3',3,'Se-K',12658,'Se-L1',1652,'Se-L2',1474.3,
    'Se-L3',1433.9,'Se-M1',229.6,'Se-M2',166.5,'Se-M3',160.7,'Se-M4',55.5,'Se-M5',54.6,
    'Se-N1',12,'Se-N2',3,'Se-N3',3,'Br-K',13474,'Br-L1',1782,'Br-L2',1596,
    'Br-L3',1550,'Br-M1',257,'Br-M2',189,'Br-M3',182,'Br-M4',70,'Br-M5',69,
    'Br-N1',27,'Br-N2',3,'Br-N3',3,'Kr-K',14326,'Kr-L1',1921,'Kr-L2',1730.9,
    'Kr-L3',1678.4,'Kr-M1',292.8,'Kr-M2',222.2,'Kr-M3',214.4,'Kr-M4',95,'Kr-M5',93.8,
    'Kr-N1',27.5,'Kr-N2',14.1,'Kr-N3',14.1,'Rb-K',15200,'Rb-L1',2065,'Rb-L2',1864,
    'Rb-L3',1804,'Rb-M1',326.7,'Rb-M2',248.7,'Rb-M3',239.1,'Rb-M4',113,'Rb-M5',112,
    'Rb-N1',30.5,'Rb-N2',16.3,'Rb-N3',15.3,'Sr-K',16105,'Sr-L1',2216,'Sr-L2',2007,
    'Sr-L3',1940,'Sr-M1',358.7,'Sr-M2',280.3,'Sr-M3',270,'Sr-M4',136,'Sr-M5',134.2,
    'Sr-N1',38.9,'Sr-N2',21.6,'Sr-N3',20.1,'Y-K',17038,'Y-L1',2373,'Y-L2',2156,
    'Y-L3',2080,'Y-M1',392,'Y-M2',310.6,'Y-M3',298.8,'Y-M4',157.7,'Y-M5',155.8,
    'Y-N1',43.8,'Y-N2',24.4,'Y-N3',23.1,'Zr-K',17998,'Zr-L1',2532,'Zr-L2',2307,
    'Zr-L3',2223,'Zr-M1',430.3,'Zr-M2',343.5,'Zr-M3',329.8,'Zr-M4',181.1,'Zr-M5',178.8,
    'Zr-N1',50.6,'Zr-N2',28.5,'Zr-N3',27.1,'Nb-K',18986,'Nb-L1',2698,'Nb-L2',2465,
    'Nb-L3',2371,'Nb-M1',466.6,'Nb-M2',376.1,'Nb-M3',360.6,'Nb-M4',205,'Nb-M5',202.3,
    'Nb-N1',56.4,'Nb-N2',32.6,'Nb-N3',30.8,'Mo-K',20000,'Mo-L1',2866,'Mo-L2',2625,
    'Mo-L3',2520,'Mo-M1',506.3,'Mo-M2',411.6,'Mo-M3',394,'Mo-M4',231.1,'Mo-M5',227.9,
    'Mo-N1',63.2,'Mo-N2',37.6,'Mo-N3',35.5,'Tc-K',21044,'Tc-L1',3043,'Tc-L2',2793,
    'Tc-L3',2677,'Tc-M1',544,'Tc-M2',447.6,'Tc-M3',417.7,'Tc-M4',257.6,'Tc-M5',253.9,
    'Tc-N1',69.5,'Tc-N2',42.3,'Tc-N3',39.9,'Ru-K',22117,'Ru-L1',3224,'Ru-L2',2967,
    'Ru-L3',2838,'Ru-M1',586.1,'Ru-M2',483.3,'Ru-M3',461.5,'Ru-M4',284.2,'Ru-M5',280,
    'Ru-N1',75,'Ru-N2',46.3,'Ru-N3',43.2,'Rh-K',23220,'Rh-L1',3412,'Rh-L2',3146,
    'Rh-L3',3004,'Rh-M1',628.1,'Rh-M2',521.3,'Rh-M3',496.5,'Rh-M4',311.9,'Rh-M5',307.2,
    'Rh-N1',81.4,'Rh-N2',50.5,'Rh-N3',47.3,'Rh-N4',2,'Rh-N5',2,'Pd-K',24350,
    'Pd-L1',3604,'Pd-L2',3330,'Pd-L3',3173,'Pd-M1',671.6,'Pd-M2',559.9,'Pd-M3',532.3,
    'Pd-M4',340.5,'Pd-M5',335.2,'Pd-N1',87.1,'Pd-N2',55.7,'Pd-N3',50.9,'Pd-N4',2,
    'Pd-N5',2,'Ag-K',25514,'Ag-L1',3806,'Ag-L2',3524,'Ag-L3',3351,'Ag-M1',719,
    'Ag-M2',603.8,'Ag-M3',573,'Ag-M4',374,'Ag-M5',368.3,'Ag-N1',97,'Ag-N2',63.7,
    'Ag-N3',58.3,'Ag-N4',4,'Ag-N5',4,'Cd-K',26711,'Cd-L1',4018,'Cd-L2',3727,
    'Cd-L3',3538,'Cd-M1',772,'Cd-M2',652.6,'Cd-M3',618.4,'Cd-M4',411.9,'Cd-M5',405.2,
    'Cd-N1',109.8,'Cd-N2',63.9,'Cd-N3',63.9,'Cd-N4',11.7,'Cd-N5',10.7,'In-K',27940,
    'In-L1',4238,'In-L2',3938,'In-L3',3730,'In-M1',827.2,'In-M2',703.2,'In-M3',665.3,
    'In-M4',451.4,'In-M5',443.9,'In-N1',122.9,'In-N2',73.5,'In-N3',73.5,'In-N4',17.7,
    'In-N5',16.9,'Sn-K',29200,'Sn-L1',4465,'Sn-L2',4156,'Sn-L3',3929,'Sn-M1',884.7,
    'Sn-M2',756.5,'Sn-M3',714.6,'Sn-M4',493.2,'Sn-M5',484.9,'Sn-N1',137.1,'Sn-N2',83.6,
    'Sn-N3',83.6,'Sn-N4',24.9,'Sn-N5',23.9,'Sb-K',30491,'Sb-L1',4698,'Sb-L2',4380,
    'Sb-L3',4132,'Sb-M1',940,'Sb-M2',812.7,'Sb-M3',766.4,'Sb-M4',537.5,'Sb-M5',528.2,
    'Sb-N1',153.2,'Sb-N2',95.6,'Sb-N3',95.6,'Sb-N4',33.3,'Sb-N5',32.1,'Sb-O1',7,
    'Sb-O2',2,'Sb-O3',2,'Te-K',31814,'Te-L1',4939,'Te-L2',4612,'Te-L3',4341,
    'Te-M1',1006,'Te-M2',870.8,'Te-M3',820.8,'Te-M4',583.4,'Te-M5',573,'Te-N1',169.4,
    'Te-N2',103.3,'Te-N3',103.3,'Te-N4',41.9,'Te-N5',40.4,'Te-O1',12,'Te-O2',2,
    'Te-O3',2,'I-K',33169,'I-L1',5188,'I-L2',4852,'I-L3',4557,'I-M1',1072,
    'I-M2',931,'I-M3',875,'I-M4',630.8,'I-M5',619.3,'I-N1',186,'I-N2',123,
    'I-N3',123,'I-N4',50.6,'I-N5',48.9,'I-O1',14,'I-O2',3,'I-O3',3,
    'Xe-K',34561,'Xe-L1',5453,'Xe-L2',5107,'Xe-L3',4786,'Xe-M1',1148.7,'Xe-M2',1002.1,
    'Xe-M3',940.6,'Xe-M4',689,'Xe-M5',676.4,'Xe-N1',213.2,'Xe-N2',146.7,'Xe-N3',145.5,
    'Xe-N4',69.5,'Xe-N5',67.5,'Xe-O1',23.3,'Xe-O2',13.4,'Xe-O3',12.1,'Cs-K',35985,
    'Cs-L1',5714,'Cs-L2',5359,'Cs-L3',5012,'Cs-M1',1211,'Cs-M2',1071,'Cs-M3',1003,
    'Cs-M4',740.5,'Cs-M5',726.6,'Cs-N1',232.3,'Cs-N2',172.4,'Cs-N3',161.3,'Cs-N4',79.8,
    'Cs-N5',77.5,'Cs-O1',22.7,'Cs-O2',14.2,'Cs-O3',12.1,'Ba-K',37441,'Ba-L1',5989,
    'Ba-L2',5624,'Ba-L3',5247,'Ba-M1',1293,'Ba-M2',1137,'Ba-M3',1063,'Ba-M4',795.7,
    'Ba-M5',780.5,'Ba-N1',253.5,'Ba-N2',192,'Ba-N3',178.6,'Ba-N4',92.6,'Ba-N5',89.9,
    'Ba-O1',30.3,'Ba-O2',17,'Ba-O3',14.8,'La-K',38925,'La-L1',6266,'La-L2',5891,
    'La-L3',5483,'La-M1',1362,'La-M2',1209,'La-M3',1128,'La-M4',853,'La-M5',836,
    'La-N1',274.7,'La-N2',205.8,'La-N3',196,'La-N4',105.3,'La-N5',102.5,'La-O1',34.3,
    'La-O2',19.3,'La-O3',16.8,'Ce-K',40443,'Ce-L1',6548,'Ce-L2',6164,'Ce-L3',5723,
    'Ce-M1',1436,'Ce-M2',1274,'Ce-M3',1187,'Ce-M4',902.4,'Ce-M5',883.8,'Ce-N1',291,
    'Ce-N2',223.2,'Ce-N3',206.5,'Ce-N4',109,'Ce-N5',109,'Ce-N6',1,'Ce-N7',1,
    'Ce-O1',37.8,'Ce-O2',19.8,'Ce-O3',17,'Pr-K',41991,'Pr-L1',6835,'Pr-L2',6440,
    'Pr-L3',5964,'Pr-M1',1511,'Pr-M2',1337,'Pr-M3',1242,'Pr-M4',948.3,'Pr-M5',928.8,
    'Pr-N1',304.5,'Pr-N2',236.3,'Pr-N3',217.6,'Pr-N4',115.1,'Pr-N5',115.1,'Pr-N6',2,
    'Pr-N7',2,'Pr-O1',37.4,'Pr-O2',22.3,'Pr-O3',22.3,'Nd-K',43569,'Nd-L1',7126,
    'Nd-L2',6722,'Nd-L3',6208,'Nd-M1',1575,'Nd-M2',1403,'Nd-M3',1297,'Nd-M4',1003.3,
    'Nd-M5',980.4,'Nd-N1',319.2,'Nd-N2',243.3,'Nd-N3',224.6,'Nd-N4',120.5,'Nd-N5',120.5,
    'Nd-N6',1.5,'Nd-N7',1.5,'Nd-O1',37.5,'Nd-O2',21.1,'Nd-O3',21.1,'Pm-K',45184,
    'Pm-L1',7428,'Pm-L2',7013,'Pm-L3',6459,'Pm-M1',1650,'Pm-M2',1471.4,'Pm-M3',1357,
    'Pm-M4',1052,'Pm-M5',1027,'Pm-N1',331,'Pm-N2',242,'Pm-N3',242,'Pm-N4',120,
    'Pm-N5',120,'Pm-N6',4,'Pm-N7',4,'Pm-O1',38,'Pm-O2',22,'Pm-O3',22,
    'Sm-K',46834,'Sm-L1',7737,'Sm-L2',7312,'Sm-L3',6716,'Sm-M1',1723,'Sm-M2',1541,
    'Sm-M3',1419.8,'Sm-M4',1110.9,'Sm-M5',1083.4,'Sm-N1',347.2,'Sm-N2',265.6,'Sm-N3',247.4,
    'Sm-N4',129,'Sm-N5',129,'Sm-N6',5.2,'Sm-N7',5.2,'Sm-O1',37.4,'Sm-O2',21.3,
    'Sm-O3',21.3,'Eu-K',48519,'Eu-L1',8052,'Eu-L2',7617,'Eu-L3',6977,'Eu-M1',1800,
    'Eu-M2',1614,'Eu-M3',1481,'Eu-M4',1158.6,'Eu-M5',1127.5,'Eu-N1',360,'Eu-N2',284,
    'Eu-N3',257,'Eu-N4',133,'Eu-N5',127.7,'Eu-N6',6,'Eu-N7',6,'Eu-O1',32,
    'Eu-O2',22,'Eu-O3',22,'Gd-K',50239,'Gd-L1',8376,'Gd-L2',7930,'Gd-L3',7243,
    'Gd-M1',1881,'Gd-M2',1688,'Gd-M3',1544,'Gd-M4',1221.9,'Gd-M5',1189.6,'Gd-N1',378.6,
    'Gd-N2',286,'Gd-N3',271,'Gd-N4',142.6,'Gd-N5',142.6,'Gd-N6',8.6,'Gd-N7',8.6,
    'Gd-O1',36,'Gd-O2',20,'Gd-O3',20,'Tb-K',51996,'Tb-L1',8708,'Tb-L2',8252,
    'Tb-L3',7514,'Tb-M1',1968,'Tb-M2',1768,'Tb-M3',1611,'Tb-M4',1276.9,'Tb-M5',1241.1,
    'Tb-N1',396,'Tb-N2',322.4,'Tb-N3',284.1,'Tb-N4',150.5,'Tb-N5',150.5,'Tb-N6',7.7,
    'Tb-N7',2.4,'Tb-O1',45.6,'Tb-O2',28.7,'Tb-O3',22.6,'Dy-K',53789,'Dy-L1',9046,
    'Dy-L2',8581,'Dy-L3',7790,'Dy-M1',2047,'Dy-M2',1842,'Dy-M3',1676,'Dy-M4',1333,
    'Dy-M5',1292,'Dy-N1',414.2,'Dy-N2',333.5,'Dy-N3',293.2,'Dy-N4',153.6,'Dy-N5',153.6,
    'Dy-N6',8,'Dy-N7',4.3,'Dy-O1',49.9,'Dy-O2',26.3,'Dy-O3',26.3,'Ho-K',55618,
    'Ho-L1',9394,'Ho-L2',8918,'Ho-L3',8071,'Ho-M1',2128,'Ho-M2',1923,'Ho-M3',1741,
    'Ho-M4',1392,'Ho-M5',1351,'Ho-N1',432.4,'Ho-N2',343.5,'Ho-N3',308.2,'Ho-N4',160,
    'Ho-N5',160,'Ho-N6',8.6,'Ho-N7',5.2,'Ho-O1',49.3,'Ho-O2',30.8,'Ho-O3',24.1,
    'Er-K',57486,'Er-L1',9751,'Er-L2',9264,'Er-L3',8358,'Er-M1',2206,'Er-M2',2006,
    'Er-M3',1812,'Er-M4',1453,'Er-M5',1409,'Er-N1',449.8,'Er-N2',366.2,'Er-N3',320.2,
    'Er-N4',167.6,'Er-N5',167.6,'Er-N6',4.7,'Er-N7',4.7,'Er-O1',50.6,'Er-O2',31.4,
    'Er-O3',24.7,'Tm-K',59390,'Tm-L1',10116,'Tm-L2',9617,'Tm-L3',8648,'Tm-M1',2307,
    'Tm-M2',2090,'Tm-M3',1885,'Tm-M4',1515,'Tm-M5',1468,'Tm-N1',470.9,'Tm-N2',385.9,
    'Tm-N3',332.6,'Tm-N4',175.5,'Tm-N5',175.5,'Tm-N6',4.6,'Tm-N7',4.6,'Tm-O1',54.7,
    'Tm-O2',31.8,'Tm-O3',25,'Yb-K',61332,'Yb-L1',10486,'Yb-L2',9978,'Yb-L3',8944,
    'Yb-M1',2398,'Yb-M2',2173,'Yb-M3',1950,'Yb-M4',1576,'Yb-M5',1528,'Yb-N1',480.5,
    'Yb-N2',388.7,'Yb-N3',339.7,'Yb-N4',191.2,'Yb-N5',182.4,'Yb-N6',2.5,'Yb-N7',1.3,
    'Yb-O1',52,'Yb-O2',30.3,'Yb-O3',24.1,'Lu-K',63314,'Lu-L1',10870,'Lu-L2',10349,
    'Lu-L3',9244,'Lu-M1',2491,'Lu-M2',2264,'Lu-M3',2024,'Lu-M4',1639,'Lu-M5',1589,
    'Lu-N1',506.8,'Lu-N2',412.4,'Lu-N3',359.2,'Lu-N4',206.1,'Lu-N5',196.3,'Lu-N6',8.9,
    'Lu-N7',7.5,'Lu-O1',57.3,'Lu-O2',33.6,'Lu-O3',26.7,'Hf-K',65351,'Hf-L1',11271,
    'Hf-L2',10739,'Hf-L3',9561,'Hf-M1',2601,'Hf-M2',2365,'Hf-M3',2107,'Hf-M4',1716,
    'Hf-M5',1662,'Hf-N1',538,'Hf-N2',438.2,'Hf-N3',380.7,'Hf-N4',220,'Hf-N5',211.5,
    'Hf-N6',15.9,'Hf-N7',14.2,'Hf-O1',64.2,'Hf-O2',38,'Hf-O3',29.9,'Ta-K',67416,
    'Ta-L1',11682,'Ta-L2',11136,'Ta-L3',9881,'Ta-M1',2708,'Ta-M2',2469,'Ta-M3',2194,
    'Ta-M4',1793,'Ta-M5',1735,'Ta-N1',563.4,'Ta-N2',463.4,'Ta-N3',400.9,'Ta-N4',237.9,
    'Ta-N5',226.4,'Ta-N6',23.5,'Ta-N7',21.6,'Ta-O1',69.7,'Ta-O2',42.2,'Ta-O3',32.7,
    'W-K',69525,'W-L1',12100,'W-L2',11544,'W-L3',10207,'W-M1',2820,'W-M2',2575,
    'W-M3',2281,'W-M4',1872,'W-M5',1809,'W-N1',594.1,'W-N2',490.4,'W-N3',423.61,
    'W-N4',255.9,'W-N5',243.5,'W-N6',33.6,'W-N7',31.4,'W-O1',75.6,'W-O2',45.3,
    'W-O3',36.8,'Re-K',71676,'Re-L1',12527,'Re-L2',11959,'Re-L3',10535,'Re-M1',2932,
    'Re-M2',2682,'Re-M3',2367,'Re-M4',1949,'Re-M5',1883,'Re-N1',625.4,'Re-N2',518.7,
    'Re-N3',446.8,'Re-N4',273.9,'Re-N5',260.5,'Re-N6',42.9,'Re-N7',40.5,'Re-O1',83,
    'Re-O2',45.6,'Re-O3',34.6,'Os-K',73871,'Os-L1',12968,'Os-L2',12385,'Os-L3',10871,
    'Os-M1',3049,'Os-M2',2792,'Os-M3',2457,'Os-M4',2031,'Os-M5',1960,'Os-N1',658.2,
    'Os-N2',549.1,'Os-N3',470.7,'Os-N4',293.1,'Os-N5',278.5,'Os-N6',53.4,'Os-N7',50.7,
    'Os-O1',84,'Os-O2',58,'Os-O3',44.5,'Ir-K',76111,'Ir-L1',13419,'Ir-L2',12824,
    'Ir-L3',11215,'Ir-M1',3174,'Ir-M2',2909,'Ir-M3',2551,'Ir-M4',2116,'Ir-M5',2040,
    'Ir-N1',691.1,'Ir-N2',577.8,'Ir-N3',495.8,'Ir-N4',311.9,'Ir-N5',296.3,'Ir-N6',63.8,
    'Ir-N7',60.8,'Ir-O1',95.2,'Ir-O2',63,'Ir-O3',48,'Pt-K',78395,'Pt-L1',13880,
    'Pt-L2',13273,'Pt-L3',11564,'Pt-M1',3296,'Pt-M2',3027,'Pt-M3',2645,'Pt-M4',2202,
    'Pt-M5',2122,'Pt-N1',725.4,'Pt-N2',609.1,'Pt-N3',519.4,'Pt-N4',331.6,'Pt-N5',314.6,
    'Pt-N6',74.5,'Pt-N7',71.2,'Pt-O1',101.7,'Pt-O2',65.3,'Pt-O3',51.7,'Au-K',80725,
    'Au-L1',14353,'Au-L2',13734,'Au-L3',11919,'Au-M1',3425,'Au-M2',3148,'Au-M3',2743,
    'Au-M4',2291,'Au-M5',2206,'Au-N1',762.1,'Au-N2',642.7,'Au-N3',546.3,'Au-N4',353.2,
    'Au-N5',335.1,'Au-N6',87.6,'Au-N7',83.9,'Au-O1',107.2,'Au-O2',74.2,'Au-O3',57.2,
    'Au-O4',5,'Au-O5',5,'Hg-K',83102,'Hg-L1',14839,'Hg-L2',14209,'Hg-L3',12284,
    'Hg-M1',3562,'Hg-M2',3279,'Hg-M3',2847,'Hg-M4',2385,'Hg-M5',2295,'Hg-N1',802.2,
    'Hg-N2',680.2,'Hg-N3',576.6,'Hg-N4',378.2,'Hg-N5',358.8,'Hg-N6',104,'Hg-N7',99.9,
    'Hg-O1',127,'Hg-O2',83.1,'Hg-O3',64.5,'Hg-O4',9.6,'Hg-O5',7.8,'Tl-K',85530,
    'Tl-L1',15347,'Tl-L2',14698,'Tl-L3',12658,'Tl-M1',3704,'Tl-M2',3416,'Tl-M3',2957,
    'Tl-M4',2485,'Tl-M5',2389,'Tl-N1',846.2,'Tl-N2',720.5,'Tl-N3',609.5,'Tl-N4',405.7,
    'Tl-N5',385,'Tl-N6',122.2,'Tl-N7',117.8,'Tl-O1',136,'Tl-O2',94.6,'Tl-O3',73.5,
    'Tl-O4',14.7,'Tl-O5',12.5,'Pb-K',88005,'Pb-L1',15861,'Pb-L2',15200,'Pb-L3',13035,
    'Pb-M1',3851,'Pb-M2',3554,'Pb-M3',3066,'Pb-M4',2586,'Pb-M5',2484,'Pb-N1',891.8,
    'Pb-N2',761.9,'Pb-N3',643.5,'Pb-N4',434.3,'Pb-N5',412.2,'Pb-N6',141.7,'Pb-N7',136.9,
    'Pb-O1',147,'Pb-O2',106.4,'Pb-O3',83.3,'Pb-O4',20.7,'Pb-O5',18.1,'Pb-P1',3,
    'Pb-P2',1,'Pb-P3',1,'Bi-K',90526,'Bi-L1',16388,'Bi-L2',15711,'Bi-L3',13419,
    'Bi-M1',3999,'Bi-M2',3696,'Bi-M3',3177,'Bi-M4',2688,'Bi-M5',2580,'Bi-N1',939,
    'Bi-N2',805.2,'Bi-N3',678.8,'Bi-N4',464,'Bi-N5',440.1,'Bi-N6',162.3,'Bi-N7',157,
    'Bi-O1',159.3,'Bi-O2',119,'Bi-O3',92.6,'Bi-O4',26.9,'Bi-O5',23.8,'Bi-P1',8,
    'Bi-P2',3,'Bi-P3',3,'Po-K',93105,'Po-L1',16939,'Po-L2',16244,'Po-L3',13814,
    'Po-M1',4149,'Po-M2',3854,'Po-M3',3302,'Po-M4',2798,'Po-M5',2683,'Po-N1',995,
    'Po-N2',851,'Po-N3',705,'Po-N4',500,'Po-N5',473,'Po-N6',184,'Po-N7',184,
    'Po-O1',177,'Po-O2',132,'Po-O3',104,'Po-O4',31,'Po-O5',31,'Po-P1',9,
    'Po-P2',4,'Po-P3',1,'At-K',95730,'At-L1',17493,'At-L2',16785,'At-L3',14214,
    'At-M1',4317,'At-M2',4008,'At-M3',3426,'At-M4',2909,'At-M5',2787,'At-N1',1042,
    'At-N2',886,'At-N3',740,'At-N4',533,'At-N5',507,'At-N6',210,'At-N7',210,
    'At-O1',195,'At-O2',148,'At-O3',115,'At-O4',40,'At-O5',40,'At-P1',13,
    'At-P2',6,'At-P3',1,'Rn-K',98404,'Rn-L1',18049,'Rn-L2',17337,'Rn-L3',14619,
    'Rn-M1',4482,'Rn-M2',4159,'Rn-M3',3538,'Rn-M4',3022,'Rn-M5',2892,'Rn-N1',1097,
    'Rn-N2',929,'Rn-N3',768,'Rn-N4',567,'Rn-N5',541,'Rn-N6',238,'Rn-N7',238,
    'Rn-O1',214,'Rn-O2',164,'Rn-O3',127,'Rn-O4',48,'Rn-O5',48,'Rn-P1',16,
    'Rn-P2',8,'Rn-P3',2,'Fr-K',101137,'Fr-L1',18639,'Fr-L2',17907,'Fr-L3',15031,
    'Fr-M1',4652,'Fr-M2',4327,'Fr-M3',3663,'Fr-M4',3136,'Fr-M5',3000,'Fr-N1',1153,
    'Fr-N2',980,'Fr-N3',810,'Fr-N4',603,'Fr-N5',577,'Fr-N6',268,'Fr-N7',268,
    'Fr-O1',234,'Fr-O2',182,'Fr-O3',140,'Fr-O4',58,'Fr-O5',58,'Fr-P1',24,
    'Fr-P2',14,'Fr-P3',7,'Ra-K',103922,'Ra-L1',19237,'Ra-L2',18484,'Ra-L3',15444,
    'Ra-M1',4822,'Ra-M2',4490,'Ra-M3',3792,'Ra-M4',3248,'Ra-M5',3105,'Ra-N1',1208,
    'Ra-N2',1058,'Ra-N3',879,'Ra-N4',636,'Ra-N5',603,'Ra-N6',299,'Ra-N7',299,
    'Ra-O1',254,'Ra-O2',200,'Ra-O3',153,'Ra-O4',68,'Ra-O5',68,'Ra-P1',31,
    'Ra-P2',20,'Ra-P3',12,'Ac-K',106755,'Ac-L1',19840,'Ac-L2',19083,'Ac-L3',15871,
    'Ac-M1',5002,'Ac-M2',4656,'Ac-M3',3909,'Ac-M4',3370,'Ac-M5',3219,'Ac-N1',1269,
    'Ac-N2',1080,'Ac-N3',890,'Ac-N4',675,'Ac-N5',639,'Ac-N6',319,'Ac-N7',319,
    'Ac-O1',272,'Ac-O2',215,'Ac-O3',167,'Ac-O4',80,'Ac-O5',80,'Ac-P1',37,
    'Ac-P2',24,'Ac-P3',15,'Th-K',109651,'Th-L1',20472,'Th-L2',19693,'Th-L3',16300,
    'Th-M1',5182,'Th-M2',4830,'Th-M3',4046,'Th-M4',3491,'Th-M5',3332,'Th-N1',1330,
    'Th-N2',1168,'Th-N3',966.4,'Th-N4',712.1,'Th-N5',675.2,'Th-N6',342.4,'Th-N7',333.1,
    'Th-O1',290,'Th-O2',229,'Th-O3',182,'Th-O4',92.5,'Th-O5',85.4,'Th-P1',41.4,
    'Th-P2',24.5,'Th-P3',16.6,'Pa-K',112601,'Pa-L1',21105,'Pa-L2',20314,'Pa-L3',16733,
    'Pa-M1',5367,'Pa-M2',5001,'Pa-M3',4174,'Pa-M4',3611,'Pa-M5',3442,'Pa-N1',1387,
    'Pa-N2',1224,'Pa-N3',1007,'Pa-N4',743,'Pa-N5',708,'Pa-N6',371,'Pa-N7',360,
    'Pa-O1',310,'Pa-O2',232,'Pa-O3',187,'Pa-O4',94,'Pa-O5',94,'Pa-P1',43,
    'Pa-P2',27,'Pa-P3',17,'U-K',115606,'U-L1',21757,'U-L2',20948,'U-L3',17166,
    'U-M1',5548,'U-M2',5182,'U-M3',4303,'U-M4',3728,'U-M5',3552,'U-N1',1439,
    'U-N2',1271,'U-N3',1043,'U-N4',778.3,'U-N5',736.2,'U-N6',388.2,'U-N7',377.4,
    'U-O1',321,'U-O2',257,'U-O3',192,'U-O4',102.8,'U-O5',94.2,'U-P1',43.9,
    'U-P2',26.8,'U-P3',16.8,'Np-K',118669,'Np-L1',22427,'Np-L2',21600,'Np-L3',17610,
    'Np-M1',5739,'Np-M2',5366,'Np-M3',4435,'Np-M4',3849,'Np-M5',3664,'Np-N1',1501,
    'Np-N2',1328,'Np-N3',1085,'Np-N4',816,'Np-N5',771,'Np-N6',414,'Np-N7',403,
    'Np-O1',338,'Np-O2',274,'Np-O3',206,'Np-O4',109,'Np-O5',101,'Np-P1',47,
    'Np-P2',29,'Np-P3',18,'Pu-K',121791,'Pu-L1',23104,'Pu-L2',22266,'Pu-L3',18057,
    'Pu-M1',5933,'Pu-M2',5547,'Pu-M3',4563,'Pu-M4',3970,'Pu-M5',3775,'Pu-N1',1559,
    'Pu-N2',1380,'Pu-N3',1123,'Pu-N4',846,'Pu-N5',798,'Pu-N6',436,'Pu-N7',424,
    'Pu-O1',350,'Pu-O2',283,'Pu-O3',213,'Pu-O4',113,'Pu-O5',102,'Pu-P1',46,
    'Pu-P2',29,'Pu-P3',16,'Am-K',124982,'Am-L1',23808,'Am-L2',22952,'Am-L3',18510,
    'Am-M1',6133,'Am-M2',5739,'Am-M3',4698,'Am-M4',4096,'Am-M5',3890,'Am-N1',1620,
    'Am-N2',1438,'Am-N3',1165,'Am-N4',880,'Am-N5',829,'Am-N6',461,'Am-N7',446,
    'Am-O1',365,'Am-O2',298,'Am-O3',219,'Am-O4',116,'Am-O5',106,'Am-P1',48,
    'Am-P2',29,'Am-P3',16,'Cm-K',128241,'Cm-L1',24526,'Cm-L2',23651,'Cm-L3',18970,
    'Cm-M1',6337,'Cm-M2',5937,'Cm-M3',4838,'Cm-M4',4224,'Cm-M5',4009,'Cm-N1',1684,
    'Cm-N2',1498,'Cm-N3',1207,'Cm-N4',916,'Cm-N5',862,'Cm-N6',484,'Cm-N7',470,
    'Cm-O1',383,'Cm-O2',313,'Cm-O3',229,'Cm-O4',124,'Cm-O5',110,'Cm-P1',50,
    'Cm-P2',30,'Cm-P3',16,'Bk-K',131556,'Bk-L1',25256,'Bk-L2',24371,'Bk-L3',19435,
    'Bk-M1',6545,'Bk-M2',6138,'Bk-M3',4976,'Bk-M4',4353,'Bk-M5',4127,'Bk-N1',1748,
    'Bk-N2',1558,'Bk-N3',1249,'Bk-N4',955,'Bk-N5',898,'Bk-N6',511,'Bk-N7',495,
    'Bk-O1',399,'Bk-O2',326,'Bk-O3',237,'Bk-O4',130,'Bk-O5',117,'Bk-P1',52,
    'Bk-P2',32,'Bk-P3',16,'Cf-K',134939,'Cf-L1',26010,'Cf-L2',25108,'Cf-L3',19907,
    'Cf-M1',6761,'Cf-M2',6345,'Cf-M3',5116,'Cf-M4',4484,'Cf-M5',4247,'Cf-N1',1813,
    'Cf-N2',1620,'Cf-N3',1292,'Cf-N4',991,'Cf-N5',930,'Cf-N6',538,'Cf-N7',520,
    'Cf-O1',416,'Cf-O2',341,'Cf-O3',245,'Cf-O4',137,'Cf-O5',122,'Cf-P1',54,
    'Cf-P2',33,'Cf-P3',17);

my %Kedges = ('H', 13.6, 'He', 24.6, 'Li', 54.7, 'Be', 111.5, 'B', 188, 'C', 284.2, 
              'N', 409.9, 'O', 543.1, 'F', 696.7, 'Ne', 870.2, 'Na', 1070.8, 'Mg', 1303, 
              'Al', 1559, 'Si', 1839, 'P', 2145.5, 'S', 2472, 'Cl', 2822.4, 
              'Ar', 3205.9, 'K', 3608.4, 'Ca', 4038.5, 'Sc', 4492, 'Ti', 4966, 
              'V', 5465, 'Cr', 5989, 'Mn', 6539, 'Fe', 7112, 'Co', 7709, 
              'Ni', 8333, 'Cu', 8979, 'Zn', 9659, 'Ga', 10367, 'Ge', 11103,
              'As', 11867, 'Se', 12658, 'Br', 13474, 'Kr', 14326, 'Rb', 15200,
              'Sr', 16105, 'Y', 17038, 'Zr', 17998, 'Nb', 18986, 'Mo', 20000, 
              'Tc', 21044, 'Ru', 22117, 'Rh', 23220, 'Pd', 24350, 'Ag', 25514, 
              'Cd', 26711, 'In', 27940, 'Sn', 29200, 'Sb', 30491, 'Te', 31814, 
              'I', 33169, 'Xe', 34561, 'Cs', 35985, 'Ba', 37441, 'La', 38925, 
              'Ce', 40443, 'Pr', 41991, 'Nd', 43569, 'Pm', 45184, 'Sm', 46834, 
              'Eu', 48519, 'Gd', 50239, 'Tb', 51996, 'Dy', 53789, 'Ho', 55618, 
              'Er', 57486, 'Tm', 59390, 'Yb', 61332, 'Lu', 63314, 'Hf', 65351, 
              'Ta', 67416, 'W', 69525, 'Re', 71676, 'Os', 73871, 'Ir', 76111, 
              'Pt', 78395, 'Au', 80725, 'Hg', 83102, 'Tl', 85530, 'Pb', 88005, 
              'Bi', 90526, 'Po', 93105, 'At', 95730, 'Rn', 98404, 'Fr',101137, 
              'Ra',103922, 'Ac',106755, 'Th',109651, 'Pa',112601, 'U',115606 );

my %L1edges =('N',  37.3, 'O',  41.6, 'F',    0., 'Ne', 48.5, 'Na', 63.5, 'Mg', 88.6,
              'Al',117.8, 'Si',149.7, 'P', 189, 'S', 230.9, 'Cl',270.0, 'Ar',326.3,
              'K', 378.6, 'Ca',438.4, 'Sc',498, 'Ti',560.9, 'V', 626.7,
              'Cr',696, 'Mn',769.1, 'Fe',844.6, 'Co',925.1, 'Ni',1008.6,
              'Cu',1096.7, 'Zn',1196.2, 'Ga', 1299, 'Ge', 1414.6, 'As', 1527, 'Se', 1652,
              'Br', 1782, 'Kr', 1921, 'Rb', 2065, 'Sr', 2216, 'Y',  2373, 'Zr', 2532,
              'Nb', 2698, 'Mo', 2866, 'Tc', 3043, 'Ru', 3224, 'Rh', 3412,
              'Pd', 3604, 'Ag', 3806, 'Cd', 4018, 'In', 4238, 'Sn', 4465,
              'Sb', 4698, 'Te', 4939, 'I',  5188, 'Xe', 5453, 'Cs', 5714,
              'Ba', 5989, 'La', 6266, 'Ce', 6548, 'Pr', 6835, 'Nd', 7126,
              'Pm', 7428, 'Sm', 7737, 'Eu', 8052, 'Gd', 8376, 'Tb', 8708,
              'Dy', 9046, 'Ho', 9394, 'Er', 9751, 'Tm', 10116, 'Yb', 10486,
              'Lu', 10870, 'Hf', 11271, 'Ta', 11682, 'W',  12100, 'Re', 12527,
              'Os', 12968, 'Ir', 13419, 'Pt', 13880, 'Au', 14353, 'Hg', 14839,
              'Tl', 15347, 'Pb', 15861, 'Bi', 16388, 'Po', 16939, 'At', 17493,
              'Rn', 18049, 'Fr', 18639, 'Ra', 19237, 'Ac', 19840, 'Th', 20472,
              'Pa', 21105, 'U' , 21757 );

my %L2edges = ('F', 48.5, 'Ne',63.5, 'Na',88.6, 'Mg',17.8, 'Al',49.7, 'Si',89,
               'P', 136, 'S', 163.6, 'Cl',202, 'Ar',250.6, 'K', 297.3, 'Ca',349.7,
               'Sc',403.6, 'Ti',460.2, 'V', 519.8, 'Cr',583.8, 'Mn',649.9, 'Fe',719.9,
               'Co',793.2, 'Ni',870, 'Cu',952.3, 'Zn', 1044.9, 'Ga', 1143.2, 'Ge', 1248.1,
               'As', 1359.1, 'Se', 1474.3, 'Br', 1596, 'Kr', 1730.9, 'Rb', 1864,
               'Sr', 2007, 'Y',  2156, 'Zr', 2307, 'Nb', 2465, 'Mo', 2625,
               'Tc', 2793, 'Ru', 2967, 'Rh', 3146, 'Pd', 3330, 'Ag', 3524,
               'Cd', 3727, 'In', 3938, 'Sn', 4156, 'Sb', 4380, 'Te', 4612,
               'I',  4852, 'Xe', 5107, 'Cs', 5359, 'Ba', 5624, 'La', 5891,
               'Ce', 6164, 'Pr', 6440, 'Nd', 6722, 'Pm', 7013, 'Sm', 7312,
               'Eu', 7617, 'Gd', 7930, 'Tb', 8252, 'Dy', 8581, 'Ho', 8918,
               'Er', 9264, 'Tm', 9617, 'Yb', 9978, 'Lu', 10349, 'Hf', 10739,
               'Ta', 11136, 'W',  11544, 'Re', 11959, 'Os', 12385, 'Ir', 12824,
               'Pt', 13273, 'Au', 13734, 'Hg', 14209, 'Tl', 14698, 'Pb', 15200,
               'Bi', 15711, 'Po', 16244, 'At', 16785, 'Rn', 17337, 'Fr', 17907,
               'Ra', 18484, 'Ac', 19083, 'Th', 19693, 'Pa', 20314, 'U' , 20948);

my %L3edges = ('Ne',21.6, 'Na',30.5, 'Mg',49.2, 'Al',72.5, 'Si',99.2, 'P', 135, 'S', 162.5, 
              'Cl',200, 'Ar',248.4, 'K', 294.6, 'Ca',346.2, 'Sc',398.7, 'Ti',453.8,
              'V', 512.1, 'Cr',574.1, 'Mn',638.7, 'Fe',706.8, 'Co',778.1, 'Ni',852.7,
              'Cu',932.7, 'Zn', 1021.8, 'Ga', 1116.4, 'Ge', 1217, 'As', 1323.6,
              'Se', 1433.9, 'Br', 1550, 'Kr', 1678.4, 'Rb', 1804, 'Sr', 1940,
              'Y',  2080, 'Zr', 2223, 'Nb', 2371, 'Mo', 2520, 'Tc', 2677,
              'Ru', 2838, 'Rh', 3004, 'Pd', 3173, 'Ag', 3351, 'Cd', 3538,
              'In', 3730, 'Sn', 3929, 'Sb', 4132, 'Te', 4341, 'I',  4557,
              'Xe', 4786, 'Cs', 5012, 'Ba', 5247, 'La', 5483, 'Ce', 5723,
              'Pr', 5964, 'Nd', 6208, 'Pm', 6459, 'Sm', 6716, 'Eu', 6977,
              'Gd', 7243, 'Tb', 7510, 'Dy', 7790, 'Ho', 8071, 'Er', 8358,
              'Tm', 8648, 'Yb', 8944, 'Lu', 9244, 'Hf', 9561, 'Ta', 9881,
              'W',  10207, 'Re', 10535, 'Os', 10871, 'Ir', 11215, 'Pt', 11564,
              'Au', 11919, 'Hg', 12284, 'Tl', 12658, 'Pb', 13035, 'Bi', 13419,
              'Po', 13814, 'At', 14214, 'Rn', 14619, 'Fr', 15031, 'Ra', 15444,
              'Ac', 15871, 'Th', 16300, 'Pa', 16733, 'U' , 17166) ;




  my @z2symb =          ( '', 'H' , 'He', 'Li', 'Be', 'B' , 'C' , 'N' ,
      'O' , 'F' , 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P' , 'S' , 'Cl', 'Ar',
      'K' , 'Ca', 'Sc', 'Ti', 'V' , 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu',
      'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y' , 'Zr',
      'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb',
      'Te', 'I' , 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm',
      'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf',
      'Ta', 'W' , 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi',
      'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U' , 'Np', 'Pu',
      'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db',
      'Sg', 'Bh', 'Hs', 'Mt' );


my $oceanData;
my $json = JSON::PP->new;
my $dataFile = "postDefaultsOceanDatafile";
if( open my $in, "<", $dataFile ) {
  local $/ = undef;
  $oceanData = $json->decode(<$in>);
  close($in);
} else {
  die "Failed to open $dataFile\n$!";
}


exit 0 if( scalar @{$oceanData->{'calc'}->{'edges'}}< 1 );


$oceanData->{'calc'}->{'edges'}[0] =~ m/(\d+)\s+(\d+)\s+(\d+)/ or die "Malformed calc->edges\n";
my $site = $1;
my $n = $2;
my $l = $3;

my $typat = $oceanData->{'structure'}->{'typat'}[$site-1];
my $znucl = $oceanData->{'structure'}->{'znucl'}[$typat-1];

my $zsymb = $z2symb[$znucl];
#print "$site $typat $znucl $zsymb\n";

my $energy = 0;
if( $n == 1 ) {
#  $energy = $Kedges{ $zsymb };
  $energy = $allEdges{ $zsymb .'-K' };
} elsif( $n == 2 ) {
  if( $l == 1 ) {
#    $energy = $L2edges{ $zsymb } + 2*$L3edges{ $zsymb };
#    $energy/= 3;
#    print "$energy   ";
    $energy = $allEdges{ $zsymb .'-L2' } + 2*$allEdges{ $zsymb .'-L3' };
    $energy /= 3;
#    print "$energy\n";
  } elsif( $l == 0 ) {
#    $energy = $L1edges{ $zsymb };
    $energy = $allEdges{ $zsymb .'-L1' };
  }
} elsif( $n == 3 ) {
  if( $l == 2 ) {
    $energy = 2*$allEdges{ $zsymb .'-M4' } + 3*$allEdges{ $zsymb .'-M5' };
    $energy /= 5;
  } elsif( $l==1) {
    $energy = $allEdges{ $zsymb .'-M2' } + 2*$allEdges{ $zsymb .'-M3' };
    $energy /= 3;
  } elsif ( $l==0 ) {
    $energy = $allEdges{ $zsymb .'-M1' };
  }
} elsif( $n==4 ) {
  if ($l == 3 ) {
    $energy = 3*$allEdges{ $zsymb .'-N6' } + 4*$allEdges{ $zsymb .'-N7' };
    $energy /= 7;
  } elsif( $l == 2 ) {
    $energy = 2*$allEdges{ $zsymb .'-N4' } + 3*$allEdges{ $zsymb .'-N5' };
    $energy /= 5; 
  } elsif( $l==1) {
    $energy = $allEdges{ $zsymb .'-N2' } + 2*$allEdges{ $zsymb .'-N3' };
    $energy /= 3;
  } elsif ( $l==0 ) {
    $energy = $allEdges{ $zsymb .'-N1' };
  }
} elsif( $n==5 ) {
  if ($l == 3 ) {
    $energy = 3*$allEdges{ $zsymb .'-O6' } + 4*$allEdges{ $zsymb .'-O7' };
    $energy /= 7;
  } elsif( $l == 2 ) {
    $energy = 2*$allEdges{ $zsymb .'-O4' } + 3*$allEdges{ $zsymb .'-O5' };
    $energy /= 5;
  } elsif( $l==1) {
    $energy = $allEdges{ $zsymb .'-O2' } + 2*$allEdges{ $zsymb .'-O3' };
    $energy /= 3;
  } elsif ( $l==0 ) {
    $energy = $allEdges{ $zsymb .'-O1' };
  }
} elsif( $n==6 ) {
  if ($l == 3 ) {
    $energy = 3*$allEdges{ $zsymb .'-P6' } + 4*$allEdges{ $zsymb .'-P7' };
    $energy /= 7;
  } elsif( $l == 2 ) {
    $energy = 2*$allEdges{ $zsymb .'-P4' } + 3*$allEdges{ $zsymb .'-P5' };
    $energy /= 5;
  } elsif( $l==1) {
    $energy = $allEdges{ $zsymb .'-P2' } + 2*$allEdges{ $zsymb .'-P3' };
    $energy /= 3;
  } elsif ( $l==0 ) {
    $energy = $allEdges{ $zsymb .'-P1' };
  }
}


print "WARNING Unsupported edge!\n" if( $energy == 0 );
#print "$energy\n";

my @dir = ( [ 1, 0, 0 ], [ 0, 1, 0 ], [0, 0, 1] );

my $quad = 0;
$quad = 1 if( $energy > 4000 ) ;

my $nphoton = 0;
for( my $i = 0; $i < 3; $i++ ) {
  for( my $j = 0; $j<3; $j++ ) {
    next if( $j == $i );
    $nphoton ++;
    my $fileName = sprintf "default_photon%i", $nphoton;
    open OUT, ">", $fileName or die;
    if( $quad == 1 ) {
      print OUT "quad\n";
    } else {
      print OUT "dipole\n";
    }
    printf OUT "cartesian %i %i %i\nend\n", $dir[$i][0], $dir[$i][1], $dir[$i][2];
    printf OUT "cartesian %i %i %i\nend\n", $dir[$j][0], $dir[$j][1], $dir[$j][2];
    printf OUT "%.1f\n", $energy;
    close OUT;
    last unless( $quad == 1 );
  }
}