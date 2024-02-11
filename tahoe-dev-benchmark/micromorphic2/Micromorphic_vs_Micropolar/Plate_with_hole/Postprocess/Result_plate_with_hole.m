clc

clear

close all

Node938_cls=load('Along_diameter_classic_0938.txt');
Node939_cls=load('Along_diameter_classic_0939.txt');
Node941_cls=load('Along_diameter_classic_0941.txt');
Node943_cls=load('Along_diameter_classic_0943.txt');
Node945_cls=load('Along_diameter_classic_0945.txt');
Node947_cls=load('Along_diameter_classic_0947.txt');
Node949_cls=load('Along_diameter_classic_0949.txt');
Node951_cls=load('Along_diameter_classic_0951.txt');
Node953_cls=load('Along_diameter_classic_0953.txt');
Node955_cls=load('Along_diameter_classic_0955.txt');
Node957_cls=load('Along_diameter_classic_0957.txt');
Node959_cls=load('Along_diameter_classic_0959.txt');
Node961_cls=load('Along_diameter_classic_0961.txt');
Node963_cls=load('Along_diameter_classic_0963.txt');
Node965_cls=load('Along_diameter_classic_0965.txt');
Node967_cls=load('Along_diameter_classic_0967.txt');
Node969_cls=load('Along_diameter_classic_0969.txt');
Node971_cls=load('Along_diameter_classic_0971.txt');
Node973_cls=load('Along_diameter_classic_0973.txt');
Node975_cls=load('Along_diameter_classic_0975.txt');
Node977_cls=load('Along_diameter_classic_0977.txt');
Node979_cls=load('Along_diameter_classic_0979.txt');
Node981_cls=load('Along_diameter_classic_0981.txt');
Node983_cls=load('Along_diameter_classic_0983.txt');
Node985_cls=load('Along_diameter_classic_0985.txt');
Node987_cls=load('Along_diameter_classic_0987.txt');
Node989_cls=load('Along_diameter_classic_0989.txt');
Node991_cls=load('Along_diameter_classic_0991.txt');
Node993_cls=load('Along_diameter_classic_0993.txt');
Node995_cls=load('Along_diameter_classic_0995.txt');
Node997_cls=load('Along_diameter_classic_0997.txt');
Node999_cls=load('Along_diameter_classic_0999.txt');
Node1001_cls=load('Along_diameter_classic_1001.txt');
Node1003_cls=load('Along_diameter_classic_1003.txt');
Node1005_cls=load('Along_diameter_classic_1005.txt');
Node1007_cls=load('Along_diameter_classic_1007.txt');


Node938_MBC_com=load('Along_diameter_Micro_BC_combine_0938.txt');
Node939_MBC_com=load('Along_diameter_Micro_BC_combine_0939.txt');
Node941_MBC_com=load('Along_diameter_Micro_BC_combine_0941.txt');
Node943_MBC_com=load('Along_diameter_Micro_BC_combine_0943.txt');
Node945_MBC_com=load('Along_diameter_Micro_BC_combine_0945.txt');
Node947_MBC_com=load('Along_diameter_Micro_BC_combine_0947.txt');
Node949_MBC_com=load('Along_diameter_Micro_BC_combine_0949.txt');
Node951_MBC_com=load('Along_diameter_Micro_BC_combine_0951.txt');
Node953_MBC_com=load('Along_diameter_Micro_BC_combine_0953.txt');
Node955_MBC_com=load('Along_diameter_Micro_BC_combine_0955.txt');
Node957_MBC_com=load('Along_diameter_Micro_BC_combine_0957.txt');
Node959_MBC_com=load('Along_diameter_Micro_BC_combine_0959.txt');
Node961_MBC_com=load('Along_diameter_Micro_BC_combine_0961.txt');
Node963_MBC_com=load('Along_diameter_Micro_BC_combine_0963.txt');
Node965_MBC_com=load('Along_diameter_Micro_BC_combine_0965.txt');
Node967_MBC_com=load('Along_diameter_Micro_BC_combine_0967.txt');
Node969_MBC_com=load('Along_diameter_Micro_BC_combine_0969.txt');
Node971_MBC_com=load('Along_diameter_Micro_BC_combine_0971.txt');
Node973_MBC_com=load('Along_diameter_Micro_BC_combine_0973.txt');
Node975_MBC_com=load('Along_diameter_Micro_BC_combine_0975.txt');
Node977_MBC_com=load('Along_diameter_Micro_BC_combine_0977.txt');
Node979_MBC_com=load('Along_diameter_Micro_BC_combine_0979.txt');
Node981_MBC_com=load('Along_diameter_Micro_BC_combine_0981.txt');
Node983_MBC_com=load('Along_diameter_Micro_BC_combine_0983.txt');
Node985_MBC_com=load('Along_diameter_Micro_BC_combine_0985.txt');
Node987_MBC_com=load('Along_diameter_Micro_BC_combine_0987.txt');
Node989_MBC_com=load('Along_diameter_Micro_BC_combine_0989.txt');
Node991_MBC_com=load('Along_diameter_Micro_BC_combine_0991.txt');
Node993_MBC_com=load('Along_diameter_Micro_BC_combine_0993.txt');
Node995_MBC_com=load('Along_diameter_Micro_BC_combine_0995.txt');
Node997_MBC_com=load('Along_diameter_Micro_BC_combine_0997.txt');
Node999_MBC_com=load('Along_diameter_Micro_BC_combine_0999.txt');
Node1001_MBC_com=load('Along_diameter_Micro_BC_combine_1001.txt');
Node1003_MBC_com=load('Along_diameter_Micro_BC_combine_1003.txt');
Node1005_MBC_com=load('Along_diameter_Micro_BC_combine_1005.txt');
Node1007_MBC_com=load('Along_diameter_Micro_BC_combine_1007.txt');


Node938_MBC_rot=load('Along_diameter_Micro_BC_rot_0938.txt');
Node939_MBC_rot=load('Along_diameter_Micro_BC_rot_0939.txt');
Node941_MBC_rot=load('Along_diameter_Micro_BC_rot_0941.txt');
Node943_MBC_rot=load('Along_diameter_Micro_BC_rot_0943.txt');
Node945_MBC_rot=load('Along_diameter_Micro_BC_rot_0945.txt');
Node947_MBC_rot=load('Along_diameter_Micro_BC_rot_0947.txt');
Node949_MBC_rot=load('Along_diameter_Micro_BC_rot_0949.txt');
Node951_MBC_rot=load('Along_diameter_Micro_BC_rot_0951.txt');
Node953_MBC_rot=load('Along_diameter_Micro_BC_rot_0953.txt');
Node955_MBC_rot=load('Along_diameter_Micro_BC_rot_0955.txt');
Node957_MBC_rot=load('Along_diameter_Micro_BC_rot_0957.txt');
Node959_MBC_rot=load('Along_diameter_Micro_BC_rot_0959.txt');
Node961_MBC_rot=load('Along_diameter_Micro_BC_rot_0961.txt');
Node963_MBC_rot=load('Along_diameter_Micro_BC_rot_0963.txt');
Node965_MBC_rot=load('Along_diameter_Micro_BC_rot_0965.txt');
Node967_MBC_rot=load('Along_diameter_Micro_BC_rot_0967.txt');
Node969_MBC_rot=load('Along_diameter_Micro_BC_rot_0969.txt');
Node971_MBC_rot=load('Along_diameter_Micro_BC_rot_0971.txt');
Node973_MBC_rot=load('Along_diameter_Micro_BC_rot_0973.txt');
Node975_MBC_rot=load('Along_diameter_Micro_BC_rot_0975.txt');
Node977_MBC_rot=load('Along_diameter_Micro_BC_rot_0977.txt');
Node979_MBC_rot=load('Along_diameter_Micro_BC_rot_0979.txt');
Node981_MBC_rot=load('Along_diameter_Micro_BC_rot_0981.txt');
Node983_MBC_rot=load('Along_diameter_Micro_BC_rot_0983.txt');
Node985_MBC_rot=load('Along_diameter_Micro_BC_rot_0985.txt');
Node987_MBC_rot=load('Along_diameter_Micro_BC_rot_0987.txt');
Node989_MBC_rot=load('Along_diameter_Micro_BC_rot_0989.txt');
Node991_MBC_rot=load('Along_diameter_Micro_BC_rot_0991.txt');
Node993_MBC_rot=load('Along_diameter_Micro_BC_rot_0993.txt');
Node995_MBC_rot=load('Along_diameter_Micro_BC_rot_0995.txt');
Node997_MBC_rot=load('Along_diameter_Micro_BC_rot_0997.txt');
Node999_MBC_rot=load('Along_diameter_Micro_BC_rot_0999.txt');
Node1001_MBC_rot=load('Along_diameter_Micro_BC_rot_1001.txt');
Node1003_MBC_rot=load('Along_diameter_Micro_BC_rot_1003.txt');
Node1005_MBC_rot=load('Along_diameter_Micro_BC_rot_1005.txt');
Node1007_MBC_rot=load('Along_diameter_Micro_BC_rot_1007.txt');



Node938_MBC_str=load('Along_diameter_Micro_BC_str_0938.txt');
Node939_MBC_str=load('Along_diameter_Micro_BC_str_0939.txt');
Node941_MBC_str=load('Along_diameter_Micro_BC_str_0941.txt');
Node943_MBC_str=load('Along_diameter_Micro_BC_str_0943.txt');
Node945_MBC_str=load('Along_diameter_Micro_BC_str_0945.txt');
Node947_MBC_str=load('Along_diameter_Micro_BC_str_0947.txt');
Node949_MBC_str=load('Along_diameter_Micro_BC_str_0949.txt');
Node951_MBC_str=load('Along_diameter_Micro_BC_str_0951.txt');
Node953_MBC_str=load('Along_diameter_Micro_BC_str_0953.txt');
Node955_MBC_str=load('Along_diameter_Micro_BC_str_0955.txt');
Node957_MBC_str=load('Along_diameter_Micro_BC_str_0957.txt');
Node959_MBC_str=load('Along_diameter_Micro_BC_str_0959.txt');
Node961_MBC_str=load('Along_diameter_Micro_BC_str_0961.txt');
Node963_MBC_str=load('Along_diameter_Micro_BC_str_0963.txt');
Node965_MBC_str=load('Along_diameter_Micro_BC_str_0965.txt');
Node967_MBC_str=load('Along_diameter_Micro_BC_str_0967.txt');
Node969_MBC_str=load('Along_diameter_Micro_BC_str_0969.txt');
Node971_MBC_str=load('Along_diameter_Micro_BC_str_0971.txt');
Node973_MBC_str=load('Along_diameter_Micro_BC_str_0973.txt');
Node975_MBC_str=load('Along_diameter_Micro_BC_str_0975.txt');
Node977_MBC_str=load('Along_diameter_Micro_BC_str_0977.txt');
Node979_MBC_str=load('Along_diameter_Micro_BC_str_0979.txt');
Node981_MBC_str=load('Along_diameter_Micro_BC_str_0981.txt');
Node983_MBC_str=load('Along_diameter_Micro_BC_str_0983.txt');
Node985_MBC_str=load('Along_diameter_Micro_BC_str_0985.txt');
Node987_MBC_str=load('Along_diameter_Micro_BC_str_0987.txt');
Node989_MBC_str=load('Along_diameter_Micro_BC_str_0989.txt');
Node991_MBC_str=load('Along_diameter_Micro_BC_str_0991.txt');
Node993_MBC_str=load('Along_diameter_Micro_BC_str_0993.txt');
Node995_MBC_str=load('Along_diameter_Micro_BC_str_0995.txt');
Node997_MBC_str=load('Along_diameter_Micro_BC_str_0997.txt');
Node999_MBC_str=load('Along_diameter_Micro_BC_str_0999.txt');
Node1001_MBC_str=load('Along_diameter_Micro_BC_str_1001.txt');
Node1003_MBC_str=load('Along_diameter_Micro_BC_str_1003.txt');
Node1005_MBC_str=load('Along_diameter_Micro_BC_str_1005.txt');
Node1007_MBC_str=load('Along_diameter_Micro_BC_str_1007.txt');


Node938_NMBC_com=load('Along_diameter_NOMicro_BC_combine_0938.txt');
Node939_NMBC_com=load('Along_diameter_NOMicro_BC_combine_0939.txt');
Node941_NMBC_com=load('Along_diameter_NOMicro_BC_combine_0941.txt');
Node943_NMBC_com=load('Along_diameter_NOMicro_BC_combine_0943.txt');
Node945_NMBC_com=load('Along_diameter_NOMicro_BC_combine_0945.txt');
Node947_NMBC_com=load('Along_diameter_NOMicro_BC_combine_0947.txt');
Node949_NMBC_com=load('Along_diameter_NOMicro_BC_combine_0949.txt');
Node951_NMBC_com=load('Along_diameter_NOMicro_BC_combine_0951.txt');
Node953_NMBC_com=load('Along_diameter_NOMicro_BC_combine_0953.txt');
Node955_NMBC_com=load('Along_diameter_NOMicro_BC_combine_0955.txt');
Node957_NMBC_com=load('Along_diameter_NOMicro_BC_combine_0957.txt');
Node959_NMBC_com=load('Along_diameter_NOMicro_BC_combine_0959.txt');
Node961_NMBC_com=load('Along_diameter_NOMicro_BC_combine_0961.txt');
Node963_NMBC_com=load('Along_diameter_NOMicro_BC_combine_0963.txt');
Node965_NMBC_com=load('Along_diameter_NOMicro_BC_combine_0965.txt');
Node967_NMBC_com=load('Along_diameter_NOMicro_BC_combine_0967.txt');
Node969_NMBC_com=load('Along_diameter_NOMicro_BC_combine_0969.txt');
Node971_NMBC_com=load('Along_diameter_NOMicro_BC_combine_0971.txt');
Node973_NMBC_com=load('Along_diameter_NOMicro_BC_combine_0973.txt');
Node975_NMBC_com=load('Along_diameter_NOMicro_BC_combine_0975.txt');
Node977_NMBC_com=load('Along_diameter_NOMicro_BC_combine_0977.txt');
Node979_NMBC_com=load('Along_diameter_NOMicro_BC_combine_0979.txt');
Node981_NMBC_com=load('Along_diameter_NOMicro_BC_combine_0981.txt');
Node983_NMBC_com=load('Along_diameter_NOMicro_BC_combine_0983.txt');
Node985_NMBC_com=load('Along_diameter_NOMicro_BC_combine_0985.txt');
Node987_NMBC_com=load('Along_diameter_NOMicro_BC_combine_0987.txt');
Node989_NMBC_com=load('Along_diameter_NOMicro_BC_combine_0989.txt');
Node991_NMBC_com=load('Along_diameter_NOMicro_BC_combine_0991.txt');
Node993_NMBC_com=load('Along_diameter_NOMicro_BC_combine_0993.txt');
Node995_NMBC_com=load('Along_diameter_NOMicro_BC_combine_0995.txt');
Node997_NMBC_com=load('Along_diameter_NOMicro_BC_combine_0997.txt');
Node999_NMBC_com=load('Along_diameter_NOMicro_BC_combine_0999.txt');
Node1001_NMBC_com=load('Along_diameter_NOMicro_BC_combine_1001.txt');
Node1003_NMBC_com=load('Along_diameter_NOMicro_BC_combine_1003.txt');
Node1005_NMBC_com=load('Along_diameter_NOMicro_BC_combine_1005.txt');
Node1007_NMBC_com=load('Along_diameter_NOMicro_BC_combine_1007.txt');


Node938_NMBC_rot=load('Along_diameter_NOMicro_BC_rot_0938.txt');
Node939_NMBC_rot=load('Along_diameter_NOMicro_BC_rot_0939.txt');
Node941_NMBC_rot=load('Along_diameter_NOMicro_BC_rot_0941.txt');
Node943_NMBC_rot=load('Along_diameter_NOMicro_BC_rot_0943.txt');
Node945_NMBC_rot=load('Along_diameter_NOMicro_BC_rot_0945.txt');
Node947_NMBC_rot=load('Along_diameter_NOMicro_BC_rot_0947.txt');
Node949_NMBC_rot=load('Along_diameter_NOMicro_BC_rot_0949.txt');
Node951_NMBC_rot=load('Along_diameter_NOMicro_BC_rot_0951.txt');
Node953_NMBC_rot=load('Along_diameter_NOMicro_BC_rot_0953.txt');
Node955_NMBC_rot=load('Along_diameter_NOMicro_BC_rot_0955.txt');
Node957_NMBC_rot=load('Along_diameter_NOMicro_BC_rot_0957.txt');
Node959_NMBC_rot=load('Along_diameter_NOMicro_BC_rot_0959.txt');
Node961_NMBC_rot=load('Along_diameter_NOMicro_BC_rot_0961.txt');
Node963_NMBC_rot=load('Along_diameter_NOMicro_BC_rot_0963.txt');
Node965_NMBC_rot=load('Along_diameter_NOMicro_BC_rot_0965.txt');
Node967_NMBC_rot=load('Along_diameter_NOMicro_BC_rot_0967.txt');
Node969_NMBC_rot=load('Along_diameter_NOMicro_BC_rot_0969.txt');
Node971_NMBC_rot=load('Along_diameter_NOMicro_BC_rot_0971.txt');
Node973_NMBC_rot=load('Along_diameter_NOMicro_BC_rot_0973.txt');
Node975_NMBC_rot=load('Along_diameter_NOMicro_BC_rot_0975.txt');
Node977_NMBC_rot=load('Along_diameter_NOMicro_BC_rot_0977.txt');
Node979_NMBC_rot=load('Along_diameter_NOMicro_BC_rot_0979.txt');
Node981_NMBC_rot=load('Along_diameter_NOMicro_BC_rot_0981.txt');
Node983_NMBC_rot=load('Along_diameter_NOMicro_BC_rot_0983.txt');
Node985_NMBC_rot=load('Along_diameter_NOMicro_BC_rot_0985.txt');
Node987_NMBC_rot=load('Along_diameter_NOMicro_BC_rot_0987.txt');
Node989_NMBC_rot=load('Along_diameter_NOMicro_BC_rot_0989.txt');
Node991_NMBC_rot=load('Along_diameter_NOMicro_BC_rot_0991.txt');
Node993_NMBC_rot=load('Along_diameter_NOMicro_BC_rot_0993.txt');
Node995_NMBC_rot=load('Along_diameter_NOMicro_BC_rot_0995.txt');
Node997_NMBC_rot=load('Along_diameter_NOMicro_BC_rot_0997.txt');
Node999_NMBC_rot=load('Along_diameter_NOMicro_BC_rot_0999.txt');
Node1001_NMBC_rot=load('Along_diameter_NOMicro_BC_rot_1001.txt');
Node1003_NMBC_rot=load('Along_diameter_NOMicro_BC_rot_1003.txt');
Node1005_NMBC_rot=load('Along_diameter_NOMicro_BC_rot_1005.txt');
Node1007_NMBC_rot=load('Along_diameter_NOMicro_BC_rot_1007.txt');



Node938_NMBC_str=load('Along_diameter_NOMicro_BC_str_0938.txt');
Node939_NMBC_str=load('Along_diameter_NOMicro_BC_str_0939.txt');
Node941_NMBC_str=load('Along_diameter_NOMicro_BC_str_0941.txt');
Node943_NMBC_str=load('Along_diameter_NOMicro_BC_str_0943.txt');
Node945_NMBC_str=load('Along_diameter_NOMicro_BC_str_0945.txt');
Node947_NMBC_str=load('Along_diameter_NOMicro_BC_str_0947.txt');
Node949_NMBC_str=load('Along_diameter_NOMicro_BC_str_0949.txt');
Node951_NMBC_str=load('Along_diameter_NOMicro_BC_str_0951.txt');
Node953_NMBC_str=load('Along_diameter_NOMicro_BC_str_0953.txt');
Node955_NMBC_str=load('Along_diameter_NOMicro_BC_str_0955.txt');
Node957_NMBC_str=load('Along_diameter_NOMicro_BC_str_0957.txt');
Node959_NMBC_str=load('Along_diameter_NOMicro_BC_str_0959.txt');
Node961_NMBC_str=load('Along_diameter_NOMicro_BC_str_0961.txt');
Node963_NMBC_str=load('Along_diameter_NOMicro_BC_str_0963.txt');
Node965_NMBC_str=load('Along_diameter_NOMicro_BC_str_0965.txt');
Node967_NMBC_str=load('Along_diameter_NOMicro_BC_str_0967.txt');
Node969_NMBC_str=load('Along_diameter_NOMicro_BC_str_0969.txt');
Node971_NMBC_str=load('Along_diameter_NOMicro_BC_str_0971.txt');
Node973_NMBC_str=load('Along_diameter_NOMicro_BC_str_0973.txt');
Node975_NMBC_str=load('Along_diameter_NOMicro_BC_str_0975.txt');
Node977_NMBC_str=load('Along_diameter_NOMicro_BC_str_0977.txt');
Node979_NMBC_str=load('Along_diameter_NOMicro_BC_str_0979.txt');
Node981_NMBC_str=load('Along_diameter_NOMicro_BC_str_0981.txt');
Node983_NMBC_str=load('Along_diameter_NOMicro_BC_str_0983.txt');
Node985_NMBC_str=load('Along_diameter_NOMicro_BC_str_0985.txt');
Node987_NMBC_str=load('Along_diameter_NOMicro_BC_str_0987.txt');
Node989_NMBC_str=load('Along_diameter_NOMicro_BC_str_0989.txt');
Node991_NMBC_str=load('Along_diameter_NOMicro_BC_str_0991.txt');
Node993_NMBC_str=load('Along_diameter_NOMicro_BC_str_0993.txt');
Node995_NMBC_str=load('Along_diameter_NOMicro_BC_str_0995.txt');
Node997_NMBC_str=load('Along_diameter_NOMicro_BC_str_0997.txt');
Node999_NMBC_str=load('Along_diameter_NOMicro_BC_str_0999.txt');
Node1001_NMBC_str=load('Along_diameter_NOMicro_BC_str_1001.txt');
Node1003_NMBC_str=load('Along_diameter_NOMicro_BC_str_1003.txt');
Node1005_NMBC_str=load('Along_diameter_NOMicro_BC_str_1005.txt');
Node1007_NMBC_str=load('Along_diameter_NOMicro_BC_str_1007.txt');


Applied_load_141=load('Applied_load_0141.txt');
Applied_load_143=load('Applied_load_0143.txt');
Applied_load_215=load('Applied_load_0215.txt');
Applied_load_287=load('Applied_load_0287.txt');
Applied_load_359=load('Applied_load_0359.txt');
Applied_load_431=load('Applied_load_0431.txt');
Applied_load_503=load('Applied_load_0503.txt');
Applied_load_575=load('Applied_load_0575.txt');
Applied_load_647=load('Applied_load_0647.txt');
Applied_load_719=load('Applied_load_0719.txt');
Applied_load_791=load('Applied_load_0791.txt');
Applied_load_863=load('Applied_load_0863.txt');
Applied_load_935=load('Applied_load_0935.txt');
Applied_load_1007=load('Applied_load_1007.txt');


Node2_AH_cls=load('Around_hole_classic_0002.txt');
Node6_AH_cls=load('Around_hole_classic_0006.txt');
Node146_AH_cls=load('Around_hole_classic_0146.txt');
Node218_AH_cls=load('Around_hole_classic_0218.txt');
Node290_AH_cls=load('Around_hole_classic_0290.txt');
Node362_AH_cls=load('Around_hole_classic_0362.txt');
Node434_AH_cls=load('Around_hole_classic_0434.txt');
Node506_AH_cls=load('Around_hole_classic_0506.txt');
Node578_AH_cls=load('Around_hole_classic_0578.txt');
Node650_AH_cls=load('Around_hole_classic_0650.txt');
Node722_AH_cls=load('Around_hole_classic_0722.txt');
Node794_AH_cls=load('Around_hole_classic_0794.txt');
Node866_AH_cls=load('Around_hole_classic_0866.txt');
Node938_AH_cls=load('Around_hole_classic_0938.txt');
Node1012_AH_cls=load('Around_hole_classic_1012.txt');
Node1016_AH_cls=load('Around_hole_classic_1016.txt');
Node1156_AH_cls=load('Around_hole_classic_1156.txt');
Node1228_AH_cls=load('Around_hole_classic_1228.txt');
Node1300_AH_cls=load('Around_hole_classic_1300.txt');
Node1372_AH_cls=load('Around_hole_classic_1372.txt');
Node1444_AH_cls=load('Around_hole_classic_1444.txt');
Node1516_AH_cls=load('Around_hole_classic_1516.txt');
Node1588_AH_cls=load('Around_hole_classic_1588.txt');
Node1660_AH_cls=load('Around_hole_classic_1660.txt');
Node1732_AH_cls=load('Around_hole_classic_1732.txt');
Node1804_AH_cls=load('Around_hole_classic_1804.txt');
Node1876_AH_cls=load('Around_hole_classic_1876.txt');


Node2_AH_MBC_com=load('Around_hole_Micro_BC_combine_0002.txt');
Node6_AH_MBC_com=load('Around_hole_Micro_BC_combine_0006.txt');
Node146_AH_MBC_com=load('Around_hole_Micro_BC_combine_0146.txt');
Node218_AH_MBC_com=load('Around_hole_Micro_BC_combine_0218.txt');
Node290_AH_MBC_com=load('Around_hole_Micro_BC_combine_0290.txt');
Node362_AH_MBC_com=load('Around_hole_Micro_BC_combine_0362.txt');
Node434_AH_MBC_com=load('Around_hole_Micro_BC_combine_0434.txt');
Node506_AH_MBC_com=load('Around_hole_Micro_BC_combine_0506.txt');
Node578_AH_MBC_com=load('Around_hole_Micro_BC_combine_0578.txt');
Node650_AH_MBC_com=load('Around_hole_Micro_BC_combine_0650.txt');
Node722_AH_MBC_com=load('Around_hole_Micro_BC_combine_0722.txt');
Node794_AH_MBC_com=load('Around_hole_Micro_BC_combine_0794.txt');
Node866_AH_MBC_com=load('Around_hole_Micro_BC_combine_0866.txt');
Node938_AH_MBC_com=load('Around_hole_Micro_BC_combine_0938.txt');
Node1012_AH_MBC_com=load('Around_hole_Micro_BC_combine_1012.txt');
Node1016_AH_MBC_com=load('Around_hole_Micro_BC_combine_1016.txt');
Node1156_AH_MBC_com=load('Around_hole_Micro_BC_combine_1156.txt');
Node1228_AH_MBC_com=load('Around_hole_Micro_BC_combine_1228.txt');
Node1300_AH_MBC_com=load('Around_hole_Micro_BC_combine_1300.txt');
Node1372_AH_MBC_com=load('Around_hole_Micro_BC_combine_1372.txt');
Node1444_AH_MBC_com=load('Around_hole_Micro_BC_combine_1444.txt');
Node1516_AH_MBC_com=load('Around_hole_Micro_BC_combine_1516.txt');
Node1588_AH_MBC_com=load('Around_hole_Micro_BC_combine_1588.txt');
Node1660_AH_MBC_com=load('Around_hole_Micro_BC_combine_1660.txt');
Node1732_AH_MBC_com=load('Around_hole_Micro_BC_combine_1732.txt');
Node1804_AH_MBC_com=load('Around_hole_Micro_BC_combine_1804.txt');
Node1876_AH_MBC_com=load('Around_hole_Micro_BC_combine_1876.txt');



Node2_AH_MBC_rot=load('Around_hole_Micro_BC_rot_0002.txt');
Node6_AH_MBC_rot=load('Around_hole_Micro_BC_rot_0006.txt');
Node146_AH_MBC_rot=load('Around_hole_Micro_BC_rot_0146.txt');
Node218_AH_MBC_rot=load('Around_hole_Micro_BC_rot_0218.txt');
Node290_AH_MBC_rot=load('Around_hole_Micro_BC_rot_0290.txt');
Node362_AH_MBC_rot=load('Around_hole_Micro_BC_rot_0362.txt');
Node434_AH_MBC_rot=load('Around_hole_Micro_BC_rot_0434.txt');
Node506_AH_MBC_rot=load('Around_hole_Micro_BC_rot_0506.txt');
Node578_AH_MBC_rot=load('Around_hole_Micro_BC_rot_0578.txt');
Node650_AH_MBC_rot=load('Around_hole_Micro_BC_rot_0650.txt');
Node722_AH_MBC_rot=load('Around_hole_Micro_BC_rot_0722.txt');
Node794_AH_MBC_rot=load('Around_hole_Micro_BC_rot_0794.txt');
Node866_AH_MBC_rot=load('Around_hole_Micro_BC_rot_0866.txt');
Node938_AH_MBC_rot=load('Around_hole_Micro_BC_rot_0938.txt');
Node1012_AH_MBC_rot=load('Around_hole_Micro_BC_rot_1012.txt');
Node1016_AH_MBC_rot=load('Around_hole_Micro_BC_rot_1016.txt');
Node1156_AH_MBC_rot=load('Around_hole_Micro_BC_rot_1156.txt');
Node1228_AH_MBC_rot=load('Around_hole_Micro_BC_rot_1228.txt');
Node1300_AH_MBC_rot=load('Around_hole_Micro_BC_rot_1300.txt');
Node1372_AH_MBC_rot=load('Around_hole_Micro_BC_rot_1372.txt');
Node1444_AH_MBC_rot=load('Around_hole_Micro_BC_rot_1444.txt');
Node1516_AH_MBC_rot=load('Around_hole_Micro_BC_rot_1516.txt');
Node1588_AH_MBC_rot=load('Around_hole_Micro_BC_rot_1588.txt');
Node1660_AH_MBC_rot=load('Around_hole_Micro_BC_rot_1660.txt');
Node1732_AH_MBC_rot=load('Around_hole_Micro_BC_rot_1732.txt');
Node1804_AH_MBC_rot=load('Around_hole_Micro_BC_rot_1804.txt');
Node1876_AH_MBC_rot=load('Around_hole_Micro_BC_rot_1876.txt');


Node2_AH_MBC_str=load('Around_hole_Micro_BC_str_0002.txt');
Node6_AH_MBC_str=load('Around_hole_Micro_BC_str_0006.txt');
Node146_AH_MBC_str=load('Around_hole_Micro_BC_str_0146.txt');
Node218_AH_MBC_str=load('Around_hole_Micro_BC_str_0218.txt');
Node290_AH_MBC_str=load('Around_hole_Micro_BC_str_0290.txt');
Node362_AH_MBC_str=load('Around_hole_Micro_BC_str_0362.txt');
Node434_AH_MBC_str=load('Around_hole_Micro_BC_str_0434.txt');
Node506_AH_MBC_str=load('Around_hole_Micro_BC_str_0506.txt');
Node578_AH_MBC_str=load('Around_hole_Micro_BC_str_0578.txt');
Node650_AH_MBC_str=load('Around_hole_Micro_BC_str_0650.txt');
Node722_AH_MBC_str=load('Around_hole_Micro_BC_str_0722.txt');
Node794_AH_MBC_str=load('Around_hole_Micro_BC_str_0794.txt');
Node866_AH_MBC_str=load('Around_hole_Micro_BC_str_0866.txt');
Node938_AH_MBC_str=load('Around_hole_Micro_BC_str_0938.txt');
Node1012_AH_MBC_str=load('Around_hole_Micro_BC_str_1012.txt');
Node1016_AH_MBC_str=load('Around_hole_Micro_BC_str_1016.txt');
Node1156_AH_MBC_str=load('Around_hole_Micro_BC_str_1156.txt');
Node1228_AH_MBC_str=load('Around_hole_Micro_BC_str_1228.txt');
Node1300_AH_MBC_str=load('Around_hole_Micro_BC_str_1300.txt');
Node1372_AH_MBC_str=load('Around_hole_Micro_BC_str_1372.txt');
Node1444_AH_MBC_str=load('Around_hole_Micro_BC_str_1444.txt');
Node1516_AH_MBC_str=load('Around_hole_Micro_BC_str_1516.txt');
Node1588_AH_MBC_str=load('Around_hole_Micro_BC_str_1588.txt');
Node1660_AH_MBC_str=load('Around_hole_Micro_BC_str_1660.txt');
Node1732_AH_MBC_str=load('Around_hole_Micro_BC_str_1732.txt');
Node1804_AH_MBC_str=load('Around_hole_Micro_BC_str_1804.txt');
Node1876_AH_MBC_str=load('Around_hole_Micro_BC_str_1876.txt');



Node2_AH_NMBC_com=load('Around_hole_NOMicro_BC_combine_0002.txt');
Node6_AH_NMBC_com=load('Around_hole_NOMicro_BC_combine_0006.txt');
Node146_AH_NMBC_com=load('Around_hole_NOMicro_BC_combine_0146.txt');
Node218_AH_NMBC_com=load('Around_hole_NOMicro_BC_combine_0218.txt');
Node290_AH_NMBC_com=load('Around_hole_NOMicro_BC_combine_0290.txt');
Node362_AH_NMBC_com=load('Around_hole_NOMicro_BC_combine_0362.txt');
Node434_AH_NMBC_com=load('Around_hole_NOMicro_BC_combine_0434.txt');
Node506_AH_NMBC_com=load('Around_hole_NOMicro_BC_combine_0506.txt');
Node578_AH_NMBC_com=load('Around_hole_NOMicro_BC_combine_0578.txt');
Node650_AH_NMBC_com=load('Around_hole_NOMicro_BC_combine_0650.txt');
Node722_AH_NMBC_com=load('Around_hole_NOMicro_BC_combine_0722.txt');
Node794_AH_NMBC_com=load('Around_hole_NOMicro_BC_combine_0794.txt');
Node866_AH_NMBC_com=load('Around_hole_NOMicro_BC_combine_0866.txt');
Node938_AH_NMBC_com=load('Around_hole_NOMicro_BC_combine_0938.txt');
Node1012_AH_NMBC_com=load('Around_hole_NOMicro_BC_combine_1012.txt');
Node1016_AH_NMBC_com=load('Around_hole_NOMicro_BC_combine_1016.txt');
Node1156_AH_NMBC_com=load('Around_hole_NOMicro_BC_combine_1156.txt');
Node1228_AH_NMBC_com=load('Around_hole_NOMicro_BC_combine_1228.txt');
Node1300_AH_NMBC_com=load('Around_hole_NOMicro_BC_combine_1300.txt');
Node1372_AH_NMBC_com=load('Around_hole_NOMicro_BC_combine_1372.txt');
Node1444_AH_NMBC_com=load('Around_hole_NOMicro_BC_combine_1444.txt');
Node1516_AH_NMBC_com=load('Around_hole_NOMicro_BC_combine_1516.txt');
Node1588_AH_NMBC_com=load('Around_hole_NOMicro_BC_combine_1588.txt');
Node1660_AH_NMBC_com=load('Around_hole_NOMicro_BC_combine_1660.txt');
Node1732_AH_NMBC_com=load('Around_hole_NOMicro_BC_combine_1732.txt');
Node1804_AH_NMBC_com=load('Around_hole_NOMicro_BC_combine_1804.txt');
Node1876_AH_NMBC_com=load('Around_hole_NOMicro_BC_combine_1876.txt');



Node2_AH_NMBC_rot=load('Around_hole_NOMicro_BC_rot_0002.txt');
Node6_AH_NMBC_rot=load('Around_hole_NOMicro_BC_rot_0006.txt');
Node146_AH_NMBC_rot=load('Around_hole_NOMicro_BC_rot_0146.txt');
Node218_AH_NMBC_rot=load('Around_hole_NOMicro_BC_rot_0218.txt');
Node290_AH_NMBC_rot=load('Around_hole_NOMicro_BC_rot_0290.txt');
Node362_AH_NMBC_rot=load('Around_hole_NOMicro_BC_rot_0362.txt');
Node434_AH_NMBC_rot=load('Around_hole_NOMicro_BC_rot_0434.txt');
Node506_AH_NMBC_rot=load('Around_hole_NOMicro_BC_rot_0506.txt');
Node578_AH_NMBC_rot=load('Around_hole_NOMicro_BC_rot_0578.txt');
Node650_AH_NMBC_rot=load('Around_hole_NOMicro_BC_rot_0650.txt');
Node722_AH_NMBC_rot=load('Around_hole_NOMicro_BC_rot_0722.txt');
Node794_AH_NMBC_rot=load('Around_hole_NOMicro_BC_rot_0794.txt');
Node866_AH_NMBC_rot=load('Around_hole_NOMicro_BC_rot_0866.txt');
Node938_AH_NMBC_rot=load('Around_hole_NOMicro_BC_rot_0938.txt');
Node1012_AH_NMBC_rot=load('Around_hole_NOMicro_BC_rot_1012.txt');
Node1016_AH_NMBC_rot=load('Around_hole_NOMicro_BC_rot_1016.txt');
Node1156_AH_NMBC_rot=load('Around_hole_NOMicro_BC_rot_1156.txt');
Node1228_AH_NMBC_rot=load('Around_hole_NOMicro_BC_rot_1228.txt');
Node1300_AH_NMBC_rot=load('Around_hole_NOMicro_BC_rot_1300.txt');
Node1372_AH_NMBC_rot=load('Around_hole_NOMicro_BC_rot_1372.txt');
Node1444_AH_NMBC_rot=load('Around_hole_NOMicro_BC_rot_1444.txt');
Node1516_AH_NMBC_rot=load('Around_hole_NOMicro_BC_rot_1516.txt');
Node1588_AH_NMBC_rot=load('Around_hole_NOMicro_BC_rot_1588.txt');
Node1660_AH_NMBC_rot=load('Around_hole_NOMicro_BC_rot_1660.txt');
Node1732_AH_NMBC_rot=load('Around_hole_NOMicro_BC_rot_1732.txt');
Node1804_AH_NMBC_rot=load('Around_hole_NOMicro_BC_rot_1804.txt');
Node1876_AH_NMBC_rot=load('Around_hole_NOMicro_BC_rot_1876.txt');


Node2_AH_NMBC_str=load('Around_hole_NOMicro_BC_str_0002.txt');
Node6_AH_NMBC_str=load('Around_hole_NOMicro_BC_str_0006.txt');
Node146_AH_NMBC_str=load('Around_hole_NOMicro_BC_str_0146.txt');
Node218_AH_NMBC_str=load('Around_hole_NOMicro_BC_str_0218.txt');
Node290_AH_NMBC_str=load('Around_hole_NOMicro_BC_str_0290.txt');
Node362_AH_NMBC_str=load('Around_hole_NOMicro_BC_str_0362.txt');
Node434_AH_NMBC_str=load('Around_hole_NOMicro_BC_str_0434.txt');
Node506_AH_NMBC_str=load('Around_hole_NOMicro_BC_str_0506.txt');
Node578_AH_NMBC_str=load('Around_hole_NOMicro_BC_str_0578.txt');
Node650_AH_NMBC_str=load('Around_hole_NOMicro_BC_str_0650.txt');
Node722_AH_NMBC_str=load('Around_hole_NOMicro_BC_str_0722.txt');
Node794_AH_NMBC_str=load('Around_hole_NOMicro_BC_str_0794.txt');
Node866_AH_NMBC_str=load('Around_hole_NOMicro_BC_str_0866.txt');
Node938_AH_NMBC_str=load('Around_hole_NOMicro_BC_str_0938.txt');
Node1012_AH_NMBC_str=load('Around_hole_NOMicro_BC_str_1012.txt');
Node1016_AH_NMBC_str=load('Around_hole_NOMicro_BC_str_1016.txt');
Node1156_AH_NMBC_str=load('Around_hole_NOMicro_BC_str_1156.txt');
Node1228_AH_NMBC_str=load('Around_hole_NOMicro_BC_str_1228.txt');
Node1300_AH_NMBC_str=load('Around_hole_NOMicro_BC_str_1300.txt');
Node1372_AH_NMBC_str=load('Around_hole_NOMicro_BC_str_1372.txt');
Node1444_AH_NMBC_str=load('Around_hole_NOMicro_BC_str_1444.txt');
Node1516_AH_NMBC_str=load('Around_hole_NOMicro_BC_str_1516.txt');
Node1588_AH_NMBC_str=load('Around_hole_NOMicro_BC_str_1588.txt');
Node1660_AH_NMBC_str=load('Around_hole_NOMicro_BC_str_1660.txt');
Node1732_AH_NMBC_str=load('Around_hole_NOMicro_BC_str_1732.txt');
Node1804_AH_NMBC_str=load('Around_hole_NOMicro_BC_str_1804.txt');
Node1876_AH_NMBC_str=load('Around_hole_NOMicro_BC_str_1876.txt');


Node1011_SI_cls=load('Stress_intensity_classic_1011.txt');
Node1012_SI_cls=load('Stress_intensity_classic_1012.txt');
Node1018_SI_cls=load('Stress_intensity_classic_1018.txt');
Node1022_SI_cls=load('Stress_intensity_classic_1022.txt');
Node1026_SI_cls=load('Stress_intensity_classic_1026.txt');
Node1030_SI_cls=load('Stress_intensity_classic_1030.txt');
Node1034_SI_cls=load('Stress_intensity_classic_1034.txt');
Node1038_SI_cls=load('Stress_intensity_classic_1038.txt');
Node1042_SI_cls=load('Stress_intensity_classic_1042.txt');
Node1046_SI_cls=load('Stress_intensity_classic_1046.txt');
Node1050_SI_cls=load('Stress_intensity_classic_1050.txt');
Node1054_SI_cls=load('Stress_intensity_classic_1054.txt');
Node1058_SI_cls=load('Stress_intensity_classic_1058.txt');
Node1062_SI_cls=load('Stress_intensity_classic_1062.txt');
Node1066_SI_cls=load('Stress_intensity_classic_1066.txt');
Node1070_SI_cls=load('Stress_intensity_classic_1070.txt');
Node1074_SI_cls=load('Stress_intensity_classic_1074.txt');
Node1078_SI_cls=load('Stress_intensity_classic_1078.txt');
Node1082_SI_cls=load('Stress_intensity_classic_1082.txt');
Node1086_SI_cls=load('Stress_intensity_classic_1086.txt');
Node1090_SI_cls=load('Stress_intensity_classic_1090.txt');
Node1094_SI_cls=load('Stress_intensity_classic_1094.txt');
Node1098_SI_cls=load('Stress_intensity_classic_1098.txt');
Node1102_SI_cls=load('Stress_intensity_classic_1102.txt');
Node1106_SI_cls=load('Stress_intensity_classic_1106.txt');
Node1110_SI_cls=load('Stress_intensity_classic_1110.txt');
Node1114_SI_cls=load('Stress_intensity_classic_1114.txt');
Node1118_SI_cls=load('Stress_intensity_classic_1118.txt');
Node1122_SI_cls=load('Stress_intensity_classic_1122.txt');
Node1126_SI_cls=load('Stress_intensity_classic_1126.txt');
Node1130_SI_cls=load('Stress_intensity_classic_1130.txt');
Node1134_SI_cls=load('Stress_intensity_classic_1134.txt');
Node1138_SI_cls=load('Stress_intensity_classic_1138.txt');
Node1142_SI_cls=load('Stress_intensity_classic_1142.txt');
Node1146_SI_cls=load('Stress_intensity_classic_1146.txt');
Node1150_SI_cls=load('Stress_intensity_classic_1150.txt');


Node1011_SI_MBC_com=load('Stress_intensity_Micro_BC_combine_1011.txt');
Node1012_SI_MBC_com=load('Stress_intensity_Micro_BC_combine_1012.txt');
Node1018_SI_MBC_com=load('Stress_intensity_Micro_BC_combine_1018.txt');
Node1022_SI_MBC_com=load('Stress_intensity_Micro_BC_combine_1022.txt');
Node1026_SI_MBC_com=load('Stress_intensity_Micro_BC_combine_1026.txt');
Node1030_SI_MBC_com=load('Stress_intensity_Micro_BC_combine_1030.txt');
Node1034_SI_MBC_com=load('Stress_intensity_Micro_BC_combine_1034.txt');
Node1038_SI_MBC_com=load('Stress_intensity_Micro_BC_combine_1038.txt');
Node1042_SI_MBC_com=load('Stress_intensity_Micro_BC_combine_1042.txt');
Node1046_SI_MBC_com=load('Stress_intensity_Micro_BC_combine_1046.txt');
Node1050_SI_MBC_com=load('Stress_intensity_Micro_BC_combine_1050.txt');
Node1054_SI_MBC_com=load('Stress_intensity_Micro_BC_combine_1054.txt');
Node1058_SI_MBC_com=load('Stress_intensity_Micro_BC_combine_1058.txt');
Node1062_SI_MBC_com=load('Stress_intensity_Micro_BC_combine_1062.txt');
Node1066_SI_MBC_com=load('Stress_intensity_Micro_BC_combine_1066.txt');
Node1070_SI_MBC_com=load('Stress_intensity_Micro_BC_combine_1070.txt');
Node1074_SI_MBC_com=load('Stress_intensity_Micro_BC_combine_1074.txt');
Node1078_SI_MBC_com=load('Stress_intensity_Micro_BC_combine_1078.txt');
Node1082_SI_MBC_com=load('Stress_intensity_Micro_BC_combine_1082.txt');
Node1086_SI_MBC_com=load('Stress_intensity_Micro_BC_combine_1086.txt');
Node1090_SI_MBC_com=load('Stress_intensity_Micro_BC_combine_1090.txt');
Node1094_SI_MBC_com=load('Stress_intensity_Micro_BC_combine_1094.txt');
Node1098_SI_MBC_com=load('Stress_intensity_Micro_BC_combine_1098.txt');
Node1102_SI_MBC_com=load('Stress_intensity_Micro_BC_combine_1102.txt');
Node1106_SI_MBC_com=load('Stress_intensity_Micro_BC_combine_1106.txt');
Node1110_SI_MBC_com=load('Stress_intensity_Micro_BC_combine_1110.txt');
Node1114_SI_MBC_com=load('Stress_intensity_Micro_BC_combine_1114.txt');
Node1118_SI_MBC_com=load('Stress_intensity_Micro_BC_combine_1118.txt');
Node1122_SI_MBC_com=load('Stress_intensity_Micro_BC_combine_1122.txt');
Node1126_SI_MBC_com=load('Stress_intensity_Micro_BC_combine_1126.txt');
Node1130_SI_MBC_com=load('Stress_intensity_Micro_BC_combine_1130.txt');
Node1134_SI_MBC_com=load('Stress_intensity_Micro_BC_combine_1134.txt');
Node1138_SI_MBC_com=load('Stress_intensity_Micro_BC_combine_1138.txt');
Node1142_SI_MBC_com=load('Stress_intensity_Micro_BC_combine_1142.txt');
Node1146_SI_MBC_com=load('Stress_intensity_Micro_BC_combine_1146.txt');
Node1150_SI_MBC_com=load('Stress_intensity_Micro_BC_combine_1150.txt');


Node1011_SI_MBC_rot=load('Stress_intensity_Micro_BC_rot_1011.txt');
Node1012_SI_MBC_rot=load('Stress_intensity_Micro_BC_rot_1012.txt');
Node1018_SI_MBC_rot=load('Stress_intensity_Micro_BC_rot_1018.txt');
Node1022_SI_MBC_rot=load('Stress_intensity_Micro_BC_rot_1022.txt');
Node1026_SI_MBC_rot=load('Stress_intensity_Micro_BC_rot_1026.txt');
Node1030_SI_MBC_rot=load('Stress_intensity_Micro_BC_rot_1030.txt');
Node1034_SI_MBC_rot=load('Stress_intensity_Micro_BC_rot_1034.txt');
Node1038_SI_MBC_rot=load('Stress_intensity_Micro_BC_rot_1038.txt');
Node1042_SI_MBC_rot=load('Stress_intensity_Micro_BC_rot_1042.txt');
Node1046_SI_MBC_rot=load('Stress_intensity_Micro_BC_rot_1046.txt');
Node1050_SI_MBC_rot=load('Stress_intensity_Micro_BC_rot_1050.txt');
Node1054_SI_MBC_rot=load('Stress_intensity_Micro_BC_rot_1054.txt');
Node1058_SI_MBC_rot=load('Stress_intensity_Micro_BC_rot_1058.txt');
Node1062_SI_MBC_rot=load('Stress_intensity_Micro_BC_rot_1062.txt');
Node1066_SI_MBC_rot=load('Stress_intensity_Micro_BC_rot_1066.txt');
Node1070_SI_MBC_rot=load('Stress_intensity_Micro_BC_rot_1070.txt');
Node1074_SI_MBC_rot=load('Stress_intensity_Micro_BC_rot_1074.txt');
Node1078_SI_MBC_rot=load('Stress_intensity_Micro_BC_rot_1078.txt');
Node1082_SI_MBC_rot=load('Stress_intensity_Micro_BC_rot_1082.txt');
Node1086_SI_MBC_rot=load('Stress_intensity_Micro_BC_rot_1086.txt');
Node1090_SI_MBC_rot=load('Stress_intensity_Micro_BC_rot_1090.txt');
Node1094_SI_MBC_rot=load('Stress_intensity_Micro_BC_rot_1094.txt');
Node1098_SI_MBC_rot=load('Stress_intensity_Micro_BC_rot_1098.txt');
Node1102_SI_MBC_rot=load('Stress_intensity_Micro_BC_rot_1102.txt');
Node1106_SI_MBC_rot=load('Stress_intensity_Micro_BC_rot_1106.txt');
Node1110_SI_MBC_rot=load('Stress_intensity_Micro_BC_rot_1110.txt');
Node1114_SI_MBC_rot=load('Stress_intensity_Micro_BC_rot_1114.txt');
Node1118_SI_MBC_rot=load('Stress_intensity_Micro_BC_rot_1118.txt');
Node1122_SI_MBC_rot=load('Stress_intensity_Micro_BC_rot_1122.txt');
Node1126_SI_MBC_rot=load('Stress_intensity_Micro_BC_rot_1126.txt');
Node1130_SI_MBC_rot=load('Stress_intensity_Micro_BC_rot_1130.txt');
Node1134_SI_MBC_rot=load('Stress_intensity_Micro_BC_rot_1134.txt');
Node1138_SI_MBC_rot=load('Stress_intensity_Micro_BC_rot_1138.txt');
Node1142_SI_MBC_rot=load('Stress_intensity_Micro_BC_rot_1142.txt');
Node1146_SI_MBC_rot=load('Stress_intensity_Micro_BC_rot_1146.txt');
Node1150_SI_MBC_rot=load('Stress_intensity_Micro_BC_rot_1150.txt');


Node1011_SI_MBC_str=load('Stress_intensity_Micro_BC_str_1011.txt');
Node1012_SI_MBC_str=load('Stress_intensity_Micro_BC_str_1012.txt');
Node1018_SI_MBC_str=load('Stress_intensity_Micro_BC_str_1018.txt');
Node1022_SI_MBC_str=load('Stress_intensity_Micro_BC_str_1022.txt');
Node1026_SI_MBC_str=load('Stress_intensity_Micro_BC_str_1026.txt');
Node1030_SI_MBC_str=load('Stress_intensity_Micro_BC_str_1030.txt');
Node1034_SI_MBC_str=load('Stress_intensity_Micro_BC_str_1034.txt');
Node1038_SI_MBC_str=load('Stress_intensity_Micro_BC_str_1038.txt');
Node1042_SI_MBC_str=load('Stress_intensity_Micro_BC_str_1042.txt');
Node1046_SI_MBC_str=load('Stress_intensity_Micro_BC_str_1046.txt');
Node1050_SI_MBC_str=load('Stress_intensity_Micro_BC_str_1050.txt');
Node1054_SI_MBC_str=load('Stress_intensity_Micro_BC_str_1054.txt');
Node1058_SI_MBC_str=load('Stress_intensity_Micro_BC_str_1058.txt');
Node1062_SI_MBC_str=load('Stress_intensity_Micro_BC_str_1062.txt');
Node1066_SI_MBC_str=load('Stress_intensity_Micro_BC_str_1066.txt');
Node1070_SI_MBC_str=load('Stress_intensity_Micro_BC_str_1070.txt');
Node1074_SI_MBC_str=load('Stress_intensity_Micro_BC_str_1074.txt');
Node1078_SI_MBC_str=load('Stress_intensity_Micro_BC_str_1078.txt');
Node1082_SI_MBC_str=load('Stress_intensity_Micro_BC_str_1082.txt');
Node1086_SI_MBC_str=load('Stress_intensity_Micro_BC_str_1086.txt');
Node1090_SI_MBC_str=load('Stress_intensity_Micro_BC_str_1090.txt');
Node1094_SI_MBC_str=load('Stress_intensity_Micro_BC_str_1094.txt');
Node1098_SI_MBC_str=load('Stress_intensity_Micro_BC_str_1098.txt');
Node1102_SI_MBC_str=load('Stress_intensity_Micro_BC_str_1102.txt');
Node1106_SI_MBC_str=load('Stress_intensity_Micro_BC_str_1106.txt');
Node1110_SI_MBC_str=load('Stress_intensity_Micro_BC_str_1110.txt');
Node1114_SI_MBC_str=load('Stress_intensity_Micro_BC_str_1114.txt');
Node1118_SI_MBC_str=load('Stress_intensity_Micro_BC_str_1118.txt');
Node1122_SI_MBC_str=load('Stress_intensity_Micro_BC_str_1122.txt');
Node1126_SI_MBC_str=load('Stress_intensity_Micro_BC_str_1126.txt');
Node1130_SI_MBC_str=load('Stress_intensity_Micro_BC_str_1130.txt');
Node1134_SI_MBC_str=load('Stress_intensity_Micro_BC_str_1134.txt');
Node1138_SI_MBC_str=load('Stress_intensity_Micro_BC_str_1138.txt');
Node1142_SI_MBC_str=load('Stress_intensity_Micro_BC_str_1142.txt');
Node1146_SI_MBC_str=load('Stress_intensity_Micro_BC_str_1146.txt');
Node1150_SI_MBC_str=load('Stress_intensity_Micro_BC_str_1150.txt');



Node1011_SI_NMBC_com=load('Stress_intensity_NOMicro_BC_combine_1011.txt');
Node1012_SI_NMBC_com=load('Stress_intensity_NOMicro_BC_combine_1012.txt');
Node1018_SI_NMBC_com=load('Stress_intensity_NOMicro_BC_combine_1018.txt');
Node1022_SI_NMBC_com=load('Stress_intensity_NOMicro_BC_combine_1022.txt');
Node1026_SI_NMBC_com=load('Stress_intensity_NOMicro_BC_combine_1026.txt');
Node1030_SI_NMBC_com=load('Stress_intensity_NOMicro_BC_combine_1030.txt');
Node1034_SI_NMBC_com=load('Stress_intensity_NOMicro_BC_combine_1034.txt');
Node1038_SI_NMBC_com=load('Stress_intensity_NOMicro_BC_combine_1038.txt');
Node1042_SI_NMBC_com=load('Stress_intensity_NOMicro_BC_combine_1042.txt');
Node1046_SI_NMBC_com=load('Stress_intensity_NOMicro_BC_combine_1046.txt');
Node1050_SI_NMBC_com=load('Stress_intensity_NOMicro_BC_combine_1050.txt');
Node1054_SI_NMBC_com=load('Stress_intensity_NOMicro_BC_combine_1054.txt');
Node1058_SI_NMBC_com=load('Stress_intensity_NOMicro_BC_combine_1058.txt');
Node1062_SI_NMBC_com=load('Stress_intensity_NOMicro_BC_combine_1062.txt');
Node1066_SI_NMBC_com=load('Stress_intensity_NOMicro_BC_combine_1066.txt');
Node1070_SI_NMBC_com=load('Stress_intensity_NOMicro_BC_combine_1070.txt');
Node1074_SI_NMBC_com=load('Stress_intensity_NOMicro_BC_combine_1074.txt');
Node1078_SI_NMBC_com=load('Stress_intensity_NOMicro_BC_combine_1078.txt');
Node1082_SI_NMBC_com=load('Stress_intensity_NOMicro_BC_combine_1082.txt');
Node1086_SI_NMBC_com=load('Stress_intensity_NOMicro_BC_combine_1086.txt');
Node1090_SI_NMBC_com=load('Stress_intensity_NOMicro_BC_combine_1090.txt');
Node1094_SI_NMBC_com=load('Stress_intensity_NOMicro_BC_combine_1094.txt');
Node1098_SI_NMBC_com=load('Stress_intensity_NOMicro_BC_combine_1098.txt');
Node1102_SI_NMBC_com=load('Stress_intensity_NOMicro_BC_combine_1102.txt');
Node1106_SI_NMBC_com=load('Stress_intensity_NOMicro_BC_combine_1106.txt');
Node1110_SI_NMBC_com=load('Stress_intensity_NOMicro_BC_combine_1110.txt');
Node1114_SI_NMBC_com=load('Stress_intensity_NOMicro_BC_combine_1114.txt');
Node1118_SI_NMBC_com=load('Stress_intensity_NOMicro_BC_combine_1118.txt');
Node1122_SI_NMBC_com=load('Stress_intensity_NOMicro_BC_combine_1122.txt');
Node1126_SI_NMBC_com=load('Stress_intensity_NOMicro_BC_combine_1126.txt');
Node1130_SI_NMBC_com=load('Stress_intensity_NOMicro_BC_combine_1130.txt');
Node1134_SI_NMBC_com=load('Stress_intensity_NOMicro_BC_combine_1134.txt');
Node1138_SI_NMBC_com=load('Stress_intensity_NOMicro_BC_combine_1138.txt');
Node1142_SI_NMBC_com=load('Stress_intensity_NOMicro_BC_combine_1142.txt');
Node1146_SI_NMBC_com=load('Stress_intensity_NOMicro_BC_combine_1146.txt');
Node1150_SI_NMBC_com=load('Stress_intensity_NOMicro_BC_combine_1150.txt');


Node1011_SI_NMBC_rot=load('Stress_intensity_NOMicro_BC_rot_1011.txt');
Node1012_SI_NMBC_rot=load('Stress_intensity_NOMicro_BC_rot_1012.txt');
Node1018_SI_NMBC_rot=load('Stress_intensity_NOMicro_BC_rot_1018.txt');
Node1022_SI_NMBC_rot=load('Stress_intensity_NOMicro_BC_rot_1022.txt');
Node1026_SI_NMBC_rot=load('Stress_intensity_NOMicro_BC_rot_1026.txt');
Node1030_SI_NMBC_rot=load('Stress_intensity_NOMicro_BC_rot_1030.txt');
Node1034_SI_NMBC_rot=load('Stress_intensity_NOMicro_BC_rot_1034.txt');
Node1038_SI_NMBC_rot=load('Stress_intensity_NOMicro_BC_rot_1038.txt');
Node1042_SI_NMBC_rot=load('Stress_intensity_NOMicro_BC_rot_1042.txt');
Node1046_SI_NMBC_rot=load('Stress_intensity_NOMicro_BC_rot_1046.txt');
Node1050_SI_NMBC_rot=load('Stress_intensity_NOMicro_BC_rot_1050.txt');
Node1054_SI_NMBC_rot=load('Stress_intensity_NOMicro_BC_rot_1054.txt');
Node1058_SI_NMBC_rot=load('Stress_intensity_NOMicro_BC_rot_1058.txt');
Node1062_SI_NMBC_rot=load('Stress_intensity_NOMicro_BC_rot_1062.txt');
Node1066_SI_NMBC_rot=load('Stress_intensity_NOMicro_BC_rot_1066.txt');
Node1070_SI_NMBC_rot=load('Stress_intensity_NOMicro_BC_rot_1070.txt');
Node1074_SI_NMBC_rot=load('Stress_intensity_NOMicro_BC_rot_1074.txt');
Node1078_SI_NMBC_rot=load('Stress_intensity_NOMicro_BC_rot_1078.txt');
Node1082_SI_NMBC_rot=load('Stress_intensity_NOMicro_BC_rot_1082.txt');
Node1086_SI_NMBC_rot=load('Stress_intensity_NOMicro_BC_rot_1086.txt');
Node1090_SI_NMBC_rot=load('Stress_intensity_NOMicro_BC_rot_1090.txt');
Node1094_SI_NMBC_rot=load('Stress_intensity_NOMicro_BC_rot_1094.txt');
Node1098_SI_NMBC_rot=load('Stress_intensity_NOMicro_BC_rot_1098.txt');
Node1102_SI_NMBC_rot=load('Stress_intensity_NOMicro_BC_rot_1102.txt');
Node1106_SI_NMBC_rot=load('Stress_intensity_NOMicro_BC_rot_1106.txt');
Node1110_SI_NMBC_rot=load('Stress_intensity_NOMicro_BC_rot_1110.txt');
Node1114_SI_NMBC_rot=load('Stress_intensity_NOMicro_BC_rot_1114.txt');
Node1118_SI_NMBC_rot=load('Stress_intensity_NOMicro_BC_rot_1118.txt');
Node1122_SI_NMBC_rot=load('Stress_intensity_NOMicro_BC_rot_1122.txt');
Node1126_SI_NMBC_rot=load('Stress_intensity_NOMicro_BC_rot_1126.txt');
Node1130_SI_NMBC_rot=load('Stress_intensity_NOMicro_BC_rot_1130.txt');
Node1134_SI_NMBC_rot=load('Stress_intensity_NOMicro_BC_rot_1134.txt');
Node1138_SI_NMBC_rot=load('Stress_intensity_NOMicro_BC_rot_1138.txt');
Node1142_SI_NMBC_rot=load('Stress_intensity_NOMicro_BC_rot_1142.txt');
Node1146_SI_NMBC_rot=load('Stress_intensity_NOMicro_BC_rot_1146.txt');
Node1150_SI_NMBC_rot=load('Stress_intensity_NOMicro_BC_rot_1150.txt');


Node1011_SI_NMBC_str=load('Stress_intensity_NOMicro_BC_str_1011.txt');
Node1012_SI_NMBC_str=load('Stress_intensity_NOMicro_BC_str_1012.txt');
Node1018_SI_NMBC_str=load('Stress_intensity_NOMicro_BC_str_1018.txt');
Node1022_SI_NMBC_str=load('Stress_intensity_NOMicro_BC_str_1022.txt');
Node1026_SI_NMBC_str=load('Stress_intensity_NOMicro_BC_str_1026.txt');
Node1030_SI_NMBC_str=load('Stress_intensity_NOMicro_BC_str_1030.txt');
Node1034_SI_NMBC_str=load('Stress_intensity_NOMicro_BC_str_1034.txt');
Node1038_SI_NMBC_str=load('Stress_intensity_NOMicro_BC_str_1038.txt');
Node1042_SI_NMBC_str=load('Stress_intensity_NOMicro_BC_str_1042.txt');
Node1046_SI_NMBC_str=load('Stress_intensity_NOMicro_BC_str_1046.txt');
Node1050_SI_NMBC_str=load('Stress_intensity_NOMicro_BC_str_1050.txt');
Node1054_SI_NMBC_str=load('Stress_intensity_NOMicro_BC_str_1054.txt');
Node1058_SI_NMBC_str=load('Stress_intensity_NOMicro_BC_str_1058.txt');
Node1062_SI_NMBC_str=load('Stress_intensity_NOMicro_BC_str_1062.txt');
Node1066_SI_NMBC_str=load('Stress_intensity_NOMicro_BC_str_1066.txt');
Node1070_SI_NMBC_str=load('Stress_intensity_NOMicro_BC_str_1070.txt');
Node1074_SI_NMBC_str=load('Stress_intensity_NOMicro_BC_str_1074.txt');
Node1078_SI_NMBC_str=load('Stress_intensity_NOMicro_BC_str_1078.txt');
Node1082_SI_NMBC_str=load('Stress_intensity_NOMicro_BC_str_1082.txt');
Node1086_SI_NMBC_str=load('Stress_intensity_NOMicro_BC_str_1086.txt');
Node1090_SI_NMBC_str=load('Stress_intensity_NOMicro_BC_str_1090.txt');
Node1094_SI_NMBC_str=load('Stress_intensity_NOMicro_BC_str_1094.txt');
Node1098_SI_NMBC_str=load('Stress_intensity_NOMicro_BC_str_1098.txt');
Node1102_SI_NMBC_str=load('Stress_intensity_NOMicro_BC_str_1102.txt');
Node1106_SI_NMBC_str=load('Stress_intensity_NOMicro_BC_str_1106.txt');
Node1110_SI_NMBC_str=load('Stress_intensity_NOMicro_BC_str_1110.txt');
Node1114_SI_NMBC_str=load('Stress_intensity_NOMicro_BC_str_1114.txt');
Node1118_SI_NMBC_str=load('Stress_intensity_NOMicro_BC_str_1118.txt');
Node1122_SI_NMBC_str=load('Stress_intensity_NOMicro_BC_str_1122.txt');
Node1126_SI_NMBC_str=load('Stress_intensity_NOMicro_BC_str_1126.txt');
Node1130_SI_NMBC_str=load('Stress_intensity_NOMicro_BC_str_1130.txt');
Node1134_SI_NMBC_str=load('Stress_intensity_NOMicro_BC_str_1134.txt');
Node1138_SI_NMBC_str=load('Stress_intensity_NOMicro_BC_str_1138.txt');
Node1142_SI_NMBC_str=load('Stress_intensity_NOMicro_BC_str_1142.txt');
Node1146_SI_NMBC_str=load('Stress_intensity_NOMicro_BC_str_1146.txt');
Node1150_SI_NMBC_str=load('Stress_intensity_NOMicro_BC_str_1150.txt');


P_applied=(Applied_load_141(1,5)+Applied_load_143(1,5)+Applied_load_215(1,5)+...
    Applied_load_287(1,5)+Applied_load_359(1,5)+Applied_load_431(1,5)+...
    Applied_load_503(1,5)+Applied_load_575(1,5)+Applied_load_647(1,5)+...
    Applied_load_719(1,5)+Applied_load_791(1,5)+Applied_load_863(1,5)+...
    Applied_load_935(1,5)+Applied_load_1007(1,5))/14;

P_applied_Sigma=(Applied_load_141(1,6)+Applied_load_143(1,6)+Applied_load_215(1,6)+...
    Applied_load_287(1,6)+Applied_load_359(1,6)+Applied_load_431(1,6)+...
    Applied_load_503(1,6)+Applied_load_575(1,6)+Applied_load_647(1,6)+...
    Applied_load_719(1,6)+Applied_load_791(1,6)+Applied_load_863(1,6)+...
    Applied_load_935(1,6)+Applied_load_1007(1,6))/14;


Result_S11_intensity_cls = [Node1012_SI_cls(1,3)+0.25 Node1012_SI_cls(1,7)/P_applied;Node1011_SI_cls(1,3)+0.25 Node1011_SI_cls(1,7)/P_applied;...
    Node1018_SI_cls(1,3)+0.25 Node1018_SI_cls(1,7)/P_applied;Node1022_SI_cls(1,3)+0.25 Node1022_SI_cls(1,7)/P_applied;...
    Node1022_SI_cls(1,3)+0.25 Node1022_SI_cls(1,7)/P_applied;Node1026_SI_cls(1,3)+0.25 Node1026_SI_cls(1,7)/P_applied;...
    Node1030_SI_cls(1,3)+0.25 Node1030_SI_cls(1,7)/P_applied;Node1034_SI_cls(1,3)+0.25 Node1034_SI_cls(1,7)/P_applied;...
    Node1038_SI_cls(1,3)+0.25 Node1038_SI_cls(1,7)/P_applied;Node1042_SI_cls(1,3)+0.25 Node1042_SI_cls(1,7)/P_applied;...
    Node1046_SI_cls(1,3)+0.25 Node1046_SI_cls(1,7)/P_applied;Node1050_SI_cls(1,3)+0.25 Node1050_SI_cls(1,7)/P_applied;...
    Node1054_SI_cls(1,3)+0.25 Node1054_SI_cls(1,7)/P_applied;Node1058_SI_cls(1,3)+0.25 Node1058_SI_cls(1,7)/P_applied;...
    Node1062_SI_cls(1,3)+0.25 Node1062_SI_cls(1,7)/P_applied;Node1066_SI_cls(1,3)+0.25 Node1066_SI_cls(1,7)/P_applied;...
    Node1070_SI_cls(1,3)+0.25 Node1070_SI_cls(1,7)/P_applied;Node1074_SI_cls(1,3)+0.25 Node1074_SI_cls(1,7)/P_applied;...
    Node1078_SI_cls(1,3)+0.25 Node1078_SI_cls(1,7)/P_applied;Node1082_SI_cls(1,3)+0.25 Node1082_SI_cls(1,7)/P_applied;...
    Node1086_SI_cls(1,3)+0.25 Node1086_SI_cls(1,7)/P_applied;Node1090_SI_cls(1,3)+0.25 Node1090_SI_cls(1,7)/P_applied;...
    Node1094_SI_cls(1,3)+0.25 Node1094_SI_cls(1,7)/P_applied;Node1098_SI_cls(1,3)+0.25 Node1098_SI_cls(1,7)/P_applied;...
    Node1102_SI_cls(1,3)+0.25 Node1102_SI_cls(1,7)/P_applied;Node1106_SI_cls(1,3)+0.25 Node1106_SI_cls(1,7)/P_applied;...
    Node1110_SI_cls(1,3)+0.25 Node1110_SI_cls(1,7)/P_applied;Node1114_SI_cls(1,3)+0.25 Node1114_SI_cls(1,7)/P_applied;...
    Node1118_SI_cls(1,3)+0.25 Node1118_SI_cls(1,7)/P_applied;Node1122_SI_cls(1,3)+0.25 Node1122_SI_cls(1,7)/P_applied;...
    Node1126_SI_cls(1,3)+0.25 Node1126_SI_cls(1,7)/P_applied;Node1130_SI_cls(1,3)+0.25 Node1130_SI_cls(1,7)/P_applied;...
    Node1134_SI_cls(1,3)+0.25 Node1134_SI_cls(1,7)/P_applied;Node1138_SI_cls(1,3)+0.25 Node1138_SI_cls(1,7)/P_applied;...
    Node1142_SI_cls(1,3)+0.25 Node1142_SI_cls(1,7)/P_applied;Node1146_SI_cls(1,3)+0.25 Node1146_SI_cls(1,7)/P_applied;...
    Node1150_SI_cls(1,3)+0.25 Node1150_SI_cls(1,7)/P_applied];


Result_S11_intensity_MBC_com = [Node1012_SI_MBC_com(1,3)+0.25 Node1012_SI_MBC_com(1,7)/P_applied;Node1011_SI_MBC_com(1,3)+0.25 Node1011_SI_MBC_com(1,7)/P_applied;...
    Node1018_SI_MBC_com(1,3)+0.25 Node1018_SI_MBC_com(1,7)/P_applied;Node1022_SI_MBC_com(1,3)+0.25 Node1022_SI_MBC_com(1,7)/P_applied;...
    Node1022_SI_MBC_com(1,3)+0.25 Node1022_SI_MBC_com(1,7)/P_applied;Node1026_SI_MBC_com(1,3)+0.25 Node1026_SI_MBC_com(1,7)/P_applied;...
    Node1030_SI_MBC_com(1,3)+0.25 Node1030_SI_MBC_com(1,7)/P_applied;Node1034_SI_MBC_com(1,3)+0.25 Node1034_SI_MBC_com(1,7)/P_applied;...
    Node1038_SI_MBC_com(1,3)+0.25 Node1038_SI_MBC_com(1,7)/P_applied;Node1042_SI_MBC_com(1,3)+0.25 Node1042_SI_MBC_com(1,7)/P_applied;...
    Node1046_SI_MBC_com(1,3)+0.25 Node1046_SI_MBC_com(1,7)/P_applied;Node1050_SI_MBC_com(1,3)+0.25 Node1050_SI_MBC_com(1,7)/P_applied;...
    Node1054_SI_MBC_com(1,3)+0.25 Node1054_SI_MBC_com(1,7)/P_applied;Node1058_SI_MBC_com(1,3)+0.25 Node1058_SI_MBC_com(1,7)/P_applied;...
    Node1062_SI_MBC_com(1,3)+0.25 Node1062_SI_MBC_com(1,7)/P_applied;Node1066_SI_MBC_com(1,3)+0.25 Node1066_SI_MBC_com(1,7)/P_applied;...
    Node1070_SI_MBC_com(1,3)+0.25 Node1070_SI_MBC_com(1,7)/P_applied;Node1074_SI_MBC_com(1,3)+0.25 Node1074_SI_MBC_com(1,7)/P_applied;...
    Node1078_SI_MBC_com(1,3)+0.25 Node1078_SI_MBC_com(1,7)/P_applied;Node1082_SI_MBC_com(1,3)+0.25 Node1082_SI_MBC_com(1,7)/P_applied;...
    Node1086_SI_MBC_com(1,3)+0.25 Node1086_SI_MBC_com(1,7)/P_applied;Node1090_SI_MBC_com(1,3)+0.25 Node1090_SI_MBC_com(1,7)/P_applied;...
    Node1094_SI_MBC_com(1,3)+0.25 Node1094_SI_MBC_com(1,7)/P_applied;Node1098_SI_MBC_com(1,3)+0.25 Node1098_SI_MBC_com(1,7)/P_applied;...
    Node1102_SI_MBC_com(1,3)+0.25 Node1102_SI_MBC_com(1,7)/P_applied;Node1106_SI_MBC_com(1,3)+0.25 Node1106_SI_MBC_com(1,7)/P_applied;...
    Node1110_SI_MBC_com(1,3)+0.25 Node1110_SI_MBC_com(1,7)/P_applied;Node1114_SI_MBC_com(1,3)+0.25 Node1114_SI_MBC_com(1,7)/P_applied;...
    Node1118_SI_MBC_com(1,3)+0.25 Node1118_SI_MBC_com(1,7)/P_applied;Node1122_SI_MBC_com(1,3)+0.25 Node1122_SI_MBC_com(1,7)/P_applied;...
    Node1126_SI_MBC_com(1,3)+0.25 Node1126_SI_MBC_com(1,7)/P_applied;Node1130_SI_MBC_com(1,3)+0.25 Node1130_SI_MBC_com(1,7)/P_applied;...
    Node1134_SI_MBC_com(1,3)+0.25 Node1134_SI_MBC_com(1,7)/P_applied;Node1138_SI_MBC_com(1,3)+0.25 Node1138_SI_MBC_com(1,7)/P_applied;...
    Node1142_SI_MBC_com(1,3)+0.25 Node1142_SI_MBC_com(1,7)/P_applied;Node1146_SI_MBC_com(1,3)+0.25 Node1146_SI_MBC_com(1,7)/P_applied;...
    Node1150_SI_MBC_com(1,3)+0.25 Node1150_SI_MBC_com(1,7)/P_applied];


Result_S11_intensity_MBC_rot = [Node1012_SI_MBC_rot(1,3)+0.25 Node1012_SI_MBC_rot(1,7)/P_applied;Node1011_SI_MBC_rot(1,3)+0.25 Node1011_SI_MBC_rot(1,7)/P_applied;...
    Node1018_SI_MBC_rot(1,3)+0.25 Node1018_SI_MBC_rot(1,7)/P_applied;Node1022_SI_MBC_rot(1,3)+0.25 Node1022_SI_MBC_rot(1,7)/P_applied;...
    Node1022_SI_MBC_rot(1,3)+0.25 Node1022_SI_MBC_rot(1,7)/P_applied;Node1026_SI_MBC_rot(1,3)+0.25 Node1026_SI_MBC_rot(1,7)/P_applied;...
    Node1030_SI_MBC_rot(1,3)+0.25 Node1030_SI_MBC_rot(1,7)/P_applied;Node1034_SI_MBC_rot(1,3)+0.25 Node1034_SI_MBC_rot(1,7)/P_applied;...
    Node1038_SI_MBC_rot(1,3)+0.25 Node1038_SI_MBC_rot(1,7)/P_applied;Node1042_SI_MBC_rot(1,3)+0.25 Node1042_SI_MBC_rot(1,7)/P_applied;...
    Node1046_SI_MBC_rot(1,3)+0.25 Node1046_SI_MBC_rot(1,7)/P_applied;Node1050_SI_MBC_rot(1,3)+0.25 Node1050_SI_MBC_rot(1,7)/P_applied;...
    Node1054_SI_MBC_rot(1,3)+0.25 Node1054_SI_MBC_rot(1,7)/P_applied;Node1058_SI_MBC_rot(1,3)+0.25 Node1058_SI_MBC_rot(1,7)/P_applied;...
    Node1062_SI_MBC_rot(1,3)+0.25 Node1062_SI_MBC_rot(1,7)/P_applied;Node1066_SI_MBC_rot(1,3)+0.25 Node1066_SI_MBC_rot(1,7)/P_applied;...
    Node1070_SI_MBC_rot(1,3)+0.25 Node1070_SI_MBC_rot(1,7)/P_applied;Node1074_SI_MBC_rot(1,3)+0.25 Node1074_SI_MBC_rot(1,7)/P_applied;...
    Node1078_SI_MBC_rot(1,3)+0.25 Node1078_SI_MBC_rot(1,7)/P_applied;Node1082_SI_MBC_rot(1,3)+0.25 Node1082_SI_MBC_rot(1,7)/P_applied;...
    Node1086_SI_MBC_rot(1,3)+0.25 Node1086_SI_MBC_rot(1,7)/P_applied;Node1090_SI_MBC_rot(1,3)+0.25 Node1090_SI_MBC_rot(1,7)/P_applied;...
    Node1094_SI_MBC_rot(1,3)+0.25 Node1094_SI_MBC_rot(1,7)/P_applied;Node1098_SI_MBC_rot(1,3)+0.25 Node1098_SI_MBC_rot(1,7)/P_applied;...
    Node1102_SI_MBC_rot(1,3)+0.25 Node1102_SI_MBC_rot(1,7)/P_applied;Node1106_SI_MBC_rot(1,3)+0.25 Node1106_SI_MBC_rot(1,7)/P_applied;...
    Node1110_SI_MBC_rot(1,3)+0.25 Node1110_SI_MBC_rot(1,7)/P_applied;Node1114_SI_MBC_rot(1,3)+0.25 Node1114_SI_MBC_rot(1,7)/P_applied;...
    Node1118_SI_MBC_rot(1,3)+0.25 Node1118_SI_MBC_rot(1,7)/P_applied;Node1122_SI_MBC_rot(1,3)+0.25 Node1122_SI_MBC_rot(1,7)/P_applied;...
    Node1126_SI_MBC_rot(1,3)+0.25 Node1126_SI_MBC_rot(1,7)/P_applied;Node1130_SI_MBC_rot(1,3)+0.25 Node1130_SI_MBC_rot(1,7)/P_applied;...
    Node1134_SI_MBC_rot(1,3)+0.25 Node1134_SI_MBC_rot(1,7)/P_applied;Node1138_SI_MBC_rot(1,3)+0.25 Node1138_SI_MBC_rot(1,7)/P_applied;...
    Node1142_SI_MBC_rot(1,3)+0.25 Node1142_SI_MBC_rot(1,7)/P_applied;Node1146_SI_MBC_rot(1,3)+0.25 Node1146_SI_MBC_rot(1,7)/P_applied;...
    Node1150_SI_MBC_rot(1,3)+0.25 Node1150_SI_MBC_rot(1,7)/P_applied];


Result_S11_intensity_MBC_str = [Node1012_SI_MBC_str(1,3)+0.25 Node1012_SI_MBC_str(1,7)/P_applied;Node1011_SI_MBC_str(1,3)+0.25 Node1011_SI_MBC_str(1,7)/P_applied;...
    Node1018_SI_MBC_str(1,3)+0.25 Node1018_SI_MBC_str(1,7)/P_applied;Node1022_SI_MBC_str(1,3)+0.25 Node1022_SI_MBC_str(1,7)/P_applied;...
    Node1022_SI_MBC_str(1,3)+0.25 Node1022_SI_MBC_str(1,7)/P_applied;Node1026_SI_MBC_str(1,3)+0.25 Node1026_SI_MBC_str(1,7)/P_applied;...
    Node1030_SI_MBC_str(1,3)+0.25 Node1030_SI_MBC_str(1,7)/P_applied;Node1034_SI_MBC_str(1,3)+0.25 Node1034_SI_MBC_str(1,7)/P_applied;...
    Node1038_SI_MBC_str(1,3)+0.25 Node1038_SI_MBC_str(1,7)/P_applied;Node1042_SI_MBC_str(1,3)+0.25 Node1042_SI_MBC_str(1,7)/P_applied;...
    Node1046_SI_MBC_str(1,3)+0.25 Node1046_SI_MBC_str(1,7)/P_applied;Node1050_SI_MBC_str(1,3)+0.25 Node1050_SI_MBC_str(1,7)/P_applied;...
    Node1054_SI_MBC_str(1,3)+0.25 Node1054_SI_MBC_str(1,7)/P_applied;Node1058_SI_MBC_str(1,3)+0.25 Node1058_SI_MBC_str(1,7)/P_applied;...
    Node1062_SI_MBC_str(1,3)+0.25 Node1062_SI_MBC_str(1,7)/P_applied;Node1066_SI_MBC_str(1,3)+0.25 Node1066_SI_MBC_str(1,7)/P_applied;...
    Node1070_SI_MBC_str(1,3)+0.25 Node1070_SI_MBC_str(1,7)/P_applied;Node1074_SI_MBC_str(1,3)+0.25 Node1074_SI_MBC_str(1,7)/P_applied;...
    Node1078_SI_MBC_str(1,3)+0.25 Node1078_SI_MBC_str(1,7)/P_applied;Node1082_SI_MBC_str(1,3)+0.25 Node1082_SI_MBC_str(1,7)/P_applied;...
    Node1086_SI_MBC_str(1,3)+0.25 Node1086_SI_MBC_str(1,7)/P_applied;Node1090_SI_MBC_str(1,3)+0.25 Node1090_SI_MBC_str(1,7)/P_applied;...
    Node1094_SI_MBC_str(1,3)+0.25 Node1094_SI_MBC_str(1,7)/P_applied;Node1098_SI_MBC_str(1,3)+0.25 Node1098_SI_MBC_str(1,7)/P_applied;...
    Node1102_SI_MBC_str(1,3)+0.25 Node1102_SI_MBC_str(1,7)/P_applied;Node1106_SI_MBC_str(1,3)+0.25 Node1106_SI_MBC_str(1,7)/P_applied;...
    Node1110_SI_MBC_str(1,3)+0.25 Node1110_SI_MBC_str(1,7)/P_applied;Node1114_SI_MBC_str(1,3)+0.25 Node1114_SI_MBC_str(1,7)/P_applied;...
    Node1118_SI_MBC_str(1,3)+0.25 Node1118_SI_MBC_str(1,7)/P_applied;Node1122_SI_MBC_str(1,3)+0.25 Node1122_SI_MBC_str(1,7)/P_applied;...
    Node1126_SI_MBC_str(1,3)+0.25 Node1126_SI_MBC_str(1,7)/P_applied;Node1130_SI_MBC_str(1,3)+0.25 Node1130_SI_MBC_str(1,7)/P_applied;...
    Node1134_SI_MBC_str(1,3)+0.25 Node1134_SI_MBC_str(1,7)/P_applied;Node1138_SI_MBC_str(1,3)+0.25 Node1138_SI_MBC_str(1,7)/P_applied;...
    Node1142_SI_MBC_str(1,3)+0.25 Node1142_SI_MBC_str(1,7)/P_applied;Node1146_SI_MBC_str(1,3)+0.25 Node1146_SI_MBC_str(1,7)/P_applied;...
    Node1150_SI_MBC_str(1,3)+0.25 Node1150_SI_MBC_str(1,7)/P_applied];


Result_S11_intensity_NMBC_com = [Node1012_SI_NMBC_com(1,3)+0.25 Node1012_SI_NMBC_com(1,7)/P_applied;Node1011_SI_NMBC_com(1,3)+0.25 Node1011_SI_NMBC_com(1,7)/P_applied;...
    Node1018_SI_NMBC_com(1,3)+0.25 Node1018_SI_NMBC_com(1,7)/P_applied;Node1022_SI_NMBC_com(1,3)+0.25 Node1022_SI_NMBC_com(1,7)/P_applied;...
    Node1022_SI_NMBC_com(1,3)+0.25 Node1022_SI_NMBC_com(1,7)/P_applied;Node1026_SI_NMBC_com(1,3)+0.25 Node1026_SI_NMBC_com(1,7)/P_applied;...
    Node1030_SI_NMBC_com(1,3)+0.25 Node1030_SI_NMBC_com(1,7)/P_applied;Node1034_SI_NMBC_com(1,3)+0.25 Node1034_SI_NMBC_com(1,7)/P_applied;...
    Node1038_SI_NMBC_com(1,3)+0.25 Node1038_SI_NMBC_com(1,7)/P_applied;Node1042_SI_NMBC_com(1,3)+0.25 Node1042_SI_NMBC_com(1,7)/P_applied;...
    Node1046_SI_NMBC_com(1,3)+0.25 Node1046_SI_NMBC_com(1,7)/P_applied;Node1050_SI_NMBC_com(1,3)+0.25 Node1050_SI_NMBC_com(1,7)/P_applied;...
    Node1054_SI_NMBC_com(1,3)+0.25 Node1054_SI_NMBC_com(1,7)/P_applied;Node1058_SI_NMBC_com(1,3)+0.25 Node1058_SI_NMBC_com(1,7)/P_applied;...
    Node1062_SI_NMBC_com(1,3)+0.25 Node1062_SI_NMBC_com(1,7)/P_applied;Node1066_SI_NMBC_com(1,3)+0.25 Node1066_SI_NMBC_com(1,7)/P_applied;...
    Node1070_SI_NMBC_com(1,3)+0.25 Node1070_SI_NMBC_com(1,7)/P_applied;Node1074_SI_NMBC_com(1,3)+0.25 Node1074_SI_NMBC_com(1,7)/P_applied;...
    Node1078_SI_NMBC_com(1,3)+0.25 Node1078_SI_NMBC_com(1,7)/P_applied;Node1082_SI_NMBC_com(1,3)+0.25 Node1082_SI_NMBC_com(1,7)/P_applied;...
    Node1086_SI_NMBC_com(1,3)+0.25 Node1086_SI_NMBC_com(1,7)/P_applied;Node1090_SI_NMBC_com(1,3)+0.25 Node1090_SI_NMBC_com(1,7)/P_applied;...
    Node1094_SI_NMBC_com(1,3)+0.25 Node1094_SI_NMBC_com(1,7)/P_applied;Node1098_SI_NMBC_com(1,3)+0.25 Node1098_SI_NMBC_com(1,7)/P_applied;...
    Node1102_SI_NMBC_com(1,3)+0.25 Node1102_SI_NMBC_com(1,7)/P_applied;Node1106_SI_NMBC_com(1,3)+0.25 Node1106_SI_NMBC_com(1,7)/P_applied;...
    Node1110_SI_NMBC_com(1,3)+0.25 Node1110_SI_NMBC_com(1,7)/P_applied;Node1114_SI_NMBC_com(1,3)+0.25 Node1114_SI_NMBC_com(1,7)/P_applied;...
    Node1118_SI_NMBC_com(1,3)+0.25 Node1118_SI_NMBC_com(1,7)/P_applied;Node1122_SI_NMBC_com(1,3)+0.25 Node1122_SI_NMBC_com(1,7)/P_applied;...
    Node1126_SI_NMBC_com(1,3)+0.25 Node1126_SI_NMBC_com(1,7)/P_applied;Node1130_SI_NMBC_com(1,3)+0.25 Node1130_SI_NMBC_com(1,7)/P_applied;...
    Node1134_SI_NMBC_com(1,3)+0.25 Node1134_SI_NMBC_com(1,7)/P_applied;Node1138_SI_NMBC_com(1,3)+0.25 Node1138_SI_NMBC_com(1,7)/P_applied;...
    Node1142_SI_NMBC_com(1,3)+0.25 Node1142_SI_NMBC_com(1,7)/P_applied;Node1146_SI_NMBC_com(1,3)+0.25 Node1146_SI_NMBC_com(1,7)/P_applied;...
    Node1150_SI_NMBC_com(1,3)+0.25 Node1150_SI_NMBC_com(1,7)/P_applied];


Result_S11_intensity_NMBC_rot = [Node1012_SI_NMBC_rot(1,3)+0.25 Node1012_SI_NMBC_rot(1,7)/P_applied;Node1011_SI_NMBC_rot(1,3)+0.25 Node1011_SI_NMBC_rot(1,7)/P_applied;...
    Node1018_SI_NMBC_rot(1,3)+0.25 Node1018_SI_NMBC_rot(1,7)/P_applied;Node1022_SI_NMBC_rot(1,3)+0.25 Node1022_SI_NMBC_rot(1,7)/P_applied;...
    Node1022_SI_NMBC_rot(1,3)+0.25 Node1022_SI_NMBC_rot(1,7)/P_applied;Node1026_SI_NMBC_rot(1,3)+0.25 Node1026_SI_NMBC_rot(1,7)/P_applied;...
    Node1030_SI_NMBC_rot(1,3)+0.25 Node1030_SI_NMBC_rot(1,7)/P_applied;Node1034_SI_NMBC_rot(1,3)+0.25 Node1034_SI_NMBC_rot(1,7)/P_applied;...
    Node1038_SI_NMBC_rot(1,3)+0.25 Node1038_SI_NMBC_rot(1,7)/P_applied;Node1042_SI_NMBC_rot(1,3)+0.25 Node1042_SI_NMBC_rot(1,7)/P_applied;...
    Node1046_SI_NMBC_rot(1,3)+0.25 Node1046_SI_NMBC_rot(1,7)/P_applied;Node1050_SI_NMBC_rot(1,3)+0.25 Node1050_SI_NMBC_rot(1,7)/P_applied;...
    Node1054_SI_NMBC_rot(1,3)+0.25 Node1054_SI_NMBC_rot(1,7)/P_applied;Node1058_SI_NMBC_rot(1,3)+0.25 Node1058_SI_NMBC_rot(1,7)/P_applied;...
    Node1062_SI_NMBC_rot(1,3)+0.25 Node1062_SI_NMBC_rot(1,7)/P_applied;Node1066_SI_NMBC_rot(1,3)+0.25 Node1066_SI_NMBC_rot(1,7)/P_applied;...
    Node1070_SI_NMBC_rot(1,3)+0.25 Node1070_SI_NMBC_rot(1,7)/P_applied;Node1074_SI_NMBC_rot(1,3)+0.25 Node1074_SI_NMBC_rot(1,7)/P_applied;...
    Node1078_SI_NMBC_rot(1,3)+0.25 Node1078_SI_NMBC_rot(1,7)/P_applied;Node1082_SI_NMBC_rot(1,3)+0.25 Node1082_SI_NMBC_rot(1,7)/P_applied;...
    Node1086_SI_NMBC_rot(1,3)+0.25 Node1086_SI_NMBC_rot(1,7)/P_applied;Node1090_SI_NMBC_rot(1,3)+0.25 Node1090_SI_NMBC_rot(1,7)/P_applied;...
    Node1094_SI_NMBC_rot(1,3)+0.25 Node1094_SI_NMBC_rot(1,7)/P_applied;Node1098_SI_NMBC_rot(1,3)+0.25 Node1098_SI_NMBC_rot(1,7)/P_applied;...
    Node1102_SI_NMBC_rot(1,3)+0.25 Node1102_SI_NMBC_rot(1,7)/P_applied;Node1106_SI_NMBC_rot(1,3)+0.25 Node1106_SI_NMBC_rot(1,7)/P_applied;...
    Node1110_SI_NMBC_rot(1,3)+0.25 Node1110_SI_NMBC_rot(1,7)/P_applied;Node1114_SI_NMBC_rot(1,3)+0.25 Node1114_SI_NMBC_rot(1,7)/P_applied;...
    Node1118_SI_NMBC_rot(1,3)+0.25 Node1118_SI_NMBC_rot(1,7)/P_applied;Node1122_SI_NMBC_rot(1,3)+0.25 Node1122_SI_NMBC_rot(1,7)/P_applied;...
    Node1126_SI_NMBC_rot(1,3)+0.25 Node1126_SI_NMBC_rot(1,7)/P_applied;Node1130_SI_NMBC_rot(1,3)+0.25 Node1130_SI_NMBC_rot(1,7)/P_applied;...
    Node1134_SI_NMBC_rot(1,3)+0.25 Node1134_SI_NMBC_rot(1,7)/P_applied;Node1138_SI_NMBC_rot(1,3)+0.25 Node1138_SI_NMBC_rot(1,7)/P_applied;...
    Node1142_SI_NMBC_rot(1,3)+0.25 Node1142_SI_NMBC_rot(1,7)/P_applied;Node1146_SI_NMBC_rot(1,3)+0.25 Node1146_SI_NMBC_rot(1,7)/P_applied;...
    Node1150_SI_NMBC_rot(1,3)+0.25 Node1150_SI_NMBC_rot(1,7)/P_applied];


Result_S11_intensity_NMBC_str = [Node1012_SI_NMBC_str(1,3)+0.25 Node1012_SI_NMBC_str(1,7)/P_applied;Node1011_SI_NMBC_str(1,3)+0.25 Node1011_SI_NMBC_str(1,7)/P_applied;...
    Node1018_SI_NMBC_str(1,3)+0.25 Node1018_SI_NMBC_str(1,7)/P_applied;Node1022_SI_NMBC_str(1,3)+0.25 Node1022_SI_NMBC_str(1,7)/P_applied;...
    Node1022_SI_NMBC_str(1,3)+0.25 Node1022_SI_NMBC_str(1,7)/P_applied;Node1026_SI_NMBC_str(1,3)+0.25 Node1026_SI_NMBC_str(1,7)/P_applied;...
    Node1030_SI_NMBC_str(1,3)+0.25 Node1030_SI_NMBC_str(1,7)/P_applied;Node1034_SI_NMBC_str(1,3)+0.25 Node1034_SI_NMBC_str(1,7)/P_applied;...
    Node1038_SI_NMBC_str(1,3)+0.25 Node1038_SI_NMBC_str(1,7)/P_applied;Node1042_SI_NMBC_str(1,3)+0.25 Node1042_SI_NMBC_str(1,7)/P_applied;...
    Node1046_SI_NMBC_str(1,3)+0.25 Node1046_SI_NMBC_str(1,7)/P_applied;Node1050_SI_NMBC_str(1,3)+0.25 Node1050_SI_NMBC_str(1,7)/P_applied;...
    Node1054_SI_NMBC_str(1,3)+0.25 Node1054_SI_NMBC_str(1,7)/P_applied;Node1058_SI_NMBC_str(1,3)+0.25 Node1058_SI_NMBC_str(1,7)/P_applied;...
    Node1062_SI_NMBC_str(1,3)+0.25 Node1062_SI_NMBC_str(1,7)/P_applied;Node1066_SI_NMBC_str(1,3)+0.25 Node1066_SI_NMBC_str(1,7)/P_applied;...
    Node1070_SI_NMBC_str(1,3)+0.25 Node1070_SI_NMBC_str(1,7)/P_applied;Node1074_SI_NMBC_str(1,3)+0.25 Node1074_SI_NMBC_str(1,7)/P_applied;...
    Node1078_SI_NMBC_str(1,3)+0.25 Node1078_SI_NMBC_str(1,7)/P_applied;Node1082_SI_NMBC_str(1,3)+0.25 Node1082_SI_NMBC_str(1,7)/P_applied;...
    Node1086_SI_NMBC_str(1,3)+0.25 Node1086_SI_NMBC_str(1,7)/P_applied;Node1090_SI_NMBC_str(1,3)+0.25 Node1090_SI_NMBC_str(1,7)/P_applied;...
    Node1094_SI_NMBC_str(1,3)+0.25 Node1094_SI_NMBC_str(1,7)/P_applied;Node1098_SI_NMBC_str(1,3)+0.25 Node1098_SI_NMBC_str(1,7)/P_applied;...
    Node1102_SI_NMBC_str(1,3)+0.25 Node1102_SI_NMBC_str(1,7)/P_applied;Node1106_SI_NMBC_str(1,3)+0.25 Node1106_SI_NMBC_str(1,7)/P_applied;...
    Node1110_SI_NMBC_str(1,3)+0.25 Node1110_SI_NMBC_str(1,7)/P_applied;Node1114_SI_NMBC_str(1,3)+0.25 Node1114_SI_NMBC_str(1,7)/P_applied;...
    Node1118_SI_NMBC_str(1,3)+0.25 Node1118_SI_NMBC_str(1,7)/P_applied;Node1122_SI_NMBC_str(1,3)+0.25 Node1122_SI_NMBC_str(1,7)/P_applied;...
    Node1126_SI_NMBC_str(1,3)+0.25 Node1126_SI_NMBC_str(1,7)/P_applied;Node1130_SI_NMBC_str(1,3)+0.25 Node1130_SI_NMBC_str(1,7)/P_applied;...
    Node1134_SI_NMBC_str(1,3)+0.25 Node1134_SI_NMBC_str(1,7)/P_applied;Node1138_SI_NMBC_str(1,3)+0.25 Node1138_SI_NMBC_str(1,7)/P_applied;...
    Node1142_SI_NMBC_str(1,3)+0.25 Node1142_SI_NMBC_str(1,7)/P_applied;Node1146_SI_NMBC_str(1,3)+0.25 Node1146_SI_NMBC_str(1,7)/P_applied;...
    Node1150_SI_NMBC_str(1,3)+0.25 Node1150_SI_NMBC_str(1,7)/P_applied];





Result_S_2_cls = [Node2_AH_cls(1,5) Node2_AH_cls(1,6) Node2_AH_cls(1,7);Node2_AH_cls(1,8) Node2_AH_cls(1,9) Node2_AH_cls(1,10);...
    Node2_AH_cls(1,11) Node2_AH_cls(1,12) Node2_AH_cls(1,13)];
X2 = sqrt(0.05^2-(Node2_AH_cls(1,3)+0.25)^2);
teta2 = asin((Node2_AH_cls(1,3)+0.25)/0.05);
R = [X2/0.05 -(Node2_AH_cls(1,3)+0.25)/0.05 0;(Node2_AH_cls(1,3)+0.25)/0.05 X2/0.05 0;...
    0 0 1];
Result_S_2_cls_trans = R'*Result_S_2_cls*R;


Result_S_6_cls = [Node6_AH_cls(1,5) Node6_AH_cls(1,6) Node6_AH_cls(1,7);Node6_AH_cls(1,8) Node6_AH_cls(1,9) Node6_AH_cls(1,10);...
    Node6_AH_cls(1,11) Node6_AH_cls(1,12) Node6_AH_cls(1,13)];
X6 = sqrt(0.05^2-(Node6_AH_cls(1,3)+0.25)^2);
teta6 = asin((Node6_AH_cls(1,3)+0.25)/0.05);
R = [X6/0.05 -(Node6_AH_cls(1,3)+0.25)/0.05 0;(Node6_AH_cls(1,3)+0.25)/0.05 X6/0.05 0;...
    0 0 1];
Result_S_6_cls_trans = R'*Result_S_6_cls*R;


Result_S_146_cls = [Node146_AH_cls(1,5) Node146_AH_cls(1,6) Node146_AH_cls(1,7);Node146_AH_cls(1,8) Node146_AH_cls(1,9) Node146_AH_cls(1,10);...
    Node146_AH_cls(1,11) Node146_AH_cls(1,12) Node146_AH_cls(1,13)];
X146 = sqrt(0.05^2-(Node146_AH_cls(1,3)+0.25)^2);
teta146 = asin((Node146_AH_cls(1,3)+0.25)/0.05);
R = [X146/0.05 -(Node146_AH_cls(1,3)+0.25)/0.05 0;(Node146_AH_cls(1,3)+0.25)/0.05 X146/0.05 0;...
    0 0 1];
Result_S_146_cls_trans = R'*Result_S_146_cls*R;


Result_S_218_cls = [Node218_AH_cls(1,5) Node218_AH_cls(1,6) Node218_AH_cls(1,7);Node218_AH_cls(1,8) Node218_AH_cls(1,9) Node218_AH_cls(1,10);...
    Node218_AH_cls(1,11) Node218_AH_cls(1,12) Node218_AH_cls(1,13)];
X218 = sqrt(0.05^2-(Node218_AH_cls(1,3)+0.25)^2);
teta218 = asin((Node218_AH_cls(1,3)+0.25)/0.05);
R = [X218/0.05 -(Node218_AH_cls(1,3)+0.25)/0.05 0;(Node218_AH_cls(1,3)+0.25)/0.05 X218/0.05 0;...
    0 0 1];
Result_S_218_cls_trans = R'*Result_S_218_cls*R;


Result_S_290_cls = [Node290_AH_cls(1,5) Node290_AH_cls(1,6) Node290_AH_cls(1,7);Node290_AH_cls(1,8) Node290_AH_cls(1,9) Node290_AH_cls(1,10);...
    Node290_AH_cls(1,11) Node290_AH_cls(1,12) Node290_AH_cls(1,13)];
X290 = sqrt(0.05^2-(Node290_AH_cls(1,3)+0.25)^2);
teta290 = asin((Node290_AH_cls(1,3)+0.25)/0.05);
R = [X290/0.05 -(Node290_AH_cls(1,3)+0.25)/0.05 0;(Node290_AH_cls(1,3)+0.25)/0.05 X290/0.05 0;...
    0 0 1];
Result_S_290_cls_trans = R'*Result_S_290_cls*R;


Result_S_362_cls = [Node362_AH_cls(1,5) Node362_AH_cls(1,6) Node362_AH_cls(1,7);Node362_AH_cls(1,8) Node362_AH_cls(1,9) Node362_AH_cls(1,10);...
    Node362_AH_cls(1,11) Node362_AH_cls(1,12) Node362_AH_cls(1,13)];
X362 = sqrt(0.05^2-(Node362_AH_cls(1,3)+0.25)^2);
teta362 = asin((Node362_AH_cls(1,3)+0.25)/0.05);
R = [X362/0.05 -(Node362_AH_cls(1,3)+0.25)/0.05 0;(Node362_AH_cls(1,3)+0.25)/0.05 X362/0.05 0;...
    0 0 1];
Result_S_362_cls_trans = R'*Result_S_362_cls*R;


Result_S_434_cls = [Node434_AH_cls(1,5) Node434_AH_cls(1,6) Node434_AH_cls(1,7);Node434_AH_cls(1,8) Node434_AH_cls(1,9) Node434_AH_cls(1,10);...
    Node434_AH_cls(1,11) Node434_AH_cls(1,12) Node434_AH_cls(1,13)];
X434 = sqrt(0.05^2-(Node434_AH_cls(1,3)+0.25)^2);
teta434 = asin((Node434_AH_cls(1,3)+0.25)/0.05);
R = [X434/0.05 -(Node434_AH_cls(1,3)+0.25)/0.05 0;(Node434_AH_cls(1,3)+0.25)/0.05 X434/0.05 0;...
    0 0 1];
Result_S_434_cls_trans = R'*Result_S_434_cls*R;


Result_S_506_cls = [Node506_AH_cls(1,5) Node506_AH_cls(1,6) Node506_AH_cls(1,7);Node506_AH_cls(1,8) Node506_AH_cls(1,9) Node506_AH_cls(1,10);...
    Node506_AH_cls(1,11) Node506_AH_cls(1,12) Node506_AH_cls(1,13)];
X506 = sqrt(0.05^2-(Node506_AH_cls(1,3)+0.25)^2);
teta506 = asin((Node506_AH_cls(1,3)+0.25)/0.05);
R = [X506/0.05 -(Node506_AH_cls(1,3)+0.25)/0.05 0;(Node506_AH_cls(1,3)+0.25)/0.05 X506/0.05 0;...
    0 0 1];
Result_S_506_cls_trans = R'*Result_S_506_cls*R;


Result_S_578_cls = [Node578_AH_cls(1,5) Node578_AH_cls(1,6) Node578_AH_cls(1,7);Node578_AH_cls(1,8) Node578_AH_cls(1,9) Node578_AH_cls(1,10);...
    Node578_AH_cls(1,11) Node578_AH_cls(1,12) Node578_AH_cls(1,13)];
X578 = sqrt(0.05^2-(Node578_AH_cls(1,3)+0.25)^2);
teta578 = asin((Node578_AH_cls(1,3)+0.25)/0.05);
R = [X578/0.05 -(Node578_AH_cls(1,3)+0.25)/0.05 0;(Node578_AH_cls(1,3)+0.25)/0.05 X578/0.05 0;...
    0 0 1];
Result_S_578_cls_trans = R'*Result_S_578_cls*R;


Result_S_650_cls = [Node650_AH_cls(1,5) Node650_AH_cls(1,6) Node650_AH_cls(1,7);Node650_AH_cls(1,8) Node650_AH_cls(1,9) Node650_AH_cls(1,10);...
    Node650_AH_cls(1,11) Node650_AH_cls(1,12) Node650_AH_cls(1,13)];
X650 = sqrt(0.05^2-(Node650_AH_cls(1,3)+0.25)^2);
teta650 = asin((Node650_AH_cls(1,3)+0.25)/0.05);
R = [X650/0.05 -(Node650_AH_cls(1,3)+0.25)/0.05 0;(Node650_AH_cls(1,3)+0.25)/0.05 X650/0.05 0;...
    0 0 1];
Result_S_650_cls_trans = R'*Result_S_650_cls*R;


Result_S_722_cls = [Node722_AH_cls(1,5) Node722_AH_cls(1,6) Node722_AH_cls(1,7);Node722_AH_cls(1,8) Node722_AH_cls(1,9) Node722_AH_cls(1,10);...
    Node722_AH_cls(1,11) Node722_AH_cls(1,12) Node722_AH_cls(1,13)];
X722 = sqrt(0.05^2-(Node722_AH_cls(1,3)+0.25)^2);
teta722 = asin((Node722_AH_cls(1,3)+0.25)/0.05);
R = [X722/0.05 -(Node722_AH_cls(1,3)+0.25)/0.05 0;(Node722_AH_cls(1,3)+0.25)/0.05 X722/0.05 0;...
    0 0 1];
Result_S_722_cls_trans = R'*Result_S_722_cls*R;


Result_S_794_cls = [Node794_AH_cls(1,5) Node794_AH_cls(1,6) Node794_AH_cls(1,7);Node794_AH_cls(1,8) Node794_AH_cls(1,9) Node794_AH_cls(1,10);...
    Node794_AH_cls(1,11) Node794_AH_cls(1,12) Node794_AH_cls(1,13)];
X794 = sqrt(0.05^2-(Node794_AH_cls(1,3)+0.25)^2);
teta794 = asin((Node794_AH_cls(1,3)+0.25)/0.05);
R = [X794/0.05 -(Node794_AH_cls(1,3)+0.25)/0.05 0;(Node794_AH_cls(1,3)+0.25)/0.05 X794/0.05 0;...
    0 0 1];
Result_S_794_cls_trans = R'*Result_S_794_cls*R;


Result_S_866_cls = [Node866_AH_cls(1,5) Node866_AH_cls(1,6) Node866_AH_cls(1,7);Node866_AH_cls(1,8) Node866_AH_cls(1,9) Node866_AH_cls(1,10);...
    Node866_AH_cls(1,11) Node866_AH_cls(1,12) Node866_AH_cls(1,13)];
X866 = sqrt(0.05^2-(Node866_AH_cls(1,3)+0.25)^2);
teta866 = asin((Node866_AH_cls(1,3)+0.25)/0.05);
R = [X866/0.05 -(Node866_AH_cls(1,3)+0.25)/0.05 0;(Node866_AH_cls(1,3)+0.25)/0.05 X866/0.05 0;...
    0 0 1];
Result_S_866_cls_trans = R'*Result_S_866_cls*R;


Result_S_938_cls = [Node938_AH_cls(1,5) Node938_AH_cls(1,6) Node938_AH_cls(1,7);Node938_AH_cls(1,8) Node938_AH_cls(1,9) Node938_AH_cls(1,10);...
    Node938_AH_cls(1,11) Node938_AH_cls(1,12) Node938_AH_cls(1,13)];
X938 = sqrt(0.05^2-(Node938_AH_cls(1,3)+0.25)^2);
teta938 = asin((Node938_AH_cls(1,3)+0.25)/0.05);
R = [X938/0.05 -(Node938_AH_cls(1,3)+0.25)/0.05 0;(Node938_AH_cls(1,3)+0.25)/0.05 X938/0.05 0;...
    0 0 1];
Result_S_938_cls_trans = R'*Result_S_938_cls*R;


Result_S_1876_cls = [Node1876_AH_cls(1,5) Node1876_AH_cls(1,6) Node1876_AH_cls(1,7);Node1876_AH_cls(1,8) Node1876_AH_cls(1,9) Node1876_AH_cls(1,10);...
    Node1876_AH_cls(1,11) Node1876_AH_cls(1,12) Node1876_AH_cls(1,13)];
X1876 = sqrt(0.05^2-(Node1876_AH_cls(1,3)+0.25)^2);
teta1876 = asin((Node1876_AH_cls(1,3)+0.25)/0.05);
R = [X1876/0.05 -(Node1876_AH_cls(1,3)+0.25)/0.05 0;(Node1876_AH_cls(1,3)+0.25)/0.05 X1876/0.05 0;...
    0 0 1];
Result_S_1876_cls_trans = R'*Result_S_1876_cls*R;


Result_S_1804_cls = [Node1804_AH_cls(1,5) Node1804_AH_cls(1,6) Node1804_AH_cls(1,7);Node1804_AH_cls(1,8) Node1804_AH_cls(1,9) Node1804_AH_cls(1,10);...
    Node1804_AH_cls(1,11) Node1804_AH_cls(1,12) Node1804_AH_cls(1,13)];
X1804 = sqrt(0.05^2-(Node1804_AH_cls(1,3)+0.25)^2);
teta1804 = asin((Node1804_AH_cls(1,3)+0.25)/0.05);
R = [X1804/0.05 -(Node1804_AH_cls(1,3)+0.25)/0.05 0;(Node1804_AH_cls(1,3)+0.25)/0.05 X1804/0.05 0;...
    0 0 1];
Result_S_1804_cls_trans = R'*Result_S_1804_cls*R;


Result_S_1732_cls = [Node1732_AH_cls(1,5) Node1732_AH_cls(1,6) Node1732_AH_cls(1,7);Node1732_AH_cls(1,8) Node1732_AH_cls(1,9) Node1732_AH_cls(1,10);...
    Node1732_AH_cls(1,11) Node1732_AH_cls(1,12) Node1732_AH_cls(1,13)];
X1732 = sqrt(0.05^2-(Node1732_AH_cls(1,3)+0.25)^2);
teta1732 = asin((Node1732_AH_cls(1,3)+0.25)/0.05);
R = [X1732/0.05 -(Node1732_AH_cls(1,3)+0.25)/0.05 0;(Node1732_AH_cls(1,3)+0.25)/0.05 X1732/0.05 0;...
    0 0 1];
Result_S_1732_cls_trans = R'*Result_S_1732_cls*R;


Result_S_1660_cls = [Node1660_AH_cls(1,5) Node1660_AH_cls(1,6) Node1660_AH_cls(1,7);Node1660_AH_cls(1,8) Node1660_AH_cls(1,9) Node1660_AH_cls(1,10);...
    Node1660_AH_cls(1,11) Node1660_AH_cls(1,12) Node1660_AH_cls(1,13)];
X1660 = sqrt(0.05^2-(Node1660_AH_cls(1,3)+0.25)^2);
teta1660 = asin((Node1660_AH_cls(1,3)+0.25)/0.05);
R = [X1660/0.05 -(Node1660_AH_cls(1,3)+0.25)/0.05 0;(Node1660_AH_cls(1,3)+0.25)/0.05 X1660/0.05 0;...
    0 0 1];
Result_S_1660_cls_trans = R'*Result_S_1660_cls*R;


Result_S_1588_cls = [Node1588_AH_cls(1,5) Node1588_AH_cls(1,6) Node1588_AH_cls(1,7);Node1588_AH_cls(1,8) Node1588_AH_cls(1,9) Node1588_AH_cls(1,10);...
    Node1588_AH_cls(1,11) Node1588_AH_cls(1,12) Node1588_AH_cls(1,13)];
X1588 = sqrt(0.05^2-(Node1588_AH_cls(1,3)+0.25)^2);
teta1588 = asin((Node1588_AH_cls(1,3)+0.25)/0.05);
R = [X1588/0.05 -(Node1588_AH_cls(1,3)+0.25)/0.05 0;(Node1588_AH_cls(1,3)+0.25)/0.05 X1588/0.05 0;...
    0 0 1];
Result_S_1588_cls_trans = R'*Result_S_1588_cls*R;


Result_S_1516_cls = [Node1516_AH_cls(1,5) Node1516_AH_cls(1,6) Node1516_AH_cls(1,7);Node1516_AH_cls(1,8) Node1516_AH_cls(1,9) Node1516_AH_cls(1,10);...
    Node1516_AH_cls(1,11) Node1516_AH_cls(1,12) Node1516_AH_cls(1,13)];
X1516 = sqrt(0.05^2-(Node1516_AH_cls(1,3)+0.25)^2);
teta1516 = asin((Node1516_AH_cls(1,3)+0.25)/0.05);
R = [X1516/0.05 -(Node1516_AH_cls(1,3)+0.25)/0.05 0;(Node1516_AH_cls(1,3)+0.25)/0.05 X1516/0.05 0;...
    0 0 1];
Result_S_1516_cls_trans = R'*Result_S_1516_cls*R;


Result_S_1444_cls = [Node1444_AH_cls(1,5) Node1444_AH_cls(1,6) Node1444_AH_cls(1,7);Node1444_AH_cls(1,8) Node1444_AH_cls(1,9) Node1444_AH_cls(1,10);...
    Node1444_AH_cls(1,11) Node1444_AH_cls(1,12) Node1444_AH_cls(1,13)];
X1444 = sqrt(0.05^2-(Node1444_AH_cls(1,3)+0.25)^2);
teta1444 = asin((Node1444_AH_cls(1,3)+0.25)/0.05);
R = [X1444/0.05 -(Node1444_AH_cls(1,3)+0.25)/0.05 0;(Node1444_AH_cls(1,3)+0.25)/0.05 X1444/0.05 0;...
    0 0 1];
Result_S_1444_cls_trans = R'*Result_S_1444_cls*R;


Result_S_1372_cls = [Node1372_AH_cls(1,5) Node1372_AH_cls(1,6) Node1372_AH_cls(1,7);Node1372_AH_cls(1,8) Node1372_AH_cls(1,9) Node1372_AH_cls(1,10);...
    Node1372_AH_cls(1,11) Node1372_AH_cls(1,12) Node1372_AH_cls(1,13)];
X1372 = sqrt(0.05^2-(Node1372_AH_cls(1,3)+0.25)^2);
teta1372 = asin((Node1372_AH_cls(1,3)+0.25)/0.05);
R = [X1372/0.05 -(Node1372_AH_cls(1,3)+0.25)/0.05 0;(Node1372_AH_cls(1,3)+0.25)/0.05 X1372/0.05 0;...
    0 0 1];
Result_S_1372_cls_trans = R'*Result_S_1372_cls*R;


Result_S_1300_cls = [Node1300_AH_cls(1,5) Node1300_AH_cls(1,6) Node1300_AH_cls(1,7);Node1300_AH_cls(1,8) Node1300_AH_cls(1,9) Node1300_AH_cls(1,10);...
    Node1300_AH_cls(1,11) Node1300_AH_cls(1,12) Node1300_AH_cls(1,13)];
X1300 = sqrt(0.05^2-(Node1300_AH_cls(1,3)+0.25)^2);
teta1300 = asin((Node1300_AH_cls(1,3)+0.25)/0.05);
R = [X1300/0.05 -(Node1300_AH_cls(1,3)+0.25)/0.05 0;(Node1300_AH_cls(1,3)+0.25)/0.05 X1300/0.05 0;...
    0 0 1];
Result_S_1300_cls_trans = R'*Result_S_1300_cls*R;


Result_S_1228_cls = [Node1228_AH_cls(1,5) Node1228_AH_cls(1,6) Node1228_AH_cls(1,7);Node1228_AH_cls(1,8) Node1228_AH_cls(1,9) Node1228_AH_cls(1,10);...
    Node1228_AH_cls(1,11) Node1228_AH_cls(1,12) Node1228_AH_cls(1,13)];
X1228 = sqrt(0.05^2-(Node1228_AH_cls(1,3)+0.25)^2);
teta1228 = asin((Node1228_AH_cls(1,3)+0.25)/0.05);
R = [X1228/0.05 -(Node1228_AH_cls(1,3)+0.25)/0.05 0;(Node1228_AH_cls(1,3)+0.25)/0.05 X1228/0.05 0;...
    0 0 1];
Result_S_1228_cls_trans = R'*Result_S_1228_cls*R;


Result_S_1156_cls = [Node1156_AH_cls(1,5) Node1156_AH_cls(1,6) Node1156_AH_cls(1,7);Node1156_AH_cls(1,8) Node1156_AH_cls(1,9) Node1156_AH_cls(1,10);...
    Node1156_AH_cls(1,11) Node1156_AH_cls(1,12) Node1156_AH_cls(1,13)];
X1156 = sqrt(0.05^2-(Node1156_AH_cls(1,3)+0.25)^2);
teta1156 = asin((Node1156_AH_cls(1,3)+0.25)/0.05);
R = [X1156/0.05 -(Node1156_AH_cls(1,3)+0.25)/0.05 0;(Node1156_AH_cls(1,3)+0.25)/0.05 X1156/0.05 0;...
    0 0 1];
Result_S_1156_cls_trans = R'*Result_S_1156_cls*R;


Result_S_1016_cls = [Node1016_AH_cls(1,5) Node1016_AH_cls(1,6) Node1016_AH_cls(1,7);Node1016_AH_cls(1,8) Node1016_AH_cls(1,9) Node1016_AH_cls(1,10);...
    Node1016_AH_cls(1,11) Node1016_AH_cls(1,12) Node1016_AH_cls(1,13)];
X1016 = sqrt(0.05^2-(Node1016_AH_cls(1,3)+0.25)^2);
teta1016 = asin((Node1016_AH_cls(1,3)+0.25)/0.05);
R = [X1016/0.05 -(Node1016_AH_cls(1,3)+0.25)/0.05 0;(Node1016_AH_cls(1,3)+0.25)/0.05 X1016/0.05 0;...
    0 0 1];
Result_S_1016_cls_trans = R'*Result_S_1016_cls*R;


Result_S_1012_cls = [Node1012_AH_cls(1,5) Node1012_AH_cls(1,6) Node1012_AH_cls(1,7);Node1012_AH_cls(1,8) Node1012_AH_cls(1,9) Node1012_AH_cls(1,10);...
    Node1012_AH_cls(1,11) Node1012_AH_cls(1,12) Node1012_AH_cls(1,13)];
X1012 = sqrt(0.05^2-(Node1012_AH_cls(1,3)+0.25)^2);
teta1012 = asin((Node1012_AH_cls(1,3)+0.25)/0.05);
R = [X1012/0.05 -(Node1012_AH_cls(1,3)+0.25)/0.05 0;(Node1012_AH_cls(1,3)+0.25)/0.05 X1012/0.05 0;...
    0 0 1];
Result_S_1012_cls_trans = R'*Result_S_1012_cls*R;


Result_S_teta_teta_cls = [teta2 Result_S_2_cls_trans(2,2)/P_applied;teta6 Result_S_6_cls_trans(2,2)/P_applied;teta146 Result_S_146_cls_trans(2,2)/P_applied;...
    teta218 Result_S_218_cls_trans(2,2)/P_applied;teta290 Result_S_290_cls_trans(2,2)/P_applied;teta362 Result_S_362_cls_trans(2,2)/P_applied;...
    teta434 Result_S_434_cls_trans(2,2)/P_applied;teta506 Result_S_506_cls_trans(2,2)/P_applied;teta578 Result_S_578_cls_trans(2,2)/P_applied;...
    teta650 Result_S_650_cls_trans(2,2)/P_applied;teta722 Result_S_722_cls_trans(2,2)/P_applied;teta794 Result_S_794_cls_trans(2,2)/P_applied;...
    teta866 Result_S_866_cls_trans(2,2)/P_applied;teta938 Result_S_938_cls_trans(2,2)/P_applied;teta1876 Result_S_1876_cls_trans(2,2)/P_applied;...
    teta1804 Result_S_1804_cls_trans(2,2)/P_applied;teta1732 Result_S_1732_cls_trans(2,2)/P_applied;teta1660 Result_S_1660_cls_trans(2,2)/P_applied;...
    teta1588 Result_S_1588_cls_trans(2,2)/P_applied;teta1516 Result_S_1516_cls_trans(2,2)/P_applied;teta1444 Result_S_1444_cls_trans(2,2)/P_applied;...
    teta1372 Result_S_1372_cls_trans(2,2)/P_applied;teta1300 Result_S_1300_cls_trans(2,2)/P_applied;teta1228 Result_S_1228_cls_trans(2,2)/P_applied;...
    teta1156 Result_S_1156_cls_trans(2,2)/P_applied;teta1016 Result_S_1016_cls_trans(2,2)/P_applied;teta1012 Result_S_1012_cls_trans(2,2)/P_applied];



Result_S_2_MBC_com = [Node2_AH_MBC_com(1,5) Node2_AH_MBC_com(1,6) Node2_AH_MBC_com(1,7);Node2_AH_MBC_com(1,8) Node2_AH_MBC_com(1,9) Node2_AH_MBC_com(1,10);...
    Node2_AH_MBC_com(1,11) Node2_AH_MBC_com(1,12) Node2_AH_MBC_com(1,13)];
X2 = sqrt(0.05^2-(Node2_AH_MBC_com(1,3)+0.25)^2);
teta2 = asin((Node2_AH_MBC_com(1,3)+0.25)/0.05);
R = [X2/0.05 -(Node2_AH_MBC_com(1,3)+0.25)/0.05 0;(Node2_AH_MBC_com(1,3)+0.25)/0.05 X2/0.05 0;...
    0 0 1];
Result_S_2_MBC_com_trans = R'*Result_S_2_MBC_com*R;


Result_S_6_MBC_com = [Node6_AH_MBC_com(1,5) Node6_AH_MBC_com(1,6) Node6_AH_MBC_com(1,7);Node6_AH_MBC_com(1,8) Node6_AH_MBC_com(1,9) Node6_AH_MBC_com(1,10);...
    Node6_AH_MBC_com(1,11) Node6_AH_MBC_com(1,12) Node6_AH_MBC_com(1,13)];
X6 = sqrt(0.05^2-(Node6_AH_MBC_com(1,3)+0.25)^2);
teta6 = asin((Node6_AH_MBC_com(1,3)+0.25)/0.05);
R = [X6/0.05 -(Node6_AH_MBC_com(1,3)+0.25)/0.05 0;(Node6_AH_MBC_com(1,3)+0.25)/0.05 X6/0.05 0;...
    0 0 1];
Result_S_6_MBC_com_trans = R'*Result_S_6_MBC_com*R;


Result_S_146_MBC_com = [Node146_AH_MBC_com(1,5) Node146_AH_MBC_com(1,6) Node146_AH_MBC_com(1,7);Node146_AH_MBC_com(1,8) Node146_AH_MBC_com(1,9) Node146_AH_MBC_com(1,10);...
    Node146_AH_MBC_com(1,11) Node146_AH_MBC_com(1,12) Node146_AH_MBC_com(1,13)];
X146 = sqrt(0.05^2-(Node146_AH_MBC_com(1,3)+0.25)^2);
teta146 = asin((Node146_AH_MBC_com(1,3)+0.25)/0.05);
R = [X146/0.05 -(Node146_AH_MBC_com(1,3)+0.25)/0.05 0;(Node146_AH_MBC_com(1,3)+0.25)/0.05 X146/0.05 0;...
    0 0 1];
Result_S_146_MBC_com_trans = R'*Result_S_146_MBC_com*R;


Result_S_218_MBC_com = [Node218_AH_MBC_com(1,5) Node218_AH_MBC_com(1,6) Node218_AH_MBC_com(1,7);Node218_AH_MBC_com(1,8) Node218_AH_MBC_com(1,9) Node218_AH_MBC_com(1,10);...
    Node218_AH_MBC_com(1,11) Node218_AH_MBC_com(1,12) Node218_AH_MBC_com(1,13)];
X218 = sqrt(0.05^2-(Node218_AH_MBC_com(1,3)+0.25)^2);
teta218 = asin((Node218_AH_MBC_com(1,3)+0.25)/0.05);
R = [X218/0.05 -(Node218_AH_MBC_com(1,3)+0.25)/0.05 0;(Node218_AH_MBC_com(1,3)+0.25)/0.05 X218/0.05 0;...
    0 0 1];
Result_S_218_MBC_com_trans = R'*Result_S_218_MBC_com*R;


Result_S_290_MBC_com = [Node290_AH_MBC_com(1,5) Node290_AH_MBC_com(1,6) Node290_AH_MBC_com(1,7);Node290_AH_MBC_com(1,8) Node290_AH_MBC_com(1,9) Node290_AH_MBC_com(1,10);...
    Node290_AH_MBC_com(1,11) Node290_AH_MBC_com(1,12) Node290_AH_MBC_com(1,13)];
X290 = sqrt(0.05^2-(Node290_AH_MBC_com(1,3)+0.25)^2);
teta290 = asin((Node290_AH_MBC_com(1,3)+0.25)/0.05);
R = [X290/0.05 -(Node290_AH_MBC_com(1,3)+0.25)/0.05 0;(Node290_AH_MBC_com(1,3)+0.25)/0.05 X290/0.05 0;...
    0 0 1];
Result_S_290_MBC_com_trans = R'*Result_S_290_MBC_com*R;


Result_S_362_MBC_com = [Node362_AH_MBC_com(1,5) Node362_AH_MBC_com(1,6) Node362_AH_MBC_com(1,7);Node362_AH_MBC_com(1,8) Node362_AH_MBC_com(1,9) Node362_AH_MBC_com(1,10);...
    Node362_AH_MBC_com(1,11) Node362_AH_MBC_com(1,12) Node362_AH_MBC_com(1,13)];
X362 = sqrt(0.05^2-(Node362_AH_MBC_com(1,3)+0.25)^2);
teta362 = asin((Node362_AH_MBC_com(1,3)+0.25)/0.05);
R = [X362/0.05 -(Node362_AH_MBC_com(1,3)+0.25)/0.05 0;(Node362_AH_MBC_com(1,3)+0.25)/0.05 X362/0.05 0;...
    0 0 1];
Result_S_362_MBC_com_trans = R'*Result_S_362_MBC_com*R;


Result_S_434_MBC_com = [Node434_AH_MBC_com(1,5) Node434_AH_MBC_com(1,6) Node434_AH_MBC_com(1,7);Node434_AH_MBC_com(1,8) Node434_AH_MBC_com(1,9) Node434_AH_MBC_com(1,10);...
    Node434_AH_MBC_com(1,11) Node434_AH_MBC_com(1,12) Node434_AH_MBC_com(1,13)];
X434 = sqrt(0.05^2-(Node434_AH_MBC_com(1,3)+0.25)^2);
teta434 = asin((Node434_AH_MBC_com(1,3)+0.25)/0.05);
R = [X434/0.05 -(Node434_AH_MBC_com(1,3)+0.25)/0.05 0;(Node434_AH_MBC_com(1,3)+0.25)/0.05 X434/0.05 0;...
    0 0 1];
Result_S_434_MBC_com_trans = R'*Result_S_434_MBC_com*R;


Result_S_506_MBC_com = [Node506_AH_MBC_com(1,5) Node506_AH_MBC_com(1,6) Node506_AH_MBC_com(1,7);Node506_AH_MBC_com(1,8) Node506_AH_MBC_com(1,9) Node506_AH_MBC_com(1,10);...
    Node506_AH_MBC_com(1,11) Node506_AH_MBC_com(1,12) Node506_AH_MBC_com(1,13)];
X506 = sqrt(0.05^2-(Node506_AH_MBC_com(1,3)+0.25)^2);
teta506 = asin((Node506_AH_MBC_com(1,3)+0.25)/0.05);
R = [X506/0.05 -(Node506_AH_MBC_com(1,3)+0.25)/0.05 0;(Node506_AH_MBC_com(1,3)+0.25)/0.05 X506/0.05 0;...
    0 0 1];
Result_S_506_MBC_com_trans = R'*Result_S_506_MBC_com*R;


Result_S_578_MBC_com = [Node578_AH_MBC_com(1,5) Node578_AH_MBC_com(1,6) Node578_AH_MBC_com(1,7);Node578_AH_MBC_com(1,8) Node578_AH_MBC_com(1,9) Node578_AH_MBC_com(1,10);...
    Node578_AH_MBC_com(1,11) Node578_AH_MBC_com(1,12) Node578_AH_MBC_com(1,13)];
X578 = sqrt(0.05^2-(Node578_AH_MBC_com(1,3)+0.25)^2);
teta578 = asin((Node578_AH_MBC_com(1,3)+0.25)/0.05);
R = [X578/0.05 -(Node578_AH_MBC_com(1,3)+0.25)/0.05 0;(Node578_AH_MBC_com(1,3)+0.25)/0.05 X578/0.05 0;...
    0 0 1];
Result_S_578_MBC_com_trans = R'*Result_S_578_MBC_com*R;


Result_S_650_MBC_com = [Node650_AH_MBC_com(1,5) Node650_AH_MBC_com(1,6) Node650_AH_MBC_com(1,7);Node650_AH_MBC_com(1,8) Node650_AH_MBC_com(1,9) Node650_AH_MBC_com(1,10);...
    Node650_AH_MBC_com(1,11) Node650_AH_MBC_com(1,12) Node650_AH_MBC_com(1,13)];
X650 = sqrt(0.05^2-(Node650_AH_MBC_com(1,3)+0.25)^2);
teta650 = asin((Node650_AH_MBC_com(1,3)+0.25)/0.05);
R = [X650/0.05 -(Node650_AH_MBC_com(1,3)+0.25)/0.05 0;(Node650_AH_MBC_com(1,3)+0.25)/0.05 X650/0.05 0;...
    0 0 1];
Result_S_650_MBC_com_trans = R'*Result_S_650_MBC_com*R;


Result_S_722_MBC_com = [Node722_AH_MBC_com(1,5) Node722_AH_MBC_com(1,6) Node722_AH_MBC_com(1,7);Node722_AH_MBC_com(1,8) Node722_AH_MBC_com(1,9) Node722_AH_MBC_com(1,10);...
    Node722_AH_MBC_com(1,11) Node722_AH_MBC_com(1,12) Node722_AH_MBC_com(1,13)];
X722 = sqrt(0.05^2-(Node722_AH_MBC_com(1,3)+0.25)^2);
teta722 = asin((Node722_AH_MBC_com(1,3)+0.25)/0.05);
R = [X722/0.05 -(Node722_AH_MBC_com(1,3)+0.25)/0.05 0;(Node722_AH_MBC_com(1,3)+0.25)/0.05 X722/0.05 0;...
    0 0 1];
Result_S_722_MBC_com_trans = R'*Result_S_722_MBC_com*R;


Result_S_794_MBC_com = [Node794_AH_MBC_com(1,5) Node794_AH_MBC_com(1,6) Node794_AH_MBC_com(1,7);Node794_AH_MBC_com(1,8) Node794_AH_MBC_com(1,9) Node794_AH_MBC_com(1,10);...
    Node794_AH_MBC_com(1,11) Node794_AH_MBC_com(1,12) Node794_AH_MBC_com(1,13)];
X794 = sqrt(0.05^2-(Node794_AH_MBC_com(1,3)+0.25)^2);
teta794 = asin((Node794_AH_MBC_com(1,3)+0.25)/0.05);
R = [X794/0.05 -(Node794_AH_MBC_com(1,3)+0.25)/0.05 0;(Node794_AH_MBC_com(1,3)+0.25)/0.05 X794/0.05 0;...
    0 0 1];
Result_S_794_MBC_com_trans = R'*Result_S_794_MBC_com*R;


Result_S_866_MBC_com = [Node866_AH_MBC_com(1,5) Node866_AH_MBC_com(1,6) Node866_AH_MBC_com(1,7);Node866_AH_MBC_com(1,8) Node866_AH_MBC_com(1,9) Node866_AH_MBC_com(1,10);...
    Node866_AH_MBC_com(1,11) Node866_AH_MBC_com(1,12) Node866_AH_MBC_com(1,13)];
X866 = sqrt(0.05^2-(Node866_AH_MBC_com(1,3)+0.25)^2);
teta866 = asin((Node866_AH_MBC_com(1,3)+0.25)/0.05);
R = [X866/0.05 -(Node866_AH_MBC_com(1,3)+0.25)/0.05 0;(Node866_AH_MBC_com(1,3)+0.25)/0.05 X866/0.05 0;...
    0 0 1];
Result_S_866_MBC_com_trans = R'*Result_S_866_MBC_com*R;


Result_S_938_MBC_com = [Node938_AH_MBC_com(1,5) Node938_AH_MBC_com(1,6) Node938_AH_MBC_com(1,7);Node938_AH_MBC_com(1,8) Node938_AH_MBC_com(1,9) Node938_AH_MBC_com(1,10);...
    Node938_AH_MBC_com(1,11) Node938_AH_MBC_com(1,12) Node938_AH_MBC_com(1,13)];
X938 = sqrt(0.05^2-(Node938_AH_MBC_com(1,3)+0.25)^2);
teta938 = asin((Node938_AH_MBC_com(1,3)+0.25)/0.05);
R = [X938/0.05 -(Node938_AH_MBC_com(1,3)+0.25)/0.05 0;(Node938_AH_MBC_com(1,3)+0.25)/0.05 X938/0.05 0;...
    0 0 1];
Result_S_938_MBC_com_trans = R'*Result_S_938_MBC_com*R;


Result_S_1876_MBC_com = [Node1876_AH_MBC_com(1,5) Node1876_AH_MBC_com(1,6) Node1876_AH_MBC_com(1,7);Node1876_AH_MBC_com(1,8) Node1876_AH_MBC_com(1,9) Node1876_AH_MBC_com(1,10);...
    Node1876_AH_MBC_com(1,11) Node1876_AH_MBC_com(1,12) Node1876_AH_MBC_com(1,13)];
X1876 = sqrt(0.05^2-(Node1876_AH_MBC_com(1,3)+0.25)^2);
teta1876 = asin((Node1876_AH_MBC_com(1,3)+0.25)/0.05);
R = [X1876/0.05 -(Node1876_AH_MBC_com(1,3)+0.25)/0.05 0;(Node1876_AH_MBC_com(1,3)+0.25)/0.05 X1876/0.05 0;...
    0 0 1];
Result_S_1876_MBC_com_trans = R'*Result_S_1876_MBC_com*R;


Result_S_1804_MBC_com = [Node1804_AH_MBC_com(1,5) Node1804_AH_MBC_com(1,6) Node1804_AH_MBC_com(1,7);Node1804_AH_MBC_com(1,8) Node1804_AH_MBC_com(1,9) Node1804_AH_MBC_com(1,10);...
    Node1804_AH_MBC_com(1,11) Node1804_AH_MBC_com(1,12) Node1804_AH_MBC_com(1,13)];
X1804 = sqrt(0.05^2-(Node1804_AH_MBC_com(1,3)+0.25)^2);
teta1804 = asin((Node1804_AH_MBC_com(1,3)+0.25)/0.05);
R = [X1804/0.05 -(Node1804_AH_MBC_com(1,3)+0.25)/0.05 0;(Node1804_AH_MBC_com(1,3)+0.25)/0.05 X1804/0.05 0;...
    0 0 1];
Result_S_1804_MBC_com_trans = R'*Result_S_1804_MBC_com*R;


Result_S_1732_MBC_com = [Node1732_AH_MBC_com(1,5) Node1732_AH_MBC_com(1,6) Node1732_AH_MBC_com(1,7);Node1732_AH_MBC_com(1,8) Node1732_AH_MBC_com(1,9) Node1732_AH_MBC_com(1,10);...
    Node1732_AH_MBC_com(1,11) Node1732_AH_MBC_com(1,12) Node1732_AH_MBC_com(1,13)];
X1732 = sqrt(0.05^2-(Node1732_AH_MBC_com(1,3)+0.25)^2);
teta1732 = asin((Node1732_AH_MBC_com(1,3)+0.25)/0.05);
R = [X1732/0.05 -(Node1732_AH_MBC_com(1,3)+0.25)/0.05 0;(Node1732_AH_MBC_com(1,3)+0.25)/0.05 X1732/0.05 0;...
    0 0 1];
Result_S_1732_MBC_com_trans = R'*Result_S_1732_MBC_com*R;


Result_S_1660_MBC_com = [Node1660_AH_MBC_com(1,5) Node1660_AH_MBC_com(1,6) Node1660_AH_MBC_com(1,7);Node1660_AH_MBC_com(1,8) Node1660_AH_MBC_com(1,9) Node1660_AH_MBC_com(1,10);...
    Node1660_AH_MBC_com(1,11) Node1660_AH_MBC_com(1,12) Node1660_AH_MBC_com(1,13)];
X1660 = sqrt(0.05^2-(Node1660_AH_MBC_com(1,3)+0.25)^2);
teta1660 = asin((Node1660_AH_MBC_com(1,3)+0.25)/0.05);
R = [X1660/0.05 -(Node1660_AH_MBC_com(1,3)+0.25)/0.05 0;(Node1660_AH_MBC_com(1,3)+0.25)/0.05 X1660/0.05 0;...
    0 0 1];
Result_S_1660_MBC_com_trans = R'*Result_S_1660_MBC_com*R;


Result_S_1588_MBC_com = [Node1588_AH_MBC_com(1,5) Node1588_AH_MBC_com(1,6) Node1588_AH_MBC_com(1,7);Node1588_AH_MBC_com(1,8) Node1588_AH_MBC_com(1,9) Node1588_AH_MBC_com(1,10);...
    Node1588_AH_MBC_com(1,11) Node1588_AH_MBC_com(1,12) Node1588_AH_MBC_com(1,13)];
X1588 = sqrt(0.05^2-(Node1588_AH_MBC_com(1,3)+0.25)^2);
teta1588 = asin((Node1588_AH_MBC_com(1,3)+0.25)/0.05);
R = [X1588/0.05 -(Node1588_AH_MBC_com(1,3)+0.25)/0.05 0;(Node1588_AH_MBC_com(1,3)+0.25)/0.05 X1588/0.05 0;...
    0 0 1];
Result_S_1588_MBC_com_trans = R'*Result_S_1588_MBC_com*R;


Result_S_1516_MBC_com = [Node1516_AH_MBC_com(1,5) Node1516_AH_MBC_com(1,6) Node1516_AH_MBC_com(1,7);Node1516_AH_MBC_com(1,8) Node1516_AH_MBC_com(1,9) Node1516_AH_MBC_com(1,10);...
    Node1516_AH_MBC_com(1,11) Node1516_AH_MBC_com(1,12) Node1516_AH_MBC_com(1,13)];
X1516 = sqrt(0.05^2-(Node1516_AH_MBC_com(1,3)+0.25)^2);
teta1516 = asin((Node1516_AH_MBC_com(1,3)+0.25)/0.05);
R = [X1516/0.05 -(Node1516_AH_MBC_com(1,3)+0.25)/0.05 0;(Node1516_AH_MBC_com(1,3)+0.25)/0.05 X1516/0.05 0;...
    0 0 1];
Result_S_1516_MBC_com_trans = R'*Result_S_1516_MBC_com*R;


Result_S_1444_MBC_com = [Node1444_AH_MBC_com(1,5) Node1444_AH_MBC_com(1,6) Node1444_AH_MBC_com(1,7);Node1444_AH_MBC_com(1,8) Node1444_AH_MBC_com(1,9) Node1444_AH_MBC_com(1,10);...
    Node1444_AH_MBC_com(1,11) Node1444_AH_MBC_com(1,12) Node1444_AH_MBC_com(1,13)];
X1444 = sqrt(0.05^2-(Node1444_AH_MBC_com(1,3)+0.25)^2);
teta1444 = asin((Node1444_AH_MBC_com(1,3)+0.25)/0.05);
R = [X1444/0.05 -(Node1444_AH_MBC_com(1,3)+0.25)/0.05 0;(Node1444_AH_MBC_com(1,3)+0.25)/0.05 X1444/0.05 0;...
    0 0 1];
Result_S_1444_MBC_com_trans = R'*Result_S_1444_MBC_com*R;


Result_S_1372_MBC_com = [Node1372_AH_MBC_com(1,5) Node1372_AH_MBC_com(1,6) Node1372_AH_MBC_com(1,7);Node1372_AH_MBC_com(1,8) Node1372_AH_MBC_com(1,9) Node1372_AH_MBC_com(1,10);...
    Node1372_AH_MBC_com(1,11) Node1372_AH_MBC_com(1,12) Node1372_AH_MBC_com(1,13)];
X1372 = sqrt(0.05^2-(Node1372_AH_MBC_com(1,3)+0.25)^2);
teta1372 = asin((Node1372_AH_MBC_com(1,3)+0.25)/0.05);
R = [X1372/0.05 -(Node1372_AH_MBC_com(1,3)+0.25)/0.05 0;(Node1372_AH_MBC_com(1,3)+0.25)/0.05 X1372/0.05 0;...
    0 0 1];
Result_S_1372_MBC_com_trans = R'*Result_S_1372_MBC_com*R;


Result_S_1300_MBC_com = [Node1300_AH_MBC_com(1,5) Node1300_AH_MBC_com(1,6) Node1300_AH_MBC_com(1,7);Node1300_AH_MBC_com(1,8) Node1300_AH_MBC_com(1,9) Node1300_AH_MBC_com(1,10);...
    Node1300_AH_MBC_com(1,11) Node1300_AH_MBC_com(1,12) Node1300_AH_MBC_com(1,13)];
X1300 = sqrt(0.05^2-(Node1300_AH_MBC_com(1,3)+0.25)^2);
teta1300 = asin((Node1300_AH_MBC_com(1,3)+0.25)/0.05);
R = [X1300/0.05 -(Node1300_AH_MBC_com(1,3)+0.25)/0.05 0;(Node1300_AH_MBC_com(1,3)+0.25)/0.05 X1300/0.05 0;...
    0 0 1];
Result_S_1300_MBC_com_trans = R'*Result_S_1300_MBC_com*R;


Result_S_1228_MBC_com = [Node1228_AH_MBC_com(1,5) Node1228_AH_MBC_com(1,6) Node1228_AH_MBC_com(1,7);Node1228_AH_MBC_com(1,8) Node1228_AH_MBC_com(1,9) Node1228_AH_MBC_com(1,10);...
    Node1228_AH_MBC_com(1,11) Node1228_AH_MBC_com(1,12) Node1228_AH_MBC_com(1,13)];
X1228 = sqrt(0.05^2-(Node1228_AH_MBC_com(1,3)+0.25)^2);
teta1228 = asin((Node1228_AH_MBC_com(1,3)+0.25)/0.05);
R = [X1228/0.05 -(Node1228_AH_MBC_com(1,3)+0.25)/0.05 0;(Node1228_AH_MBC_com(1,3)+0.25)/0.05 X1228/0.05 0;...
    0 0 1];
Result_S_1228_MBC_com_trans = R'*Result_S_1228_MBC_com*R;


Result_S_1156_MBC_com = [Node1156_AH_MBC_com(1,5) Node1156_AH_MBC_com(1,6) Node1156_AH_MBC_com(1,7);Node1156_AH_MBC_com(1,8) Node1156_AH_MBC_com(1,9) Node1156_AH_MBC_com(1,10);...
    Node1156_AH_MBC_com(1,11) Node1156_AH_MBC_com(1,12) Node1156_AH_MBC_com(1,13)];
X1156 = sqrt(0.05^2-(Node1156_AH_MBC_com(1,3)+0.25)^2);
teta1156 = asin((Node1156_AH_MBC_com(1,3)+0.25)/0.05);
R = [X1156/0.05 -(Node1156_AH_MBC_com(1,3)+0.25)/0.05 0;(Node1156_AH_MBC_com(1,3)+0.25)/0.05 X1156/0.05 0;...
    0 0 1];
Result_S_1156_MBC_com_trans = R'*Result_S_1156_MBC_com*R;


Result_S_1016_MBC_com = [Node1016_AH_MBC_com(1,5) Node1016_AH_MBC_com(1,6) Node1016_AH_MBC_com(1,7);Node1016_AH_MBC_com(1,8) Node1016_AH_MBC_com(1,9) Node1016_AH_MBC_com(1,10);...
    Node1016_AH_MBC_com(1,11) Node1016_AH_MBC_com(1,12) Node1016_AH_MBC_com(1,13)];
X1016 = sqrt(0.05^2-(Node1016_AH_MBC_com(1,3)+0.25)^2);
teta1016 = asin((Node1016_AH_MBC_com(1,3)+0.25)/0.05);
R = [X1016/0.05 -(Node1016_AH_MBC_com(1,3)+0.25)/0.05 0;(Node1016_AH_MBC_com(1,3)+0.25)/0.05 X1016/0.05 0;...
    0 0 1];
Result_S_1016_MBC_com_trans = R'*Result_S_1016_MBC_com*R;


Result_S_1012_MBC_com = [Node1012_AH_MBC_com(1,5) Node1012_AH_MBC_com(1,6) Node1012_AH_MBC_com(1,7);Node1012_AH_MBC_com(1,8) Node1012_AH_MBC_com(1,9) Node1012_AH_MBC_com(1,10);...
    Node1012_AH_MBC_com(1,11) Node1012_AH_MBC_com(1,12) Node1012_AH_MBC_com(1,13)];
X1012 = sqrt(0.05^2-(Node1012_AH_MBC_com(1,3)+0.25)^2);
teta1012 = asin((Node1012_AH_MBC_com(1,3)+0.25)/0.05);
R = [X1012/0.05 -(Node1012_AH_MBC_com(1,3)+0.25)/0.05 0;(Node1012_AH_MBC_com(1,3)+0.25)/0.05 X1012/0.05 0;...
    0 0 1];
Result_S_1012_MBC_com_trans = R'*Result_S_1012_MBC_com*R;


Result_S_teta_teta_MBC_com = [teta2 Result_S_2_MBC_com_trans(2,2)/P_applied;teta6 Result_S_6_MBC_com_trans(2,2)/P_applied;teta146 Result_S_146_MBC_com_trans(2,2)/P_applied;...
    teta218 Result_S_218_MBC_com_trans(2,2)/P_applied;teta290 Result_S_290_MBC_com_trans(2,2)/P_applied;teta362 Result_S_362_MBC_com_trans(2,2)/P_applied;...
    teta434 Result_S_434_MBC_com_trans(2,2)/P_applied;teta506 Result_S_506_MBC_com_trans(2,2)/P_applied;teta578 Result_S_578_MBC_com_trans(2,2)/P_applied;...
    teta650 Result_S_650_MBC_com_trans(2,2)/P_applied;teta722 Result_S_722_MBC_com_trans(2,2)/P_applied;teta794 Result_S_794_MBC_com_trans(2,2)/P_applied;...
    teta866 Result_S_866_MBC_com_trans(2,2)/P_applied;teta938 Result_S_938_MBC_com_trans(2,2)/P_applied;teta1876 Result_S_1876_MBC_com_trans(2,2)/P_applied;...
    teta1804 Result_S_1804_MBC_com_trans(2,2)/P_applied;teta1732 Result_S_1732_MBC_com_trans(2,2)/P_applied;teta1660 Result_S_1660_MBC_com_trans(2,2)/P_applied;...
    teta1588 Result_S_1588_MBC_com_trans(2,2)/P_applied;teta1516 Result_S_1516_MBC_com_trans(2,2)/P_applied;teta1444 Result_S_1444_MBC_com_trans(2,2)/P_applied;...
    teta1372 Result_S_1372_MBC_com_trans(2,2)/P_applied;teta1300 Result_S_1300_MBC_com_trans(2,2)/P_applied;teta1228 Result_S_1228_MBC_com_trans(2,2)/P_applied;...
    teta1156 Result_S_1156_MBC_com_trans(2,2)/P_applied;teta1016 Result_S_1016_MBC_com_trans(2,2)/P_applied;teta1012 Result_S_1012_MBC_com_trans(2,2)/P_applied];



Result_S_2_MBC_rot = [Node2_AH_MBC_rot(1,5) Node2_AH_MBC_rot(1,6) Node2_AH_MBC_rot(1,7);Node2_AH_MBC_rot(1,8) Node2_AH_MBC_rot(1,9) Node2_AH_MBC_rot(1,10);...
    Node2_AH_MBC_rot(1,11) Node2_AH_MBC_rot(1,12) Node2_AH_MBC_rot(1,13)];
X2 = sqrt(0.05^2-(Node2_AH_MBC_rot(1,3)+0.25)^2);
teta2 = asin((Node2_AH_MBC_rot(1,3)+0.25)/0.05);
R = [X2/0.05 -(Node2_AH_MBC_rot(1,3)+0.25)/0.05 0;(Node2_AH_MBC_rot(1,3)+0.25)/0.05 X2/0.05 0;...
    0 0 1];
Result_S_2_MBC_rot_trans = R'*Result_S_2_MBC_rot*R;


Result_S_6_MBC_rot = [Node6_AH_MBC_rot(1,5) Node6_AH_MBC_rot(1,6) Node6_AH_MBC_rot(1,7);Node6_AH_MBC_rot(1,8) Node6_AH_MBC_rot(1,9) Node6_AH_MBC_rot(1,10);...
    Node6_AH_MBC_rot(1,11) Node6_AH_MBC_rot(1,12) Node6_AH_MBC_rot(1,13)];
X6 = sqrt(0.05^2-(Node6_AH_MBC_rot(1,3)+0.25)^2);
teta6 = asin((Node6_AH_MBC_rot(1,3)+0.25)/0.05);
R = [X6/0.05 -(Node6_AH_MBC_rot(1,3)+0.25)/0.05 0;(Node6_AH_MBC_rot(1,3)+0.25)/0.05 X6/0.05 0;...
    0 0 1];
Result_S_6_MBC_rot_trans = R'*Result_S_6_MBC_rot*R;


Result_S_146_MBC_rot = [Node146_AH_MBC_rot(1,5) Node146_AH_MBC_rot(1,6) Node146_AH_MBC_rot(1,7);Node146_AH_MBC_rot(1,8) Node146_AH_MBC_rot(1,9) Node146_AH_MBC_rot(1,10);...
    Node146_AH_MBC_rot(1,11) Node146_AH_MBC_rot(1,12) Node146_AH_MBC_rot(1,13)];
X146 = sqrt(0.05^2-(Node146_AH_MBC_rot(1,3)+0.25)^2);
teta146 = asin((Node146_AH_MBC_rot(1,3)+0.25)/0.05);
R = [X146/0.05 -(Node146_AH_MBC_rot(1,3)+0.25)/0.05 0;(Node146_AH_MBC_rot(1,3)+0.25)/0.05 X146/0.05 0;...
    0 0 1];
Result_S_146_MBC_rot_trans = R'*Result_S_146_MBC_rot*R;


Result_S_218_MBC_rot = [Node218_AH_MBC_rot(1,5) Node218_AH_MBC_rot(1,6) Node218_AH_MBC_rot(1,7);Node218_AH_MBC_rot(1,8) Node218_AH_MBC_rot(1,9) Node218_AH_MBC_rot(1,10);...
    Node218_AH_MBC_rot(1,11) Node218_AH_MBC_rot(1,12) Node218_AH_MBC_rot(1,13)];
X218 = sqrt(0.05^2-(Node218_AH_MBC_rot(1,3)+0.25)^2);
teta218 = asin((Node218_AH_MBC_rot(1,3)+0.25)/0.05);
R = [X218/0.05 -(Node218_AH_MBC_rot(1,3)+0.25)/0.05 0;(Node218_AH_MBC_rot(1,3)+0.25)/0.05 X218/0.05 0;...
    0 0 1];
Result_S_218_MBC_rot_trans = R'*Result_S_218_MBC_rot*R;


Result_S_290_MBC_rot = [Node290_AH_MBC_rot(1,5) Node290_AH_MBC_rot(1,6) Node290_AH_MBC_rot(1,7);Node290_AH_MBC_rot(1,8) Node290_AH_MBC_rot(1,9) Node290_AH_MBC_rot(1,10);...
    Node290_AH_MBC_rot(1,11) Node290_AH_MBC_rot(1,12) Node290_AH_MBC_rot(1,13)];
X290 = sqrt(0.05^2-(Node290_AH_MBC_rot(1,3)+0.25)^2);
teta290 = asin((Node290_AH_MBC_rot(1,3)+0.25)/0.05);
R = [X290/0.05 -(Node290_AH_MBC_rot(1,3)+0.25)/0.05 0;(Node290_AH_MBC_rot(1,3)+0.25)/0.05 X290/0.05 0;...
    0 0 1];
Result_S_290_MBC_rot_trans = R'*Result_S_290_MBC_rot*R;


Result_S_362_MBC_rot = [Node362_AH_MBC_rot(1,5) Node362_AH_MBC_rot(1,6) Node362_AH_MBC_rot(1,7);Node362_AH_MBC_rot(1,8) Node362_AH_MBC_rot(1,9) Node362_AH_MBC_rot(1,10);...
    Node362_AH_MBC_rot(1,11) Node362_AH_MBC_rot(1,12) Node362_AH_MBC_rot(1,13)];
X362 = sqrt(0.05^2-(Node362_AH_MBC_rot(1,3)+0.25)^2);
teta362 = asin((Node362_AH_MBC_rot(1,3)+0.25)/0.05);
R = [X362/0.05 -(Node362_AH_MBC_rot(1,3)+0.25)/0.05 0;(Node362_AH_MBC_rot(1,3)+0.25)/0.05 X362/0.05 0;...
    0 0 1];
Result_S_362_MBC_rot_trans = R'*Result_S_362_MBC_rot*R;


Result_S_434_MBC_rot = [Node434_AH_MBC_rot(1,5) Node434_AH_MBC_rot(1,6) Node434_AH_MBC_rot(1,7);Node434_AH_MBC_rot(1,8) Node434_AH_MBC_rot(1,9) Node434_AH_MBC_rot(1,10);...
    Node434_AH_MBC_rot(1,11) Node434_AH_MBC_rot(1,12) Node434_AH_MBC_rot(1,13)];
X434 = sqrt(0.05^2-(Node434_AH_MBC_rot(1,3)+0.25)^2);
teta434 = asin((Node434_AH_MBC_rot(1,3)+0.25)/0.05);
R = [X434/0.05 -(Node434_AH_MBC_rot(1,3)+0.25)/0.05 0;(Node434_AH_MBC_rot(1,3)+0.25)/0.05 X434/0.05 0;...
    0 0 1];
Result_S_434_MBC_rot_trans = R'*Result_S_434_MBC_rot*R;


Result_S_506_MBC_rot = [Node506_AH_MBC_rot(1,5) Node506_AH_MBC_rot(1,6) Node506_AH_MBC_rot(1,7);Node506_AH_MBC_rot(1,8) Node506_AH_MBC_rot(1,9) Node506_AH_MBC_rot(1,10);...
    Node506_AH_MBC_rot(1,11) Node506_AH_MBC_rot(1,12) Node506_AH_MBC_rot(1,13)];
X506 = sqrt(0.05^2-(Node506_AH_MBC_rot(1,3)+0.25)^2);
teta506 = asin((Node506_AH_MBC_rot(1,3)+0.25)/0.05);
R = [X506/0.05 -(Node506_AH_MBC_rot(1,3)+0.25)/0.05 0;(Node506_AH_MBC_rot(1,3)+0.25)/0.05 X506/0.05 0;...
    0 0 1];
Result_S_506_MBC_rot_trans = R'*Result_S_506_MBC_rot*R;


Result_S_578_MBC_rot = [Node578_AH_MBC_rot(1,5) Node578_AH_MBC_rot(1,6) Node578_AH_MBC_rot(1,7);Node578_AH_MBC_rot(1,8) Node578_AH_MBC_rot(1,9) Node578_AH_MBC_rot(1,10);...
    Node578_AH_MBC_rot(1,11) Node578_AH_MBC_rot(1,12) Node578_AH_MBC_rot(1,13)];
X578 = sqrt(0.05^2-(Node578_AH_MBC_rot(1,3)+0.25)^2);
teta578 = asin((Node578_AH_MBC_rot(1,3)+0.25)/0.05);
R = [X578/0.05 -(Node578_AH_MBC_rot(1,3)+0.25)/0.05 0;(Node578_AH_MBC_rot(1,3)+0.25)/0.05 X578/0.05 0;...
    0 0 1];
Result_S_578_MBC_rot_trans = R'*Result_S_578_MBC_rot*R;


Result_S_650_MBC_rot = [Node650_AH_MBC_rot(1,5) Node650_AH_MBC_rot(1,6) Node650_AH_MBC_rot(1,7);Node650_AH_MBC_rot(1,8) Node650_AH_MBC_rot(1,9) Node650_AH_MBC_rot(1,10);...
    Node650_AH_MBC_rot(1,11) Node650_AH_MBC_rot(1,12) Node650_AH_MBC_rot(1,13)];
X650 = sqrt(0.05^2-(Node650_AH_MBC_rot(1,3)+0.25)^2);
teta650 = asin((Node650_AH_MBC_rot(1,3)+0.25)/0.05);
R = [X650/0.05 -(Node650_AH_MBC_rot(1,3)+0.25)/0.05 0;(Node650_AH_MBC_rot(1,3)+0.25)/0.05 X650/0.05 0;...
    0 0 1];
Result_S_650_MBC_rot_trans = R'*Result_S_650_MBC_rot*R;


Result_S_722_MBC_rot = [Node722_AH_MBC_rot(1,5) Node722_AH_MBC_rot(1,6) Node722_AH_MBC_rot(1,7);Node722_AH_MBC_rot(1,8) Node722_AH_MBC_rot(1,9) Node722_AH_MBC_rot(1,10);...
    Node722_AH_MBC_rot(1,11) Node722_AH_MBC_rot(1,12) Node722_AH_MBC_rot(1,13)];
X722 = sqrt(0.05^2-(Node722_AH_MBC_rot(1,3)+0.25)^2);
teta722 = asin((Node722_AH_MBC_rot(1,3)+0.25)/0.05);
R = [X722/0.05 -(Node722_AH_MBC_rot(1,3)+0.25)/0.05 0;(Node722_AH_MBC_rot(1,3)+0.25)/0.05 X722/0.05 0;...
    0 0 1];
Result_S_722_MBC_rot_trans = R'*Result_S_722_MBC_rot*R;


Result_S_794_MBC_rot = [Node794_AH_MBC_rot(1,5) Node794_AH_MBC_rot(1,6) Node794_AH_MBC_rot(1,7);Node794_AH_MBC_rot(1,8) Node794_AH_MBC_rot(1,9) Node794_AH_MBC_rot(1,10);...
    Node794_AH_MBC_rot(1,11) Node794_AH_MBC_rot(1,12) Node794_AH_MBC_rot(1,13)];
X794 = sqrt(0.05^2-(Node794_AH_MBC_rot(1,3)+0.25)^2);
teta794 = asin((Node794_AH_MBC_rot(1,3)+0.25)/0.05);
R = [X794/0.05 -(Node794_AH_MBC_rot(1,3)+0.25)/0.05 0;(Node794_AH_MBC_rot(1,3)+0.25)/0.05 X794/0.05 0;...
    0 0 1];
Result_S_794_MBC_rot_trans = R'*Result_S_794_MBC_rot*R;


Result_S_866_MBC_rot = [Node866_AH_MBC_rot(1,5) Node866_AH_MBC_rot(1,6) Node866_AH_MBC_rot(1,7);Node866_AH_MBC_rot(1,8) Node866_AH_MBC_rot(1,9) Node866_AH_MBC_rot(1,10);...
    Node866_AH_MBC_rot(1,11) Node866_AH_MBC_rot(1,12) Node866_AH_MBC_rot(1,13)];
X866 = sqrt(0.05^2-(Node866_AH_MBC_rot(1,3)+0.25)^2);
teta866 = asin((Node866_AH_MBC_rot(1,3)+0.25)/0.05);
R = [X866/0.05 -(Node866_AH_MBC_rot(1,3)+0.25)/0.05 0;(Node866_AH_MBC_rot(1,3)+0.25)/0.05 X866/0.05 0;...
    0 0 1];
Result_S_866_MBC_rot_trans = R'*Result_S_866_MBC_rot*R;


Result_S_938_MBC_rot = [Node938_AH_MBC_rot(1,5) Node938_AH_MBC_rot(1,6) Node938_AH_MBC_rot(1,7);Node938_AH_MBC_rot(1,8) Node938_AH_MBC_rot(1,9) Node938_AH_MBC_rot(1,10);...
    Node938_AH_MBC_rot(1,11) Node938_AH_MBC_rot(1,12) Node938_AH_MBC_rot(1,13)];
X938 = sqrt(0.05^2-(Node938_AH_MBC_rot(1,3)+0.25)^2);
teta938 = asin((Node938_AH_MBC_rot(1,3)+0.25)/0.05);
R = [X938/0.05 -(Node938_AH_MBC_rot(1,3)+0.25)/0.05 0;(Node938_AH_MBC_rot(1,3)+0.25)/0.05 X938/0.05 0;...
    0 0 1];
Result_S_938_MBC_rot_trans = R'*Result_S_938_MBC_rot*R;


Result_S_1876_MBC_rot = [Node1876_AH_MBC_rot(1,5) Node1876_AH_MBC_rot(1,6) Node1876_AH_MBC_rot(1,7);Node1876_AH_MBC_rot(1,8) Node1876_AH_MBC_rot(1,9) Node1876_AH_MBC_rot(1,10);...
    Node1876_AH_MBC_rot(1,11) Node1876_AH_MBC_rot(1,12) Node1876_AH_MBC_rot(1,13)];
X1876 = sqrt(0.05^2-(Node1876_AH_MBC_rot(1,3)+0.25)^2);
teta1876 = asin((Node1876_AH_MBC_rot(1,3)+0.25)/0.05);
R = [X1876/0.05 -(Node1876_AH_MBC_rot(1,3)+0.25)/0.05 0;(Node1876_AH_MBC_rot(1,3)+0.25)/0.05 X1876/0.05 0;...
    0 0 1];
Result_S_1876_MBC_rot_trans = R'*Result_S_1876_MBC_rot*R;


Result_S_1804_MBC_rot = [Node1804_AH_MBC_rot(1,5) Node1804_AH_MBC_rot(1,6) Node1804_AH_MBC_rot(1,7);Node1804_AH_MBC_rot(1,8) Node1804_AH_MBC_rot(1,9) Node1804_AH_MBC_rot(1,10);...
    Node1804_AH_MBC_rot(1,11) Node1804_AH_MBC_rot(1,12) Node1804_AH_MBC_rot(1,13)];
X1804 = sqrt(0.05^2-(Node1804_AH_MBC_rot(1,3)+0.25)^2);
teta1804 = asin((Node1804_AH_MBC_rot(1,3)+0.25)/0.05);
R = [X1804/0.05 -(Node1804_AH_MBC_rot(1,3)+0.25)/0.05 0;(Node1804_AH_MBC_rot(1,3)+0.25)/0.05 X1804/0.05 0;...
    0 0 1];
Result_S_1804_MBC_rot_trans = R'*Result_S_1804_MBC_rot*R;


Result_S_1732_MBC_rot = [Node1732_AH_MBC_rot(1,5) Node1732_AH_MBC_rot(1,6) Node1732_AH_MBC_rot(1,7);Node1732_AH_MBC_rot(1,8) Node1732_AH_MBC_rot(1,9) Node1732_AH_MBC_rot(1,10);...
    Node1732_AH_MBC_rot(1,11) Node1732_AH_MBC_rot(1,12) Node1732_AH_MBC_rot(1,13)];
X1732 = sqrt(0.05^2-(Node1732_AH_MBC_rot(1,3)+0.25)^2);
teta1732 = asin((Node1732_AH_MBC_rot(1,3)+0.25)/0.05);
R = [X1732/0.05 -(Node1732_AH_MBC_rot(1,3)+0.25)/0.05 0;(Node1732_AH_MBC_rot(1,3)+0.25)/0.05 X1732/0.05 0;...
    0 0 1];
Result_S_1732_MBC_rot_trans = R'*Result_S_1732_MBC_rot*R;


Result_S_1660_MBC_rot = [Node1660_AH_MBC_rot(1,5) Node1660_AH_MBC_rot(1,6) Node1660_AH_MBC_rot(1,7);Node1660_AH_MBC_rot(1,8) Node1660_AH_MBC_rot(1,9) Node1660_AH_MBC_rot(1,10);...
    Node1660_AH_MBC_rot(1,11) Node1660_AH_MBC_rot(1,12) Node1660_AH_MBC_rot(1,13)];
X1660 = sqrt(0.05^2-(Node1660_AH_MBC_rot(1,3)+0.25)^2);
teta1660 = asin((Node1660_AH_MBC_rot(1,3)+0.25)/0.05);
R = [X1660/0.05 -(Node1660_AH_MBC_rot(1,3)+0.25)/0.05 0;(Node1660_AH_MBC_rot(1,3)+0.25)/0.05 X1660/0.05 0;...
    0 0 1];
Result_S_1660_MBC_rot_trans = R'*Result_S_1660_MBC_rot*R;


Result_S_1588_MBC_rot = [Node1588_AH_MBC_rot(1,5) Node1588_AH_MBC_rot(1,6) Node1588_AH_MBC_rot(1,7);Node1588_AH_MBC_rot(1,8) Node1588_AH_MBC_rot(1,9) Node1588_AH_MBC_rot(1,10);...
    Node1588_AH_MBC_rot(1,11) Node1588_AH_MBC_rot(1,12) Node1588_AH_MBC_rot(1,13)];
X1588 = sqrt(0.05^2-(Node1588_AH_MBC_rot(1,3)+0.25)^2);
teta1588 = asin((Node1588_AH_MBC_rot(1,3)+0.25)/0.05);
R = [X1588/0.05 -(Node1588_AH_MBC_rot(1,3)+0.25)/0.05 0;(Node1588_AH_MBC_rot(1,3)+0.25)/0.05 X1588/0.05 0;...
    0 0 1];
Result_S_1588_MBC_rot_trans = R'*Result_S_1588_MBC_rot*R;


Result_S_1516_MBC_rot = [Node1516_AH_MBC_rot(1,5) Node1516_AH_MBC_rot(1,6) Node1516_AH_MBC_rot(1,7);Node1516_AH_MBC_rot(1,8) Node1516_AH_MBC_rot(1,9) Node1516_AH_MBC_rot(1,10);...
    Node1516_AH_MBC_rot(1,11) Node1516_AH_MBC_rot(1,12) Node1516_AH_MBC_rot(1,13)];
X1516 = sqrt(0.05^2-(Node1516_AH_MBC_rot(1,3)+0.25)^2);
teta1516 = asin((Node1516_AH_MBC_rot(1,3)+0.25)/0.05);
R = [X1516/0.05 -(Node1516_AH_MBC_rot(1,3)+0.25)/0.05 0;(Node1516_AH_MBC_rot(1,3)+0.25)/0.05 X1516/0.05 0;...
    0 0 1];
Result_S_1516_MBC_rot_trans = R'*Result_S_1516_MBC_rot*R;


Result_S_1444_MBC_rot = [Node1444_AH_MBC_rot(1,5) Node1444_AH_MBC_rot(1,6) Node1444_AH_MBC_rot(1,7);Node1444_AH_MBC_rot(1,8) Node1444_AH_MBC_rot(1,9) Node1444_AH_MBC_rot(1,10);...
    Node1444_AH_MBC_rot(1,11) Node1444_AH_MBC_rot(1,12) Node1444_AH_MBC_rot(1,13)];
X1444 = sqrt(0.05^2-(Node1444_AH_MBC_rot(1,3)+0.25)^2);
teta1444 = asin((Node1444_AH_MBC_rot(1,3)+0.25)/0.05);
R = [X1444/0.05 -(Node1444_AH_MBC_rot(1,3)+0.25)/0.05 0;(Node1444_AH_MBC_rot(1,3)+0.25)/0.05 X1444/0.05 0;...
    0 0 1];
Result_S_1444_MBC_rot_trans = R'*Result_S_1444_MBC_rot*R;


Result_S_1372_MBC_rot = [Node1372_AH_MBC_rot(1,5) Node1372_AH_MBC_rot(1,6) Node1372_AH_MBC_rot(1,7);Node1372_AH_MBC_rot(1,8) Node1372_AH_MBC_rot(1,9) Node1372_AH_MBC_rot(1,10);...
    Node1372_AH_MBC_rot(1,11) Node1372_AH_MBC_rot(1,12) Node1372_AH_MBC_rot(1,13)];
X1372 = sqrt(0.05^2-(Node1372_AH_MBC_rot(1,3)+0.25)^2);
teta1372 = asin((Node1372_AH_MBC_rot(1,3)+0.25)/0.05);
R = [X1372/0.05 -(Node1372_AH_MBC_rot(1,3)+0.25)/0.05 0;(Node1372_AH_MBC_rot(1,3)+0.25)/0.05 X1372/0.05 0;...
    0 0 1];
Result_S_1372_MBC_rot_trans = R'*Result_S_1372_MBC_rot*R;


Result_S_1300_MBC_rot = [Node1300_AH_MBC_rot(1,5) Node1300_AH_MBC_rot(1,6) Node1300_AH_MBC_rot(1,7);Node1300_AH_MBC_rot(1,8) Node1300_AH_MBC_rot(1,9) Node1300_AH_MBC_rot(1,10);...
    Node1300_AH_MBC_rot(1,11) Node1300_AH_MBC_rot(1,12) Node1300_AH_MBC_rot(1,13)];
X1300 = sqrt(0.05^2-(Node1300_AH_MBC_rot(1,3)+0.25)^2);
teta1300 = asin((Node1300_AH_MBC_rot(1,3)+0.25)/0.05);
R = [X1300/0.05 -(Node1300_AH_MBC_rot(1,3)+0.25)/0.05 0;(Node1300_AH_MBC_rot(1,3)+0.25)/0.05 X1300/0.05 0;...
    0 0 1];
Result_S_1300_MBC_rot_trans = R'*Result_S_1300_MBC_rot*R;


Result_S_1228_MBC_rot = [Node1228_AH_MBC_rot(1,5) Node1228_AH_MBC_rot(1,6) Node1228_AH_MBC_rot(1,7);Node1228_AH_MBC_rot(1,8) Node1228_AH_MBC_rot(1,9) Node1228_AH_MBC_rot(1,10);...
    Node1228_AH_MBC_rot(1,11) Node1228_AH_MBC_rot(1,12) Node1228_AH_MBC_rot(1,13)];
X1228 = sqrt(0.05^2-(Node1228_AH_MBC_rot(1,3)+0.25)^2);
teta1228 = asin((Node1228_AH_MBC_rot(1,3)+0.25)/0.05);
R = [X1228/0.05 -(Node1228_AH_MBC_rot(1,3)+0.25)/0.05 0;(Node1228_AH_MBC_rot(1,3)+0.25)/0.05 X1228/0.05 0;...
    0 0 1];
Result_S_1228_MBC_rot_trans = R'*Result_S_1228_MBC_rot*R;


Result_S_1156_MBC_rot = [Node1156_AH_MBC_rot(1,5) Node1156_AH_MBC_rot(1,6) Node1156_AH_MBC_rot(1,7);Node1156_AH_MBC_rot(1,8) Node1156_AH_MBC_rot(1,9) Node1156_AH_MBC_rot(1,10);...
    Node1156_AH_MBC_rot(1,11) Node1156_AH_MBC_rot(1,12) Node1156_AH_MBC_rot(1,13)];
X1156 = sqrt(0.05^2-(Node1156_AH_MBC_rot(1,3)+0.25)^2);
teta1156 = asin((Node1156_AH_MBC_rot(1,3)+0.25)/0.05);
R = [X1156/0.05 -(Node1156_AH_MBC_rot(1,3)+0.25)/0.05 0;(Node1156_AH_MBC_rot(1,3)+0.25)/0.05 X1156/0.05 0;...
    0 0 1];
Result_S_1156_MBC_rot_trans = R'*Result_S_1156_MBC_rot*R;


Result_S_1016_MBC_rot = [Node1016_AH_MBC_rot(1,5) Node1016_AH_MBC_rot(1,6) Node1016_AH_MBC_rot(1,7);Node1016_AH_MBC_rot(1,8) Node1016_AH_MBC_rot(1,9) Node1016_AH_MBC_rot(1,10);...
    Node1016_AH_MBC_rot(1,11) Node1016_AH_MBC_rot(1,12) Node1016_AH_MBC_rot(1,13)];
X1016 = sqrt(0.05^2-(Node1016_AH_MBC_rot(1,3)+0.25)^2);
teta1016 = asin((Node1016_AH_MBC_rot(1,3)+0.25)/0.05);
R = [X1016/0.05 -(Node1016_AH_MBC_rot(1,3)+0.25)/0.05 0;(Node1016_AH_MBC_rot(1,3)+0.25)/0.05 X1016/0.05 0;...
    0 0 1];
Result_S_1016_MBC_rot_trans = R'*Result_S_1016_MBC_rot*R;


Result_S_1012_MBC_rot = [Node1012_AH_MBC_rot(1,5) Node1012_AH_MBC_rot(1,6) Node1012_AH_MBC_rot(1,7);Node1012_AH_MBC_rot(1,8) Node1012_AH_MBC_rot(1,9) Node1012_AH_MBC_rot(1,10);...
    Node1012_AH_MBC_rot(1,11) Node1012_AH_MBC_rot(1,12) Node1012_AH_MBC_rot(1,13)];
X1012 = sqrt(0.05^2-(Node1012_AH_MBC_rot(1,3)+0.25)^2);
teta1012 = asin((Node1012_AH_MBC_rot(1,3)+0.25)/0.05);
R = [X1012/0.05 -(Node1012_AH_MBC_rot(1,3)+0.25)/0.05 0;(Node1012_AH_MBC_rot(1,3)+0.25)/0.05 X1012/0.05 0;...
    0 0 1];
Result_S_1012_MBC_rot_trans = R'*Result_S_1012_MBC_rot*R;


Result_S_teta_teta_MBC_rot = [teta2 Result_S_2_MBC_rot_trans(2,2)/P_applied;teta6 Result_S_6_MBC_rot_trans(2,2)/P_applied;teta146 Result_S_146_MBC_rot_trans(2,2)/P_applied;...
    teta218 Result_S_218_MBC_rot_trans(2,2)/P_applied;teta290 Result_S_290_MBC_rot_trans(2,2)/P_applied;teta362 Result_S_362_MBC_rot_trans(2,2)/P_applied;...
    teta434 Result_S_434_MBC_rot_trans(2,2)/P_applied;teta506 Result_S_506_MBC_rot_trans(2,2)/P_applied;teta578 Result_S_578_MBC_rot_trans(2,2)/P_applied;...
    teta650 Result_S_650_MBC_rot_trans(2,2)/P_applied;teta722 Result_S_722_MBC_rot_trans(2,2)/P_applied;teta794 Result_S_794_MBC_rot_trans(2,2)/P_applied;...
    teta866 Result_S_866_MBC_rot_trans(2,2)/P_applied;teta938 Result_S_938_MBC_rot_trans(2,2)/P_applied;teta1876 Result_S_1876_MBC_rot_trans(2,2)/P_applied;...
    teta1804 Result_S_1804_MBC_rot_trans(2,2)/P_applied;teta1732 Result_S_1732_MBC_rot_trans(2,2)/P_applied;teta1660 Result_S_1660_MBC_rot_trans(2,2)/P_applied;...
    teta1588 Result_S_1588_MBC_rot_trans(2,2)/P_applied;teta1516 Result_S_1516_MBC_rot_trans(2,2)/P_applied;teta1444 Result_S_1444_MBC_rot_trans(2,2)/P_applied;...
    teta1372 Result_S_1372_MBC_rot_trans(2,2)/P_applied;teta1300 Result_S_1300_MBC_rot_trans(2,2)/P_applied;teta1228 Result_S_1228_MBC_rot_trans(2,2)/P_applied;...
    teta1156 Result_S_1156_MBC_rot_trans(2,2)/P_applied;teta1016 Result_S_1016_MBC_rot_trans(2,2)/P_applied;teta1012 Result_S_1012_MBC_rot_trans(2,2)/P_applied];



Result_S_2_MBC_str = [Node2_AH_MBC_str(1,5) Node2_AH_MBC_str(1,6) Node2_AH_MBC_str(1,7);Node2_AH_MBC_str(1,8) Node2_AH_MBC_str(1,9) Node2_AH_MBC_str(1,10);...
    Node2_AH_MBC_str(1,11) Node2_AH_MBC_str(1,12) Node2_AH_MBC_str(1,13)];
X2 = sqrt(0.05^2-(Node2_AH_MBC_str(1,3)+0.25)^2);
teta2 = asin((Node2_AH_MBC_str(1,3)+0.25)/0.05);
R = [X2/0.05 -(Node2_AH_MBC_str(1,3)+0.25)/0.05 0;(Node2_AH_MBC_str(1,3)+0.25)/0.05 X2/0.05 0;...
    0 0 1];
Result_S_2_MBC_str_trans = R'*Result_S_2_MBC_str*R;


Result_S_6_MBC_str = [Node6_AH_MBC_str(1,5) Node6_AH_MBC_str(1,6) Node6_AH_MBC_str(1,7);Node6_AH_MBC_str(1,8) Node6_AH_MBC_str(1,9) Node6_AH_MBC_str(1,10);...
    Node6_AH_MBC_str(1,11) Node6_AH_MBC_str(1,12) Node6_AH_MBC_str(1,13)];
X6 = sqrt(0.05^2-(Node6_AH_MBC_str(1,3)+0.25)^2);
teta6 = asin((Node6_AH_MBC_str(1,3)+0.25)/0.05);
R = [X6/0.05 -(Node6_AH_MBC_str(1,3)+0.25)/0.05 0;(Node6_AH_MBC_str(1,3)+0.25)/0.05 X6/0.05 0;...
    0 0 1];
Result_S_6_MBC_str_trans = R'*Result_S_6_MBC_str*R;


Result_S_146_MBC_str = [Node146_AH_MBC_str(1,5) Node146_AH_MBC_str(1,6) Node146_AH_MBC_str(1,7);Node146_AH_MBC_str(1,8) Node146_AH_MBC_str(1,9) Node146_AH_MBC_str(1,10);...
    Node146_AH_MBC_str(1,11) Node146_AH_MBC_str(1,12) Node146_AH_MBC_str(1,13)];
X146 = sqrt(0.05^2-(Node146_AH_MBC_str(1,3)+0.25)^2);
teta146 = asin((Node146_AH_MBC_str(1,3)+0.25)/0.05);
R = [X146/0.05 -(Node146_AH_MBC_str(1,3)+0.25)/0.05 0;(Node146_AH_MBC_str(1,3)+0.25)/0.05 X146/0.05 0;...
    0 0 1];
Result_S_146_MBC_str_trans = R'*Result_S_146_MBC_str*R;


Result_S_218_MBC_str = [Node218_AH_MBC_str(1,5) Node218_AH_MBC_str(1,6) Node218_AH_MBC_str(1,7);Node218_AH_MBC_str(1,8) Node218_AH_MBC_str(1,9) Node218_AH_MBC_str(1,10);...
    Node218_AH_MBC_str(1,11) Node218_AH_MBC_str(1,12) Node218_AH_MBC_str(1,13)];
X218 = sqrt(0.05^2-(Node218_AH_MBC_str(1,3)+0.25)^2);
teta218 = asin((Node218_AH_MBC_str(1,3)+0.25)/0.05);
R = [X218/0.05 -(Node218_AH_MBC_str(1,3)+0.25)/0.05 0;(Node218_AH_MBC_str(1,3)+0.25)/0.05 X218/0.05 0;...
    0 0 1];
Result_S_218_MBC_str_trans = R'*Result_S_218_MBC_str*R;


Result_S_290_MBC_str = [Node290_AH_MBC_str(1,5) Node290_AH_MBC_str(1,6) Node290_AH_MBC_str(1,7);Node290_AH_MBC_str(1,8) Node290_AH_MBC_str(1,9) Node290_AH_MBC_str(1,10);...
    Node290_AH_MBC_str(1,11) Node290_AH_MBC_str(1,12) Node290_AH_MBC_str(1,13)];
X290 = sqrt(0.05^2-(Node290_AH_MBC_str(1,3)+0.25)^2);
teta290 = asin((Node290_AH_MBC_str(1,3)+0.25)/0.05);
R = [X290/0.05 -(Node290_AH_MBC_str(1,3)+0.25)/0.05 0;(Node290_AH_MBC_str(1,3)+0.25)/0.05 X290/0.05 0;...
    0 0 1];
Result_S_290_MBC_str_trans = R'*Result_S_290_MBC_str*R;


Result_S_362_MBC_str = [Node362_AH_MBC_str(1,5) Node362_AH_MBC_str(1,6) Node362_AH_MBC_str(1,7);Node362_AH_MBC_str(1,8) Node362_AH_MBC_str(1,9) Node362_AH_MBC_str(1,10);...
    Node362_AH_MBC_str(1,11) Node362_AH_MBC_str(1,12) Node362_AH_MBC_str(1,13)];
X362 = sqrt(0.05^2-(Node362_AH_MBC_str(1,3)+0.25)^2);
teta362 = asin((Node362_AH_MBC_str(1,3)+0.25)/0.05);
R = [X362/0.05 -(Node362_AH_MBC_str(1,3)+0.25)/0.05 0;(Node362_AH_MBC_str(1,3)+0.25)/0.05 X362/0.05 0;...
    0 0 1];
Result_S_362_MBC_str_trans = R'*Result_S_362_MBC_str*R;


Result_S_434_MBC_str = [Node434_AH_MBC_str(1,5) Node434_AH_MBC_str(1,6) Node434_AH_MBC_str(1,7);Node434_AH_MBC_str(1,8) Node434_AH_MBC_str(1,9) Node434_AH_MBC_str(1,10);...
    Node434_AH_MBC_str(1,11) Node434_AH_MBC_str(1,12) Node434_AH_MBC_str(1,13)];
X434 = sqrt(0.05^2-(Node434_AH_MBC_str(1,3)+0.25)^2);
teta434 = asin((Node434_AH_MBC_str(1,3)+0.25)/0.05);
R = [X434/0.05 -(Node434_AH_MBC_str(1,3)+0.25)/0.05 0;(Node434_AH_MBC_str(1,3)+0.25)/0.05 X434/0.05 0;...
    0 0 1];
Result_S_434_MBC_str_trans = R'*Result_S_434_MBC_str*R;


Result_S_506_MBC_str = [Node506_AH_MBC_str(1,5) Node506_AH_MBC_str(1,6) Node506_AH_MBC_str(1,7);Node506_AH_MBC_str(1,8) Node506_AH_MBC_str(1,9) Node506_AH_MBC_str(1,10);...
    Node506_AH_MBC_str(1,11) Node506_AH_MBC_str(1,12) Node506_AH_MBC_str(1,13)];
X506 = sqrt(0.05^2-(Node506_AH_MBC_str(1,3)+0.25)^2);
teta506 = asin((Node506_AH_MBC_str(1,3)+0.25)/0.05);
R = [X506/0.05 -(Node506_AH_MBC_str(1,3)+0.25)/0.05 0;(Node506_AH_MBC_str(1,3)+0.25)/0.05 X506/0.05 0;...
    0 0 1];
Result_S_506_MBC_str_trans = R'*Result_S_506_MBC_str*R;


Result_S_578_MBC_str = [Node578_AH_MBC_str(1,5) Node578_AH_MBC_str(1,6) Node578_AH_MBC_str(1,7);Node578_AH_MBC_str(1,8) Node578_AH_MBC_str(1,9) Node578_AH_MBC_str(1,10);...
    Node578_AH_MBC_str(1,11) Node578_AH_MBC_str(1,12) Node578_AH_MBC_str(1,13)];
X578 = sqrt(0.05^2-(Node578_AH_MBC_str(1,3)+0.25)^2);
teta578 = asin((Node578_AH_MBC_str(1,3)+0.25)/0.05);
R = [X578/0.05 -(Node578_AH_MBC_str(1,3)+0.25)/0.05 0;(Node578_AH_MBC_str(1,3)+0.25)/0.05 X578/0.05 0;...
    0 0 1];
Result_S_578_MBC_str_trans = R'*Result_S_578_MBC_str*R;


Result_S_650_MBC_str = [Node650_AH_MBC_str(1,5) Node650_AH_MBC_str(1,6) Node650_AH_MBC_str(1,7);Node650_AH_MBC_str(1,8) Node650_AH_MBC_str(1,9) Node650_AH_MBC_str(1,10);...
    Node650_AH_MBC_str(1,11) Node650_AH_MBC_str(1,12) Node650_AH_MBC_str(1,13)];
X650 = sqrt(0.05^2-(Node650_AH_MBC_str(1,3)+0.25)^2);
teta650 = asin((Node650_AH_MBC_str(1,3)+0.25)/0.05);
R = [X650/0.05 -(Node650_AH_MBC_str(1,3)+0.25)/0.05 0;(Node650_AH_MBC_str(1,3)+0.25)/0.05 X650/0.05 0;...
    0 0 1];
Result_S_650_MBC_str_trans = R'*Result_S_650_MBC_str*R;


Result_S_722_MBC_str = [Node722_AH_MBC_str(1,5) Node722_AH_MBC_str(1,6) Node722_AH_MBC_str(1,7);Node722_AH_MBC_str(1,8) Node722_AH_MBC_str(1,9) Node722_AH_MBC_str(1,10);...
    Node722_AH_MBC_str(1,11) Node722_AH_MBC_str(1,12) Node722_AH_MBC_str(1,13)];
X722 = sqrt(0.05^2-(Node722_AH_MBC_str(1,3)+0.25)^2);
teta722 = asin((Node722_AH_MBC_str(1,3)+0.25)/0.05);
R = [X722/0.05 -(Node722_AH_MBC_str(1,3)+0.25)/0.05 0;(Node722_AH_MBC_str(1,3)+0.25)/0.05 X722/0.05 0;...
    0 0 1];
Result_S_722_MBC_str_trans = R'*Result_S_722_MBC_str*R;


Result_S_794_MBC_str = [Node794_AH_MBC_str(1,5) Node794_AH_MBC_str(1,6) Node794_AH_MBC_str(1,7);Node794_AH_MBC_str(1,8) Node794_AH_MBC_str(1,9) Node794_AH_MBC_str(1,10);...
    Node794_AH_MBC_str(1,11) Node794_AH_MBC_str(1,12) Node794_AH_MBC_str(1,13)];
X794 = sqrt(0.05^2-(Node794_AH_MBC_str(1,3)+0.25)^2);
teta794 = asin((Node794_AH_MBC_str(1,3)+0.25)/0.05);
R = [X794/0.05 -(Node794_AH_MBC_str(1,3)+0.25)/0.05 0;(Node794_AH_MBC_str(1,3)+0.25)/0.05 X794/0.05 0;...
    0 0 1];
Result_S_794_MBC_str_trans = R'*Result_S_794_MBC_str*R;


Result_S_866_MBC_str = [Node866_AH_MBC_str(1,5) Node866_AH_MBC_str(1,6) Node866_AH_MBC_str(1,7);Node866_AH_MBC_str(1,8) Node866_AH_MBC_str(1,9) Node866_AH_MBC_str(1,10);...
    Node866_AH_MBC_str(1,11) Node866_AH_MBC_str(1,12) Node866_AH_MBC_str(1,13)];
X866 = sqrt(0.05^2-(Node866_AH_MBC_str(1,3)+0.25)^2);
teta866 = asin((Node866_AH_MBC_str(1,3)+0.25)/0.05);
R = [X866/0.05 -(Node866_AH_MBC_str(1,3)+0.25)/0.05 0;(Node866_AH_MBC_str(1,3)+0.25)/0.05 X866/0.05 0;...
    0 0 1];
Result_S_866_MBC_str_trans = R'*Result_S_866_MBC_str*R;


Result_S_938_MBC_str = [Node938_AH_MBC_str(1,5) Node938_AH_MBC_str(1,6) Node938_AH_MBC_str(1,7);Node938_AH_MBC_str(1,8) Node938_AH_MBC_str(1,9) Node938_AH_MBC_str(1,10);...
    Node938_AH_MBC_str(1,11) Node938_AH_MBC_str(1,12) Node938_AH_MBC_str(1,13)];
X938 = sqrt(0.05^2-(Node938_AH_MBC_str(1,3)+0.25)^2);
teta938 = asin((Node938_AH_MBC_str(1,3)+0.25)/0.05);
R = [X938/0.05 -(Node938_AH_MBC_str(1,3)+0.25)/0.05 0;(Node938_AH_MBC_str(1,3)+0.25)/0.05 X938/0.05 0;...
    0 0 1];
Result_S_938_MBC_str_trans = R'*Result_S_938_MBC_str*R;


Result_S_1876_MBC_str = [Node1876_AH_MBC_str(1,5) Node1876_AH_MBC_str(1,6) Node1876_AH_MBC_str(1,7);Node1876_AH_MBC_str(1,8) Node1876_AH_MBC_str(1,9) Node1876_AH_MBC_str(1,10);...
    Node1876_AH_MBC_str(1,11) Node1876_AH_MBC_str(1,12) Node1876_AH_MBC_str(1,13)];
X1876 = sqrt(0.05^2-(Node1876_AH_MBC_str(1,3)+0.25)^2);
teta1876 = asin((Node1876_AH_MBC_str(1,3)+0.25)/0.05);
R = [X1876/0.05 -(Node1876_AH_MBC_str(1,3)+0.25)/0.05 0;(Node1876_AH_MBC_str(1,3)+0.25)/0.05 X1876/0.05 0;...
    0 0 1];
Result_S_1876_MBC_str_trans = R'*Result_S_1876_MBC_str*R;


Result_S_1804_MBC_str = [Node1804_AH_MBC_str(1,5) Node1804_AH_MBC_str(1,6) Node1804_AH_MBC_str(1,7);Node1804_AH_MBC_str(1,8) Node1804_AH_MBC_str(1,9) Node1804_AH_MBC_str(1,10);...
    Node1804_AH_MBC_str(1,11) Node1804_AH_MBC_str(1,12) Node1804_AH_MBC_str(1,13)];
X1804 = sqrt(0.05^2-(Node1804_AH_MBC_str(1,3)+0.25)^2);
teta1804 = asin((Node1804_AH_MBC_str(1,3)+0.25)/0.05);
R = [X1804/0.05 -(Node1804_AH_MBC_str(1,3)+0.25)/0.05 0;(Node1804_AH_MBC_str(1,3)+0.25)/0.05 X1804/0.05 0;...
    0 0 1];
Result_S_1804_MBC_str_trans = R'*Result_S_1804_MBC_str*R;


Result_S_1732_MBC_str = [Node1732_AH_MBC_str(1,5) Node1732_AH_MBC_str(1,6) Node1732_AH_MBC_str(1,7);Node1732_AH_MBC_str(1,8) Node1732_AH_MBC_str(1,9) Node1732_AH_MBC_str(1,10);...
    Node1732_AH_MBC_str(1,11) Node1732_AH_MBC_str(1,12) Node1732_AH_MBC_str(1,13)];
X1732 = sqrt(0.05^2-(Node1732_AH_MBC_str(1,3)+0.25)^2);
teta1732 = asin((Node1732_AH_MBC_str(1,3)+0.25)/0.05);
R = [X1732/0.05 -(Node1732_AH_MBC_str(1,3)+0.25)/0.05 0;(Node1732_AH_MBC_str(1,3)+0.25)/0.05 X1732/0.05 0;...
    0 0 1];
Result_S_1732_MBC_str_trans = R'*Result_S_1732_MBC_str*R;


Result_S_1660_MBC_str = [Node1660_AH_MBC_str(1,5) Node1660_AH_MBC_str(1,6) Node1660_AH_MBC_str(1,7);Node1660_AH_MBC_str(1,8) Node1660_AH_MBC_str(1,9) Node1660_AH_MBC_str(1,10);...
    Node1660_AH_MBC_str(1,11) Node1660_AH_MBC_str(1,12) Node1660_AH_MBC_str(1,13)];
X1660 = sqrt(0.05^2-(Node1660_AH_MBC_str(1,3)+0.25)^2);
teta1660 = asin((Node1660_AH_MBC_str(1,3)+0.25)/0.05);
R = [X1660/0.05 -(Node1660_AH_MBC_str(1,3)+0.25)/0.05 0;(Node1660_AH_MBC_str(1,3)+0.25)/0.05 X1660/0.05 0;...
    0 0 1];
Result_S_1660_MBC_str_trans = R'*Result_S_1660_MBC_str*R;


Result_S_1588_MBC_str = [Node1588_AH_MBC_str(1,5) Node1588_AH_MBC_str(1,6) Node1588_AH_MBC_str(1,7);Node1588_AH_MBC_str(1,8) Node1588_AH_MBC_str(1,9) Node1588_AH_MBC_str(1,10);...
    Node1588_AH_MBC_str(1,11) Node1588_AH_MBC_str(1,12) Node1588_AH_MBC_str(1,13)];
X1588 = sqrt(0.05^2-(Node1588_AH_MBC_str(1,3)+0.25)^2);
teta1588 = asin((Node1588_AH_MBC_str(1,3)+0.25)/0.05);
R = [X1588/0.05 -(Node1588_AH_MBC_str(1,3)+0.25)/0.05 0;(Node1588_AH_MBC_str(1,3)+0.25)/0.05 X1588/0.05 0;...
    0 0 1];
Result_S_1588_MBC_str_trans = R'*Result_S_1588_MBC_str*R;


Result_S_1516_MBC_str = [Node1516_AH_MBC_str(1,5) Node1516_AH_MBC_str(1,6) Node1516_AH_MBC_str(1,7);Node1516_AH_MBC_str(1,8) Node1516_AH_MBC_str(1,9) Node1516_AH_MBC_str(1,10);...
    Node1516_AH_MBC_str(1,11) Node1516_AH_MBC_str(1,12) Node1516_AH_MBC_str(1,13)];
X1516 = sqrt(0.05^2-(Node1516_AH_MBC_str(1,3)+0.25)^2);
teta1516 = asin((Node1516_AH_MBC_str(1,3)+0.25)/0.05);
R = [X1516/0.05 -(Node1516_AH_MBC_str(1,3)+0.25)/0.05 0;(Node1516_AH_MBC_str(1,3)+0.25)/0.05 X1516/0.05 0;...
    0 0 1];
Result_S_1516_MBC_str_trans = R'*Result_S_1516_MBC_str*R;


Result_S_1444_MBC_str = [Node1444_AH_MBC_str(1,5) Node1444_AH_MBC_str(1,6) Node1444_AH_MBC_str(1,7);Node1444_AH_MBC_str(1,8) Node1444_AH_MBC_str(1,9) Node1444_AH_MBC_str(1,10);...
    Node1444_AH_MBC_str(1,11) Node1444_AH_MBC_str(1,12) Node1444_AH_MBC_str(1,13)];
X1444 = sqrt(0.05^2-(Node1444_AH_MBC_str(1,3)+0.25)^2);
teta1444 = asin((Node1444_AH_MBC_str(1,3)+0.25)/0.05);
R = [X1444/0.05 -(Node1444_AH_MBC_str(1,3)+0.25)/0.05 0;(Node1444_AH_MBC_str(1,3)+0.25)/0.05 X1444/0.05 0;...
    0 0 1];
Result_S_1444_MBC_str_trans = R'*Result_S_1444_MBC_str*R;


Result_S_1372_MBC_str = [Node1372_AH_MBC_str(1,5) Node1372_AH_MBC_str(1,6) Node1372_AH_MBC_str(1,7);Node1372_AH_MBC_str(1,8) Node1372_AH_MBC_str(1,9) Node1372_AH_MBC_str(1,10);...
    Node1372_AH_MBC_str(1,11) Node1372_AH_MBC_str(1,12) Node1372_AH_MBC_str(1,13)];
X1372 = sqrt(0.05^2-(Node1372_AH_MBC_str(1,3)+0.25)^2);
teta1372 = asin((Node1372_AH_MBC_str(1,3)+0.25)/0.05);
R = [X1372/0.05 -(Node1372_AH_MBC_str(1,3)+0.25)/0.05 0;(Node1372_AH_MBC_str(1,3)+0.25)/0.05 X1372/0.05 0;...
    0 0 1];
Result_S_1372_MBC_str_trans = R'*Result_S_1372_MBC_str*R;


Result_S_1300_MBC_str = [Node1300_AH_MBC_str(1,5) Node1300_AH_MBC_str(1,6) Node1300_AH_MBC_str(1,7);Node1300_AH_MBC_str(1,8) Node1300_AH_MBC_str(1,9) Node1300_AH_MBC_str(1,10);...
    Node1300_AH_MBC_str(1,11) Node1300_AH_MBC_str(1,12) Node1300_AH_MBC_str(1,13)];
X1300 = sqrt(0.05^2-(Node1300_AH_MBC_str(1,3)+0.25)^2);
teta1300 = asin((Node1300_AH_MBC_str(1,3)+0.25)/0.05);
R = [X1300/0.05 -(Node1300_AH_MBC_str(1,3)+0.25)/0.05 0;(Node1300_AH_MBC_str(1,3)+0.25)/0.05 X1300/0.05 0;...
    0 0 1];
Result_S_1300_MBC_str_trans = R'*Result_S_1300_MBC_str*R;


Result_S_1228_MBC_str = [Node1228_AH_MBC_str(1,5) Node1228_AH_MBC_str(1,6) Node1228_AH_MBC_str(1,7);Node1228_AH_MBC_str(1,8) Node1228_AH_MBC_str(1,9) Node1228_AH_MBC_str(1,10);...
    Node1228_AH_MBC_str(1,11) Node1228_AH_MBC_str(1,12) Node1228_AH_MBC_str(1,13)];
X1228 = sqrt(0.05^2-(Node1228_AH_MBC_str(1,3)+0.25)^2);
teta1228 = asin((Node1228_AH_MBC_str(1,3)+0.25)/0.05);
R = [X1228/0.05 -(Node1228_AH_MBC_str(1,3)+0.25)/0.05 0;(Node1228_AH_MBC_str(1,3)+0.25)/0.05 X1228/0.05 0;...
    0 0 1];
Result_S_1228_MBC_str_trans = R'*Result_S_1228_MBC_str*R;


Result_S_1156_MBC_str = [Node1156_AH_MBC_str(1,5) Node1156_AH_MBC_str(1,6) Node1156_AH_MBC_str(1,7);Node1156_AH_MBC_str(1,8) Node1156_AH_MBC_str(1,9) Node1156_AH_MBC_str(1,10);...
    Node1156_AH_MBC_str(1,11) Node1156_AH_MBC_str(1,12) Node1156_AH_MBC_str(1,13)];
X1156 = sqrt(0.05^2-(Node1156_AH_MBC_str(1,3)+0.25)^2);
teta1156 = asin((Node1156_AH_MBC_str(1,3)+0.25)/0.05);
R = [X1156/0.05 -(Node1156_AH_MBC_str(1,3)+0.25)/0.05 0;(Node1156_AH_MBC_str(1,3)+0.25)/0.05 X1156/0.05 0;...
    0 0 1];
Result_S_1156_MBC_str_trans = R'*Result_S_1156_MBC_str*R;


Result_S_1016_MBC_str = [Node1016_AH_MBC_str(1,5) Node1016_AH_MBC_str(1,6) Node1016_AH_MBC_str(1,7);Node1016_AH_MBC_str(1,8) Node1016_AH_MBC_str(1,9) Node1016_AH_MBC_str(1,10);...
    Node1016_AH_MBC_str(1,11) Node1016_AH_MBC_str(1,12) Node1016_AH_MBC_str(1,13)];
X1016 = sqrt(0.05^2-(Node1016_AH_MBC_str(1,3)+0.25)^2);
teta1016 = asin((Node1016_AH_MBC_str(1,3)+0.25)/0.05);
R = [X1016/0.05 -(Node1016_AH_MBC_str(1,3)+0.25)/0.05 0;(Node1016_AH_MBC_str(1,3)+0.25)/0.05 X1016/0.05 0;...
    0 0 1];
Result_S_1016_MBC_str_trans = R'*Result_S_1016_MBC_str*R;


Result_S_1012_MBC_str = [Node1012_AH_MBC_str(1,5) Node1012_AH_MBC_str(1,6) Node1012_AH_MBC_str(1,7);Node1012_AH_MBC_str(1,8) Node1012_AH_MBC_str(1,9) Node1012_AH_MBC_str(1,10);...
    Node1012_AH_MBC_str(1,11) Node1012_AH_MBC_str(1,12) Node1012_AH_MBC_str(1,13)];
X1012 = sqrt(0.05^2-(Node1012_AH_MBC_str(1,3)+0.25)^2);
teta1012 = asin((Node1012_AH_MBC_str(1,3)+0.25)/0.05);
R = [X1012/0.05 -(Node1012_AH_MBC_str(1,3)+0.25)/0.05 0;(Node1012_AH_MBC_str(1,3)+0.25)/0.05 X1012/0.05 0;...
    0 0 1];
Result_S_1012_MBC_str_trans = R'*Result_S_1012_MBC_str*R;


Result_S_teta_teta_MBC_str = [teta2 Result_S_2_MBC_str_trans(2,2)/P_applied;teta6 Result_S_6_MBC_str_trans(2,2)/P_applied;teta146 Result_S_146_MBC_str_trans(2,2)/P_applied;...
    teta218 Result_S_218_MBC_str_trans(2,2)/P_applied;teta290 Result_S_290_MBC_str_trans(2,2)/P_applied;teta362 Result_S_362_MBC_str_trans(2,2)/P_applied;...
    teta434 Result_S_434_MBC_str_trans(2,2)/P_applied;teta506 Result_S_506_MBC_str_trans(2,2)/P_applied;teta578 Result_S_578_MBC_str_trans(2,2)/P_applied;...
    teta650 Result_S_650_MBC_str_trans(2,2)/P_applied;teta722 Result_S_722_MBC_str_trans(2,2)/P_applied;teta794 Result_S_794_MBC_str_trans(2,2)/P_applied;...
    teta866 Result_S_866_MBC_str_trans(2,2)/P_applied;teta938 Result_S_938_MBC_str_trans(2,2)/P_applied;teta1876 Result_S_1876_MBC_str_trans(2,2)/P_applied;...
    teta1804 Result_S_1804_MBC_str_trans(2,2)/P_applied;teta1732 Result_S_1732_MBC_str_trans(2,2)/P_applied;teta1660 Result_S_1660_MBC_str_trans(2,2)/P_applied;...
    teta1588 Result_S_1588_MBC_str_trans(2,2)/P_applied;teta1516 Result_S_1516_MBC_str_trans(2,2)/P_applied;teta1444 Result_S_1444_MBC_str_trans(2,2)/P_applied;...
    teta1372 Result_S_1372_MBC_str_trans(2,2)/P_applied;teta1300 Result_S_1300_MBC_str_trans(2,2)/P_applied;teta1228 Result_S_1228_MBC_str_trans(2,2)/P_applied;...
    teta1156 Result_S_1156_MBC_str_trans(2,2)/P_applied;teta1016 Result_S_1016_MBC_str_trans(2,2)/P_applied;teta1012 Result_S_1012_MBC_str_trans(2,2)/P_applied];



Result_S_2_NMBC_com = [Node2_AH_NMBC_com(1,5) Node2_AH_NMBC_com(1,6) Node2_AH_NMBC_com(1,7);Node2_AH_NMBC_com(1,8) Node2_AH_NMBC_com(1,9) Node2_AH_NMBC_com(1,10);...
    Node2_AH_NMBC_com(1,11) Node2_AH_NMBC_com(1,12) Node2_AH_NMBC_com(1,13)];
X2 = sqrt(0.05^2-(Node2_AH_NMBC_com(1,3)+0.25)^2);
teta2 = asin((Node2_AH_NMBC_com(1,3)+0.25)/0.05);
R = [X2/0.05 -(Node2_AH_NMBC_com(1,3)+0.25)/0.05 0;(Node2_AH_NMBC_com(1,3)+0.25)/0.05 X2/0.05 0;...
    0 0 1];
Result_S_2_NMBC_com_trans = R'*Result_S_2_NMBC_com*R;


Result_S_6_NMBC_com = [Node6_AH_NMBC_com(1,5) Node6_AH_NMBC_com(1,6) Node6_AH_NMBC_com(1,7);Node6_AH_NMBC_com(1,8) Node6_AH_NMBC_com(1,9) Node6_AH_NMBC_com(1,10);...
    Node6_AH_NMBC_com(1,11) Node6_AH_NMBC_com(1,12) Node6_AH_NMBC_com(1,13)];
X6 = sqrt(0.05^2-(Node6_AH_NMBC_com(1,3)+0.25)^2);
teta6 = asin((Node6_AH_NMBC_com(1,3)+0.25)/0.05);
R = [X6/0.05 -(Node6_AH_NMBC_com(1,3)+0.25)/0.05 0;(Node6_AH_NMBC_com(1,3)+0.25)/0.05 X6/0.05 0;...
    0 0 1];
Result_S_6_NMBC_com_trans = R'*Result_S_6_NMBC_com*R;


Result_S_146_NMBC_com = [Node146_AH_NMBC_com(1,5) Node146_AH_NMBC_com(1,6) Node146_AH_NMBC_com(1,7);Node146_AH_NMBC_com(1,8) Node146_AH_NMBC_com(1,9) Node146_AH_NMBC_com(1,10);...
    Node146_AH_NMBC_com(1,11) Node146_AH_NMBC_com(1,12) Node146_AH_NMBC_com(1,13)];
X146 = sqrt(0.05^2-(Node146_AH_NMBC_com(1,3)+0.25)^2);
teta146 = asin((Node146_AH_NMBC_com(1,3)+0.25)/0.05);
R = [X146/0.05 -(Node146_AH_NMBC_com(1,3)+0.25)/0.05 0;(Node146_AH_NMBC_com(1,3)+0.25)/0.05 X146/0.05 0;...
    0 0 1];
Result_S_146_NMBC_com_trans = R'*Result_S_146_NMBC_com*R;


Result_S_218_NMBC_com = [Node218_AH_NMBC_com(1,5) Node218_AH_NMBC_com(1,6) Node218_AH_NMBC_com(1,7);Node218_AH_NMBC_com(1,8) Node218_AH_NMBC_com(1,9) Node218_AH_NMBC_com(1,10);...
    Node218_AH_NMBC_com(1,11) Node218_AH_NMBC_com(1,12) Node218_AH_NMBC_com(1,13)];
X218 = sqrt(0.05^2-(Node218_AH_NMBC_com(1,3)+0.25)^2);
teta218 = asin((Node218_AH_NMBC_com(1,3)+0.25)/0.05);
R = [X218/0.05 -(Node218_AH_NMBC_com(1,3)+0.25)/0.05 0;(Node218_AH_NMBC_com(1,3)+0.25)/0.05 X218/0.05 0;...
    0 0 1];
Result_S_218_NMBC_com_trans = R'*Result_S_218_NMBC_com*R;


Result_S_290_NMBC_com = [Node290_AH_NMBC_com(1,5) Node290_AH_NMBC_com(1,6) Node290_AH_NMBC_com(1,7);Node290_AH_NMBC_com(1,8) Node290_AH_NMBC_com(1,9) Node290_AH_NMBC_com(1,10);...
    Node290_AH_NMBC_com(1,11) Node290_AH_NMBC_com(1,12) Node290_AH_NMBC_com(1,13)];
X290 = sqrt(0.05^2-(Node290_AH_NMBC_com(1,3)+0.25)^2);
teta290 = asin((Node290_AH_NMBC_com(1,3)+0.25)/0.05);
R = [X290/0.05 -(Node290_AH_NMBC_com(1,3)+0.25)/0.05 0;(Node290_AH_NMBC_com(1,3)+0.25)/0.05 X290/0.05 0;...
    0 0 1];
Result_S_290_NMBC_com_trans = R'*Result_S_290_NMBC_com*R;


Result_S_362_NMBC_com = [Node362_AH_NMBC_com(1,5) Node362_AH_NMBC_com(1,6) Node362_AH_NMBC_com(1,7);Node362_AH_NMBC_com(1,8) Node362_AH_NMBC_com(1,9) Node362_AH_NMBC_com(1,10);...
    Node362_AH_NMBC_com(1,11) Node362_AH_NMBC_com(1,12) Node362_AH_NMBC_com(1,13)];
X362 = sqrt(0.05^2-(Node362_AH_NMBC_com(1,3)+0.25)^2);
teta362 = asin((Node362_AH_NMBC_com(1,3)+0.25)/0.05);
R = [X362/0.05 -(Node362_AH_NMBC_com(1,3)+0.25)/0.05 0;(Node362_AH_NMBC_com(1,3)+0.25)/0.05 X362/0.05 0;...
    0 0 1];
Result_S_362_NMBC_com_trans = R'*Result_S_362_NMBC_com*R;


Result_S_434_NMBC_com = [Node434_AH_NMBC_com(1,5) Node434_AH_NMBC_com(1,6) Node434_AH_NMBC_com(1,7);Node434_AH_NMBC_com(1,8) Node434_AH_NMBC_com(1,9) Node434_AH_NMBC_com(1,10);...
    Node434_AH_NMBC_com(1,11) Node434_AH_NMBC_com(1,12) Node434_AH_NMBC_com(1,13)];
X434 = sqrt(0.05^2-(Node434_AH_NMBC_com(1,3)+0.25)^2);
teta434 = asin((Node434_AH_NMBC_com(1,3)+0.25)/0.05);
R = [X434/0.05 -(Node434_AH_NMBC_com(1,3)+0.25)/0.05 0;(Node434_AH_NMBC_com(1,3)+0.25)/0.05 X434/0.05 0;...
    0 0 1];
Result_S_434_NMBC_com_trans = R'*Result_S_434_NMBC_com*R;


Result_S_506_NMBC_com = [Node506_AH_NMBC_com(1,5) Node506_AH_NMBC_com(1,6) Node506_AH_NMBC_com(1,7);Node506_AH_NMBC_com(1,8) Node506_AH_NMBC_com(1,9) Node506_AH_NMBC_com(1,10);...
    Node506_AH_NMBC_com(1,11) Node506_AH_NMBC_com(1,12) Node506_AH_NMBC_com(1,13)];
X506 = sqrt(0.05^2-(Node506_AH_NMBC_com(1,3)+0.25)^2);
teta506 = asin((Node506_AH_NMBC_com(1,3)+0.25)/0.05);
R = [X506/0.05 -(Node506_AH_NMBC_com(1,3)+0.25)/0.05 0;(Node506_AH_NMBC_com(1,3)+0.25)/0.05 X506/0.05 0;...
    0 0 1];
Result_S_506_NMBC_com_trans = R'*Result_S_506_NMBC_com*R;


Result_S_578_NMBC_com = [Node578_AH_NMBC_com(1,5) Node578_AH_NMBC_com(1,6) Node578_AH_NMBC_com(1,7);Node578_AH_NMBC_com(1,8) Node578_AH_NMBC_com(1,9) Node578_AH_NMBC_com(1,10);...
    Node578_AH_NMBC_com(1,11) Node578_AH_NMBC_com(1,12) Node578_AH_NMBC_com(1,13)];
X578 = sqrt(0.05^2-(Node578_AH_NMBC_com(1,3)+0.25)^2);
teta578 = asin((Node578_AH_NMBC_com(1,3)+0.25)/0.05);
R = [X578/0.05 -(Node578_AH_NMBC_com(1,3)+0.25)/0.05 0;(Node578_AH_NMBC_com(1,3)+0.25)/0.05 X578/0.05 0;...
    0 0 1];
Result_S_578_NMBC_com_trans = R'*Result_S_578_NMBC_com*R;


Result_S_650_NMBC_com = [Node650_AH_NMBC_com(1,5) Node650_AH_NMBC_com(1,6) Node650_AH_NMBC_com(1,7);Node650_AH_NMBC_com(1,8) Node650_AH_NMBC_com(1,9) Node650_AH_NMBC_com(1,10);...
    Node650_AH_NMBC_com(1,11) Node650_AH_NMBC_com(1,12) Node650_AH_NMBC_com(1,13)];
X650 = sqrt(0.05^2-(Node650_AH_NMBC_com(1,3)+0.25)^2);
teta650 = asin((Node650_AH_NMBC_com(1,3)+0.25)/0.05);
R = [X650/0.05 -(Node650_AH_NMBC_com(1,3)+0.25)/0.05 0;(Node650_AH_NMBC_com(1,3)+0.25)/0.05 X650/0.05 0;...
    0 0 1];
Result_S_650_NMBC_com_trans = R'*Result_S_650_NMBC_com*R;


Result_S_722_NMBC_com = [Node722_AH_NMBC_com(1,5) Node722_AH_NMBC_com(1,6) Node722_AH_NMBC_com(1,7);Node722_AH_NMBC_com(1,8) Node722_AH_NMBC_com(1,9) Node722_AH_NMBC_com(1,10);...
    Node722_AH_NMBC_com(1,11) Node722_AH_NMBC_com(1,12) Node722_AH_NMBC_com(1,13)];
X722 = sqrt(0.05^2-(Node722_AH_NMBC_com(1,3)+0.25)^2);
teta722 = asin((Node722_AH_NMBC_com(1,3)+0.25)/0.05);
R = [X722/0.05 -(Node722_AH_NMBC_com(1,3)+0.25)/0.05 0;(Node722_AH_NMBC_com(1,3)+0.25)/0.05 X722/0.05 0;...
    0 0 1];
Result_S_722_NMBC_com_trans = R'*Result_S_722_NMBC_com*R;


Result_S_794_NMBC_com = [Node794_AH_NMBC_com(1,5) Node794_AH_NMBC_com(1,6) Node794_AH_NMBC_com(1,7);Node794_AH_NMBC_com(1,8) Node794_AH_NMBC_com(1,9) Node794_AH_NMBC_com(1,10);...
    Node794_AH_NMBC_com(1,11) Node794_AH_NMBC_com(1,12) Node794_AH_NMBC_com(1,13)];
X794 = sqrt(0.05^2-(Node794_AH_NMBC_com(1,3)+0.25)^2);
teta794 = asin((Node794_AH_NMBC_com(1,3)+0.25)/0.05);
R = [X794/0.05 -(Node794_AH_NMBC_com(1,3)+0.25)/0.05 0;(Node794_AH_NMBC_com(1,3)+0.25)/0.05 X794/0.05 0;...
    0 0 1];
Result_S_794_NMBC_com_trans = R'*Result_S_794_NMBC_com*R;


Result_S_866_NMBC_com = [Node866_AH_NMBC_com(1,5) Node866_AH_NMBC_com(1,6) Node866_AH_NMBC_com(1,7);Node866_AH_NMBC_com(1,8) Node866_AH_NMBC_com(1,9) Node866_AH_NMBC_com(1,10);...
    Node866_AH_NMBC_com(1,11) Node866_AH_NMBC_com(1,12) Node866_AH_NMBC_com(1,13)];
X866 = sqrt(0.05^2-(Node866_AH_NMBC_com(1,3)+0.25)^2);
teta866 = asin((Node866_AH_NMBC_com(1,3)+0.25)/0.05);
R = [X866/0.05 -(Node866_AH_NMBC_com(1,3)+0.25)/0.05 0;(Node866_AH_NMBC_com(1,3)+0.25)/0.05 X866/0.05 0;...
    0 0 1];
Result_S_866_NMBC_com_trans = R'*Result_S_866_NMBC_com*R;


Result_S_938_NMBC_com = [Node938_AH_NMBC_com(1,5) Node938_AH_NMBC_com(1,6) Node938_AH_NMBC_com(1,7);Node938_AH_NMBC_com(1,8) Node938_AH_NMBC_com(1,9) Node938_AH_NMBC_com(1,10);...
    Node938_AH_NMBC_com(1,11) Node938_AH_NMBC_com(1,12) Node938_AH_NMBC_com(1,13)];
X938 = sqrt(0.05^2-(Node938_AH_NMBC_com(1,3)+0.25)^2);
teta938 = asin((Node938_AH_NMBC_com(1,3)+0.25)/0.05);
R = [X938/0.05 -(Node938_AH_NMBC_com(1,3)+0.25)/0.05 0;(Node938_AH_NMBC_com(1,3)+0.25)/0.05 X938/0.05 0;...
    0 0 1];
Result_S_938_NMBC_com_trans = R'*Result_S_938_NMBC_com*R;


Result_S_1876_NMBC_com = [Node1876_AH_NMBC_com(1,5) Node1876_AH_NMBC_com(1,6) Node1876_AH_NMBC_com(1,7);Node1876_AH_NMBC_com(1,8) Node1876_AH_NMBC_com(1,9) Node1876_AH_NMBC_com(1,10);...
    Node1876_AH_NMBC_com(1,11) Node1876_AH_NMBC_com(1,12) Node1876_AH_NMBC_com(1,13)];
X1876 = sqrt(0.05^2-(Node1876_AH_NMBC_com(1,3)+0.25)^2);
teta1876 = asin((Node1876_AH_NMBC_com(1,3)+0.25)/0.05);
R = [X1876/0.05 -(Node1876_AH_NMBC_com(1,3)+0.25)/0.05 0;(Node1876_AH_NMBC_com(1,3)+0.25)/0.05 X1876/0.05 0;...
    0 0 1];
Result_S_1876_NMBC_com_trans = R'*Result_S_1876_NMBC_com*R;


Result_S_1804_NMBC_com = [Node1804_AH_NMBC_com(1,5) Node1804_AH_NMBC_com(1,6) Node1804_AH_NMBC_com(1,7);Node1804_AH_NMBC_com(1,8) Node1804_AH_NMBC_com(1,9) Node1804_AH_NMBC_com(1,10);...
    Node1804_AH_NMBC_com(1,11) Node1804_AH_NMBC_com(1,12) Node1804_AH_NMBC_com(1,13)];
X1804 = sqrt(0.05^2-(Node1804_AH_NMBC_com(1,3)+0.25)^2);
teta1804 = asin((Node1804_AH_NMBC_com(1,3)+0.25)/0.05);
R = [X1804/0.05 -(Node1804_AH_NMBC_com(1,3)+0.25)/0.05 0;(Node1804_AH_NMBC_com(1,3)+0.25)/0.05 X1804/0.05 0;...
    0 0 1];
Result_S_1804_NMBC_com_trans = R'*Result_S_1804_NMBC_com*R;


Result_S_1732_NMBC_com = [Node1732_AH_NMBC_com(1,5) Node1732_AH_NMBC_com(1,6) Node1732_AH_NMBC_com(1,7);Node1732_AH_NMBC_com(1,8) Node1732_AH_NMBC_com(1,9) Node1732_AH_NMBC_com(1,10);...
    Node1732_AH_NMBC_com(1,11) Node1732_AH_NMBC_com(1,12) Node1732_AH_NMBC_com(1,13)];
X1732 = sqrt(0.05^2-(Node1732_AH_NMBC_com(1,3)+0.25)^2);
teta1732 = asin((Node1732_AH_NMBC_com(1,3)+0.25)/0.05);
R = [X1732/0.05 -(Node1732_AH_NMBC_com(1,3)+0.25)/0.05 0;(Node1732_AH_NMBC_com(1,3)+0.25)/0.05 X1732/0.05 0;...
    0 0 1];
Result_S_1732_NMBC_com_trans = R'*Result_S_1732_NMBC_com*R;


Result_S_1660_NMBC_com = [Node1660_AH_NMBC_com(1,5) Node1660_AH_NMBC_com(1,6) Node1660_AH_NMBC_com(1,7);Node1660_AH_NMBC_com(1,8) Node1660_AH_NMBC_com(1,9) Node1660_AH_NMBC_com(1,10);...
    Node1660_AH_NMBC_com(1,11) Node1660_AH_NMBC_com(1,12) Node1660_AH_NMBC_com(1,13)];
X1660 = sqrt(0.05^2-(Node1660_AH_NMBC_com(1,3)+0.25)^2);
teta1660 = asin((Node1660_AH_NMBC_com(1,3)+0.25)/0.05);
R = [X1660/0.05 -(Node1660_AH_NMBC_com(1,3)+0.25)/0.05 0;(Node1660_AH_NMBC_com(1,3)+0.25)/0.05 X1660/0.05 0;...
    0 0 1];
Result_S_1660_NMBC_com_trans = R'*Result_S_1660_NMBC_com*R;


Result_S_1588_NMBC_com = [Node1588_AH_NMBC_com(1,5) Node1588_AH_NMBC_com(1,6) Node1588_AH_NMBC_com(1,7);Node1588_AH_NMBC_com(1,8) Node1588_AH_NMBC_com(1,9) Node1588_AH_NMBC_com(1,10);...
    Node1588_AH_NMBC_com(1,11) Node1588_AH_NMBC_com(1,12) Node1588_AH_NMBC_com(1,13)];
X1588 = sqrt(0.05^2-(Node1588_AH_NMBC_com(1,3)+0.25)^2);
teta1588 = asin((Node1588_AH_NMBC_com(1,3)+0.25)/0.05);
R = [X1588/0.05 -(Node1588_AH_NMBC_com(1,3)+0.25)/0.05 0;(Node1588_AH_NMBC_com(1,3)+0.25)/0.05 X1588/0.05 0;...
    0 0 1];
Result_S_1588_NMBC_com_trans = R'*Result_S_1588_NMBC_com*R;


Result_S_1516_NMBC_com = [Node1516_AH_NMBC_com(1,5) Node1516_AH_NMBC_com(1,6) Node1516_AH_NMBC_com(1,7);Node1516_AH_NMBC_com(1,8) Node1516_AH_NMBC_com(1,9) Node1516_AH_NMBC_com(1,10);...
    Node1516_AH_NMBC_com(1,11) Node1516_AH_NMBC_com(1,12) Node1516_AH_NMBC_com(1,13)];
X1516 = sqrt(0.05^2-(Node1516_AH_NMBC_com(1,3)+0.25)^2);
teta1516 = asin((Node1516_AH_NMBC_com(1,3)+0.25)/0.05);
R = [X1516/0.05 -(Node1516_AH_NMBC_com(1,3)+0.25)/0.05 0;(Node1516_AH_NMBC_com(1,3)+0.25)/0.05 X1516/0.05 0;...
    0 0 1];
Result_S_1516_NMBC_com_trans = R'*Result_S_1516_NMBC_com*R;


Result_S_1444_NMBC_com = [Node1444_AH_NMBC_com(1,5) Node1444_AH_NMBC_com(1,6) Node1444_AH_NMBC_com(1,7);Node1444_AH_NMBC_com(1,8) Node1444_AH_NMBC_com(1,9) Node1444_AH_NMBC_com(1,10);...
    Node1444_AH_NMBC_com(1,11) Node1444_AH_NMBC_com(1,12) Node1444_AH_NMBC_com(1,13)];
X1444 = sqrt(0.05^2-(Node1444_AH_NMBC_com(1,3)+0.25)^2);
teta1444 = asin((Node1444_AH_NMBC_com(1,3)+0.25)/0.05);
R = [X1444/0.05 -(Node1444_AH_NMBC_com(1,3)+0.25)/0.05 0;(Node1444_AH_NMBC_com(1,3)+0.25)/0.05 X1444/0.05 0;...
    0 0 1];
Result_S_1444_NMBC_com_trans = R'*Result_S_1444_NMBC_com*R;


Result_S_1372_NMBC_com = [Node1372_AH_NMBC_com(1,5) Node1372_AH_NMBC_com(1,6) Node1372_AH_NMBC_com(1,7);Node1372_AH_NMBC_com(1,8) Node1372_AH_NMBC_com(1,9) Node1372_AH_NMBC_com(1,10);...
    Node1372_AH_NMBC_com(1,11) Node1372_AH_NMBC_com(1,12) Node1372_AH_NMBC_com(1,13)];
X1372 = sqrt(0.05^2-(Node1372_AH_NMBC_com(1,3)+0.25)^2);
teta1372 = asin((Node1372_AH_NMBC_com(1,3)+0.25)/0.05);
R = [X1372/0.05 -(Node1372_AH_NMBC_com(1,3)+0.25)/0.05 0;(Node1372_AH_NMBC_com(1,3)+0.25)/0.05 X1372/0.05 0;...
    0 0 1];
Result_S_1372_NMBC_com_trans = R'*Result_S_1372_NMBC_com*R;


Result_S_1300_NMBC_com = [Node1300_AH_NMBC_com(1,5) Node1300_AH_NMBC_com(1,6) Node1300_AH_NMBC_com(1,7);Node1300_AH_NMBC_com(1,8) Node1300_AH_NMBC_com(1,9) Node1300_AH_NMBC_com(1,10);...
    Node1300_AH_NMBC_com(1,11) Node1300_AH_NMBC_com(1,12) Node1300_AH_NMBC_com(1,13)];
X1300 = sqrt(0.05^2-(Node1300_AH_NMBC_com(1,3)+0.25)^2);
teta1300 = asin((Node1300_AH_NMBC_com(1,3)+0.25)/0.05);
R = [X1300/0.05 -(Node1300_AH_NMBC_com(1,3)+0.25)/0.05 0;(Node1300_AH_NMBC_com(1,3)+0.25)/0.05 X1300/0.05 0;...
    0 0 1];
Result_S_1300_NMBC_com_trans = R'*Result_S_1300_NMBC_com*R;


Result_S_1228_NMBC_com = [Node1228_AH_NMBC_com(1,5) Node1228_AH_NMBC_com(1,6) Node1228_AH_NMBC_com(1,7);Node1228_AH_NMBC_com(1,8) Node1228_AH_NMBC_com(1,9) Node1228_AH_NMBC_com(1,10);...
    Node1228_AH_NMBC_com(1,11) Node1228_AH_NMBC_com(1,12) Node1228_AH_NMBC_com(1,13)];
X1228 = sqrt(0.05^2-(Node1228_AH_NMBC_com(1,3)+0.25)^2);
teta1228 = asin((Node1228_AH_NMBC_com(1,3)+0.25)/0.05);
R = [X1228/0.05 -(Node1228_AH_NMBC_com(1,3)+0.25)/0.05 0;(Node1228_AH_NMBC_com(1,3)+0.25)/0.05 X1228/0.05 0;...
    0 0 1];
Result_S_1228_NMBC_com_trans = R'*Result_S_1228_NMBC_com*R;


Result_S_1156_NMBC_com = [Node1156_AH_NMBC_com(1,5) Node1156_AH_NMBC_com(1,6) Node1156_AH_NMBC_com(1,7);Node1156_AH_NMBC_com(1,8) Node1156_AH_NMBC_com(1,9) Node1156_AH_NMBC_com(1,10);...
    Node1156_AH_NMBC_com(1,11) Node1156_AH_NMBC_com(1,12) Node1156_AH_NMBC_com(1,13)];
X1156 = sqrt(0.05^2-(Node1156_AH_NMBC_com(1,3)+0.25)^2);
teta1156 = asin((Node1156_AH_NMBC_com(1,3)+0.25)/0.05);
R = [X1156/0.05 -(Node1156_AH_NMBC_com(1,3)+0.25)/0.05 0;(Node1156_AH_NMBC_com(1,3)+0.25)/0.05 X1156/0.05 0;...
    0 0 1];
Result_S_1156_NMBC_com_trans = R'*Result_S_1156_NMBC_com*R;


Result_S_1016_NMBC_com = [Node1016_AH_NMBC_com(1,5) Node1016_AH_NMBC_com(1,6) Node1016_AH_NMBC_com(1,7);Node1016_AH_NMBC_com(1,8) Node1016_AH_NMBC_com(1,9) Node1016_AH_NMBC_com(1,10);...
    Node1016_AH_NMBC_com(1,11) Node1016_AH_NMBC_com(1,12) Node1016_AH_NMBC_com(1,13)];
X1016 = sqrt(0.05^2-(Node1016_AH_NMBC_com(1,3)+0.25)^2);
teta1016 = asin((Node1016_AH_NMBC_com(1,3)+0.25)/0.05);
R = [X1016/0.05 -(Node1016_AH_NMBC_com(1,3)+0.25)/0.05 0;(Node1016_AH_NMBC_com(1,3)+0.25)/0.05 X1016/0.05 0;...
    0 0 1];
Result_S_1016_NMBC_com_trans = R'*Result_S_1016_NMBC_com*R;


Result_S_1012_NMBC_com = [Node1012_AH_NMBC_com(1,5) Node1012_AH_NMBC_com(1,6) Node1012_AH_NMBC_com(1,7);Node1012_AH_NMBC_com(1,8) Node1012_AH_NMBC_com(1,9) Node1012_AH_NMBC_com(1,10);...
    Node1012_AH_NMBC_com(1,11) Node1012_AH_NMBC_com(1,12) Node1012_AH_NMBC_com(1,13)];
X1012 = sqrt(0.05^2-(Node1012_AH_NMBC_com(1,3)+0.25)^2);
teta1012 = asin((Node1012_AH_NMBC_com(1,3)+0.25)/0.05);
R = [X1012/0.05 -(Node1012_AH_NMBC_com(1,3)+0.25)/0.05 0;(Node1012_AH_NMBC_com(1,3)+0.25)/0.05 X1012/0.05 0;...
    0 0 1];
Result_S_1012_NMBC_com_trans = R'*Result_S_1012_NMBC_com*R;


Result_S_teta_teta_NMBC_com = [teta2 Result_S_2_NMBC_com_trans(2,2)/P_applied;teta6 Result_S_6_NMBC_com_trans(2,2)/P_applied;teta146 Result_S_146_NMBC_com_trans(2,2)/P_applied;...
    teta218 Result_S_218_NMBC_com_trans(2,2)/P_applied;teta290 Result_S_290_NMBC_com_trans(2,2)/P_applied;teta362 Result_S_362_NMBC_com_trans(2,2)/P_applied;...
    teta434 Result_S_434_NMBC_com_trans(2,2)/P_applied;teta506 Result_S_506_NMBC_com_trans(2,2)/P_applied;teta578 Result_S_578_NMBC_com_trans(2,2)/P_applied;...
    teta650 Result_S_650_NMBC_com_trans(2,2)/P_applied;teta722 Result_S_722_NMBC_com_trans(2,2)/P_applied;teta794 Result_S_794_NMBC_com_trans(2,2)/P_applied;...
    teta866 Result_S_866_NMBC_com_trans(2,2)/P_applied;teta938 Result_S_938_NMBC_com_trans(2,2)/P_applied;teta1876 Result_S_1876_NMBC_com_trans(2,2)/P_applied;...
    teta1804 Result_S_1804_NMBC_com_trans(2,2)/P_applied;teta1732 Result_S_1732_NMBC_com_trans(2,2)/P_applied;teta1660 Result_S_1660_NMBC_com_trans(2,2)/P_applied;...
    teta1588 Result_S_1588_NMBC_com_trans(2,2)/P_applied;teta1516 Result_S_1516_NMBC_com_trans(2,2)/P_applied;teta1444 Result_S_1444_NMBC_com_trans(2,2)/P_applied;...
    teta1372 Result_S_1372_NMBC_com_trans(2,2)/P_applied;teta1300 Result_S_1300_NMBC_com_trans(2,2)/P_applied;teta1228 Result_S_1228_NMBC_com_trans(2,2)/P_applied;...
    teta1156 Result_S_1156_NMBC_com_trans(2,2)/P_applied;teta1016 Result_S_1016_NMBC_com_trans(2,2)/P_applied;teta1012 Result_S_1012_NMBC_com_trans(2,2)/P_applied];



Result_S_2_NMBC_rot = [Node2_AH_NMBC_rot(1,5) Node2_AH_NMBC_rot(1,6) Node2_AH_NMBC_rot(1,7);Node2_AH_NMBC_rot(1,8) Node2_AH_NMBC_rot(1,9) Node2_AH_NMBC_rot(1,10);...
    Node2_AH_NMBC_rot(1,11) Node2_AH_NMBC_rot(1,12) Node2_AH_NMBC_rot(1,13)];
X2 = sqrt(0.05^2-(Node2_AH_NMBC_rot(1,3)+0.25)^2);
teta2 = asin((Node2_AH_NMBC_rot(1,3)+0.25)/0.05);
R = [X2/0.05 -(Node2_AH_NMBC_rot(1,3)+0.25)/0.05 0;(Node2_AH_NMBC_rot(1,3)+0.25)/0.05 X2/0.05 0;...
    0 0 1];
Result_S_2_NMBC_rot_trans = R'*Result_S_2_NMBC_rot*R;


Result_S_6_NMBC_rot = [Node6_AH_NMBC_rot(1,5) Node6_AH_NMBC_rot(1,6) Node6_AH_NMBC_rot(1,7);Node6_AH_NMBC_rot(1,8) Node6_AH_NMBC_rot(1,9) Node6_AH_NMBC_rot(1,10);...
    Node6_AH_NMBC_rot(1,11) Node6_AH_NMBC_rot(1,12) Node6_AH_NMBC_rot(1,13)];
X6 = sqrt(0.05^2-(Node6_AH_NMBC_rot(1,3)+0.25)^2);
teta6 = asin((Node6_AH_NMBC_rot(1,3)+0.25)/0.05);
R = [X6/0.05 -(Node6_AH_NMBC_rot(1,3)+0.25)/0.05 0;(Node6_AH_NMBC_rot(1,3)+0.25)/0.05 X6/0.05 0;...
    0 0 1];
Result_S_6_NMBC_rot_trans = R'*Result_S_6_NMBC_rot*R;


Result_S_146_NMBC_rot = [Node146_AH_NMBC_rot(1,5) Node146_AH_NMBC_rot(1,6) Node146_AH_NMBC_rot(1,7);Node146_AH_NMBC_rot(1,8) Node146_AH_NMBC_rot(1,9) Node146_AH_NMBC_rot(1,10);...
    Node146_AH_NMBC_rot(1,11) Node146_AH_NMBC_rot(1,12) Node146_AH_NMBC_rot(1,13)];
X146 = sqrt(0.05^2-(Node146_AH_NMBC_rot(1,3)+0.25)^2);
teta146 = asin((Node146_AH_NMBC_rot(1,3)+0.25)/0.05);
R = [X146/0.05 -(Node146_AH_NMBC_rot(1,3)+0.25)/0.05 0;(Node146_AH_NMBC_rot(1,3)+0.25)/0.05 X146/0.05 0;...
    0 0 1];
Result_S_146_NMBC_rot_trans = R'*Result_S_146_NMBC_rot*R;


Result_S_218_NMBC_rot = [Node218_AH_NMBC_rot(1,5) Node218_AH_NMBC_rot(1,6) Node218_AH_NMBC_rot(1,7);Node218_AH_NMBC_rot(1,8) Node218_AH_NMBC_rot(1,9) Node218_AH_NMBC_rot(1,10);...
    Node218_AH_NMBC_rot(1,11) Node218_AH_NMBC_rot(1,12) Node218_AH_NMBC_rot(1,13)];
X218 = sqrt(0.05^2-(Node218_AH_NMBC_rot(1,3)+0.25)^2);
teta218 = asin((Node218_AH_NMBC_rot(1,3)+0.25)/0.05);
R = [X218/0.05 -(Node218_AH_NMBC_rot(1,3)+0.25)/0.05 0;(Node218_AH_NMBC_rot(1,3)+0.25)/0.05 X218/0.05 0;...
    0 0 1];
Result_S_218_NMBC_rot_trans = R'*Result_S_218_NMBC_rot*R;


Result_S_290_NMBC_rot = [Node290_AH_NMBC_rot(1,5) Node290_AH_NMBC_rot(1,6) Node290_AH_NMBC_rot(1,7);Node290_AH_NMBC_rot(1,8) Node290_AH_NMBC_rot(1,9) Node290_AH_NMBC_rot(1,10);...
    Node290_AH_NMBC_rot(1,11) Node290_AH_NMBC_rot(1,12) Node290_AH_NMBC_rot(1,13)];
X290 = sqrt(0.05^2-(Node290_AH_NMBC_rot(1,3)+0.25)^2);
teta290 = asin((Node290_AH_NMBC_rot(1,3)+0.25)/0.05);
R = [X290/0.05 -(Node290_AH_NMBC_rot(1,3)+0.25)/0.05 0;(Node290_AH_NMBC_rot(1,3)+0.25)/0.05 X290/0.05 0;...
    0 0 1];
Result_S_290_NMBC_rot_trans = R'*Result_S_290_NMBC_rot*R;


Result_S_362_NMBC_rot = [Node362_AH_NMBC_rot(1,5) Node362_AH_NMBC_rot(1,6) Node362_AH_NMBC_rot(1,7);Node362_AH_NMBC_rot(1,8) Node362_AH_NMBC_rot(1,9) Node362_AH_NMBC_rot(1,10);...
    Node362_AH_NMBC_rot(1,11) Node362_AH_NMBC_rot(1,12) Node362_AH_NMBC_rot(1,13)];
X362 = sqrt(0.05^2-(Node362_AH_NMBC_rot(1,3)+0.25)^2);
teta362 = asin((Node362_AH_NMBC_rot(1,3)+0.25)/0.05);
R = [X362/0.05 -(Node362_AH_NMBC_rot(1,3)+0.25)/0.05 0;(Node362_AH_NMBC_rot(1,3)+0.25)/0.05 X362/0.05 0;...
    0 0 1];
Result_S_362_NMBC_rot_trans = R'*Result_S_362_NMBC_rot*R;


Result_S_434_NMBC_rot = [Node434_AH_NMBC_rot(1,5) Node434_AH_NMBC_rot(1,6) Node434_AH_NMBC_rot(1,7);Node434_AH_NMBC_rot(1,8) Node434_AH_NMBC_rot(1,9) Node434_AH_NMBC_rot(1,10);...
    Node434_AH_NMBC_rot(1,11) Node434_AH_NMBC_rot(1,12) Node434_AH_NMBC_rot(1,13)];
X434 = sqrt(0.05^2-(Node434_AH_NMBC_rot(1,3)+0.25)^2);
teta434 = asin((Node434_AH_NMBC_rot(1,3)+0.25)/0.05);
R = [X434/0.05 -(Node434_AH_NMBC_rot(1,3)+0.25)/0.05 0;(Node434_AH_NMBC_rot(1,3)+0.25)/0.05 X434/0.05 0;...
    0 0 1];
Result_S_434_NMBC_rot_trans = R'*Result_S_434_NMBC_rot*R;


Result_S_506_NMBC_rot = [Node506_AH_NMBC_rot(1,5) Node506_AH_NMBC_rot(1,6) Node506_AH_NMBC_rot(1,7);Node506_AH_NMBC_rot(1,8) Node506_AH_NMBC_rot(1,9) Node506_AH_NMBC_rot(1,10);...
    Node506_AH_NMBC_rot(1,11) Node506_AH_NMBC_rot(1,12) Node506_AH_NMBC_rot(1,13)];
X506 = sqrt(0.05^2-(Node506_AH_NMBC_rot(1,3)+0.25)^2);
teta506 = asin((Node506_AH_NMBC_rot(1,3)+0.25)/0.05);
R = [X506/0.05 -(Node506_AH_NMBC_rot(1,3)+0.25)/0.05 0;(Node506_AH_NMBC_rot(1,3)+0.25)/0.05 X506/0.05 0;...
    0 0 1];
Result_S_506_NMBC_rot_trans = R'*Result_S_506_NMBC_rot*R;


Result_S_578_NMBC_rot = [Node578_AH_NMBC_rot(1,5) Node578_AH_NMBC_rot(1,6) Node578_AH_NMBC_rot(1,7);Node578_AH_NMBC_rot(1,8) Node578_AH_NMBC_rot(1,9) Node578_AH_NMBC_rot(1,10);...
    Node578_AH_NMBC_rot(1,11) Node578_AH_NMBC_rot(1,12) Node578_AH_NMBC_rot(1,13)];
X578 = sqrt(0.05^2-(Node578_AH_NMBC_rot(1,3)+0.25)^2);
teta578 = asin((Node578_AH_NMBC_rot(1,3)+0.25)/0.05);
R = [X578/0.05 -(Node578_AH_NMBC_rot(1,3)+0.25)/0.05 0;(Node578_AH_NMBC_rot(1,3)+0.25)/0.05 X578/0.05 0;...
    0 0 1];
Result_S_578_NMBC_rot_trans = R'*Result_S_578_NMBC_rot*R;


Result_S_650_NMBC_rot = [Node650_AH_NMBC_rot(1,5) Node650_AH_NMBC_rot(1,6) Node650_AH_NMBC_rot(1,7);Node650_AH_NMBC_rot(1,8) Node650_AH_NMBC_rot(1,9) Node650_AH_NMBC_rot(1,10);...
    Node650_AH_NMBC_rot(1,11) Node650_AH_NMBC_rot(1,12) Node650_AH_NMBC_rot(1,13)];
X650 = sqrt(0.05^2-(Node650_AH_NMBC_rot(1,3)+0.25)^2);
teta650 = asin((Node650_AH_NMBC_rot(1,3)+0.25)/0.05);
R = [X650/0.05 -(Node650_AH_NMBC_rot(1,3)+0.25)/0.05 0;(Node650_AH_NMBC_rot(1,3)+0.25)/0.05 X650/0.05 0;...
    0 0 1];
Result_S_650_NMBC_rot_trans = R'*Result_S_650_NMBC_rot*R;


Result_S_722_NMBC_rot = [Node722_AH_NMBC_rot(1,5) Node722_AH_NMBC_rot(1,6) Node722_AH_NMBC_rot(1,7);Node722_AH_NMBC_rot(1,8) Node722_AH_NMBC_rot(1,9) Node722_AH_NMBC_rot(1,10);...
    Node722_AH_NMBC_rot(1,11) Node722_AH_NMBC_rot(1,12) Node722_AH_NMBC_rot(1,13)];
X722 = sqrt(0.05^2-(Node722_AH_NMBC_rot(1,3)+0.25)^2);
teta722 = asin((Node722_AH_NMBC_rot(1,3)+0.25)/0.05);
R = [X722/0.05 -(Node722_AH_NMBC_rot(1,3)+0.25)/0.05 0;(Node722_AH_NMBC_rot(1,3)+0.25)/0.05 X722/0.05 0;...
    0 0 1];
Result_S_722_NMBC_rot_trans = R'*Result_S_722_NMBC_rot*R;


Result_S_794_NMBC_rot = [Node794_AH_NMBC_rot(1,5) Node794_AH_NMBC_rot(1,6) Node794_AH_NMBC_rot(1,7);Node794_AH_NMBC_rot(1,8) Node794_AH_NMBC_rot(1,9) Node794_AH_NMBC_rot(1,10);...
    Node794_AH_NMBC_rot(1,11) Node794_AH_NMBC_rot(1,12) Node794_AH_NMBC_rot(1,13)];
X794 = sqrt(0.05^2-(Node794_AH_NMBC_rot(1,3)+0.25)^2);
teta794 = asin((Node794_AH_NMBC_rot(1,3)+0.25)/0.05);
R = [X794/0.05 -(Node794_AH_NMBC_rot(1,3)+0.25)/0.05 0;(Node794_AH_NMBC_rot(1,3)+0.25)/0.05 X794/0.05 0;...
    0 0 1];
Result_S_794_NMBC_rot_trans = R'*Result_S_794_NMBC_rot*R;


Result_S_866_NMBC_rot = [Node866_AH_NMBC_rot(1,5) Node866_AH_NMBC_rot(1,6) Node866_AH_NMBC_rot(1,7);Node866_AH_NMBC_rot(1,8) Node866_AH_NMBC_rot(1,9) Node866_AH_NMBC_rot(1,10);...
    Node866_AH_NMBC_rot(1,11) Node866_AH_NMBC_rot(1,12) Node866_AH_NMBC_rot(1,13)];
X866 = sqrt(0.05^2-(Node866_AH_NMBC_rot(1,3)+0.25)^2);
teta866 = asin((Node866_AH_NMBC_rot(1,3)+0.25)/0.05);
R = [X866/0.05 -(Node866_AH_NMBC_rot(1,3)+0.25)/0.05 0;(Node866_AH_NMBC_rot(1,3)+0.25)/0.05 X866/0.05 0;...
    0 0 1];
Result_S_866_NMBC_rot_trans = R'*Result_S_866_NMBC_rot*R;


Result_S_938_NMBC_rot = [Node938_AH_NMBC_rot(1,5) Node938_AH_NMBC_rot(1,6) Node938_AH_NMBC_rot(1,7);Node938_AH_NMBC_rot(1,8) Node938_AH_NMBC_rot(1,9) Node938_AH_NMBC_rot(1,10);...
    Node938_AH_NMBC_rot(1,11) Node938_AH_NMBC_rot(1,12) Node938_AH_NMBC_rot(1,13)];
X938 = sqrt(0.05^2-(Node938_AH_NMBC_rot(1,3)+0.25)^2);
teta938 = asin((Node938_AH_NMBC_rot(1,3)+0.25)/0.05);
R = [X938/0.05 -(Node938_AH_NMBC_rot(1,3)+0.25)/0.05 0;(Node938_AH_NMBC_rot(1,3)+0.25)/0.05 X938/0.05 0;...
    0 0 1];
Result_S_938_NMBC_rot_trans = R'*Result_S_938_NMBC_rot*R;


Result_S_1876_NMBC_rot = [Node1876_AH_NMBC_rot(1,5) Node1876_AH_NMBC_rot(1,6) Node1876_AH_NMBC_rot(1,7);Node1876_AH_NMBC_rot(1,8) Node1876_AH_NMBC_rot(1,9) Node1876_AH_NMBC_rot(1,10);...
    Node1876_AH_NMBC_rot(1,11) Node1876_AH_NMBC_rot(1,12) Node1876_AH_NMBC_rot(1,13)];
X1876 = sqrt(0.05^2-(Node1876_AH_NMBC_rot(1,3)+0.25)^2);
teta1876 = asin((Node1876_AH_NMBC_rot(1,3)+0.25)/0.05);
R = [X1876/0.05 -(Node1876_AH_NMBC_rot(1,3)+0.25)/0.05 0;(Node1876_AH_NMBC_rot(1,3)+0.25)/0.05 X1876/0.05 0;...
    0 0 1];
Result_S_1876_NMBC_rot_trans = R'*Result_S_1876_NMBC_rot*R;


Result_S_1804_NMBC_rot = [Node1804_AH_NMBC_rot(1,5) Node1804_AH_NMBC_rot(1,6) Node1804_AH_NMBC_rot(1,7);Node1804_AH_NMBC_rot(1,8) Node1804_AH_NMBC_rot(1,9) Node1804_AH_NMBC_rot(1,10);...
    Node1804_AH_NMBC_rot(1,11) Node1804_AH_NMBC_rot(1,12) Node1804_AH_NMBC_rot(1,13)];
X1804 = sqrt(0.05^2-(Node1804_AH_NMBC_rot(1,3)+0.25)^2);
teta1804 = asin((Node1804_AH_NMBC_rot(1,3)+0.25)/0.05);
R = [X1804/0.05 -(Node1804_AH_NMBC_rot(1,3)+0.25)/0.05 0;(Node1804_AH_NMBC_rot(1,3)+0.25)/0.05 X1804/0.05 0;...
    0 0 1];
Result_S_1804_NMBC_rot_trans = R'*Result_S_1804_NMBC_rot*R;


Result_S_1732_NMBC_rot = [Node1732_AH_NMBC_rot(1,5) Node1732_AH_NMBC_rot(1,6) Node1732_AH_NMBC_rot(1,7);Node1732_AH_NMBC_rot(1,8) Node1732_AH_NMBC_rot(1,9) Node1732_AH_NMBC_rot(1,10);...
    Node1732_AH_NMBC_rot(1,11) Node1732_AH_NMBC_rot(1,12) Node1732_AH_NMBC_rot(1,13)];
X1732 = sqrt(0.05^2-(Node1732_AH_NMBC_rot(1,3)+0.25)^2);
teta1732 = asin((Node1732_AH_NMBC_rot(1,3)+0.25)/0.05);
R = [X1732/0.05 -(Node1732_AH_NMBC_rot(1,3)+0.25)/0.05 0;(Node1732_AH_NMBC_rot(1,3)+0.25)/0.05 X1732/0.05 0;...
    0 0 1];
Result_S_1732_NMBC_rot_trans = R'*Result_S_1732_NMBC_rot*R;


Result_S_1660_NMBC_rot = [Node1660_AH_NMBC_rot(1,5) Node1660_AH_NMBC_rot(1,6) Node1660_AH_NMBC_rot(1,7);Node1660_AH_NMBC_rot(1,8) Node1660_AH_NMBC_rot(1,9) Node1660_AH_NMBC_rot(1,10);...
    Node1660_AH_NMBC_rot(1,11) Node1660_AH_NMBC_rot(1,12) Node1660_AH_NMBC_rot(1,13)];
X1660 = sqrt(0.05^2-(Node1660_AH_NMBC_rot(1,3)+0.25)^2);
teta1660 = asin((Node1660_AH_NMBC_rot(1,3)+0.25)/0.05);
R = [X1660/0.05 -(Node1660_AH_NMBC_rot(1,3)+0.25)/0.05 0;(Node1660_AH_NMBC_rot(1,3)+0.25)/0.05 X1660/0.05 0;...
    0 0 1];
Result_S_1660_NMBC_rot_trans = R'*Result_S_1660_NMBC_rot*R;


Result_S_1588_NMBC_rot = [Node1588_AH_NMBC_rot(1,5) Node1588_AH_NMBC_rot(1,6) Node1588_AH_NMBC_rot(1,7);Node1588_AH_NMBC_rot(1,8) Node1588_AH_NMBC_rot(1,9) Node1588_AH_NMBC_rot(1,10);...
    Node1588_AH_NMBC_rot(1,11) Node1588_AH_NMBC_rot(1,12) Node1588_AH_NMBC_rot(1,13)];
X1588 = sqrt(0.05^2-(Node1588_AH_NMBC_rot(1,3)+0.25)^2);
teta1588 = asin((Node1588_AH_NMBC_rot(1,3)+0.25)/0.05);
R = [X1588/0.05 -(Node1588_AH_NMBC_rot(1,3)+0.25)/0.05 0;(Node1588_AH_NMBC_rot(1,3)+0.25)/0.05 X1588/0.05 0;...
    0 0 1];
Result_S_1588_NMBC_rot_trans = R'*Result_S_1588_NMBC_rot*R;


Result_S_1516_NMBC_rot = [Node1516_AH_NMBC_rot(1,5) Node1516_AH_NMBC_rot(1,6) Node1516_AH_NMBC_rot(1,7);Node1516_AH_NMBC_rot(1,8) Node1516_AH_NMBC_rot(1,9) Node1516_AH_NMBC_rot(1,10);...
    Node1516_AH_NMBC_rot(1,11) Node1516_AH_NMBC_rot(1,12) Node1516_AH_NMBC_rot(1,13)];
X1516 = sqrt(0.05^2-(Node1516_AH_NMBC_rot(1,3)+0.25)^2);
teta1516 = asin((Node1516_AH_NMBC_rot(1,3)+0.25)/0.05);
R = [X1516/0.05 -(Node1516_AH_NMBC_rot(1,3)+0.25)/0.05 0;(Node1516_AH_NMBC_rot(1,3)+0.25)/0.05 X1516/0.05 0;...
    0 0 1];
Result_S_1516_NMBC_rot_trans = R'*Result_S_1516_NMBC_rot*R;


Result_S_1444_NMBC_rot = [Node1444_AH_NMBC_rot(1,5) Node1444_AH_NMBC_rot(1,6) Node1444_AH_NMBC_rot(1,7);Node1444_AH_NMBC_rot(1,8) Node1444_AH_NMBC_rot(1,9) Node1444_AH_NMBC_rot(1,10);...
    Node1444_AH_NMBC_rot(1,11) Node1444_AH_NMBC_rot(1,12) Node1444_AH_NMBC_rot(1,13)];
X1444 = sqrt(0.05^2-(Node1444_AH_NMBC_rot(1,3)+0.25)^2);
teta1444 = asin((Node1444_AH_NMBC_rot(1,3)+0.25)/0.05);
R = [X1444/0.05 -(Node1444_AH_NMBC_rot(1,3)+0.25)/0.05 0;(Node1444_AH_NMBC_rot(1,3)+0.25)/0.05 X1444/0.05 0;...
    0 0 1];
Result_S_1444_NMBC_rot_trans = R'*Result_S_1444_NMBC_rot*R;


Result_S_1372_NMBC_rot = [Node1372_AH_NMBC_rot(1,5) Node1372_AH_NMBC_rot(1,6) Node1372_AH_NMBC_rot(1,7);Node1372_AH_NMBC_rot(1,8) Node1372_AH_NMBC_rot(1,9) Node1372_AH_NMBC_rot(1,10);...
    Node1372_AH_NMBC_rot(1,11) Node1372_AH_NMBC_rot(1,12) Node1372_AH_NMBC_rot(1,13)];
X1372 = sqrt(0.05^2-(Node1372_AH_NMBC_rot(1,3)+0.25)^2);
teta1372 = asin((Node1372_AH_NMBC_rot(1,3)+0.25)/0.05);
R = [X1372/0.05 -(Node1372_AH_NMBC_rot(1,3)+0.25)/0.05 0;(Node1372_AH_NMBC_rot(1,3)+0.25)/0.05 X1372/0.05 0;...
    0 0 1];
Result_S_1372_NMBC_rot_trans = R'*Result_S_1372_NMBC_rot*R;


Result_S_1300_NMBC_rot = [Node1300_AH_NMBC_rot(1,5) Node1300_AH_NMBC_rot(1,6) Node1300_AH_NMBC_rot(1,7);Node1300_AH_NMBC_rot(1,8) Node1300_AH_NMBC_rot(1,9) Node1300_AH_NMBC_rot(1,10);...
    Node1300_AH_NMBC_rot(1,11) Node1300_AH_NMBC_rot(1,12) Node1300_AH_NMBC_rot(1,13)];
X1300 = sqrt(0.05^2-(Node1300_AH_NMBC_rot(1,3)+0.25)^2);
teta1300 = asin((Node1300_AH_NMBC_rot(1,3)+0.25)/0.05);
R = [X1300/0.05 -(Node1300_AH_NMBC_rot(1,3)+0.25)/0.05 0;(Node1300_AH_NMBC_rot(1,3)+0.25)/0.05 X1300/0.05 0;...
    0 0 1];
Result_S_1300_NMBC_rot_trans = R'*Result_S_1300_NMBC_rot*R;


Result_S_1228_NMBC_rot = [Node1228_AH_NMBC_rot(1,5) Node1228_AH_NMBC_rot(1,6) Node1228_AH_NMBC_rot(1,7);Node1228_AH_NMBC_rot(1,8) Node1228_AH_NMBC_rot(1,9) Node1228_AH_NMBC_rot(1,10);...
    Node1228_AH_NMBC_rot(1,11) Node1228_AH_NMBC_rot(1,12) Node1228_AH_NMBC_rot(1,13)];
X1228 = sqrt(0.05^2-(Node1228_AH_NMBC_rot(1,3)+0.25)^2);
teta1228 = asin((Node1228_AH_NMBC_rot(1,3)+0.25)/0.05);
R = [X1228/0.05 -(Node1228_AH_NMBC_rot(1,3)+0.25)/0.05 0;(Node1228_AH_NMBC_rot(1,3)+0.25)/0.05 X1228/0.05 0;...
    0 0 1];
Result_S_1228_NMBC_rot_trans = R'*Result_S_1228_NMBC_rot*R;


Result_S_1156_NMBC_rot = [Node1156_AH_NMBC_rot(1,5) Node1156_AH_NMBC_rot(1,6) Node1156_AH_NMBC_rot(1,7);Node1156_AH_NMBC_rot(1,8) Node1156_AH_NMBC_rot(1,9) Node1156_AH_NMBC_rot(1,10);...
    Node1156_AH_NMBC_rot(1,11) Node1156_AH_NMBC_rot(1,12) Node1156_AH_NMBC_rot(1,13)];
X1156 = sqrt(0.05^2-(Node1156_AH_NMBC_rot(1,3)+0.25)^2);
teta1156 = asin((Node1156_AH_NMBC_rot(1,3)+0.25)/0.05);
R = [X1156/0.05 -(Node1156_AH_NMBC_rot(1,3)+0.25)/0.05 0;(Node1156_AH_NMBC_rot(1,3)+0.25)/0.05 X1156/0.05 0;...
    0 0 1];
Result_S_1156_NMBC_rot_trans = R'*Result_S_1156_NMBC_rot*R;


Result_S_1016_NMBC_rot = [Node1016_AH_NMBC_rot(1,5) Node1016_AH_NMBC_rot(1,6) Node1016_AH_NMBC_rot(1,7);Node1016_AH_NMBC_rot(1,8) Node1016_AH_NMBC_rot(1,9) Node1016_AH_NMBC_rot(1,10);...
    Node1016_AH_NMBC_rot(1,11) Node1016_AH_NMBC_rot(1,12) Node1016_AH_NMBC_rot(1,13)];
X1016 = sqrt(0.05^2-(Node1016_AH_NMBC_rot(1,3)+0.25)^2);
teta1016 = asin((Node1016_AH_NMBC_rot(1,3)+0.25)/0.05);
R = [X1016/0.05 -(Node1016_AH_NMBC_rot(1,3)+0.25)/0.05 0;(Node1016_AH_NMBC_rot(1,3)+0.25)/0.05 X1016/0.05 0;...
    0 0 1];
Result_S_1016_NMBC_rot_trans = R'*Result_S_1016_NMBC_rot*R;


Result_S_1012_NMBC_rot = [Node1012_AH_NMBC_rot(1,5) Node1012_AH_NMBC_rot(1,6) Node1012_AH_NMBC_rot(1,7);Node1012_AH_NMBC_rot(1,8) Node1012_AH_NMBC_rot(1,9) Node1012_AH_NMBC_rot(1,10);...
    Node1012_AH_NMBC_rot(1,11) Node1012_AH_NMBC_rot(1,12) Node1012_AH_NMBC_rot(1,13)];
X1012 = sqrt(0.05^2-(Node1012_AH_NMBC_rot(1,3)+0.25)^2);
teta1012 = asin((Node1012_AH_NMBC_rot(1,3)+0.25)/0.05);
R = [X1012/0.05 -(Node1012_AH_NMBC_rot(1,3)+0.25)/0.05 0;(Node1012_AH_NMBC_rot(1,3)+0.25)/0.05 X1012/0.05 0;...
    0 0 1];
Result_S_1012_NMBC_rot_trans = R'*Result_S_1012_NMBC_rot*R;


Result_S_teta_teta_NMBC_rot = [teta2 Result_S_2_NMBC_rot_trans(2,2)/P_applied;teta6 Result_S_6_NMBC_rot_trans(2,2)/P_applied;teta146 Result_S_146_NMBC_rot_trans(2,2)/P_applied;...
    teta218 Result_S_218_NMBC_rot_trans(2,2)/P_applied;teta290 Result_S_290_NMBC_rot_trans(2,2)/P_applied;teta362 Result_S_362_NMBC_rot_trans(2,2)/P_applied;...
    teta434 Result_S_434_NMBC_rot_trans(2,2)/P_applied;teta506 Result_S_506_NMBC_rot_trans(2,2)/P_applied;teta578 Result_S_578_NMBC_rot_trans(2,2)/P_applied;...
    teta650 Result_S_650_NMBC_rot_trans(2,2)/P_applied;teta722 Result_S_722_NMBC_rot_trans(2,2)/P_applied;teta794 Result_S_794_NMBC_rot_trans(2,2)/P_applied;...
    teta866 Result_S_866_NMBC_rot_trans(2,2)/P_applied;teta938 Result_S_938_NMBC_rot_trans(2,2)/P_applied;teta1876 Result_S_1876_NMBC_rot_trans(2,2)/P_applied;...
    teta1804 Result_S_1804_NMBC_rot_trans(2,2)/P_applied;teta1732 Result_S_1732_NMBC_rot_trans(2,2)/P_applied;teta1660 Result_S_1660_NMBC_rot_trans(2,2)/P_applied;...
    teta1588 Result_S_1588_NMBC_rot_trans(2,2)/P_applied;teta1516 Result_S_1516_NMBC_rot_trans(2,2)/P_applied;teta1444 Result_S_1444_NMBC_rot_trans(2,2)/P_applied;...
    teta1372 Result_S_1372_NMBC_rot_trans(2,2)/P_applied;teta1300 Result_S_1300_NMBC_rot_trans(2,2)/P_applied;teta1228 Result_S_1228_NMBC_rot_trans(2,2)/P_applied;...
    teta1156 Result_S_1156_NMBC_rot_trans(2,2)/P_applied;teta1016 Result_S_1016_NMBC_rot_trans(2,2)/P_applied;teta1012 Result_S_1012_NMBC_rot_trans(2,2)/P_applied];



Result_S_2_NMBC_str = [Node2_AH_NMBC_str(1,5) Node2_AH_NMBC_str(1,6) Node2_AH_NMBC_str(1,7);Node2_AH_NMBC_str(1,8) Node2_AH_NMBC_str(1,9) Node2_AH_NMBC_str(1,10);...
    Node2_AH_NMBC_str(1,11) Node2_AH_NMBC_str(1,12) Node2_AH_NMBC_str(1,13)];
X2 = sqrt(0.05^2-(Node2_AH_NMBC_str(1,3)+0.25)^2);
teta2 = asin((Node2_AH_NMBC_str(1,3)+0.25)/0.05);
R = [X2/0.05 -(Node2_AH_NMBC_str(1,3)+0.25)/0.05 0;(Node2_AH_NMBC_str(1,3)+0.25)/0.05 X2/0.05 0;...
    0 0 1];
Result_S_2_NMBC_str_trans = R'*Result_S_2_NMBC_str*R;


Result_S_6_NMBC_str = [Node6_AH_NMBC_str(1,5) Node6_AH_NMBC_str(1,6) Node6_AH_NMBC_str(1,7);Node6_AH_NMBC_str(1,8) Node6_AH_NMBC_str(1,9) Node6_AH_NMBC_str(1,10);...
    Node6_AH_NMBC_str(1,11) Node6_AH_NMBC_str(1,12) Node6_AH_NMBC_str(1,13)];
X6 = sqrt(0.05^2-(Node6_AH_NMBC_str(1,3)+0.25)^2);
teta6 = asin((Node6_AH_NMBC_str(1,3)+0.25)/0.05);
R = [X6/0.05 -(Node6_AH_NMBC_str(1,3)+0.25)/0.05 0;(Node6_AH_NMBC_str(1,3)+0.25)/0.05 X6/0.05 0;...
    0 0 1];
Result_S_6_NMBC_str_trans = R'*Result_S_6_NMBC_str*R;


Result_S_146_NMBC_str = [Node146_AH_NMBC_str(1,5) Node146_AH_NMBC_str(1,6) Node146_AH_NMBC_str(1,7);Node146_AH_NMBC_str(1,8) Node146_AH_NMBC_str(1,9) Node146_AH_NMBC_str(1,10);...
    Node146_AH_NMBC_str(1,11) Node146_AH_NMBC_str(1,12) Node146_AH_NMBC_str(1,13)];
X146 = sqrt(0.05^2-(Node146_AH_NMBC_str(1,3)+0.25)^2);
teta146 = asin((Node146_AH_NMBC_str(1,3)+0.25)/0.05);
R = [X146/0.05 -(Node146_AH_NMBC_str(1,3)+0.25)/0.05 0;(Node146_AH_NMBC_str(1,3)+0.25)/0.05 X146/0.05 0;...
    0 0 1];
Result_S_146_NMBC_str_trans = R'*Result_S_146_NMBC_str*R;


Result_S_218_NMBC_str = [Node218_AH_NMBC_str(1,5) Node218_AH_NMBC_str(1,6) Node218_AH_NMBC_str(1,7);Node218_AH_NMBC_str(1,8) Node218_AH_NMBC_str(1,9) Node218_AH_NMBC_str(1,10);...
    Node218_AH_NMBC_str(1,11) Node218_AH_NMBC_str(1,12) Node218_AH_NMBC_str(1,13)];
X218 = sqrt(0.05^2-(Node218_AH_NMBC_str(1,3)+0.25)^2);
teta218 = asin((Node218_AH_NMBC_str(1,3)+0.25)/0.05);
R = [X218/0.05 -(Node218_AH_NMBC_str(1,3)+0.25)/0.05 0;(Node218_AH_NMBC_str(1,3)+0.25)/0.05 X218/0.05 0;...
    0 0 1];
Result_S_218_NMBC_str_trans = R'*Result_S_218_NMBC_str*R;


Result_S_290_NMBC_str = [Node290_AH_NMBC_str(1,5) Node290_AH_NMBC_str(1,6) Node290_AH_NMBC_str(1,7);Node290_AH_NMBC_str(1,8) Node290_AH_NMBC_str(1,9) Node290_AH_NMBC_str(1,10);...
    Node290_AH_NMBC_str(1,11) Node290_AH_NMBC_str(1,12) Node290_AH_NMBC_str(1,13)];
X290 = sqrt(0.05^2-(Node290_AH_NMBC_str(1,3)+0.25)^2);
teta290 = asin((Node290_AH_NMBC_str(1,3)+0.25)/0.05);
R = [X290/0.05 -(Node290_AH_NMBC_str(1,3)+0.25)/0.05 0;(Node290_AH_NMBC_str(1,3)+0.25)/0.05 X290/0.05 0;...
    0 0 1];
Result_S_290_NMBC_str_trans = R'*Result_S_290_NMBC_str*R;


Result_S_362_NMBC_str = [Node362_AH_NMBC_str(1,5) Node362_AH_NMBC_str(1,6) Node362_AH_NMBC_str(1,7);Node362_AH_NMBC_str(1,8) Node362_AH_NMBC_str(1,9) Node362_AH_NMBC_str(1,10);...
    Node362_AH_NMBC_str(1,11) Node362_AH_NMBC_str(1,12) Node362_AH_NMBC_str(1,13)];
X362 = sqrt(0.05^2-(Node362_AH_NMBC_str(1,3)+0.25)^2);
teta362 = asin((Node362_AH_NMBC_str(1,3)+0.25)/0.05);
R = [X362/0.05 -(Node362_AH_NMBC_str(1,3)+0.25)/0.05 0;(Node362_AH_NMBC_str(1,3)+0.25)/0.05 X362/0.05 0;...
    0 0 1];
Result_S_362_NMBC_str_trans = R'*Result_S_362_NMBC_str*R;


Result_S_434_NMBC_str = [Node434_AH_NMBC_str(1,5) Node434_AH_NMBC_str(1,6) Node434_AH_NMBC_str(1,7);Node434_AH_NMBC_str(1,8) Node434_AH_NMBC_str(1,9) Node434_AH_NMBC_str(1,10);...
    Node434_AH_NMBC_str(1,11) Node434_AH_NMBC_str(1,12) Node434_AH_NMBC_str(1,13)];
X434 = sqrt(0.05^2-(Node434_AH_NMBC_str(1,3)+0.25)^2);
teta434 = asin((Node434_AH_NMBC_str(1,3)+0.25)/0.05);
R = [X434/0.05 -(Node434_AH_NMBC_str(1,3)+0.25)/0.05 0;(Node434_AH_NMBC_str(1,3)+0.25)/0.05 X434/0.05 0;...
    0 0 1];
Result_S_434_NMBC_str_trans = R'*Result_S_434_NMBC_str*R;


Result_S_506_NMBC_str = [Node506_AH_NMBC_str(1,5) Node506_AH_NMBC_str(1,6) Node506_AH_NMBC_str(1,7);Node506_AH_NMBC_str(1,8) Node506_AH_NMBC_str(1,9) Node506_AH_NMBC_str(1,10);...
    Node506_AH_NMBC_str(1,11) Node506_AH_NMBC_str(1,12) Node506_AH_NMBC_str(1,13)];
X506 = sqrt(0.05^2-(Node506_AH_NMBC_str(1,3)+0.25)^2);
teta506 = asin((Node506_AH_NMBC_str(1,3)+0.25)/0.05);
R = [X506/0.05 -(Node506_AH_NMBC_str(1,3)+0.25)/0.05 0;(Node506_AH_NMBC_str(1,3)+0.25)/0.05 X506/0.05 0;...
    0 0 1];
Result_S_506_NMBC_str_trans = R'*Result_S_506_NMBC_str*R;


Result_S_578_NMBC_str = [Node578_AH_NMBC_str(1,5) Node578_AH_NMBC_str(1,6) Node578_AH_NMBC_str(1,7);Node578_AH_NMBC_str(1,8) Node578_AH_NMBC_str(1,9) Node578_AH_NMBC_str(1,10);...
    Node578_AH_NMBC_str(1,11) Node578_AH_NMBC_str(1,12) Node578_AH_NMBC_str(1,13)];
X578 = sqrt(0.05^2-(Node578_AH_NMBC_str(1,3)+0.25)^2);
teta578 = asin((Node578_AH_NMBC_str(1,3)+0.25)/0.05);
R = [X578/0.05 -(Node578_AH_NMBC_str(1,3)+0.25)/0.05 0;(Node578_AH_NMBC_str(1,3)+0.25)/0.05 X578/0.05 0;...
    0 0 1];
Result_S_578_NMBC_str_trans = R'*Result_S_578_NMBC_str*R;


Result_S_650_NMBC_str = [Node650_AH_NMBC_str(1,5) Node650_AH_NMBC_str(1,6) Node650_AH_NMBC_str(1,7);Node650_AH_NMBC_str(1,8) Node650_AH_NMBC_str(1,9) Node650_AH_NMBC_str(1,10);...
    Node650_AH_NMBC_str(1,11) Node650_AH_NMBC_str(1,12) Node650_AH_NMBC_str(1,13)];
X650 = sqrt(0.05^2-(Node650_AH_NMBC_str(1,3)+0.25)^2);
teta650 = asin((Node650_AH_NMBC_str(1,3)+0.25)/0.05);
R = [X650/0.05 -(Node650_AH_NMBC_str(1,3)+0.25)/0.05 0;(Node650_AH_NMBC_str(1,3)+0.25)/0.05 X650/0.05 0;...
    0 0 1];
Result_S_650_NMBC_str_trans = R'*Result_S_650_NMBC_str*R;


Result_S_722_NMBC_str = [Node722_AH_NMBC_str(1,5) Node722_AH_NMBC_str(1,6) Node722_AH_NMBC_str(1,7);Node722_AH_NMBC_str(1,8) Node722_AH_NMBC_str(1,9) Node722_AH_NMBC_str(1,10);...
    Node722_AH_NMBC_str(1,11) Node722_AH_NMBC_str(1,12) Node722_AH_NMBC_str(1,13)];
X722 = sqrt(0.05^2-(Node722_AH_NMBC_str(1,3)+0.25)^2);
teta722 = asin((Node722_AH_NMBC_str(1,3)+0.25)/0.05);
R = [X722/0.05 -(Node722_AH_NMBC_str(1,3)+0.25)/0.05 0;(Node722_AH_NMBC_str(1,3)+0.25)/0.05 X722/0.05 0;...
    0 0 1];
Result_S_722_NMBC_str_trans = R'*Result_S_722_NMBC_str*R;


Result_S_794_NMBC_str = [Node794_AH_NMBC_str(1,5) Node794_AH_NMBC_str(1,6) Node794_AH_NMBC_str(1,7);Node794_AH_NMBC_str(1,8) Node794_AH_NMBC_str(1,9) Node794_AH_NMBC_str(1,10);...
    Node794_AH_NMBC_str(1,11) Node794_AH_NMBC_str(1,12) Node794_AH_NMBC_str(1,13)];
X794 = sqrt(0.05^2-(Node794_AH_NMBC_str(1,3)+0.25)^2);
teta794 = asin((Node794_AH_NMBC_str(1,3)+0.25)/0.05);
R = [X794/0.05 -(Node794_AH_NMBC_str(1,3)+0.25)/0.05 0;(Node794_AH_NMBC_str(1,3)+0.25)/0.05 X794/0.05 0;...
    0 0 1];
Result_S_794_NMBC_str_trans = R'*Result_S_794_NMBC_str*R;


Result_S_866_NMBC_str = [Node866_AH_NMBC_str(1,5) Node866_AH_NMBC_str(1,6) Node866_AH_NMBC_str(1,7);Node866_AH_NMBC_str(1,8) Node866_AH_NMBC_str(1,9) Node866_AH_NMBC_str(1,10);...
    Node866_AH_NMBC_str(1,11) Node866_AH_NMBC_str(1,12) Node866_AH_NMBC_str(1,13)];
X866 = sqrt(0.05^2-(Node866_AH_NMBC_str(1,3)+0.25)^2);
teta866 = asin((Node866_AH_NMBC_str(1,3)+0.25)/0.05);
R = [X866/0.05 -(Node866_AH_NMBC_str(1,3)+0.25)/0.05 0;(Node866_AH_NMBC_str(1,3)+0.25)/0.05 X866/0.05 0;...
    0 0 1];
Result_S_866_NMBC_str_trans = R'*Result_S_866_NMBC_str*R;


Result_S_938_NMBC_str = [Node938_AH_NMBC_str(1,5) Node938_AH_NMBC_str(1,6) Node938_AH_NMBC_str(1,7);Node938_AH_NMBC_str(1,8) Node938_AH_NMBC_str(1,9) Node938_AH_NMBC_str(1,10);...
    Node938_AH_NMBC_str(1,11) Node938_AH_NMBC_str(1,12) Node938_AH_NMBC_str(1,13)];
X938 = sqrt(0.05^2-(Node938_AH_NMBC_str(1,3)+0.25)^2);
teta938 = asin((Node938_AH_NMBC_str(1,3)+0.25)/0.05);
R = [X938/0.05 -(Node938_AH_NMBC_str(1,3)+0.25)/0.05 0;(Node938_AH_NMBC_str(1,3)+0.25)/0.05 X938/0.05 0;...
    0 0 1];
Result_S_938_NMBC_str_trans = R'*Result_S_938_NMBC_str*R;


Result_S_1876_NMBC_str = [Node1876_AH_NMBC_str(1,5) Node1876_AH_NMBC_str(1,6) Node1876_AH_NMBC_str(1,7);Node1876_AH_NMBC_str(1,8) Node1876_AH_NMBC_str(1,9) Node1876_AH_NMBC_str(1,10);...
    Node1876_AH_NMBC_str(1,11) Node1876_AH_NMBC_str(1,12) Node1876_AH_NMBC_str(1,13)];
X1876 = sqrt(0.05^2-(Node1876_AH_NMBC_str(1,3)+0.25)^2);
teta1876 = asin((Node1876_AH_NMBC_str(1,3)+0.25)/0.05);
R = [X1876/0.05 -(Node1876_AH_NMBC_str(1,3)+0.25)/0.05 0;(Node1876_AH_NMBC_str(1,3)+0.25)/0.05 X1876/0.05 0;...
    0 0 1];
Result_S_1876_NMBC_str_trans = R'*Result_S_1876_NMBC_str*R;


Result_S_1804_NMBC_str = [Node1804_AH_NMBC_str(1,5) Node1804_AH_NMBC_str(1,6) Node1804_AH_NMBC_str(1,7);Node1804_AH_NMBC_str(1,8) Node1804_AH_NMBC_str(1,9) Node1804_AH_NMBC_str(1,10);...
    Node1804_AH_NMBC_str(1,11) Node1804_AH_NMBC_str(1,12) Node1804_AH_NMBC_str(1,13)];
X1804 = sqrt(0.05^2-(Node1804_AH_NMBC_str(1,3)+0.25)^2);
teta1804 = asin((Node1804_AH_NMBC_str(1,3)+0.25)/0.05);
R = [X1804/0.05 -(Node1804_AH_NMBC_str(1,3)+0.25)/0.05 0;(Node1804_AH_NMBC_str(1,3)+0.25)/0.05 X1804/0.05 0;...
    0 0 1];
Result_S_1804_NMBC_str_trans = R'*Result_S_1804_NMBC_str*R;


Result_S_1732_NMBC_str = [Node1732_AH_NMBC_str(1,5) Node1732_AH_NMBC_str(1,6) Node1732_AH_NMBC_str(1,7);Node1732_AH_NMBC_str(1,8) Node1732_AH_NMBC_str(1,9) Node1732_AH_NMBC_str(1,10);...
    Node1732_AH_NMBC_str(1,11) Node1732_AH_NMBC_str(1,12) Node1732_AH_NMBC_str(1,13)];
X1732 = sqrt(0.05^2-(Node1732_AH_NMBC_str(1,3)+0.25)^2);
teta1732 = asin((Node1732_AH_NMBC_str(1,3)+0.25)/0.05);
R = [X1732/0.05 -(Node1732_AH_NMBC_str(1,3)+0.25)/0.05 0;(Node1732_AH_NMBC_str(1,3)+0.25)/0.05 X1732/0.05 0;...
    0 0 1];
Result_S_1732_NMBC_str_trans = R'*Result_S_1732_NMBC_str*R;


Result_S_1660_NMBC_str = [Node1660_AH_NMBC_str(1,5) Node1660_AH_NMBC_str(1,6) Node1660_AH_NMBC_str(1,7);Node1660_AH_NMBC_str(1,8) Node1660_AH_NMBC_str(1,9) Node1660_AH_NMBC_str(1,10);...
    Node1660_AH_NMBC_str(1,11) Node1660_AH_NMBC_str(1,12) Node1660_AH_NMBC_str(1,13)];
X1660 = sqrt(0.05^2-(Node1660_AH_NMBC_str(1,3)+0.25)^2);
teta1660 = asin((Node1660_AH_NMBC_str(1,3)+0.25)/0.05);
R = [X1660/0.05 -(Node1660_AH_NMBC_str(1,3)+0.25)/0.05 0;(Node1660_AH_NMBC_str(1,3)+0.25)/0.05 X1660/0.05 0;...
    0 0 1];
Result_S_1660_NMBC_str_trans = R'*Result_S_1660_NMBC_str*R;


Result_S_1588_NMBC_str = [Node1588_AH_NMBC_str(1,5) Node1588_AH_NMBC_str(1,6) Node1588_AH_NMBC_str(1,7);Node1588_AH_NMBC_str(1,8) Node1588_AH_NMBC_str(1,9) Node1588_AH_NMBC_str(1,10);...
    Node1588_AH_NMBC_str(1,11) Node1588_AH_NMBC_str(1,12) Node1588_AH_NMBC_str(1,13)];
X1588 = sqrt(0.05^2-(Node1588_AH_NMBC_str(1,3)+0.25)^2);
teta1588 = asin((Node1588_AH_NMBC_str(1,3)+0.25)/0.05);
R = [X1588/0.05 -(Node1588_AH_NMBC_str(1,3)+0.25)/0.05 0;(Node1588_AH_NMBC_str(1,3)+0.25)/0.05 X1588/0.05 0;...
    0 0 1];
Result_S_1588_NMBC_str_trans = R'*Result_S_1588_NMBC_str*R;


Result_S_1516_NMBC_str = [Node1516_AH_NMBC_str(1,5) Node1516_AH_NMBC_str(1,6) Node1516_AH_NMBC_str(1,7);Node1516_AH_NMBC_str(1,8) Node1516_AH_NMBC_str(1,9) Node1516_AH_NMBC_str(1,10);...
    Node1516_AH_NMBC_str(1,11) Node1516_AH_NMBC_str(1,12) Node1516_AH_NMBC_str(1,13)];
X1516 = sqrt(0.05^2-(Node1516_AH_NMBC_str(1,3)+0.25)^2);
teta1516 = asin((Node1516_AH_NMBC_str(1,3)+0.25)/0.05);
R = [X1516/0.05 -(Node1516_AH_NMBC_str(1,3)+0.25)/0.05 0;(Node1516_AH_NMBC_str(1,3)+0.25)/0.05 X1516/0.05 0;...
    0 0 1];
Result_S_1516_NMBC_str_trans = R'*Result_S_1516_NMBC_str*R;


Result_S_1444_NMBC_str = [Node1444_AH_NMBC_str(1,5) Node1444_AH_NMBC_str(1,6) Node1444_AH_NMBC_str(1,7);Node1444_AH_NMBC_str(1,8) Node1444_AH_NMBC_str(1,9) Node1444_AH_NMBC_str(1,10);...
    Node1444_AH_NMBC_str(1,11) Node1444_AH_NMBC_str(1,12) Node1444_AH_NMBC_str(1,13)];
X1444 = sqrt(0.05^2-(Node1444_AH_NMBC_str(1,3)+0.25)^2);
teta1444 = asin((Node1444_AH_NMBC_str(1,3)+0.25)/0.05);
R = [X1444/0.05 -(Node1444_AH_NMBC_str(1,3)+0.25)/0.05 0;(Node1444_AH_NMBC_str(1,3)+0.25)/0.05 X1444/0.05 0;...
    0 0 1];
Result_S_1444_NMBC_str_trans = R'*Result_S_1444_NMBC_str*R;


Result_S_1372_NMBC_str = [Node1372_AH_NMBC_str(1,5) Node1372_AH_NMBC_str(1,6) Node1372_AH_NMBC_str(1,7);Node1372_AH_NMBC_str(1,8) Node1372_AH_NMBC_str(1,9) Node1372_AH_NMBC_str(1,10);...
    Node1372_AH_NMBC_str(1,11) Node1372_AH_NMBC_str(1,12) Node1372_AH_NMBC_str(1,13)];
X1372 = sqrt(0.05^2-(Node1372_AH_NMBC_str(1,3)+0.25)^2);
teta1372 = asin((Node1372_AH_NMBC_str(1,3)+0.25)/0.05);
R = [X1372/0.05 -(Node1372_AH_NMBC_str(1,3)+0.25)/0.05 0;(Node1372_AH_NMBC_str(1,3)+0.25)/0.05 X1372/0.05 0;...
    0 0 1];
Result_S_1372_NMBC_str_trans = R'*Result_S_1372_NMBC_str*R;


Result_S_1300_NMBC_str = [Node1300_AH_NMBC_str(1,5) Node1300_AH_NMBC_str(1,6) Node1300_AH_NMBC_str(1,7);Node1300_AH_NMBC_str(1,8) Node1300_AH_NMBC_str(1,9) Node1300_AH_NMBC_str(1,10);...
    Node1300_AH_NMBC_str(1,11) Node1300_AH_NMBC_str(1,12) Node1300_AH_NMBC_str(1,13)];
X1300 = sqrt(0.05^2-(Node1300_AH_NMBC_str(1,3)+0.25)^2);
teta1300 = asin((Node1300_AH_NMBC_str(1,3)+0.25)/0.05);
R = [X1300/0.05 -(Node1300_AH_NMBC_str(1,3)+0.25)/0.05 0;(Node1300_AH_NMBC_str(1,3)+0.25)/0.05 X1300/0.05 0;...
    0 0 1];
Result_S_1300_NMBC_str_trans = R'*Result_S_1300_NMBC_str*R;


Result_S_1228_NMBC_str = [Node1228_AH_NMBC_str(1,5) Node1228_AH_NMBC_str(1,6) Node1228_AH_NMBC_str(1,7);Node1228_AH_NMBC_str(1,8) Node1228_AH_NMBC_str(1,9) Node1228_AH_NMBC_str(1,10);...
    Node1228_AH_NMBC_str(1,11) Node1228_AH_NMBC_str(1,12) Node1228_AH_NMBC_str(1,13)];
X1228 = sqrt(0.05^2-(Node1228_AH_NMBC_str(1,3)+0.25)^2);
teta1228 = asin((Node1228_AH_NMBC_str(1,3)+0.25)/0.05);
R = [X1228/0.05 -(Node1228_AH_NMBC_str(1,3)+0.25)/0.05 0;(Node1228_AH_NMBC_str(1,3)+0.25)/0.05 X1228/0.05 0;...
    0 0 1];
Result_S_1228_NMBC_str_trans = R'*Result_S_1228_NMBC_str*R;


Result_S_1156_NMBC_str = [Node1156_AH_NMBC_str(1,5) Node1156_AH_NMBC_str(1,6) Node1156_AH_NMBC_str(1,7);Node1156_AH_NMBC_str(1,8) Node1156_AH_NMBC_str(1,9) Node1156_AH_NMBC_str(1,10);...
    Node1156_AH_NMBC_str(1,11) Node1156_AH_NMBC_str(1,12) Node1156_AH_NMBC_str(1,13)];
X1156 = sqrt(0.05^2-(Node1156_AH_NMBC_str(1,3)+0.25)^2);
teta1156 = asin((Node1156_AH_NMBC_str(1,3)+0.25)/0.05);
R = [X1156/0.05 -(Node1156_AH_NMBC_str(1,3)+0.25)/0.05 0;(Node1156_AH_NMBC_str(1,3)+0.25)/0.05 X1156/0.05 0;...
    0 0 1];
Result_S_1156_NMBC_str_trans = R'*Result_S_1156_NMBC_str*R;


Result_S_1016_NMBC_str = [Node1016_AH_NMBC_str(1,5) Node1016_AH_NMBC_str(1,6) Node1016_AH_NMBC_str(1,7);Node1016_AH_NMBC_str(1,8) Node1016_AH_NMBC_str(1,9) Node1016_AH_NMBC_str(1,10);...
    Node1016_AH_NMBC_str(1,11) Node1016_AH_NMBC_str(1,12) Node1016_AH_NMBC_str(1,13)];
X1016 = sqrt(0.05^2-(Node1016_AH_NMBC_str(1,3)+0.25)^2);
teta1016 = asin((Node1016_AH_NMBC_str(1,3)+0.25)/0.05);
R = [X1016/0.05 -(Node1016_AH_NMBC_str(1,3)+0.25)/0.05 0;(Node1016_AH_NMBC_str(1,3)+0.25)/0.05 X1016/0.05 0;...
    0 0 1];
Result_S_1016_NMBC_str_trans = R'*Result_S_1016_NMBC_str*R;


Result_S_1012_NMBC_str = [Node1012_AH_NMBC_str(1,5) Node1012_AH_NMBC_str(1,6) Node1012_AH_NMBC_str(1,7);Node1012_AH_NMBC_str(1,8) Node1012_AH_NMBC_str(1,9) Node1012_AH_NMBC_str(1,10);...
    Node1012_AH_NMBC_str(1,11) Node1012_AH_NMBC_str(1,12) Node1012_AH_NMBC_str(1,13)];
X1012 = sqrt(0.05^2-(Node1012_AH_NMBC_str(1,3)+0.25)^2);
teta1012 = asin((Node1012_AH_NMBC_str(1,3)+0.25)/0.05);
R = [X1012/0.05 -(Node1012_AH_NMBC_str(1,3)+0.25)/0.05 0;(Node1012_AH_NMBC_str(1,3)+0.25)/0.05 X1012/0.05 0;...
    0 0 1];
Result_S_1012_NMBC_str_trans = R'*Result_S_1012_NMBC_str*R;


Result_S_teta_teta_NMBC_str = [teta2 Result_S_2_NMBC_str_trans(2,2)/P_applied;teta6 Result_S_6_NMBC_str_trans(2,2)/P_applied;teta146 Result_S_146_NMBC_str_trans(2,2)/P_applied;...
    teta218 Result_S_218_NMBC_str_trans(2,2)/P_applied;teta290 Result_S_290_NMBC_str_trans(2,2)/P_applied;teta362 Result_S_362_NMBC_str_trans(2,2)/P_applied;...
    teta434 Result_S_434_NMBC_str_trans(2,2)/P_applied;teta506 Result_S_506_NMBC_str_trans(2,2)/P_applied;teta578 Result_S_578_NMBC_str_trans(2,2)/P_applied;...
    teta650 Result_S_650_NMBC_str_trans(2,2)/P_applied;teta722 Result_S_722_NMBC_str_trans(2,2)/P_applied;teta794 Result_S_794_NMBC_str_trans(2,2)/P_applied;...
    teta866 Result_S_866_NMBC_str_trans(2,2)/P_applied;teta938 Result_S_938_NMBC_str_trans(2,2)/P_applied;teta1876 Result_S_1876_NMBC_str_trans(2,2)/P_applied;...
    teta1804 Result_S_1804_NMBC_str_trans(2,2)/P_applied;teta1732 Result_S_1732_NMBC_str_trans(2,2)/P_applied;teta1660 Result_S_1660_NMBC_str_trans(2,2)/P_applied;...
    teta1588 Result_S_1588_NMBC_str_trans(2,2)/P_applied;teta1516 Result_S_1516_NMBC_str_trans(2,2)/P_applied;teta1444 Result_S_1444_NMBC_str_trans(2,2)/P_applied;...
    teta1372 Result_S_1372_NMBC_str_trans(2,2)/P_applied;teta1300 Result_S_1300_NMBC_str_trans(2,2)/P_applied;teta1228 Result_S_1228_NMBC_str_trans(2,2)/P_applied;...
    teta1156 Result_S_1156_NMBC_str_trans(2,2)/P_applied;teta1016 Result_S_1016_NMBC_str_trans(2,2)/P_applied;teta1012 Result_S_1012_NMBC_str_trans(2,2)/P_applied];


Result_S_938_AD_cls = [Node938_cls(1,7) Node938_cls(1,8) Node938_cls(1,9);Node938_cls(1,10) Node938_cls(1,11) Node938_cls(1,12);...
    Node938_cls(1,13) Node938_cls(1,14) Node938_cls(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_938_AD_cls_trans = R'*Result_S_938_AD_cls*R;


Result_S_939_AD_cls = [Node939_cls(1,7) Node939_cls(1,8) Node939_cls(1,9);Node939_cls(1,10) Node939_cls(1,11) Node939_cls(1,12);...
    Node939_cls(1,13) Node939_cls(1,14) Node939_cls(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_939_AD_cls_trans = R'*Result_S_939_AD_cls*R;


Result_S_941_AD_cls = [Node941_cls(1,7) Node941_cls(1,8) Node941_cls(1,9);Node941_cls(1,10) Node941_cls(1,11) Node941_cls(1,12);...
    Node941_cls(1,13) Node941_cls(1,14) Node941_cls(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_941_AD_cls_trans = R'*Result_S_941_AD_cls*R;


Result_S_943_AD_cls = [Node943_cls(1,7) Node943_cls(1,8) Node943_cls(1,9);Node943_cls(1,10) Node943_cls(1,11) Node943_cls(1,12);...
    Node943_cls(1,13) Node943_cls(1,14) Node943_cls(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_943_AD_cls_trans = R'*Result_S_943_AD_cls*R;


Result_S_945_AD_cls = [Node945_cls(1,7) Node945_cls(1,8) Node945_cls(1,9);Node945_cls(1,10) Node945_cls(1,11) Node945_cls(1,12);...
    Node945_cls(1,13) Node945_cls(1,14) Node945_cls(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_945_AD_cls_trans = R'*Result_S_945_AD_cls*R;


Result_S_947_AD_cls = [Node947_cls(1,7) Node947_cls(1,8) Node947_cls(1,9);Node947_cls(1,10) Node947_cls(1,11) Node947_cls(1,12);...
    Node947_cls(1,13) Node947_cls(1,14) Node947_cls(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_947_AD_cls_trans = R'*Result_S_947_AD_cls*R;


Result_S_949_AD_cls = [Node949_cls(1,7) Node949_cls(1,8) Node949_cls(1,9);Node949_cls(1,10) Node949_cls(1,11) Node949_cls(1,12);...
    Node949_cls(1,13) Node949_cls(1,14) Node949_cls(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_949_AD_cls_trans = R'*Result_S_949_AD_cls*R;


Result_S_951_AD_cls = [Node951_cls(1,7) Node951_cls(1,8) Node951_cls(1,9);Node951_cls(1,10) Node951_cls(1,11) Node951_cls(1,12);...
    Node951_cls(1,13) Node951_cls(1,14) Node951_cls(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_951_AD_cls_trans = R'*Result_S_951_AD_cls*R;


Result_S_953_AD_cls = [Node953_cls(1,7) Node953_cls(1,8) Node953_cls(1,9);Node953_cls(1,10) Node953_cls(1,11) Node953_cls(1,12);...
    Node953_cls(1,13) Node953_cls(1,14) Node953_cls(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_953_AD_cls_trans = R'*Result_S_953_AD_cls*R;


Result_S_955_AD_cls = [Node955_cls(1,7) Node955_cls(1,8) Node955_cls(1,9);Node955_cls(1,10) Node955_cls(1,11) Node955_cls(1,12);...
    Node955_cls(1,13) Node955_cls(1,14) Node955_cls(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_955_AD_cls_trans = R'*Result_S_955_AD_cls*R;


Result_S_957_AD_cls = [Node957_cls(1,7) Node957_cls(1,8) Node957_cls(1,9);Node957_cls(1,10) Node957_cls(1,11) Node957_cls(1,12);...
    Node957_cls(1,13) Node957_cls(1,14) Node957_cls(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_957_AD_cls_trans = R'*Result_S_957_AD_cls*R;


Result_S_959_AD_cls = [Node959_cls(1,7) Node959_cls(1,8) Node959_cls(1,9);Node959_cls(1,10) Node959_cls(1,11) Node959_cls(1,12);...
    Node959_cls(1,13) Node959_cls(1,14) Node959_cls(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_959_AD_cls_trans = R'*Result_S_959_AD_cls*R;


Result_S_961_AD_cls = [Node961_cls(1,7) Node961_cls(1,8) Node961_cls(1,9);Node961_cls(1,10) Node961_cls(1,11) Node961_cls(1,12);...
    Node961_cls(1,13) Node961_cls(1,14) Node961_cls(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_961_AD_cls_trans = R'*Result_S_961_AD_cls*R;

Result_S_963_AD_cls = [Node963_cls(1,7) Node963_cls(1,8) Node963_cls(1,9);Node963_cls(1,10) Node963_cls(1,11) Node963_cls(1,12);...
    Node963_cls(1,13) Node963_cls(1,14) Node963_cls(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_963_AD_cls_trans = R'*Result_S_963_AD_cls*R;


Result_S_965_AD_cls = [Node965_cls(1,7) Node965_cls(1,8) Node965_cls(1,9);Node965_cls(1,10) Node965_cls(1,11) Node965_cls(1,12);...
    Node965_cls(1,13) Node965_cls(1,14) Node965_cls(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_965_AD_cls_trans = R'*Result_S_965_AD_cls*R;

Result_S_967_AD_cls = [Node967_cls(1,7) Node967_cls(1,8) Node967_cls(1,9);Node967_cls(1,10) Node967_cls(1,11) Node967_cls(1,12);...
    Node967_cls(1,13) Node967_cls(1,14) Node967_cls(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_967_AD_cls_trans = R'*Result_S_967_AD_cls*R;

Result_S_969_AD_cls = [Node969_cls(1,7) Node969_cls(1,8) Node969_cls(1,9);Node969_cls(1,10) Node969_cls(1,11) Node969_cls(1,12);...
    Node969_cls(1,13) Node969_cls(1,14) Node969_cls(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_969_AD_cls_trans = R'*Result_S_969_AD_cls*R;

Result_S_971_AD_cls = [Node971_cls(1,7) Node971_cls(1,8) Node971_cls(1,9);Node971_cls(1,10) Node971_cls(1,11) Node971_cls(1,12);...
    Node971_cls(1,13) Node971_cls(1,14) Node971_cls(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_971_AD_cls_trans = R'*Result_S_971_AD_cls*R;

Result_S_973_AD_cls = [Node973_cls(1,7) Node973_cls(1,8) Node973_cls(1,9);Node973_cls(1,10) Node973_cls(1,11) Node973_cls(1,12);...
    Node973_cls(1,13) Node973_cls(1,14) Node973_cls(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_973_AD_cls_trans = R'*Result_S_973_AD_cls*R;

Result_S_975_AD_cls = [Node975_cls(1,7) Node975_cls(1,8) Node975_cls(1,9);Node975_cls(1,10) Node975_cls(1,11) Node975_cls(1,12);...
    Node975_cls(1,13) Node975_cls(1,14) Node975_cls(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_975_AD_cls_trans = R'*Result_S_975_AD_cls*R;

Result_S_977_AD_cls = [Node977_cls(1,7) Node977_cls(1,8) Node977_cls(1,9);Node977_cls(1,10) Node977_cls(1,11) Node977_cls(1,12);...
    Node977_cls(1,13) Node977_cls(1,14) Node977_cls(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_977_AD_cls_trans = R'*Result_S_977_AD_cls*R;

Result_S_979_AD_cls = [Node979_cls(1,7) Node979_cls(1,8) Node979_cls(1,9);Node979_cls(1,10) Node979_cls(1,11) Node979_cls(1,12);...
    Node979_cls(1,13) Node979_cls(1,14) Node979_cls(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_979_AD_cls_trans = R'*Result_S_979_AD_cls*R;

Result_S_981_AD_cls = [Node981_cls(1,7) Node981_cls(1,8) Node981_cls(1,9);Node981_cls(1,10) Node981_cls(1,11) Node981_cls(1,12);...
    Node981_cls(1,13) Node981_cls(1,14) Node981_cls(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_981_AD_cls_trans = R'*Result_S_981_AD_cls*R;

Result_S_983_AD_cls = [Node983_cls(1,7) Node983_cls(1,8) Node983_cls(1,9);Node983_cls(1,10) Node983_cls(1,11) Node983_cls(1,12);...
    Node983_cls(1,13) Node983_cls(1,14) Node983_cls(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_983_AD_cls_trans = R'*Result_S_983_AD_cls*R;

Result_S_985_AD_cls = [Node985_cls(1,7) Node985_cls(1,8) Node985_cls(1,9);Node985_cls(1,10) Node985_cls(1,11) Node985_cls(1,12);...
    Node985_cls(1,13) Node985_cls(1,14) Node985_cls(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_985_AD_cls_trans = R'*Result_S_985_AD_cls*R;

Result_S_987_AD_cls = [Node987_cls(1,7) Node987_cls(1,8) Node987_cls(1,9);Node987_cls(1,10) Node987_cls(1,11) Node987_cls(1,12);...
    Node987_cls(1,13) Node987_cls(1,14) Node987_cls(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_987_AD_cls_trans = R'*Result_S_987_AD_cls*R;

Result_S_989_AD_cls = [Node989_cls(1,7) Node989_cls(1,8) Node989_cls(1,9);Node989_cls(1,10) Node989_cls(1,11) Node989_cls(1,12);...
    Node989_cls(1,13) Node989_cls(1,14) Node989_cls(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_989_AD_cls_trans = R'*Result_S_989_AD_cls*R;

Result_S_991_AD_cls = [Node991_cls(1,7) Node991_cls(1,8) Node991_cls(1,9);Node991_cls(1,10) Node991_cls(1,11) Node991_cls(1,12);...
    Node991_cls(1,13) Node991_cls(1,14) Node991_cls(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_991_AD_cls_trans = R'*Result_S_991_AD_cls*R;

Result_S_993_AD_cls = [Node993_cls(1,7) Node993_cls(1,8) Node993_cls(1,9);Node993_cls(1,10) Node993_cls(1,11) Node993_cls(1,12);...
    Node993_cls(1,13) Node993_cls(1,14) Node993_cls(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_993_AD_cls_trans = R'*Result_S_993_AD_cls*R;

Result_S_995_AD_cls = [Node995_cls(1,7) Node995_cls(1,8) Node995_cls(1,9);Node995_cls(1,10) Node995_cls(1,11) Node995_cls(1,12);...
    Node995_cls(1,13) Node995_cls(1,14) Node995_cls(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_995_AD_cls_trans = R'*Result_S_995_AD_cls*R;

Result_S_997_AD_cls = [Node997_cls(1,7) Node997_cls(1,8) Node997_cls(1,9);Node997_cls(1,10) Node997_cls(1,11) Node997_cls(1,12);...
    Node997_cls(1,13) Node997_cls(1,14) Node997_cls(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_997_AD_cls_trans = R'*Result_S_997_AD_cls*R;

Result_S_999_AD_cls = [Node999_cls(1,7) Node999_cls(1,8) Node999_cls(1,9);Node999_cls(1,10) Node999_cls(1,11) Node999_cls(1,12);...
    Node999_cls(1,13) Node999_cls(1,14) Node999_cls(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_999_AD_cls_trans = R'*Result_S_999_AD_cls*R;

Result_S_1001_AD_cls = [Node1001_cls(1,7) Node1001_cls(1,8) Node1001_cls(1,9);Node1001_cls(1,10) Node1001_cls(1,11) Node1001_cls(1,12);...
    Node1001_cls(1,13) Node1001_cls(1,14) Node1001_cls(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_1001_AD_cls_trans = R'*Result_S_1001_AD_cls*R;

Result_S_1003_AD_cls = [Node1003_cls(1,7) Node1003_cls(1,8) Node1003_cls(1,9);Node1003_cls(1,10) Node1003_cls(1,11) Node1003_cls(1,12);...
    Node1003_cls(1,13) Node1003_cls(1,14) Node1003_cls(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_1003_AD_cls_trans = R'*Result_S_1003_AD_cls*R;

Result_S_1005_AD_cls = [Node1005_cls(1,7) Node1005_cls(1,8) Node1005_cls(1,9);Node1005_cls(1,10) Node1005_cls(1,11) Node1005_cls(1,12);...
    Node1005_cls(1,13) Node1005_cls(1,14) Node1005_cls(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_1005_AD_cls_trans = R'*Result_S_1005_AD_cls*R;

Result_S_1007_AD_cls = [Node1007_cls(1,7) Node1007_cls(1,8) Node1007_cls(1,9);Node1007_cls(1,10) Node1007_cls(1,11) Node1007_cls(1,12);...
    Node1007_cls(1,13) Node1007_cls(1,14) Node1007_cls(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_1007_AD_cls_trans = R'*Result_S_1007_AD_cls*R;

r_938 = sqrt((Node938_cls(1,2)+0.25)^2+(Node938_cls(1,3)+0.25)^2);
r_939 = sqrt((Node939_cls(1,2)+0.25)^2+(Node939_cls(1,3)+0.25)^2);
r_941 = sqrt((Node941_cls(1,2)+0.25)^2+(Node941_cls(1,3)+0.25)^2);
r_943 = sqrt((Node943_cls(1,2)+0.25)^2+(Node943_cls(1,3)+0.25)^2);
r_945 = sqrt((Node945_cls(1,2)+0.25)^2+(Node945_cls(1,3)+0.25)^2);
r_947 = sqrt((Node947_cls(1,2)+0.25)^2+(Node947_cls(1,3)+0.25)^2);
r_949 = sqrt((Node949_cls(1,2)+0.25)^2+(Node949_cls(1,3)+0.25)^2);
r_951 = sqrt((Node951_cls(1,2)+0.25)^2+(Node951_cls(1,3)+0.25)^2);
r_953 = sqrt((Node953_cls(1,2)+0.25)^2+(Node953_cls(1,3)+0.25)^2);
r_955 = sqrt((Node955_cls(1,2)+0.25)^2+(Node955_cls(1,3)+0.25)^2);
r_957 = sqrt((Node957_cls(1,2)+0.25)^2+(Node957_cls(1,3)+0.25)^2);
r_959 = sqrt((Node959_cls(1,2)+0.25)^2+(Node959_cls(1,3)+0.25)^2);
r_961 = sqrt((Node961_cls(1,2)+0.25)^2+(Node961_cls(1,3)+0.25)^2);
r_963 = sqrt((Node963_cls(1,2)+0.25)^2+(Node963_cls(1,3)+0.25)^2);
r_965 = sqrt((Node965_cls(1,2)+0.25)^2+(Node965_cls(1,3)+0.25)^2);
r_967 = sqrt((Node967_cls(1,2)+0.25)^2+(Node967_cls(1,3)+0.25)^2);
r_969 = sqrt((Node969_cls(1,2)+0.25)^2+(Node969_cls(1,3)+0.25)^2);
r_971 = sqrt((Node971_cls(1,2)+0.25)^2+(Node971_cls(1,3)+0.25)^2);
r_973 = sqrt((Node973_cls(1,2)+0.25)^2+(Node973_cls(1,3)+0.25)^2);
r_975 = sqrt((Node975_cls(1,2)+0.25)^2+(Node975_cls(1,3)+0.25)^2);
r_977 = sqrt((Node977_cls(1,2)+0.25)^2+(Node977_cls(1,3)+0.25)^2);
r_979 = sqrt((Node979_cls(1,2)+0.25)^2+(Node979_cls(1,3)+0.25)^2);
r_981 = sqrt((Node981_cls(1,2)+0.25)^2+(Node981_cls(1,3)+0.25)^2);
r_983 = sqrt((Node983_cls(1,2)+0.25)^2+(Node983_cls(1,3)+0.25)^2);
r_985 = sqrt((Node985_cls(1,2)+0.25)^2+(Node985_cls(1,3)+0.25)^2);
r_987 = sqrt((Node987_cls(1,2)+0.25)^2+(Node987_cls(1,3)+0.25)^2);
r_989 = sqrt((Node989_cls(1,2)+0.25)^2+(Node989_cls(1,3)+0.25)^2);
r_991 = sqrt((Node991_cls(1,2)+0.25)^2+(Node991_cls(1,3)+0.25)^2);
r_993 = sqrt((Node993_cls(1,2)+0.25)^2+(Node993_cls(1,3)+0.25)^2);
r_995 = sqrt((Node995_cls(1,2)+0.25)^2+(Node995_cls(1,3)+0.25)^2);
r_997 = sqrt((Node997_cls(1,2)+0.25)^2+(Node997_cls(1,3)+0.25)^2);
r_999 = sqrt((Node999_cls(1,2)+0.25)^2+(Node999_cls(1,3)+0.25)^2);
r_1001 = sqrt((Node1001_cls(1,2)+0.25)^2+(Node1001_cls(1,3)+0.25)^2);
r_1003 = sqrt((Node1003_cls(1,2)+0.25)^2+(Node1003_cls(1,3)+0.25)^2);
r_1005 = sqrt((Node1005_cls(1,2)+0.25)^2+(Node1005_cls(1,3)+0.25)^2);
r_1007 = sqrt((Node1007_cls(1,2)+0.25)^2+(Node1007_cls(1,3)+0.25)^2);

Result_S_r_teta_cls = [r_938/0.05 Result_S_938_AD_cls_trans(1,2)/P_applied;r_939/0.05  Result_S_939_AD_cls_trans(1,2)/P_applied;...
    r_941/0.05  Result_S_941_AD_cls_trans(1,2)/P_applied;r_943/0.05  Result_S_943_AD_cls_trans(1,2)/P_applied;r_945/0.05  Result_S_945_AD_cls_trans(1,2)/P_applied;...
    r_947/0.05  Result_S_947_AD_cls_trans(1,2)/P_applied;r_949/0.05  Result_S_949_AD_cls_trans(1,2)/P_applied;r_951/0.05  Result_S_951_AD_cls_trans(1,2)/P_applied;...
    r_953/0.05  Result_S_953_AD_cls_trans(1,2)/P_applied;r_955/0.05  Result_S_955_AD_cls_trans(1,2)/P_applied;r_957/0.05  Result_S_957_AD_cls_trans(1,2)/P_applied;...
    r_959/0.05  Result_S_959_AD_cls_trans(1,2)/P_applied;r_961/0.05  Result_S_961_AD_cls_trans(1,2)/P_applied;r_963/0.05  Result_S_963_AD_cls_trans(1,2)/P_applied;...
    r_965/0.05  Result_S_965_AD_cls_trans(1,2)/P_applied;r_967/0.05  Result_S_967_AD_cls_trans(1,2)/P_applied;r_969/0.05  Result_S_969_AD_cls_trans(1,2)/P_applied;...
    r_971/0.05  Result_S_971_AD_cls_trans(1,2)/P_applied;r_973/0.05  Result_S_973_AD_cls_trans(1,2)/P_applied;r_975/0.05  Result_S_975_AD_cls_trans(1,2)/P_applied;...
    r_977/0.05  Result_S_977_AD_cls_trans(1,2)/P_applied;r_979/0.05  Result_S_979_AD_cls_trans(1,2)/P_applied;r_981/0.05  Result_S_981_AD_cls_trans(1,2)/P_applied;...
    r_983/0.05  Result_S_983_AD_cls_trans(1,2)/P_applied;r_985/0.05  Result_S_985_AD_cls_trans(1,2)/P_applied;r_987/0.05  Result_S_987_AD_cls_trans(1,2)/P_applied;...
    r_989/0.05  Result_S_989_AD_cls_trans(1,2)/P_applied;r_991/0.05  Result_S_991_AD_cls_trans(1,2)/P_applied;r_993/0.05  Result_S_993_AD_cls_trans(1,2)/P_applied;...
    r_995/0.05  Result_S_995_AD_cls_trans(1,2)/P_applied;r_997/0.05  Result_S_997_AD_cls_trans(1,2)/P_applied;r_999/0.05  Result_S_999_AD_cls_trans(1,2)/P_applied;...
    r_1001/0.05  Result_S_1001_AD_cls_trans(1,2)/P_applied;r_1003/0.05  Result_S_1003_AD_cls_trans(1,2)/P_applied;r_1005/0.05  Result_S_1005_AD_cls_trans(1,2)/P_applied;...
    r_1007/0.05  Result_S_1007_AD_cls_trans(1,2)/P_applied];


Result_S_938_AD_MBC_com = [Node938_MBC_com(1,7) Node938_MBC_com(1,8) Node938_MBC_com(1,9);Node938_MBC_com(1,10) Node938_MBC_com(1,11) Node938_MBC_com(1,12);...
    Node938_MBC_com(1,13) Node938_MBC_com(1,14) Node938_MBC_com(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_938_AD_MBC_com_trans = R'*Result_S_938_AD_MBC_com*R;


Result_S_939_AD_MBC_com = [Node939_MBC_com(1,7) Node939_MBC_com(1,8) Node939_MBC_com(1,9);Node939_MBC_com(1,10) Node939_MBC_com(1,11) Node939_MBC_com(1,12);...
    Node939_MBC_com(1,13) Node939_MBC_com(1,14) Node939_MBC_com(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_939_AD_MBC_com_trans = R'*Result_S_939_AD_MBC_com*R;


Result_S_941_AD_MBC_com = [Node941_MBC_com(1,7) Node941_MBC_com(1,8) Node941_MBC_com(1,9);Node941_MBC_com(1,10) Node941_MBC_com(1,11) Node941_MBC_com(1,12);...
    Node941_MBC_com(1,13) Node941_MBC_com(1,14) Node941_MBC_com(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_941_AD_MBC_com_trans = R'*Result_S_941_AD_MBC_com*R;


Result_S_943_AD_MBC_com = [Node943_MBC_com(1,7) Node943_MBC_com(1,8) Node943_MBC_com(1,9);Node943_MBC_com(1,10) Node943_MBC_com(1,11) Node943_MBC_com(1,12);...
    Node943_MBC_com(1,13) Node943_MBC_com(1,14) Node943_MBC_com(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_943_AD_MBC_com_trans = R'*Result_S_943_AD_MBC_com*R;


Result_S_945_AD_MBC_com = [Node945_MBC_com(1,7) Node945_MBC_com(1,8) Node945_MBC_com(1,9);Node945_MBC_com(1,10) Node945_MBC_com(1,11) Node945_MBC_com(1,12);...
    Node945_MBC_com(1,13) Node945_MBC_com(1,14) Node945_MBC_com(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_945_AD_MBC_com_trans = R'*Result_S_945_AD_MBC_com*R;


Result_S_947_AD_MBC_com = [Node947_MBC_com(1,7) Node947_MBC_com(1,8) Node947_MBC_com(1,9);Node947_MBC_com(1,10) Node947_MBC_com(1,11) Node947_MBC_com(1,12);...
    Node947_MBC_com(1,13) Node947_MBC_com(1,14) Node947_MBC_com(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_947_AD_MBC_com_trans = R'*Result_S_947_AD_MBC_com*R;


Result_S_949_AD_MBC_com = [Node949_MBC_com(1,7) Node949_MBC_com(1,8) Node949_MBC_com(1,9);Node949_MBC_com(1,10) Node949_MBC_com(1,11) Node949_MBC_com(1,12);...
    Node949_MBC_com(1,13) Node949_MBC_com(1,14) Node949_MBC_com(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_949_AD_MBC_com_trans = R'*Result_S_949_AD_MBC_com*R;


Result_S_951_AD_MBC_com = [Node951_MBC_com(1,7) Node951_MBC_com(1,8) Node951_MBC_com(1,9);Node951_MBC_com(1,10) Node951_MBC_com(1,11) Node951_MBC_com(1,12);...
    Node951_MBC_com(1,13) Node951_MBC_com(1,14) Node951_MBC_com(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_951_AD_MBC_com_trans = R'*Result_S_951_AD_MBC_com*R;


Result_S_953_AD_MBC_com = [Node953_MBC_com(1,7) Node953_MBC_com(1,8) Node953_MBC_com(1,9);Node953_MBC_com(1,10) Node953_MBC_com(1,11) Node953_MBC_com(1,12);...
    Node953_MBC_com(1,13) Node953_MBC_com(1,14) Node953_MBC_com(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_953_AD_MBC_com_trans = R'*Result_S_953_AD_MBC_com*R;


Result_S_955_AD_MBC_com = [Node955_MBC_com(1,7) Node955_MBC_com(1,8) Node955_MBC_com(1,9);Node955_MBC_com(1,10) Node955_MBC_com(1,11) Node955_MBC_com(1,12);...
    Node955_MBC_com(1,13) Node955_MBC_com(1,14) Node955_MBC_com(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_955_AD_MBC_com_trans = R'*Result_S_955_AD_MBC_com*R;


Result_S_957_AD_MBC_com = [Node957_MBC_com(1,7) Node957_MBC_com(1,8) Node957_MBC_com(1,9);Node957_MBC_com(1,10) Node957_MBC_com(1,11) Node957_MBC_com(1,12);...
    Node957_MBC_com(1,13) Node957_MBC_com(1,14) Node957_MBC_com(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_957_AD_MBC_com_trans = R'*Result_S_957_AD_MBC_com*R;


Result_S_959_AD_MBC_com = [Node959_MBC_com(1,7) Node959_MBC_com(1,8) Node959_MBC_com(1,9);Node959_MBC_com(1,10) Node959_MBC_com(1,11) Node959_MBC_com(1,12);...
    Node959_MBC_com(1,13) Node959_MBC_com(1,14) Node959_MBC_com(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_959_AD_MBC_com_trans = R'*Result_S_959_AD_MBC_com*R;


Result_S_961_AD_MBC_com = [Node961_MBC_com(1,7) Node961_MBC_com(1,8) Node961_MBC_com(1,9);Node961_MBC_com(1,10) Node961_MBC_com(1,11) Node961_MBC_com(1,12);...
    Node961_MBC_com(1,13) Node961_MBC_com(1,14) Node961_MBC_com(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_961_AD_MBC_com_trans = R'*Result_S_961_AD_MBC_com*R;

Result_S_963_AD_MBC_com = [Node963_MBC_com(1,7) Node963_MBC_com(1,8) Node963_MBC_com(1,9);Node963_MBC_com(1,10) Node963_MBC_com(1,11) Node963_MBC_com(1,12);...
    Node963_MBC_com(1,13) Node963_MBC_com(1,14) Node963_MBC_com(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_963_AD_MBC_com_trans = R'*Result_S_963_AD_MBC_com*R;


Result_S_965_AD_MBC_com = [Node965_MBC_com(1,7) Node965_MBC_com(1,8) Node965_MBC_com(1,9);Node965_MBC_com(1,10) Node965_MBC_com(1,11) Node965_MBC_com(1,12);...
    Node965_MBC_com(1,13) Node965_MBC_com(1,14) Node965_MBC_com(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_965_AD_MBC_com_trans = R'*Result_S_965_AD_MBC_com*R;

Result_S_967_AD_MBC_com = [Node967_MBC_com(1,7) Node967_MBC_com(1,8) Node967_MBC_com(1,9);Node967_MBC_com(1,10) Node967_MBC_com(1,11) Node967_MBC_com(1,12);...
    Node967_MBC_com(1,13) Node967_MBC_com(1,14) Node967_MBC_com(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_967_AD_MBC_com_trans = R'*Result_S_967_AD_MBC_com*R;

Result_S_969_AD_MBC_com = [Node969_MBC_com(1,7) Node969_MBC_com(1,8) Node969_MBC_com(1,9);Node969_MBC_com(1,10) Node969_MBC_com(1,11) Node969_MBC_com(1,12);...
    Node969_MBC_com(1,13) Node969_MBC_com(1,14) Node969_MBC_com(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_969_AD_MBC_com_trans = R'*Result_S_969_AD_MBC_com*R;

Result_S_971_AD_MBC_com = [Node971_MBC_com(1,7) Node971_MBC_com(1,8) Node971_MBC_com(1,9);Node971_MBC_com(1,10) Node971_MBC_com(1,11) Node971_MBC_com(1,12);...
    Node971_MBC_com(1,13) Node971_MBC_com(1,14) Node971_MBC_com(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_971_AD_MBC_com_trans = R'*Result_S_971_AD_MBC_com*R;

Result_S_973_AD_MBC_com = [Node973_MBC_com(1,7) Node973_MBC_com(1,8) Node973_MBC_com(1,9);Node973_MBC_com(1,10) Node973_MBC_com(1,11) Node973_MBC_com(1,12);...
    Node973_MBC_com(1,13) Node973_MBC_com(1,14) Node973_MBC_com(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_973_AD_MBC_com_trans = R'*Result_S_973_AD_MBC_com*R;

Result_S_975_AD_MBC_com = [Node975_MBC_com(1,7) Node975_MBC_com(1,8) Node975_MBC_com(1,9);Node975_MBC_com(1,10) Node975_MBC_com(1,11) Node975_MBC_com(1,12);...
    Node975_MBC_com(1,13) Node975_MBC_com(1,14) Node975_MBC_com(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_975_AD_MBC_com_trans = R'*Result_S_975_AD_MBC_com*R;

Result_S_977_AD_MBC_com = [Node977_MBC_com(1,7) Node977_MBC_com(1,8) Node977_MBC_com(1,9);Node977_MBC_com(1,10) Node977_MBC_com(1,11) Node977_MBC_com(1,12);...
    Node977_MBC_com(1,13) Node977_MBC_com(1,14) Node977_MBC_com(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_977_AD_MBC_com_trans = R'*Result_S_977_AD_MBC_com*R;

Result_S_979_AD_MBC_com = [Node979_MBC_com(1,7) Node979_MBC_com(1,8) Node979_MBC_com(1,9);Node979_MBC_com(1,10) Node979_MBC_com(1,11) Node979_MBC_com(1,12);...
    Node979_MBC_com(1,13) Node979_MBC_com(1,14) Node979_MBC_com(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_979_AD_MBC_com_trans = R'*Result_S_979_AD_MBC_com*R;

Result_S_981_AD_MBC_com = [Node981_MBC_com(1,7) Node981_MBC_com(1,8) Node981_MBC_com(1,9);Node981_MBC_com(1,10) Node981_MBC_com(1,11) Node981_MBC_com(1,12);...
    Node981_MBC_com(1,13) Node981_MBC_com(1,14) Node981_MBC_com(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_981_AD_MBC_com_trans = R'*Result_S_981_AD_MBC_com*R;

Result_S_983_AD_MBC_com = [Node983_MBC_com(1,7) Node983_MBC_com(1,8) Node983_MBC_com(1,9);Node983_MBC_com(1,10) Node983_MBC_com(1,11) Node983_MBC_com(1,12);...
    Node983_MBC_com(1,13) Node983_MBC_com(1,14) Node983_MBC_com(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_983_AD_MBC_com_trans = R'*Result_S_983_AD_MBC_com*R;

Result_S_985_AD_MBC_com = [Node985_MBC_com(1,7) Node985_MBC_com(1,8) Node985_MBC_com(1,9);Node985_MBC_com(1,10) Node985_MBC_com(1,11) Node985_MBC_com(1,12);...
    Node985_MBC_com(1,13) Node985_MBC_com(1,14) Node985_MBC_com(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_985_AD_MBC_com_trans = R'*Result_S_985_AD_MBC_com*R;

Result_S_987_AD_MBC_com = [Node987_MBC_com(1,7) Node987_MBC_com(1,8) Node987_MBC_com(1,9);Node987_MBC_com(1,10) Node987_MBC_com(1,11) Node987_MBC_com(1,12);...
    Node987_MBC_com(1,13) Node987_MBC_com(1,14) Node987_MBC_com(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_987_AD_MBC_com_trans = R'*Result_S_987_AD_MBC_com*R;

Result_S_989_AD_MBC_com = [Node989_MBC_com(1,7) Node989_MBC_com(1,8) Node989_MBC_com(1,9);Node989_MBC_com(1,10) Node989_MBC_com(1,11) Node989_MBC_com(1,12);...
    Node989_MBC_com(1,13) Node989_MBC_com(1,14) Node989_MBC_com(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_989_AD_MBC_com_trans = R'*Result_S_989_AD_MBC_com*R;

Result_S_991_AD_MBC_com = [Node991_MBC_com(1,7) Node991_MBC_com(1,8) Node991_MBC_com(1,9);Node991_MBC_com(1,10) Node991_MBC_com(1,11) Node991_MBC_com(1,12);...
    Node991_MBC_com(1,13) Node991_MBC_com(1,14) Node991_MBC_com(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_991_AD_MBC_com_trans = R'*Result_S_991_AD_MBC_com*R;

Result_S_993_AD_MBC_com = [Node993_MBC_com(1,7) Node993_MBC_com(1,8) Node993_MBC_com(1,9);Node993_MBC_com(1,10) Node993_MBC_com(1,11) Node993_MBC_com(1,12);...
    Node993_MBC_com(1,13) Node993_MBC_com(1,14) Node993_MBC_com(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_993_AD_MBC_com_trans = R'*Result_S_993_AD_MBC_com*R;

Result_S_995_AD_MBC_com = [Node995_MBC_com(1,7) Node995_MBC_com(1,8) Node995_MBC_com(1,9);Node995_MBC_com(1,10) Node995_MBC_com(1,11) Node995_MBC_com(1,12);...
    Node995_MBC_com(1,13) Node995_MBC_com(1,14) Node995_MBC_com(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_995_AD_MBC_com_trans = R'*Result_S_995_AD_MBC_com*R;

Result_S_997_AD_MBC_com = [Node997_MBC_com(1,7) Node997_MBC_com(1,8) Node997_MBC_com(1,9);Node997_MBC_com(1,10) Node997_MBC_com(1,11) Node997_MBC_com(1,12);...
    Node997_MBC_com(1,13) Node997_MBC_com(1,14) Node997_MBC_com(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_997_AD_MBC_com_trans = R'*Result_S_997_AD_MBC_com*R;

Result_S_999_AD_MBC_com = [Node999_MBC_com(1,7) Node999_MBC_com(1,8) Node999_MBC_com(1,9);Node999_MBC_com(1,10) Node999_MBC_com(1,11) Node999_MBC_com(1,12);...
    Node999_MBC_com(1,13) Node999_MBC_com(1,14) Node999_MBC_com(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_999_AD_MBC_com_trans = R'*Result_S_999_AD_MBC_com*R;

Result_S_1001_AD_MBC_com = [Node1001_MBC_com(1,7) Node1001_MBC_com(1,8) Node1001_MBC_com(1,9);Node1001_MBC_com(1,10) Node1001_MBC_com(1,11) Node1001_MBC_com(1,12);...
    Node1001_MBC_com(1,13) Node1001_MBC_com(1,14) Node1001_MBC_com(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_1001_AD_MBC_com_trans = R'*Result_S_1001_AD_MBC_com*R;

Result_S_1003_AD_MBC_com = [Node1003_MBC_com(1,7) Node1003_MBC_com(1,8) Node1003_MBC_com(1,9);Node1003_MBC_com(1,10) Node1003_MBC_com(1,11) Node1003_MBC_com(1,12);...
    Node1003_MBC_com(1,13) Node1003_MBC_com(1,14) Node1003_MBC_com(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_1003_AD_MBC_com_trans = R'*Result_S_1003_AD_MBC_com*R;

Result_S_1005_AD_MBC_com = [Node1005_MBC_com(1,7) Node1005_MBC_com(1,8) Node1005_MBC_com(1,9);Node1005_MBC_com(1,10) Node1005_MBC_com(1,11) Node1005_MBC_com(1,12);...
    Node1005_MBC_com(1,13) Node1005_MBC_com(1,14) Node1005_MBC_com(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_1005_AD_MBC_com_trans = R'*Result_S_1005_AD_MBC_com*R;

Result_S_1007_AD_MBC_com = [Node1007_MBC_com(1,7) Node1007_MBC_com(1,8) Node1007_MBC_com(1,9);Node1007_MBC_com(1,10) Node1007_MBC_com(1,11) Node1007_MBC_com(1,12);...
    Node1007_MBC_com(1,13) Node1007_MBC_com(1,14) Node1007_MBC_com(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_1007_AD_MBC_com_trans = R'*Result_S_1007_AD_MBC_com*R;


Result_S_r_teta_MBC_com = [r_938/0.05 Result_S_938_AD_MBC_com_trans(1,2)/P_applied;r_939/0.05  Result_S_939_AD_MBC_com_trans(1,2)/P_applied;...
    r_941/0.05  Result_S_941_AD_MBC_com_trans(1,2)/P_applied;r_943/0.05  Result_S_943_AD_MBC_com_trans(1,2)/P_applied;r_945/0.05  Result_S_945_AD_MBC_com_trans(1,2)/P_applied;...
    r_947/0.05  Result_S_947_AD_MBC_com_trans(1,2)/P_applied;r_949/0.05  Result_S_949_AD_MBC_com_trans(1,2)/P_applied;r_951/0.05  Result_S_951_AD_MBC_com_trans(1,2)/P_applied;...
    r_953/0.05  Result_S_953_AD_MBC_com_trans(1,2)/P_applied;r_955/0.05  Result_S_955_AD_MBC_com_trans(1,2)/P_applied;r_957/0.05  Result_S_957_AD_MBC_com_trans(1,2)/P_applied;...
    r_959/0.05  Result_S_959_AD_MBC_com_trans(1,2)/P_applied;r_961/0.05  Result_S_961_AD_MBC_com_trans(1,2)/P_applied;r_963/0.05  Result_S_963_AD_MBC_com_trans(1,2)/P_applied;...
    r_965/0.05  Result_S_965_AD_MBC_com_trans(1,2)/P_applied;r_967/0.05  Result_S_967_AD_MBC_com_trans(1,2)/P_applied;r_969/0.05  Result_S_969_AD_MBC_com_trans(1,2)/P_applied;...
    r_971/0.05  Result_S_971_AD_MBC_com_trans(1,2)/P_applied;r_973/0.05  Result_S_973_AD_MBC_com_trans(1,2)/P_applied;r_975/0.05  Result_S_975_AD_MBC_com_trans(1,2)/P_applied;...
    r_977/0.05  Result_S_977_AD_MBC_com_trans(1,2)/P_applied;r_979/0.05  Result_S_979_AD_MBC_com_trans(1,2)/P_applied;r_981/0.05  Result_S_981_AD_MBC_com_trans(1,2)/P_applied;...
    r_983/0.05  Result_S_983_AD_MBC_com_trans(1,2)/P_applied;r_985/0.05  Result_S_985_AD_MBC_com_trans(1,2)/P_applied;r_987/0.05  Result_S_987_AD_MBC_com_trans(1,2)/P_applied;...
    r_989/0.05  Result_S_989_AD_MBC_com_trans(1,2)/P_applied;r_991/0.05  Result_S_991_AD_MBC_com_trans(1,2)/P_applied;r_993/0.05  Result_S_993_AD_MBC_com_trans(1,2)/P_applied;...
    r_995/0.05  Result_S_995_AD_MBC_com_trans(1,2)/P_applied;r_997/0.05  Result_S_997_AD_MBC_com_trans(1,2)/P_applied;r_999/0.05  Result_S_999_AD_MBC_com_trans(1,2)/P_applied;...
    r_1001/0.05  Result_S_1001_AD_MBC_com_trans(1,2)/P_applied;r_1003/0.05  Result_S_1003_AD_MBC_com_trans(1,2)/P_applied;r_1005/0.05  Result_S_1005_AD_MBC_com_trans(1,2)/P_applied;...
    r_1007/0.05  Result_S_1007_AD_MBC_com_trans(1,2)/P_applied];


Result_S_938_AD_MBC_rot = [Node938_MBC_rot(1,7) Node938_MBC_rot(1,8) Node938_MBC_rot(1,9);Node938_MBC_rot(1,10) Node938_MBC_rot(1,11) Node938_MBC_rot(1,12);...
    Node938_MBC_rot(1,13) Node938_MBC_rot(1,14) Node938_MBC_rot(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_938_AD_MBC_rot_trans = R'*Result_S_938_AD_MBC_rot*R;


Result_S_939_AD_MBC_rot = [Node939_MBC_rot(1,7) Node939_MBC_rot(1,8) Node939_MBC_rot(1,9);Node939_MBC_rot(1,10) Node939_MBC_rot(1,11) Node939_MBC_rot(1,12);...
    Node939_MBC_rot(1,13) Node939_MBC_rot(1,14) Node939_MBC_rot(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_939_AD_MBC_rot_trans = R'*Result_S_939_AD_MBC_rot*R;


Result_S_941_AD_MBC_rot = [Node941_MBC_rot(1,7) Node941_MBC_rot(1,8) Node941_MBC_rot(1,9);Node941_MBC_rot(1,10) Node941_MBC_rot(1,11) Node941_MBC_rot(1,12);...
    Node941_MBC_rot(1,13) Node941_MBC_rot(1,14) Node941_MBC_rot(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_941_AD_MBC_rot_trans = R'*Result_S_941_AD_MBC_rot*R;


Result_S_943_AD_MBC_rot = [Node943_MBC_rot(1,7) Node943_MBC_rot(1,8) Node943_MBC_rot(1,9);Node943_MBC_rot(1,10) Node943_MBC_rot(1,11) Node943_MBC_rot(1,12);...
    Node943_MBC_rot(1,13) Node943_MBC_rot(1,14) Node943_MBC_rot(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_943_AD_MBC_rot_trans = R'*Result_S_943_AD_MBC_rot*R;


Result_S_945_AD_MBC_rot = [Node945_MBC_rot(1,7) Node945_MBC_rot(1,8) Node945_MBC_rot(1,9);Node945_MBC_rot(1,10) Node945_MBC_rot(1,11) Node945_MBC_rot(1,12);...
    Node945_MBC_rot(1,13) Node945_MBC_rot(1,14) Node945_MBC_rot(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_945_AD_MBC_rot_trans = R'*Result_S_945_AD_MBC_rot*R;


Result_S_947_AD_MBC_rot = [Node947_MBC_rot(1,7) Node947_MBC_rot(1,8) Node947_MBC_rot(1,9);Node947_MBC_rot(1,10) Node947_MBC_rot(1,11) Node947_MBC_rot(1,12);...
    Node947_MBC_rot(1,13) Node947_MBC_rot(1,14) Node947_MBC_rot(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_947_AD_MBC_rot_trans = R'*Result_S_947_AD_MBC_rot*R;


Result_S_949_AD_MBC_rot = [Node949_MBC_rot(1,7) Node949_MBC_rot(1,8) Node949_MBC_rot(1,9);Node949_MBC_rot(1,10) Node949_MBC_rot(1,11) Node949_MBC_rot(1,12);...
    Node949_MBC_rot(1,13) Node949_MBC_rot(1,14) Node949_MBC_rot(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_949_AD_MBC_rot_trans = R'*Result_S_949_AD_MBC_rot*R;


Result_S_951_AD_MBC_rot = [Node951_MBC_rot(1,7) Node951_MBC_rot(1,8) Node951_MBC_rot(1,9);Node951_MBC_rot(1,10) Node951_MBC_rot(1,11) Node951_MBC_rot(1,12);...
    Node951_MBC_rot(1,13) Node951_MBC_rot(1,14) Node951_MBC_rot(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_951_AD_MBC_rot_trans = R'*Result_S_951_AD_MBC_rot*R;


Result_S_953_AD_MBC_rot = [Node953_MBC_rot(1,7) Node953_MBC_rot(1,8) Node953_MBC_rot(1,9);Node953_MBC_rot(1,10) Node953_MBC_rot(1,11) Node953_MBC_rot(1,12);...
    Node953_MBC_rot(1,13) Node953_MBC_rot(1,14) Node953_MBC_rot(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_953_AD_MBC_rot_trans = R'*Result_S_953_AD_MBC_rot*R;


Result_S_955_AD_MBC_rot = [Node955_MBC_rot(1,7) Node955_MBC_rot(1,8) Node955_MBC_rot(1,9);Node955_MBC_rot(1,10) Node955_MBC_rot(1,11) Node955_MBC_rot(1,12);...
    Node955_MBC_rot(1,13) Node955_MBC_rot(1,14) Node955_MBC_rot(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_955_AD_MBC_rot_trans = R'*Result_S_955_AD_MBC_rot*R;


Result_S_957_AD_MBC_rot = [Node957_MBC_rot(1,7) Node957_MBC_rot(1,8) Node957_MBC_rot(1,9);Node957_MBC_rot(1,10) Node957_MBC_rot(1,11) Node957_MBC_rot(1,12);...
    Node957_MBC_rot(1,13) Node957_MBC_rot(1,14) Node957_MBC_rot(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_957_AD_MBC_rot_trans = R'*Result_S_957_AD_MBC_rot*R;


Result_S_959_AD_MBC_rot = [Node959_MBC_rot(1,7) Node959_MBC_rot(1,8) Node959_MBC_rot(1,9);Node959_MBC_rot(1,10) Node959_MBC_rot(1,11) Node959_MBC_rot(1,12);...
    Node959_MBC_rot(1,13) Node959_MBC_rot(1,14) Node959_MBC_rot(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_959_AD_MBC_rot_trans = R'*Result_S_959_AD_MBC_rot*R;


Result_S_961_AD_MBC_rot = [Node961_MBC_rot(1,7) Node961_MBC_rot(1,8) Node961_MBC_rot(1,9);Node961_MBC_rot(1,10) Node961_MBC_rot(1,11) Node961_MBC_rot(1,12);...
    Node961_MBC_rot(1,13) Node961_MBC_rot(1,14) Node961_MBC_rot(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_961_AD_MBC_rot_trans = R'*Result_S_961_AD_MBC_rot*R;

Result_S_963_AD_MBC_rot = [Node963_MBC_rot(1,7) Node963_MBC_rot(1,8) Node963_MBC_rot(1,9);Node963_MBC_rot(1,10) Node963_MBC_rot(1,11) Node963_MBC_rot(1,12);...
    Node963_MBC_rot(1,13) Node963_MBC_rot(1,14) Node963_MBC_rot(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_963_AD_MBC_rot_trans = R'*Result_S_963_AD_MBC_rot*R;


Result_S_965_AD_MBC_rot = [Node965_MBC_rot(1,7) Node965_MBC_rot(1,8) Node965_MBC_rot(1,9);Node965_MBC_rot(1,10) Node965_MBC_rot(1,11) Node965_MBC_rot(1,12);...
    Node965_MBC_rot(1,13) Node965_MBC_rot(1,14) Node965_MBC_rot(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_965_AD_MBC_rot_trans = R'*Result_S_965_AD_MBC_rot*R;

Result_S_967_AD_MBC_rot = [Node967_MBC_rot(1,7) Node967_MBC_rot(1,8) Node967_MBC_rot(1,9);Node967_MBC_rot(1,10) Node967_MBC_rot(1,11) Node967_MBC_rot(1,12);...
    Node967_MBC_rot(1,13) Node967_MBC_rot(1,14) Node967_MBC_rot(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_967_AD_MBC_rot_trans = R'*Result_S_967_AD_MBC_rot*R;

Result_S_969_AD_MBC_rot = [Node969_MBC_rot(1,7) Node969_MBC_rot(1,8) Node969_MBC_rot(1,9);Node969_MBC_rot(1,10) Node969_MBC_rot(1,11) Node969_MBC_rot(1,12);...
    Node969_MBC_rot(1,13) Node969_MBC_rot(1,14) Node969_MBC_rot(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_969_AD_MBC_rot_trans = R'*Result_S_969_AD_MBC_rot*R;

Result_S_971_AD_MBC_rot = [Node971_MBC_rot(1,7) Node971_MBC_rot(1,8) Node971_MBC_rot(1,9);Node971_MBC_rot(1,10) Node971_MBC_rot(1,11) Node971_MBC_rot(1,12);...
    Node971_MBC_rot(1,13) Node971_MBC_rot(1,14) Node971_MBC_rot(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_971_AD_MBC_rot_trans = R'*Result_S_971_AD_MBC_rot*R;

Result_S_973_AD_MBC_rot = [Node973_MBC_rot(1,7) Node973_MBC_rot(1,8) Node973_MBC_rot(1,9);Node973_MBC_rot(1,10) Node973_MBC_rot(1,11) Node973_MBC_rot(1,12);...
    Node973_MBC_rot(1,13) Node973_MBC_rot(1,14) Node973_MBC_rot(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_973_AD_MBC_rot_trans = R'*Result_S_973_AD_MBC_rot*R;

Result_S_975_AD_MBC_rot = [Node975_MBC_rot(1,7) Node975_MBC_rot(1,8) Node975_MBC_rot(1,9);Node975_MBC_rot(1,10) Node975_MBC_rot(1,11) Node975_MBC_rot(1,12);...
    Node975_MBC_rot(1,13) Node975_MBC_rot(1,14) Node975_MBC_rot(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_975_AD_MBC_rot_trans = R'*Result_S_975_AD_MBC_rot*R;

Result_S_977_AD_MBC_rot = [Node977_MBC_rot(1,7) Node977_MBC_rot(1,8) Node977_MBC_rot(1,9);Node977_MBC_rot(1,10) Node977_MBC_rot(1,11) Node977_MBC_rot(1,12);...
    Node977_MBC_rot(1,13) Node977_MBC_rot(1,14) Node977_MBC_rot(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_977_AD_MBC_rot_trans = R'*Result_S_977_AD_MBC_rot*R;

Result_S_979_AD_MBC_rot = [Node979_MBC_rot(1,7) Node979_MBC_rot(1,8) Node979_MBC_rot(1,9);Node979_MBC_rot(1,10) Node979_MBC_rot(1,11) Node979_MBC_rot(1,12);...
    Node979_MBC_rot(1,13) Node979_MBC_rot(1,14) Node979_MBC_rot(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_979_AD_MBC_rot_trans = R'*Result_S_979_AD_MBC_rot*R;

Result_S_981_AD_MBC_rot = [Node981_MBC_rot(1,7) Node981_MBC_rot(1,8) Node981_MBC_rot(1,9);Node981_MBC_rot(1,10) Node981_MBC_rot(1,11) Node981_MBC_rot(1,12);...
    Node981_MBC_rot(1,13) Node981_MBC_rot(1,14) Node981_MBC_rot(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_981_AD_MBC_rot_trans = R'*Result_S_981_AD_MBC_rot*R;

Result_S_983_AD_MBC_rot = [Node983_MBC_rot(1,7) Node983_MBC_rot(1,8) Node983_MBC_rot(1,9);Node983_MBC_rot(1,10) Node983_MBC_rot(1,11) Node983_MBC_rot(1,12);...
    Node983_MBC_rot(1,13) Node983_MBC_rot(1,14) Node983_MBC_rot(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_983_AD_MBC_rot_trans = R'*Result_S_983_AD_MBC_rot*R;

Result_S_985_AD_MBC_rot = [Node985_MBC_rot(1,7) Node985_MBC_rot(1,8) Node985_MBC_rot(1,9);Node985_MBC_rot(1,10) Node985_MBC_rot(1,11) Node985_MBC_rot(1,12);...
    Node985_MBC_rot(1,13) Node985_MBC_rot(1,14) Node985_MBC_rot(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_985_AD_MBC_rot_trans = R'*Result_S_985_AD_MBC_rot*R;

Result_S_987_AD_MBC_rot = [Node987_MBC_rot(1,7) Node987_MBC_rot(1,8) Node987_MBC_rot(1,9);Node987_MBC_rot(1,10) Node987_MBC_rot(1,11) Node987_MBC_rot(1,12);...
    Node987_MBC_rot(1,13) Node987_MBC_rot(1,14) Node987_MBC_rot(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_987_AD_MBC_rot_trans = R'*Result_S_987_AD_MBC_rot*R;

Result_S_989_AD_MBC_rot = [Node989_MBC_rot(1,7) Node989_MBC_rot(1,8) Node989_MBC_rot(1,9);Node989_MBC_rot(1,10) Node989_MBC_rot(1,11) Node989_MBC_rot(1,12);...
    Node989_MBC_rot(1,13) Node989_MBC_rot(1,14) Node989_MBC_rot(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_989_AD_MBC_rot_trans = R'*Result_S_989_AD_MBC_rot*R;

Result_S_991_AD_MBC_rot = [Node991_MBC_rot(1,7) Node991_MBC_rot(1,8) Node991_MBC_rot(1,9);Node991_MBC_rot(1,10) Node991_MBC_rot(1,11) Node991_MBC_rot(1,12);...
    Node991_MBC_rot(1,13) Node991_MBC_rot(1,14) Node991_MBC_rot(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_991_AD_MBC_rot_trans = R'*Result_S_991_AD_MBC_rot*R;

Result_S_993_AD_MBC_rot = [Node993_MBC_rot(1,7) Node993_MBC_rot(1,8) Node993_MBC_rot(1,9);Node993_MBC_rot(1,10) Node993_MBC_rot(1,11) Node993_MBC_rot(1,12);...
    Node993_MBC_rot(1,13) Node993_MBC_rot(1,14) Node993_MBC_rot(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_993_AD_MBC_rot_trans = R'*Result_S_993_AD_MBC_rot*R;

Result_S_995_AD_MBC_rot = [Node995_MBC_rot(1,7) Node995_MBC_rot(1,8) Node995_MBC_rot(1,9);Node995_MBC_rot(1,10) Node995_MBC_rot(1,11) Node995_MBC_rot(1,12);...
    Node995_MBC_rot(1,13) Node995_MBC_rot(1,14) Node995_MBC_rot(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_995_AD_MBC_rot_trans = R'*Result_S_995_AD_MBC_rot*R;

Result_S_997_AD_MBC_rot = [Node997_MBC_rot(1,7) Node997_MBC_rot(1,8) Node997_MBC_rot(1,9);Node997_MBC_rot(1,10) Node997_MBC_rot(1,11) Node997_MBC_rot(1,12);...
    Node997_MBC_rot(1,13) Node997_MBC_rot(1,14) Node997_MBC_rot(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_997_AD_MBC_rot_trans = R'*Result_S_997_AD_MBC_rot*R;

Result_S_999_AD_MBC_rot = [Node999_MBC_rot(1,7) Node999_MBC_rot(1,8) Node999_MBC_rot(1,9);Node999_MBC_rot(1,10) Node999_MBC_rot(1,11) Node999_MBC_rot(1,12);...
    Node999_MBC_rot(1,13) Node999_MBC_rot(1,14) Node999_MBC_rot(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_999_AD_MBC_rot_trans = R'*Result_S_999_AD_MBC_rot*R;

Result_S_1001_AD_MBC_rot = [Node1001_MBC_rot(1,7) Node1001_MBC_rot(1,8) Node1001_MBC_rot(1,9);Node1001_MBC_rot(1,10) Node1001_MBC_rot(1,11) Node1001_MBC_rot(1,12);...
    Node1001_MBC_rot(1,13) Node1001_MBC_rot(1,14) Node1001_MBC_rot(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_1001_AD_MBC_rot_trans = R'*Result_S_1001_AD_MBC_rot*R;

Result_S_1003_AD_MBC_rot = [Node1003_MBC_rot(1,7) Node1003_MBC_rot(1,8) Node1003_MBC_rot(1,9);Node1003_MBC_rot(1,10) Node1003_MBC_rot(1,11) Node1003_MBC_rot(1,12);...
    Node1003_MBC_rot(1,13) Node1003_MBC_rot(1,14) Node1003_MBC_rot(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_1003_AD_MBC_rot_trans = R'*Result_S_1003_AD_MBC_rot*R;

Result_S_1005_AD_MBC_rot = [Node1005_MBC_rot(1,7) Node1005_MBC_rot(1,8) Node1005_MBC_rot(1,9);Node1005_MBC_rot(1,10) Node1005_MBC_rot(1,11) Node1005_MBC_rot(1,12);...
    Node1005_MBC_rot(1,13) Node1005_MBC_rot(1,14) Node1005_MBC_rot(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_1005_AD_MBC_rot_trans = R'*Result_S_1005_AD_MBC_rot*R;

Result_S_1007_AD_MBC_rot = [Node1007_MBC_rot(1,7) Node1007_MBC_rot(1,8) Node1007_MBC_rot(1,9);Node1007_MBC_rot(1,10) Node1007_MBC_rot(1,11) Node1007_MBC_rot(1,12);...
    Node1007_MBC_rot(1,13) Node1007_MBC_rot(1,14) Node1007_MBC_rot(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_1007_AD_MBC_rot_trans = R'*Result_S_1007_AD_MBC_rot*R;


Result_S_r_teta_MBC_rot = [r_938/0.05 Result_S_938_AD_MBC_rot_trans(1,2)/P_applied;r_939/0.05  Result_S_939_AD_MBC_rot_trans(1,2)/P_applied;...
    r_941/0.05  Result_S_941_AD_MBC_rot_trans(1,2)/P_applied;r_943/0.05  Result_S_943_AD_MBC_rot_trans(1,2)/P_applied;r_945/0.05  Result_S_945_AD_MBC_rot_trans(1,2)/P_applied;...
    r_947/0.05  Result_S_947_AD_MBC_rot_trans(1,2)/P_applied;r_949/0.05  Result_S_949_AD_MBC_rot_trans(1,2)/P_applied;r_951/0.05  Result_S_951_AD_MBC_rot_trans(1,2)/P_applied;...
    r_953/0.05  Result_S_953_AD_MBC_rot_trans(1,2)/P_applied;r_955/0.05  Result_S_955_AD_MBC_rot_trans(1,2)/P_applied;r_957/0.05  Result_S_957_AD_MBC_rot_trans(1,2)/P_applied;...
    r_959/0.05  Result_S_959_AD_MBC_rot_trans(1,2)/P_applied;r_961/0.05  Result_S_961_AD_MBC_rot_trans(1,2)/P_applied;r_963/0.05  Result_S_963_AD_MBC_rot_trans(1,2)/P_applied;...
    r_965/0.05  Result_S_965_AD_MBC_rot_trans(1,2)/P_applied;r_967/0.05  Result_S_967_AD_MBC_rot_trans(1,2)/P_applied;r_969/0.05  Result_S_969_AD_MBC_rot_trans(1,2)/P_applied;...
    r_971/0.05  Result_S_971_AD_MBC_rot_trans(1,2)/P_applied;r_973/0.05  Result_S_973_AD_MBC_rot_trans(1,2)/P_applied;r_975/0.05  Result_S_975_AD_MBC_rot_trans(1,2)/P_applied;...
    r_977/0.05  Result_S_977_AD_MBC_rot_trans(1,2)/P_applied;r_979/0.05  Result_S_979_AD_MBC_rot_trans(1,2)/P_applied;r_981/0.05  Result_S_981_AD_MBC_rot_trans(1,2)/P_applied;...
    r_983/0.05  Result_S_983_AD_MBC_rot_trans(1,2)/P_applied;r_985/0.05  Result_S_985_AD_MBC_rot_trans(1,2)/P_applied;r_987/0.05  Result_S_987_AD_MBC_rot_trans(1,2)/P_applied;...
    r_989/0.05  Result_S_989_AD_MBC_rot_trans(1,2)/P_applied;r_991/0.05  Result_S_991_AD_MBC_rot_trans(1,2)/P_applied;r_993/0.05  Result_S_993_AD_MBC_rot_trans(1,2)/P_applied;...
    r_995/0.05  Result_S_995_AD_MBC_rot_trans(1,2)/P_applied;r_997/0.05  Result_S_997_AD_MBC_rot_trans(1,2)/P_applied;r_999/0.05  Result_S_999_AD_MBC_rot_trans(1,2)/P_applied;...
    r_1001/0.05  Result_S_1001_AD_MBC_rot_trans(1,2)/P_applied;r_1003/0.05  Result_S_1003_AD_MBC_rot_trans(1,2)/P_applied;r_1005/0.05  Result_S_1005_AD_MBC_rot_trans(1,2)/P_applied;...
    r_1007/0.05  Result_S_1007_AD_MBC_rot_trans(1,2)/P_applied];


Result_S_938_AD_MBC_str = [Node938_MBC_str(1,7) Node938_MBC_str(1,8) Node938_MBC_str(1,9);Node938_MBC_str(1,10) Node938_MBC_str(1,11) Node938_MBC_str(1,12);...
    Node938_MBC_str(1,13) Node938_MBC_str(1,14) Node938_MBC_str(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_938_AD_MBC_str_trans = R'*Result_S_938_AD_MBC_str*R;


Result_S_939_AD_MBC_str = [Node939_MBC_str(1,7) Node939_MBC_str(1,8) Node939_MBC_str(1,9);Node939_MBC_str(1,10) Node939_MBC_str(1,11) Node939_MBC_str(1,12);...
    Node939_MBC_str(1,13) Node939_MBC_str(1,14) Node939_MBC_str(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_939_AD_MBC_str_trans = R'*Result_S_939_AD_MBC_str*R;


Result_S_941_AD_MBC_str = [Node941_MBC_str(1,7) Node941_MBC_str(1,8) Node941_MBC_str(1,9);Node941_MBC_str(1,10) Node941_MBC_str(1,11) Node941_MBC_str(1,12);...
    Node941_MBC_str(1,13) Node941_MBC_str(1,14) Node941_MBC_str(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_941_AD_MBC_str_trans = R'*Result_S_941_AD_MBC_str*R;


Result_S_943_AD_MBC_str = [Node943_MBC_str(1,7) Node943_MBC_str(1,8) Node943_MBC_str(1,9);Node943_MBC_str(1,10) Node943_MBC_str(1,11) Node943_MBC_str(1,12);...
    Node943_MBC_str(1,13) Node943_MBC_str(1,14) Node943_MBC_str(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_943_AD_MBC_str_trans = R'*Result_S_943_AD_MBC_str*R;


Result_S_945_AD_MBC_str = [Node945_MBC_str(1,7) Node945_MBC_str(1,8) Node945_MBC_str(1,9);Node945_MBC_str(1,10) Node945_MBC_str(1,11) Node945_MBC_str(1,12);...
    Node945_MBC_str(1,13) Node945_MBC_str(1,14) Node945_MBC_str(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_945_AD_MBC_str_trans = R'*Result_S_945_AD_MBC_str*R;


Result_S_947_AD_MBC_str = [Node947_MBC_str(1,7) Node947_MBC_str(1,8) Node947_MBC_str(1,9);Node947_MBC_str(1,10) Node947_MBC_str(1,11) Node947_MBC_str(1,12);...
    Node947_MBC_str(1,13) Node947_MBC_str(1,14) Node947_MBC_str(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_947_AD_MBC_str_trans = R'*Result_S_947_AD_MBC_str*R;


Result_S_949_AD_MBC_str = [Node949_MBC_str(1,7) Node949_MBC_str(1,8) Node949_MBC_str(1,9);Node949_MBC_str(1,10) Node949_MBC_str(1,11) Node949_MBC_str(1,12);...
    Node949_MBC_str(1,13) Node949_MBC_str(1,14) Node949_MBC_str(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_949_AD_MBC_str_trans = R'*Result_S_949_AD_MBC_str*R;


Result_S_951_AD_MBC_str = [Node951_MBC_str(1,7) Node951_MBC_str(1,8) Node951_MBC_str(1,9);Node951_MBC_str(1,10) Node951_MBC_str(1,11) Node951_MBC_str(1,12);...
    Node951_MBC_str(1,13) Node951_MBC_str(1,14) Node951_MBC_str(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_951_AD_MBC_str_trans = R'*Result_S_951_AD_MBC_str*R;


Result_S_953_AD_MBC_str = [Node953_MBC_str(1,7) Node953_MBC_str(1,8) Node953_MBC_str(1,9);Node953_MBC_str(1,10) Node953_MBC_str(1,11) Node953_MBC_str(1,12);...
    Node953_MBC_str(1,13) Node953_MBC_str(1,14) Node953_MBC_str(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_953_AD_MBC_str_trans = R'*Result_S_953_AD_MBC_str*R;


Result_S_955_AD_MBC_str = [Node955_MBC_str(1,7) Node955_MBC_str(1,8) Node955_MBC_str(1,9);Node955_MBC_str(1,10) Node955_MBC_str(1,11) Node955_MBC_str(1,12);...
    Node955_MBC_str(1,13) Node955_MBC_str(1,14) Node955_MBC_str(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_955_AD_MBC_str_trans = R'*Result_S_955_AD_MBC_str*R;


Result_S_957_AD_MBC_str = [Node957_MBC_str(1,7) Node957_MBC_str(1,8) Node957_MBC_str(1,9);Node957_MBC_str(1,10) Node957_MBC_str(1,11) Node957_MBC_str(1,12);...
    Node957_MBC_str(1,13) Node957_MBC_str(1,14) Node957_MBC_str(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_957_AD_MBC_str_trans = R'*Result_S_957_AD_MBC_str*R;


Result_S_959_AD_MBC_str = [Node959_MBC_str(1,7) Node959_MBC_str(1,8) Node959_MBC_str(1,9);Node959_MBC_str(1,10) Node959_MBC_str(1,11) Node959_MBC_str(1,12);...
    Node959_MBC_str(1,13) Node959_MBC_str(1,14) Node959_MBC_str(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_959_AD_MBC_str_trans = R'*Result_S_959_AD_MBC_str*R;


Result_S_961_AD_MBC_str = [Node961_MBC_str(1,7) Node961_MBC_str(1,8) Node961_MBC_str(1,9);Node961_MBC_str(1,10) Node961_MBC_str(1,11) Node961_MBC_str(1,12);...
    Node961_MBC_str(1,13) Node961_MBC_str(1,14) Node961_MBC_str(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_961_AD_MBC_str_trans = R'*Result_S_961_AD_MBC_str*R;

Result_S_963_AD_MBC_str = [Node963_MBC_str(1,7) Node963_MBC_str(1,8) Node963_MBC_str(1,9);Node963_MBC_str(1,10) Node963_MBC_str(1,11) Node963_MBC_str(1,12);...
    Node963_MBC_str(1,13) Node963_MBC_str(1,14) Node963_MBC_str(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_963_AD_MBC_str_trans = R'*Result_S_963_AD_MBC_str*R;


Result_S_965_AD_MBC_str = [Node965_MBC_str(1,7) Node965_MBC_str(1,8) Node965_MBC_str(1,9);Node965_MBC_str(1,10) Node965_MBC_str(1,11) Node965_MBC_str(1,12);...
    Node965_MBC_str(1,13) Node965_MBC_str(1,14) Node965_MBC_str(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_965_AD_MBC_str_trans = R'*Result_S_965_AD_MBC_str*R;

Result_S_967_AD_MBC_str = [Node967_MBC_str(1,7) Node967_MBC_str(1,8) Node967_MBC_str(1,9);Node967_MBC_str(1,10) Node967_MBC_str(1,11) Node967_MBC_str(1,12);...
    Node967_MBC_str(1,13) Node967_MBC_str(1,14) Node967_MBC_str(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_967_AD_MBC_str_trans = R'*Result_S_967_AD_MBC_str*R;

Result_S_969_AD_MBC_str = [Node969_MBC_str(1,7) Node969_MBC_str(1,8) Node969_MBC_str(1,9);Node969_MBC_str(1,10) Node969_MBC_str(1,11) Node969_MBC_str(1,12);...
    Node969_MBC_str(1,13) Node969_MBC_str(1,14) Node969_MBC_str(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_969_AD_MBC_str_trans = R'*Result_S_969_AD_MBC_str*R;

Result_S_971_AD_MBC_str = [Node971_MBC_str(1,7) Node971_MBC_str(1,8) Node971_MBC_str(1,9);Node971_MBC_str(1,10) Node971_MBC_str(1,11) Node971_MBC_str(1,12);...
    Node971_MBC_str(1,13) Node971_MBC_str(1,14) Node971_MBC_str(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_971_AD_MBC_str_trans = R'*Result_S_971_AD_MBC_str*R;

Result_S_973_AD_MBC_str = [Node973_MBC_str(1,7) Node973_MBC_str(1,8) Node973_MBC_str(1,9);Node973_MBC_str(1,10) Node973_MBC_str(1,11) Node973_MBC_str(1,12);...
    Node973_MBC_str(1,13) Node973_MBC_str(1,14) Node973_MBC_str(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_973_AD_MBC_str_trans = R'*Result_S_973_AD_MBC_str*R;

Result_S_975_AD_MBC_str = [Node975_MBC_str(1,7) Node975_MBC_str(1,8) Node975_MBC_str(1,9);Node975_MBC_str(1,10) Node975_MBC_str(1,11) Node975_MBC_str(1,12);...
    Node975_MBC_str(1,13) Node975_MBC_str(1,14) Node975_MBC_str(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_975_AD_MBC_str_trans = R'*Result_S_975_AD_MBC_str*R;

Result_S_977_AD_MBC_str = [Node977_MBC_str(1,7) Node977_MBC_str(1,8) Node977_MBC_str(1,9);Node977_MBC_str(1,10) Node977_MBC_str(1,11) Node977_MBC_str(1,12);...
    Node977_MBC_str(1,13) Node977_MBC_str(1,14) Node977_MBC_str(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_977_AD_MBC_str_trans = R'*Result_S_977_AD_MBC_str*R;

Result_S_979_AD_MBC_str = [Node979_MBC_str(1,7) Node979_MBC_str(1,8) Node979_MBC_str(1,9);Node979_MBC_str(1,10) Node979_MBC_str(1,11) Node979_MBC_str(1,12);...
    Node979_MBC_str(1,13) Node979_MBC_str(1,14) Node979_MBC_str(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_979_AD_MBC_str_trans = R'*Result_S_979_AD_MBC_str*R;

Result_S_981_AD_MBC_str = [Node981_MBC_str(1,7) Node981_MBC_str(1,8) Node981_MBC_str(1,9);Node981_MBC_str(1,10) Node981_MBC_str(1,11) Node981_MBC_str(1,12);...
    Node981_MBC_str(1,13) Node981_MBC_str(1,14) Node981_MBC_str(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_981_AD_MBC_str_trans = R'*Result_S_981_AD_MBC_str*R;

Result_S_983_AD_MBC_str = [Node983_MBC_str(1,7) Node983_MBC_str(1,8) Node983_MBC_str(1,9);Node983_MBC_str(1,10) Node983_MBC_str(1,11) Node983_MBC_str(1,12);...
    Node983_MBC_str(1,13) Node983_MBC_str(1,14) Node983_MBC_str(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_983_AD_MBC_str_trans = R'*Result_S_983_AD_MBC_str*R;

Result_S_985_AD_MBC_str = [Node985_MBC_str(1,7) Node985_MBC_str(1,8) Node985_MBC_str(1,9);Node985_MBC_str(1,10) Node985_MBC_str(1,11) Node985_MBC_str(1,12);...
    Node985_MBC_str(1,13) Node985_MBC_str(1,14) Node985_MBC_str(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_985_AD_MBC_str_trans = R'*Result_S_985_AD_MBC_str*R;

Result_S_987_AD_MBC_str = [Node987_MBC_str(1,7) Node987_MBC_str(1,8) Node987_MBC_str(1,9);Node987_MBC_str(1,10) Node987_MBC_str(1,11) Node987_MBC_str(1,12);...
    Node987_MBC_str(1,13) Node987_MBC_str(1,14) Node987_MBC_str(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_987_AD_MBC_str_trans = R'*Result_S_987_AD_MBC_str*R;

Result_S_989_AD_MBC_str = [Node989_MBC_str(1,7) Node989_MBC_str(1,8) Node989_MBC_str(1,9);Node989_MBC_str(1,10) Node989_MBC_str(1,11) Node989_MBC_str(1,12);...
    Node989_MBC_str(1,13) Node989_MBC_str(1,14) Node989_MBC_str(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_989_AD_MBC_str_trans = R'*Result_S_989_AD_MBC_str*R;

Result_S_991_AD_MBC_str = [Node991_MBC_str(1,7) Node991_MBC_str(1,8) Node991_MBC_str(1,9);Node991_MBC_str(1,10) Node991_MBC_str(1,11) Node991_MBC_str(1,12);...
    Node991_MBC_str(1,13) Node991_MBC_str(1,14) Node991_MBC_str(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_991_AD_MBC_str_trans = R'*Result_S_991_AD_MBC_str*R;

Result_S_993_AD_MBC_str = [Node993_MBC_str(1,7) Node993_MBC_str(1,8) Node993_MBC_str(1,9);Node993_MBC_str(1,10) Node993_MBC_str(1,11) Node993_MBC_str(1,12);...
    Node993_MBC_str(1,13) Node993_MBC_str(1,14) Node993_MBC_str(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_993_AD_MBC_str_trans = R'*Result_S_993_AD_MBC_str*R;

Result_S_995_AD_MBC_str = [Node995_MBC_str(1,7) Node995_MBC_str(1,8) Node995_MBC_str(1,9);Node995_MBC_str(1,10) Node995_MBC_str(1,11) Node995_MBC_str(1,12);...
    Node995_MBC_str(1,13) Node995_MBC_str(1,14) Node995_MBC_str(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_995_AD_MBC_str_trans = R'*Result_S_995_AD_MBC_str*R;

Result_S_997_AD_MBC_str = [Node997_MBC_str(1,7) Node997_MBC_str(1,8) Node997_MBC_str(1,9);Node997_MBC_str(1,10) Node997_MBC_str(1,11) Node997_MBC_str(1,12);...
    Node997_MBC_str(1,13) Node997_MBC_str(1,14) Node997_MBC_str(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_997_AD_MBC_str_trans = R'*Result_S_997_AD_MBC_str*R;

Result_S_999_AD_MBC_str = [Node999_MBC_str(1,7) Node999_MBC_str(1,8) Node999_MBC_str(1,9);Node999_MBC_str(1,10) Node999_MBC_str(1,11) Node999_MBC_str(1,12);...
    Node999_MBC_str(1,13) Node999_MBC_str(1,14) Node999_MBC_str(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_999_AD_MBC_str_trans = R'*Result_S_999_AD_MBC_str*R;

Result_S_1001_AD_MBC_str = [Node1001_MBC_str(1,7) Node1001_MBC_str(1,8) Node1001_MBC_str(1,9);Node1001_MBC_str(1,10) Node1001_MBC_str(1,11) Node1001_MBC_str(1,12);...
    Node1001_MBC_str(1,13) Node1001_MBC_str(1,14) Node1001_MBC_str(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_1001_AD_MBC_str_trans = R'*Result_S_1001_AD_MBC_str*R;

Result_S_1003_AD_MBC_str = [Node1003_MBC_str(1,7) Node1003_MBC_str(1,8) Node1003_MBC_str(1,9);Node1003_MBC_str(1,10) Node1003_MBC_str(1,11) Node1003_MBC_str(1,12);...
    Node1003_MBC_str(1,13) Node1003_MBC_str(1,14) Node1003_MBC_str(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_1003_AD_MBC_str_trans = R'*Result_S_1003_AD_MBC_str*R;

Result_S_1005_AD_MBC_str = [Node1005_MBC_str(1,7) Node1005_MBC_str(1,8) Node1005_MBC_str(1,9);Node1005_MBC_str(1,10) Node1005_MBC_str(1,11) Node1005_MBC_str(1,12);...
    Node1005_MBC_str(1,13) Node1005_MBC_str(1,14) Node1005_MBC_str(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_1005_AD_MBC_str_trans = R'*Result_S_1005_AD_MBC_str*R;

Result_S_1007_AD_MBC_str = [Node1007_MBC_str(1,7) Node1007_MBC_str(1,8) Node1007_MBC_str(1,9);Node1007_MBC_str(1,10) Node1007_MBC_str(1,11) Node1007_MBC_str(1,12);...
    Node1007_MBC_str(1,13) Node1007_MBC_str(1,14) Node1007_MBC_str(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_1007_AD_MBC_str_trans = R'*Result_S_1007_AD_MBC_str*R;


Result_S_r_teta_MBC_str = [r_938/0.05 Result_S_938_AD_MBC_str_trans(1,2)/P_applied;r_939/0.05  Result_S_939_AD_MBC_str_trans(1,2)/P_applied;...
    r_941/0.05  Result_S_941_AD_MBC_str_trans(1,2)/P_applied;r_943/0.05  Result_S_943_AD_MBC_str_trans(1,2)/P_applied;r_945/0.05  Result_S_945_AD_MBC_str_trans(1,2)/P_applied;...
    r_947/0.05  Result_S_947_AD_MBC_str_trans(1,2)/P_applied;r_949/0.05  Result_S_949_AD_MBC_str_trans(1,2)/P_applied;r_951/0.05  Result_S_951_AD_MBC_str_trans(1,2)/P_applied;...
    r_953/0.05  Result_S_953_AD_MBC_str_trans(1,2)/P_applied;r_955/0.05  Result_S_955_AD_MBC_str_trans(1,2)/P_applied;r_957/0.05  Result_S_957_AD_MBC_str_trans(1,2)/P_applied;...
    r_959/0.05  Result_S_959_AD_MBC_str_trans(1,2)/P_applied;r_961/0.05  Result_S_961_AD_MBC_str_trans(1,2)/P_applied;r_963/0.05  Result_S_963_AD_MBC_str_trans(1,2)/P_applied;...
    r_965/0.05  Result_S_965_AD_MBC_str_trans(1,2)/P_applied;r_967/0.05  Result_S_967_AD_MBC_str_trans(1,2)/P_applied;r_969/0.05  Result_S_969_AD_MBC_str_trans(1,2)/P_applied;...
    r_971/0.05  Result_S_971_AD_MBC_str_trans(1,2)/P_applied;r_973/0.05  Result_S_973_AD_MBC_str_trans(1,2)/P_applied;r_975/0.05  Result_S_975_AD_MBC_str_trans(1,2)/P_applied;...
    r_977/0.05  Result_S_977_AD_MBC_str_trans(1,2)/P_applied;r_979/0.05  Result_S_979_AD_MBC_str_trans(1,2)/P_applied;r_981/0.05  Result_S_981_AD_MBC_str_trans(1,2)/P_applied;...
    r_983/0.05  Result_S_983_AD_MBC_str_trans(1,2)/P_applied;r_985/0.05  Result_S_985_AD_MBC_str_trans(1,2)/P_applied;r_987/0.05  Result_S_987_AD_MBC_str_trans(1,2)/P_applied;...
    r_989/0.05  Result_S_989_AD_MBC_str_trans(1,2)/P_applied;r_991/0.05  Result_S_991_AD_MBC_str_trans(1,2)/P_applied;r_993/0.05  Result_S_993_AD_MBC_str_trans(1,2)/P_applied;...
    r_995/0.05  Result_S_995_AD_MBC_str_trans(1,2)/P_applied;r_997/0.05  Result_S_997_AD_MBC_str_trans(1,2)/P_applied;r_999/0.05  Result_S_999_AD_MBC_str_trans(1,2)/P_applied;...
    r_1001/0.05  Result_S_1001_AD_MBC_str_trans(1,2)/P_applied;r_1003/0.05  Result_S_1003_AD_MBC_str_trans(1,2)/P_applied;r_1005/0.05  Result_S_1005_AD_MBC_str_trans(1,2)/P_applied;...
    r_1007/0.05  Result_S_1007_AD_MBC_str_trans(1,2)/P_applied];



Result_S_938_AD_NMBC_com = [Node938_NMBC_com(1,7) Node938_NMBC_com(1,8) Node938_NMBC_com(1,9);Node938_NMBC_com(1,10) Node938_NMBC_com(1,11) Node938_NMBC_com(1,12);...
    Node938_NMBC_com(1,13) Node938_NMBC_com(1,14) Node938_NMBC_com(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_938_AD_NMBC_com_trans = R'*Result_S_938_AD_NMBC_com*R;


Result_S_939_AD_NMBC_com = [Node939_NMBC_com(1,7) Node939_NMBC_com(1,8) Node939_NMBC_com(1,9);Node939_NMBC_com(1,10) Node939_NMBC_com(1,11) Node939_NMBC_com(1,12);...
    Node939_NMBC_com(1,13) Node939_NMBC_com(1,14) Node939_NMBC_com(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_939_AD_NMBC_com_trans = R'*Result_S_939_AD_NMBC_com*R;


Result_S_941_AD_NMBC_com = [Node941_NMBC_com(1,7) Node941_NMBC_com(1,8) Node941_NMBC_com(1,9);Node941_NMBC_com(1,10) Node941_NMBC_com(1,11) Node941_NMBC_com(1,12);...
    Node941_NMBC_com(1,13) Node941_NMBC_com(1,14) Node941_NMBC_com(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_941_AD_NMBC_com_trans = R'*Result_S_941_AD_NMBC_com*R;


Result_S_943_AD_NMBC_com = [Node943_NMBC_com(1,7) Node943_NMBC_com(1,8) Node943_NMBC_com(1,9);Node943_NMBC_com(1,10) Node943_NMBC_com(1,11) Node943_NMBC_com(1,12);...
    Node943_NMBC_com(1,13) Node943_NMBC_com(1,14) Node943_NMBC_com(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_943_AD_NMBC_com_trans = R'*Result_S_943_AD_NMBC_com*R;


Result_S_945_AD_NMBC_com = [Node945_NMBC_com(1,7) Node945_NMBC_com(1,8) Node945_NMBC_com(1,9);Node945_NMBC_com(1,10) Node945_NMBC_com(1,11) Node945_NMBC_com(1,12);...
    Node945_NMBC_com(1,13) Node945_NMBC_com(1,14) Node945_NMBC_com(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_945_AD_NMBC_com_trans = R'*Result_S_945_AD_NMBC_com*R;


Result_S_947_AD_NMBC_com = [Node947_NMBC_com(1,7) Node947_NMBC_com(1,8) Node947_NMBC_com(1,9);Node947_NMBC_com(1,10) Node947_NMBC_com(1,11) Node947_NMBC_com(1,12);...
    Node947_NMBC_com(1,13) Node947_NMBC_com(1,14) Node947_NMBC_com(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_947_AD_NMBC_com_trans = R'*Result_S_947_AD_NMBC_com*R;


Result_S_949_AD_NMBC_com = [Node949_NMBC_com(1,7) Node949_NMBC_com(1,8) Node949_NMBC_com(1,9);Node949_NMBC_com(1,10) Node949_NMBC_com(1,11) Node949_NMBC_com(1,12);...
    Node949_NMBC_com(1,13) Node949_NMBC_com(1,14) Node949_NMBC_com(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_949_AD_NMBC_com_trans = R'*Result_S_949_AD_NMBC_com*R;


Result_S_951_AD_NMBC_com = [Node951_NMBC_com(1,7) Node951_NMBC_com(1,8) Node951_NMBC_com(1,9);Node951_NMBC_com(1,10) Node951_NMBC_com(1,11) Node951_NMBC_com(1,12);...
    Node951_NMBC_com(1,13) Node951_NMBC_com(1,14) Node951_NMBC_com(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_951_AD_NMBC_com_trans = R'*Result_S_951_AD_NMBC_com*R;


Result_S_953_AD_NMBC_com = [Node953_NMBC_com(1,7) Node953_NMBC_com(1,8) Node953_NMBC_com(1,9);Node953_NMBC_com(1,10) Node953_NMBC_com(1,11) Node953_NMBC_com(1,12);...
    Node953_NMBC_com(1,13) Node953_NMBC_com(1,14) Node953_NMBC_com(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_953_AD_NMBC_com_trans = R'*Result_S_953_AD_NMBC_com*R;


Result_S_955_AD_NMBC_com = [Node955_NMBC_com(1,7) Node955_NMBC_com(1,8) Node955_NMBC_com(1,9);Node955_NMBC_com(1,10) Node955_NMBC_com(1,11) Node955_NMBC_com(1,12);...
    Node955_NMBC_com(1,13) Node955_NMBC_com(1,14) Node955_NMBC_com(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_955_AD_NMBC_com_trans = R'*Result_S_955_AD_NMBC_com*R;


Result_S_957_AD_NMBC_com = [Node957_NMBC_com(1,7) Node957_NMBC_com(1,8) Node957_NMBC_com(1,9);Node957_NMBC_com(1,10) Node957_NMBC_com(1,11) Node957_NMBC_com(1,12);...
    Node957_NMBC_com(1,13) Node957_NMBC_com(1,14) Node957_NMBC_com(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_957_AD_NMBC_com_trans = R'*Result_S_957_AD_NMBC_com*R;


Result_S_959_AD_NMBC_com = [Node959_NMBC_com(1,7) Node959_NMBC_com(1,8) Node959_NMBC_com(1,9);Node959_NMBC_com(1,10) Node959_NMBC_com(1,11) Node959_NMBC_com(1,12);...
    Node959_NMBC_com(1,13) Node959_NMBC_com(1,14) Node959_NMBC_com(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_959_AD_NMBC_com_trans = R'*Result_S_959_AD_NMBC_com*R;


Result_S_961_AD_NMBC_com = [Node961_NMBC_com(1,7) Node961_NMBC_com(1,8) Node961_NMBC_com(1,9);Node961_NMBC_com(1,10) Node961_NMBC_com(1,11) Node961_NMBC_com(1,12);...
    Node961_NMBC_com(1,13) Node961_NMBC_com(1,14) Node961_NMBC_com(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_961_AD_NMBC_com_trans = R'*Result_S_961_AD_NMBC_com*R;

Result_S_963_AD_NMBC_com = [Node963_NMBC_com(1,7) Node963_NMBC_com(1,8) Node963_NMBC_com(1,9);Node963_NMBC_com(1,10) Node963_NMBC_com(1,11) Node963_NMBC_com(1,12);...
    Node963_NMBC_com(1,13) Node963_NMBC_com(1,14) Node963_NMBC_com(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_963_AD_NMBC_com_trans = R'*Result_S_963_AD_NMBC_com*R;


Result_S_965_AD_NMBC_com = [Node965_NMBC_com(1,7) Node965_NMBC_com(1,8) Node965_NMBC_com(1,9);Node965_NMBC_com(1,10) Node965_NMBC_com(1,11) Node965_NMBC_com(1,12);...
    Node965_NMBC_com(1,13) Node965_NMBC_com(1,14) Node965_NMBC_com(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_965_AD_NMBC_com_trans = R'*Result_S_965_AD_NMBC_com*R;

Result_S_967_AD_NMBC_com = [Node967_NMBC_com(1,7) Node967_NMBC_com(1,8) Node967_NMBC_com(1,9);Node967_NMBC_com(1,10) Node967_NMBC_com(1,11) Node967_NMBC_com(1,12);...
    Node967_NMBC_com(1,13) Node967_NMBC_com(1,14) Node967_NMBC_com(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_967_AD_NMBC_com_trans = R'*Result_S_967_AD_NMBC_com*R;

Result_S_969_AD_NMBC_com = [Node969_NMBC_com(1,7) Node969_NMBC_com(1,8) Node969_NMBC_com(1,9);Node969_NMBC_com(1,10) Node969_NMBC_com(1,11) Node969_NMBC_com(1,12);...
    Node969_NMBC_com(1,13) Node969_NMBC_com(1,14) Node969_NMBC_com(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_969_AD_NMBC_com_trans = R'*Result_S_969_AD_NMBC_com*R;

Result_S_971_AD_NMBC_com = [Node971_NMBC_com(1,7) Node971_NMBC_com(1,8) Node971_NMBC_com(1,9);Node971_NMBC_com(1,10) Node971_NMBC_com(1,11) Node971_NMBC_com(1,12);...
    Node971_NMBC_com(1,13) Node971_NMBC_com(1,14) Node971_NMBC_com(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_971_AD_NMBC_com_trans = R'*Result_S_971_AD_NMBC_com*R;

Result_S_973_AD_NMBC_com = [Node973_NMBC_com(1,7) Node973_NMBC_com(1,8) Node973_NMBC_com(1,9);Node973_NMBC_com(1,10) Node973_NMBC_com(1,11) Node973_NMBC_com(1,12);...
    Node973_NMBC_com(1,13) Node973_NMBC_com(1,14) Node973_NMBC_com(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_973_AD_NMBC_com_trans = R'*Result_S_973_AD_NMBC_com*R;

Result_S_975_AD_NMBC_com = [Node975_NMBC_com(1,7) Node975_NMBC_com(1,8) Node975_NMBC_com(1,9);Node975_NMBC_com(1,10) Node975_NMBC_com(1,11) Node975_NMBC_com(1,12);...
    Node975_NMBC_com(1,13) Node975_NMBC_com(1,14) Node975_NMBC_com(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_975_AD_NMBC_com_trans = R'*Result_S_975_AD_NMBC_com*R;

Result_S_977_AD_NMBC_com = [Node977_NMBC_com(1,7) Node977_NMBC_com(1,8) Node977_NMBC_com(1,9);Node977_NMBC_com(1,10) Node977_NMBC_com(1,11) Node977_NMBC_com(1,12);...
    Node977_NMBC_com(1,13) Node977_NMBC_com(1,14) Node977_NMBC_com(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_977_AD_NMBC_com_trans = R'*Result_S_977_AD_NMBC_com*R;

Result_S_979_AD_NMBC_com = [Node979_NMBC_com(1,7) Node979_NMBC_com(1,8) Node979_NMBC_com(1,9);Node979_NMBC_com(1,10) Node979_NMBC_com(1,11) Node979_NMBC_com(1,12);...
    Node979_NMBC_com(1,13) Node979_NMBC_com(1,14) Node979_NMBC_com(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_979_AD_NMBC_com_trans = R'*Result_S_979_AD_NMBC_com*R;

Result_S_981_AD_NMBC_com = [Node981_NMBC_com(1,7) Node981_NMBC_com(1,8) Node981_NMBC_com(1,9);Node981_NMBC_com(1,10) Node981_NMBC_com(1,11) Node981_NMBC_com(1,12);...
    Node981_NMBC_com(1,13) Node981_NMBC_com(1,14) Node981_NMBC_com(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_981_AD_NMBC_com_trans = R'*Result_S_981_AD_NMBC_com*R;

Result_S_983_AD_NMBC_com = [Node983_NMBC_com(1,7) Node983_NMBC_com(1,8) Node983_NMBC_com(1,9);Node983_NMBC_com(1,10) Node983_NMBC_com(1,11) Node983_NMBC_com(1,12);...
    Node983_NMBC_com(1,13) Node983_NMBC_com(1,14) Node983_NMBC_com(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_983_AD_NMBC_com_trans = R'*Result_S_983_AD_NMBC_com*R;

Result_S_985_AD_NMBC_com = [Node985_NMBC_com(1,7) Node985_NMBC_com(1,8) Node985_NMBC_com(1,9);Node985_NMBC_com(1,10) Node985_NMBC_com(1,11) Node985_NMBC_com(1,12);...
    Node985_NMBC_com(1,13) Node985_NMBC_com(1,14) Node985_NMBC_com(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_985_AD_NMBC_com_trans = R'*Result_S_985_AD_NMBC_com*R;

Result_S_987_AD_NMBC_com = [Node987_NMBC_com(1,7) Node987_NMBC_com(1,8) Node987_NMBC_com(1,9);Node987_NMBC_com(1,10) Node987_NMBC_com(1,11) Node987_NMBC_com(1,12);...
    Node987_NMBC_com(1,13) Node987_NMBC_com(1,14) Node987_NMBC_com(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_987_AD_NMBC_com_trans = R'*Result_S_987_AD_NMBC_com*R;

Result_S_989_AD_NMBC_com = [Node989_NMBC_com(1,7) Node989_NMBC_com(1,8) Node989_NMBC_com(1,9);Node989_NMBC_com(1,10) Node989_NMBC_com(1,11) Node989_NMBC_com(1,12);...
    Node989_NMBC_com(1,13) Node989_NMBC_com(1,14) Node989_NMBC_com(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_989_AD_NMBC_com_trans = R'*Result_S_989_AD_NMBC_com*R;

Result_S_991_AD_NMBC_com = [Node991_NMBC_com(1,7) Node991_NMBC_com(1,8) Node991_NMBC_com(1,9);Node991_NMBC_com(1,10) Node991_NMBC_com(1,11) Node991_NMBC_com(1,12);...
    Node991_NMBC_com(1,13) Node991_NMBC_com(1,14) Node991_NMBC_com(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_991_AD_NMBC_com_trans = R'*Result_S_991_AD_NMBC_com*R;

Result_S_993_AD_NMBC_com = [Node993_NMBC_com(1,7) Node993_NMBC_com(1,8) Node993_NMBC_com(1,9);Node993_NMBC_com(1,10) Node993_NMBC_com(1,11) Node993_NMBC_com(1,12);...
    Node993_NMBC_com(1,13) Node993_NMBC_com(1,14) Node993_NMBC_com(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_993_AD_NMBC_com_trans = R'*Result_S_993_AD_NMBC_com*R;

Result_S_995_AD_NMBC_com = [Node995_NMBC_com(1,7) Node995_NMBC_com(1,8) Node995_NMBC_com(1,9);Node995_NMBC_com(1,10) Node995_NMBC_com(1,11) Node995_NMBC_com(1,12);...
    Node995_NMBC_com(1,13) Node995_NMBC_com(1,14) Node995_NMBC_com(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_995_AD_NMBC_com_trans = R'*Result_S_995_AD_NMBC_com*R;

Result_S_997_AD_NMBC_com = [Node997_NMBC_com(1,7) Node997_NMBC_com(1,8) Node997_NMBC_com(1,9);Node997_NMBC_com(1,10) Node997_NMBC_com(1,11) Node997_NMBC_com(1,12);...
    Node997_NMBC_com(1,13) Node997_NMBC_com(1,14) Node997_NMBC_com(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_997_AD_NMBC_com_trans = R'*Result_S_997_AD_NMBC_com*R;

Result_S_999_AD_NMBC_com = [Node999_NMBC_com(1,7) Node999_NMBC_com(1,8) Node999_NMBC_com(1,9);Node999_NMBC_com(1,10) Node999_NMBC_com(1,11) Node999_NMBC_com(1,12);...
    Node999_NMBC_com(1,13) Node999_NMBC_com(1,14) Node999_NMBC_com(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_999_AD_NMBC_com_trans = R'*Result_S_999_AD_NMBC_com*R;

Result_S_1001_AD_NMBC_com = [Node1001_NMBC_com(1,7) Node1001_NMBC_com(1,8) Node1001_NMBC_com(1,9);Node1001_NMBC_com(1,10) Node1001_NMBC_com(1,11) Node1001_NMBC_com(1,12);...
    Node1001_NMBC_com(1,13) Node1001_NMBC_com(1,14) Node1001_NMBC_com(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_1001_AD_NMBC_com_trans = R'*Result_S_1001_AD_NMBC_com*R;

Result_S_1003_AD_NMBC_com = [Node1003_NMBC_com(1,7) Node1003_NMBC_com(1,8) Node1003_NMBC_com(1,9);Node1003_NMBC_com(1,10) Node1003_NMBC_com(1,11) Node1003_NMBC_com(1,12);...
    Node1003_NMBC_com(1,13) Node1003_NMBC_com(1,14) Node1003_NMBC_com(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_1003_AD_NMBC_com_trans = R'*Result_S_1003_AD_NMBC_com*R;

Result_S_1005_AD_NMBC_com = [Node1005_NMBC_com(1,7) Node1005_NMBC_com(1,8) Node1005_NMBC_com(1,9);Node1005_NMBC_com(1,10) Node1005_NMBC_com(1,11) Node1005_NMBC_com(1,12);...
    Node1005_NMBC_com(1,13) Node1005_NMBC_com(1,14) Node1005_NMBC_com(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_1005_AD_NMBC_com_trans = R'*Result_S_1005_AD_NMBC_com*R;

Result_S_1007_AD_NMBC_com = [Node1007_NMBC_com(1,7) Node1007_NMBC_com(1,8) Node1007_NMBC_com(1,9);Node1007_NMBC_com(1,10) Node1007_NMBC_com(1,11) Node1007_NMBC_com(1,12);...
    Node1007_NMBC_com(1,13) Node1007_NMBC_com(1,14) Node1007_NMBC_com(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_1007_AD_NMBC_com_trans = R'*Result_S_1007_AD_NMBC_com*R;


Result_S_r_teta_NMBC_com = [r_938/0.05 Result_S_938_AD_NMBC_com_trans(1,2)/P_applied;r_939/0.05  Result_S_939_AD_NMBC_com_trans(1,2)/P_applied;...
    r_941/0.05  Result_S_941_AD_NMBC_com_trans(1,2)/P_applied;r_943/0.05  Result_S_943_AD_NMBC_com_trans(1,2)/P_applied;r_945/0.05  Result_S_945_AD_NMBC_com_trans(1,2)/P_applied;...
    r_947/0.05  Result_S_947_AD_NMBC_com_trans(1,2)/P_applied;r_949/0.05  Result_S_949_AD_NMBC_com_trans(1,2)/P_applied;r_951/0.05  Result_S_951_AD_NMBC_com_trans(1,2)/P_applied;...
    r_953/0.05  Result_S_953_AD_NMBC_com_trans(1,2)/P_applied;r_955/0.05  Result_S_955_AD_NMBC_com_trans(1,2)/P_applied;r_957/0.05  Result_S_957_AD_NMBC_com_trans(1,2)/P_applied;...
    r_959/0.05  Result_S_959_AD_NMBC_com_trans(1,2)/P_applied;r_961/0.05  Result_S_961_AD_NMBC_com_trans(1,2)/P_applied;r_963/0.05  Result_S_963_AD_NMBC_com_trans(1,2)/P_applied;...
    r_965/0.05  Result_S_965_AD_NMBC_com_trans(1,2)/P_applied;r_967/0.05  Result_S_967_AD_NMBC_com_trans(1,2)/P_applied;r_969/0.05  Result_S_969_AD_NMBC_com_trans(1,2)/P_applied;...
    r_971/0.05  Result_S_971_AD_NMBC_com_trans(1,2)/P_applied;r_973/0.05  Result_S_973_AD_NMBC_com_trans(1,2)/P_applied;r_975/0.05  Result_S_975_AD_NMBC_com_trans(1,2)/P_applied;...
    r_977/0.05  Result_S_977_AD_NMBC_com_trans(1,2)/P_applied;r_979/0.05  Result_S_979_AD_NMBC_com_trans(1,2)/P_applied;r_981/0.05  Result_S_981_AD_NMBC_com_trans(1,2)/P_applied;...
    r_983/0.05  Result_S_983_AD_NMBC_com_trans(1,2)/P_applied;r_985/0.05  Result_S_985_AD_NMBC_com_trans(1,2)/P_applied;r_987/0.05  Result_S_987_AD_NMBC_com_trans(1,2)/P_applied;...
    r_989/0.05  Result_S_989_AD_NMBC_com_trans(1,2)/P_applied;r_991/0.05  Result_S_991_AD_NMBC_com_trans(1,2)/P_applied;r_993/0.05  Result_S_993_AD_NMBC_com_trans(1,2)/P_applied;...
    r_995/0.05  Result_S_995_AD_NMBC_com_trans(1,2)/P_applied;r_997/0.05  Result_S_997_AD_NMBC_com_trans(1,2)/P_applied;r_999/0.05  Result_S_999_AD_NMBC_com_trans(1,2)/P_applied;...
    r_1001/0.05  Result_S_1001_AD_NMBC_com_trans(1,2)/P_applied;r_1003/0.05  Result_S_1003_AD_NMBC_com_trans(1,2)/P_applied;r_1005/0.05  Result_S_1005_AD_NMBC_com_trans(1,2)/P_applied;...
    r_1007/0.05  Result_S_1007_AD_NMBC_com_trans(1,2)/P_applied];


Result_S_938_AD_NMBC_rot = [Node938_NMBC_rot(1,7) Node938_NMBC_rot(1,8) Node938_NMBC_rot(1,9);Node938_NMBC_rot(1,10) Node938_NMBC_rot(1,11) Node938_NMBC_rot(1,12);...
    Node938_NMBC_rot(1,13) Node938_NMBC_rot(1,14) Node938_NMBC_rot(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_938_AD_NMBC_rot_trans = R'*Result_S_938_AD_NMBC_rot*R;


Result_S_939_AD_NMBC_rot = [Node939_NMBC_rot(1,7) Node939_NMBC_rot(1,8) Node939_NMBC_rot(1,9);Node939_NMBC_rot(1,10) Node939_NMBC_rot(1,11) Node939_NMBC_rot(1,12);...
    Node939_NMBC_rot(1,13) Node939_NMBC_rot(1,14) Node939_NMBC_rot(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_939_AD_NMBC_rot_trans = R'*Result_S_939_AD_NMBC_rot*R;


Result_S_941_AD_NMBC_rot = [Node941_NMBC_rot(1,7) Node941_NMBC_rot(1,8) Node941_NMBC_rot(1,9);Node941_NMBC_rot(1,10) Node941_NMBC_rot(1,11) Node941_NMBC_rot(1,12);...
    Node941_NMBC_rot(1,13) Node941_NMBC_rot(1,14) Node941_NMBC_rot(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_941_AD_NMBC_rot_trans = R'*Result_S_941_AD_NMBC_rot*R;


Result_S_943_AD_NMBC_rot = [Node943_NMBC_rot(1,7) Node943_NMBC_rot(1,8) Node943_NMBC_rot(1,9);Node943_NMBC_rot(1,10) Node943_NMBC_rot(1,11) Node943_NMBC_rot(1,12);...
    Node943_NMBC_rot(1,13) Node943_NMBC_rot(1,14) Node943_NMBC_rot(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_943_AD_NMBC_rot_trans = R'*Result_S_943_AD_NMBC_rot*R;


Result_S_945_AD_NMBC_rot = [Node945_NMBC_rot(1,7) Node945_NMBC_rot(1,8) Node945_NMBC_rot(1,9);Node945_NMBC_rot(1,10) Node945_NMBC_rot(1,11) Node945_NMBC_rot(1,12);...
    Node945_NMBC_rot(1,13) Node945_NMBC_rot(1,14) Node945_NMBC_rot(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_945_AD_NMBC_rot_trans = R'*Result_S_945_AD_NMBC_rot*R;


Result_S_947_AD_NMBC_rot = [Node947_NMBC_rot(1,7) Node947_NMBC_rot(1,8) Node947_NMBC_rot(1,9);Node947_NMBC_rot(1,10) Node947_NMBC_rot(1,11) Node947_NMBC_rot(1,12);...
    Node947_NMBC_rot(1,13) Node947_NMBC_rot(1,14) Node947_NMBC_rot(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_947_AD_NMBC_rot_trans = R'*Result_S_947_AD_NMBC_rot*R;


Result_S_949_AD_NMBC_rot = [Node949_NMBC_rot(1,7) Node949_NMBC_rot(1,8) Node949_NMBC_rot(1,9);Node949_NMBC_rot(1,10) Node949_NMBC_rot(1,11) Node949_NMBC_rot(1,12);...
    Node949_NMBC_rot(1,13) Node949_NMBC_rot(1,14) Node949_NMBC_rot(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_949_AD_NMBC_rot_trans = R'*Result_S_949_AD_NMBC_rot*R;


Result_S_951_AD_NMBC_rot = [Node951_NMBC_rot(1,7) Node951_NMBC_rot(1,8) Node951_NMBC_rot(1,9);Node951_NMBC_rot(1,10) Node951_NMBC_rot(1,11) Node951_NMBC_rot(1,12);...
    Node951_NMBC_rot(1,13) Node951_NMBC_rot(1,14) Node951_NMBC_rot(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_951_AD_NMBC_rot_trans = R'*Result_S_951_AD_NMBC_rot*R;


Result_S_953_AD_NMBC_rot = [Node953_NMBC_rot(1,7) Node953_NMBC_rot(1,8) Node953_NMBC_rot(1,9);Node953_NMBC_rot(1,10) Node953_NMBC_rot(1,11) Node953_NMBC_rot(1,12);...
    Node953_NMBC_rot(1,13) Node953_NMBC_rot(1,14) Node953_NMBC_rot(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_953_AD_NMBC_rot_trans = R'*Result_S_953_AD_NMBC_rot*R;


Result_S_955_AD_NMBC_rot = [Node955_NMBC_rot(1,7) Node955_NMBC_rot(1,8) Node955_NMBC_rot(1,9);Node955_NMBC_rot(1,10) Node955_NMBC_rot(1,11) Node955_NMBC_rot(1,12);...
    Node955_NMBC_rot(1,13) Node955_NMBC_rot(1,14) Node955_NMBC_rot(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_955_AD_NMBC_rot_trans = R'*Result_S_955_AD_NMBC_rot*R;


Result_S_957_AD_NMBC_rot = [Node957_NMBC_rot(1,7) Node957_NMBC_rot(1,8) Node957_NMBC_rot(1,9);Node957_NMBC_rot(1,10) Node957_NMBC_rot(1,11) Node957_NMBC_rot(1,12);...
    Node957_NMBC_rot(1,13) Node957_NMBC_rot(1,14) Node957_NMBC_rot(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_957_AD_NMBC_rot_trans = R'*Result_S_957_AD_NMBC_rot*R;


Result_S_959_AD_NMBC_rot = [Node959_NMBC_rot(1,7) Node959_NMBC_rot(1,8) Node959_NMBC_rot(1,9);Node959_NMBC_rot(1,10) Node959_NMBC_rot(1,11) Node959_NMBC_rot(1,12);...
    Node959_NMBC_rot(1,13) Node959_NMBC_rot(1,14) Node959_NMBC_rot(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_959_AD_NMBC_rot_trans = R'*Result_S_959_AD_NMBC_rot*R;


Result_S_961_AD_NMBC_rot = [Node961_NMBC_rot(1,7) Node961_NMBC_rot(1,8) Node961_NMBC_rot(1,9);Node961_NMBC_rot(1,10) Node961_NMBC_rot(1,11) Node961_NMBC_rot(1,12);...
    Node961_NMBC_rot(1,13) Node961_NMBC_rot(1,14) Node961_NMBC_rot(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_961_AD_NMBC_rot_trans = R'*Result_S_961_AD_NMBC_rot*R;

Result_S_963_AD_NMBC_rot = [Node963_NMBC_rot(1,7) Node963_NMBC_rot(1,8) Node963_NMBC_rot(1,9);Node963_NMBC_rot(1,10) Node963_NMBC_rot(1,11) Node963_NMBC_rot(1,12);...
    Node963_NMBC_rot(1,13) Node963_NMBC_rot(1,14) Node963_NMBC_rot(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_963_AD_NMBC_rot_trans = R'*Result_S_963_AD_NMBC_rot*R;


Result_S_965_AD_NMBC_rot = [Node965_NMBC_rot(1,7) Node965_NMBC_rot(1,8) Node965_NMBC_rot(1,9);Node965_NMBC_rot(1,10) Node965_NMBC_rot(1,11) Node965_NMBC_rot(1,12);...
    Node965_NMBC_rot(1,13) Node965_NMBC_rot(1,14) Node965_NMBC_rot(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_965_AD_NMBC_rot_trans = R'*Result_S_965_AD_NMBC_rot*R;

Result_S_967_AD_NMBC_rot = [Node967_NMBC_rot(1,7) Node967_NMBC_rot(1,8) Node967_NMBC_rot(1,9);Node967_NMBC_rot(1,10) Node967_NMBC_rot(1,11) Node967_NMBC_rot(1,12);...
    Node967_NMBC_rot(1,13) Node967_NMBC_rot(1,14) Node967_NMBC_rot(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_967_AD_NMBC_rot_trans = R'*Result_S_967_AD_NMBC_rot*R;

Result_S_969_AD_NMBC_rot = [Node969_NMBC_rot(1,7) Node969_NMBC_rot(1,8) Node969_NMBC_rot(1,9);Node969_NMBC_rot(1,10) Node969_NMBC_rot(1,11) Node969_NMBC_rot(1,12);...
    Node969_NMBC_rot(1,13) Node969_NMBC_rot(1,14) Node969_NMBC_rot(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_969_AD_NMBC_rot_trans = R'*Result_S_969_AD_NMBC_rot*R;

Result_S_971_AD_NMBC_rot = [Node971_NMBC_rot(1,7) Node971_NMBC_rot(1,8) Node971_NMBC_rot(1,9);Node971_NMBC_rot(1,10) Node971_NMBC_rot(1,11) Node971_NMBC_rot(1,12);...
    Node971_NMBC_rot(1,13) Node971_NMBC_rot(1,14) Node971_NMBC_rot(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_971_AD_NMBC_rot_trans = R'*Result_S_971_AD_NMBC_rot*R;

Result_S_973_AD_NMBC_rot = [Node973_NMBC_rot(1,7) Node973_NMBC_rot(1,8) Node973_NMBC_rot(1,9);Node973_NMBC_rot(1,10) Node973_NMBC_rot(1,11) Node973_NMBC_rot(1,12);...
    Node973_NMBC_rot(1,13) Node973_NMBC_rot(1,14) Node973_NMBC_rot(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_973_AD_NMBC_rot_trans = R'*Result_S_973_AD_NMBC_rot*R;

Result_S_975_AD_NMBC_rot = [Node975_NMBC_rot(1,7) Node975_NMBC_rot(1,8) Node975_NMBC_rot(1,9);Node975_NMBC_rot(1,10) Node975_NMBC_rot(1,11) Node975_NMBC_rot(1,12);...
    Node975_NMBC_rot(1,13) Node975_NMBC_rot(1,14) Node975_NMBC_rot(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_975_AD_NMBC_rot_trans = R'*Result_S_975_AD_NMBC_rot*R;

Result_S_977_AD_NMBC_rot = [Node977_NMBC_rot(1,7) Node977_NMBC_rot(1,8) Node977_NMBC_rot(1,9);Node977_NMBC_rot(1,10) Node977_NMBC_rot(1,11) Node977_NMBC_rot(1,12);...
    Node977_NMBC_rot(1,13) Node977_NMBC_rot(1,14) Node977_NMBC_rot(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_977_AD_NMBC_rot_trans = R'*Result_S_977_AD_NMBC_rot*R;

Result_S_979_AD_NMBC_rot = [Node979_NMBC_rot(1,7) Node979_NMBC_rot(1,8) Node979_NMBC_rot(1,9);Node979_NMBC_rot(1,10) Node979_NMBC_rot(1,11) Node979_NMBC_rot(1,12);...
    Node979_NMBC_rot(1,13) Node979_NMBC_rot(1,14) Node979_NMBC_rot(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_979_AD_NMBC_rot_trans = R'*Result_S_979_AD_NMBC_rot*R;

Result_S_981_AD_NMBC_rot = [Node981_NMBC_rot(1,7) Node981_NMBC_rot(1,8) Node981_NMBC_rot(1,9);Node981_NMBC_rot(1,10) Node981_NMBC_rot(1,11) Node981_NMBC_rot(1,12);...
    Node981_NMBC_rot(1,13) Node981_NMBC_rot(1,14) Node981_NMBC_rot(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_981_AD_NMBC_rot_trans = R'*Result_S_981_AD_NMBC_rot*R;

Result_S_983_AD_NMBC_rot = [Node983_NMBC_rot(1,7) Node983_NMBC_rot(1,8) Node983_NMBC_rot(1,9);Node983_NMBC_rot(1,10) Node983_NMBC_rot(1,11) Node983_NMBC_rot(1,12);...
    Node983_NMBC_rot(1,13) Node983_NMBC_rot(1,14) Node983_NMBC_rot(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_983_AD_NMBC_rot_trans = R'*Result_S_983_AD_NMBC_rot*R;

Result_S_985_AD_NMBC_rot = [Node985_NMBC_rot(1,7) Node985_NMBC_rot(1,8) Node985_NMBC_rot(1,9);Node985_NMBC_rot(1,10) Node985_NMBC_rot(1,11) Node985_NMBC_rot(1,12);...
    Node985_NMBC_rot(1,13) Node985_NMBC_rot(1,14) Node985_NMBC_rot(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_985_AD_NMBC_rot_trans = R'*Result_S_985_AD_NMBC_rot*R;

Result_S_987_AD_NMBC_rot = [Node987_NMBC_rot(1,7) Node987_NMBC_rot(1,8) Node987_NMBC_rot(1,9);Node987_NMBC_rot(1,10) Node987_NMBC_rot(1,11) Node987_NMBC_rot(1,12);...
    Node987_NMBC_rot(1,13) Node987_NMBC_rot(1,14) Node987_NMBC_rot(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_987_AD_NMBC_rot_trans = R'*Result_S_987_AD_NMBC_rot*R;

Result_S_989_AD_NMBC_rot = [Node989_NMBC_rot(1,7) Node989_NMBC_rot(1,8) Node989_NMBC_rot(1,9);Node989_NMBC_rot(1,10) Node989_NMBC_rot(1,11) Node989_NMBC_rot(1,12);...
    Node989_NMBC_rot(1,13) Node989_NMBC_rot(1,14) Node989_NMBC_rot(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_989_AD_NMBC_rot_trans = R'*Result_S_989_AD_NMBC_rot*R;

Result_S_991_AD_NMBC_rot = [Node991_NMBC_rot(1,7) Node991_NMBC_rot(1,8) Node991_NMBC_rot(1,9);Node991_NMBC_rot(1,10) Node991_NMBC_rot(1,11) Node991_NMBC_rot(1,12);...
    Node991_NMBC_rot(1,13) Node991_NMBC_rot(1,14) Node991_NMBC_rot(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_991_AD_NMBC_rot_trans = R'*Result_S_991_AD_NMBC_rot*R;

Result_S_993_AD_NMBC_rot = [Node993_NMBC_rot(1,7) Node993_NMBC_rot(1,8) Node993_NMBC_rot(1,9);Node993_NMBC_rot(1,10) Node993_NMBC_rot(1,11) Node993_NMBC_rot(1,12);...
    Node993_NMBC_rot(1,13) Node993_NMBC_rot(1,14) Node993_NMBC_rot(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_993_AD_NMBC_rot_trans = R'*Result_S_993_AD_NMBC_rot*R;

Result_S_995_AD_NMBC_rot = [Node995_NMBC_rot(1,7) Node995_NMBC_rot(1,8) Node995_NMBC_rot(1,9);Node995_NMBC_rot(1,10) Node995_NMBC_rot(1,11) Node995_NMBC_rot(1,12);...
    Node995_NMBC_rot(1,13) Node995_NMBC_rot(1,14) Node995_NMBC_rot(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_995_AD_NMBC_rot_trans = R'*Result_S_995_AD_NMBC_rot*R;

Result_S_997_AD_NMBC_rot = [Node997_NMBC_rot(1,7) Node997_NMBC_rot(1,8) Node997_NMBC_rot(1,9);Node997_NMBC_rot(1,10) Node997_NMBC_rot(1,11) Node997_NMBC_rot(1,12);...
    Node997_NMBC_rot(1,13) Node997_NMBC_rot(1,14) Node997_NMBC_rot(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_997_AD_NMBC_rot_trans = R'*Result_S_997_AD_NMBC_rot*R;

Result_S_999_AD_NMBC_rot = [Node999_NMBC_rot(1,7) Node999_NMBC_rot(1,8) Node999_NMBC_rot(1,9);Node999_NMBC_rot(1,10) Node999_NMBC_rot(1,11) Node999_NMBC_rot(1,12);...
    Node999_NMBC_rot(1,13) Node999_NMBC_rot(1,14) Node999_NMBC_rot(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_999_AD_NMBC_rot_trans = R'*Result_S_999_AD_NMBC_rot*R;

Result_S_1001_AD_NMBC_rot = [Node1001_NMBC_rot(1,7) Node1001_NMBC_rot(1,8) Node1001_NMBC_rot(1,9);Node1001_NMBC_rot(1,10) Node1001_NMBC_rot(1,11) Node1001_NMBC_rot(1,12);...
    Node1001_NMBC_rot(1,13) Node1001_NMBC_rot(1,14) Node1001_NMBC_rot(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_1001_AD_NMBC_rot_trans = R'*Result_S_1001_AD_NMBC_rot*R;

Result_S_1003_AD_NMBC_rot = [Node1003_NMBC_rot(1,7) Node1003_NMBC_rot(1,8) Node1003_NMBC_rot(1,9);Node1003_NMBC_rot(1,10) Node1003_NMBC_rot(1,11) Node1003_NMBC_rot(1,12);...
    Node1003_NMBC_rot(1,13) Node1003_NMBC_rot(1,14) Node1003_NMBC_rot(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_1003_AD_NMBC_rot_trans = R'*Result_S_1003_AD_NMBC_rot*R;

Result_S_1005_AD_NMBC_rot = [Node1005_NMBC_rot(1,7) Node1005_NMBC_rot(1,8) Node1005_NMBC_rot(1,9);Node1005_NMBC_rot(1,10) Node1005_NMBC_rot(1,11) Node1005_NMBC_rot(1,12);...
    Node1005_NMBC_rot(1,13) Node1005_NMBC_rot(1,14) Node1005_NMBC_rot(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_1005_AD_NMBC_rot_trans = R'*Result_S_1005_AD_NMBC_rot*R;

Result_S_1007_AD_NMBC_rot = [Node1007_NMBC_rot(1,7) Node1007_NMBC_rot(1,8) Node1007_NMBC_rot(1,9);Node1007_NMBC_rot(1,10) Node1007_NMBC_rot(1,11) Node1007_NMBC_rot(1,12);...
    Node1007_NMBC_rot(1,13) Node1007_NMBC_rot(1,14) Node1007_NMBC_rot(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_1007_AD_NMBC_rot_trans = R'*Result_S_1007_AD_NMBC_rot*R;


Result_S_r_teta_NMBC_rot = [r_938/0.05 Result_S_938_AD_NMBC_rot_trans(1,2)/P_applied;r_939/0.05  Result_S_939_AD_NMBC_rot_trans(1,2)/P_applied;...
    r_941/0.05  Result_S_941_AD_NMBC_rot_trans(1,2)/P_applied;r_943/0.05  Result_S_943_AD_NMBC_rot_trans(1,2)/P_applied;r_945/0.05  Result_S_945_AD_NMBC_rot_trans(1,2)/P_applied;...
    r_947/0.05  Result_S_947_AD_NMBC_rot_trans(1,2)/P_applied;r_949/0.05  Result_S_949_AD_NMBC_rot_trans(1,2)/P_applied;r_951/0.05  Result_S_951_AD_NMBC_rot_trans(1,2)/P_applied;...
    r_953/0.05  Result_S_953_AD_NMBC_rot_trans(1,2)/P_applied;r_955/0.05  Result_S_955_AD_NMBC_rot_trans(1,2)/P_applied;r_957/0.05  Result_S_957_AD_NMBC_rot_trans(1,2)/P_applied;...
    r_959/0.05  Result_S_959_AD_NMBC_rot_trans(1,2)/P_applied;r_961/0.05  Result_S_961_AD_NMBC_rot_trans(1,2)/P_applied;r_963/0.05  Result_S_963_AD_NMBC_rot_trans(1,2)/P_applied;...
    r_965/0.05  Result_S_965_AD_NMBC_rot_trans(1,2)/P_applied;r_967/0.05  Result_S_967_AD_NMBC_rot_trans(1,2)/P_applied;r_969/0.05  Result_S_969_AD_NMBC_rot_trans(1,2)/P_applied;...
    r_971/0.05  Result_S_971_AD_NMBC_rot_trans(1,2)/P_applied;r_973/0.05  Result_S_973_AD_NMBC_rot_trans(1,2)/P_applied;r_975/0.05  Result_S_975_AD_NMBC_rot_trans(1,2)/P_applied;...
    r_977/0.05  Result_S_977_AD_NMBC_rot_trans(1,2)/P_applied;r_979/0.05  Result_S_979_AD_NMBC_rot_trans(1,2)/P_applied;r_981/0.05  Result_S_981_AD_NMBC_rot_trans(1,2)/P_applied;...
    r_983/0.05  Result_S_983_AD_NMBC_rot_trans(1,2)/P_applied;r_985/0.05  Result_S_985_AD_NMBC_rot_trans(1,2)/P_applied;r_987/0.05  Result_S_987_AD_NMBC_rot_trans(1,2)/P_applied;...
    r_989/0.05  Result_S_989_AD_NMBC_rot_trans(1,2)/P_applied;r_991/0.05  Result_S_991_AD_NMBC_rot_trans(1,2)/P_applied;r_993/0.05  Result_S_993_AD_NMBC_rot_trans(1,2)/P_applied;...
    r_995/0.05  Result_S_995_AD_NMBC_rot_trans(1,2)/P_applied;r_997/0.05  Result_S_997_AD_NMBC_rot_trans(1,2)/P_applied;r_999/0.05  Result_S_999_AD_NMBC_rot_trans(1,2)/P_applied;...
    r_1001/0.05  Result_S_1001_AD_NMBC_rot_trans(1,2)/P_applied;r_1003/0.05  Result_S_1003_AD_NMBC_rot_trans(1,2)/P_applied;r_1005/0.05  Result_S_1005_AD_NMBC_rot_trans(1,2)/P_applied;...
    r_1007/0.05  Result_S_1007_AD_NMBC_rot_trans(1,2)/P_applied];


Result_S_938_AD_NMBC_str = [Node938_NMBC_str(1,7) Node938_NMBC_str(1,8) Node938_NMBC_str(1,9);Node938_NMBC_str(1,10) Node938_NMBC_str(1,11) Node938_NMBC_str(1,12);...
    Node938_NMBC_str(1,13) Node938_NMBC_str(1,14) Node938_NMBC_str(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_938_AD_NMBC_str_trans = R'*Result_S_938_AD_NMBC_str*R;


Result_S_939_AD_NMBC_str = [Node939_NMBC_str(1,7) Node939_NMBC_str(1,8) Node939_NMBC_str(1,9);Node939_NMBC_str(1,10) Node939_NMBC_str(1,11) Node939_NMBC_str(1,12);...
    Node939_NMBC_str(1,13) Node939_NMBC_str(1,14) Node939_NMBC_str(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_939_AD_NMBC_str_trans = R'*Result_S_939_AD_NMBC_str*R;


Result_S_941_AD_NMBC_str = [Node941_NMBC_str(1,7) Node941_NMBC_str(1,8) Node941_NMBC_str(1,9);Node941_NMBC_str(1,10) Node941_NMBC_str(1,11) Node941_NMBC_str(1,12);...
    Node941_NMBC_str(1,13) Node941_NMBC_str(1,14) Node941_NMBC_str(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_941_AD_NMBC_str_trans = R'*Result_S_941_AD_NMBC_str*R;


Result_S_943_AD_NMBC_str = [Node943_NMBC_str(1,7) Node943_NMBC_str(1,8) Node943_NMBC_str(1,9);Node943_NMBC_str(1,10) Node943_NMBC_str(1,11) Node943_NMBC_str(1,12);...
    Node943_NMBC_str(1,13) Node943_NMBC_str(1,14) Node943_NMBC_str(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_943_AD_NMBC_str_trans = R'*Result_S_943_AD_NMBC_str*R;


Result_S_945_AD_NMBC_str = [Node945_NMBC_str(1,7) Node945_NMBC_str(1,8) Node945_NMBC_str(1,9);Node945_NMBC_str(1,10) Node945_NMBC_str(1,11) Node945_NMBC_str(1,12);...
    Node945_NMBC_str(1,13) Node945_NMBC_str(1,14) Node945_NMBC_str(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_945_AD_NMBC_str_trans = R'*Result_S_945_AD_NMBC_str*R;


Result_S_947_AD_NMBC_str = [Node947_NMBC_str(1,7) Node947_NMBC_str(1,8) Node947_NMBC_str(1,9);Node947_NMBC_str(1,10) Node947_NMBC_str(1,11) Node947_NMBC_str(1,12);...
    Node947_NMBC_str(1,13) Node947_NMBC_str(1,14) Node947_NMBC_str(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_947_AD_NMBC_str_trans = R'*Result_S_947_AD_NMBC_str*R;


Result_S_949_AD_NMBC_str = [Node949_NMBC_str(1,7) Node949_NMBC_str(1,8) Node949_NMBC_str(1,9);Node949_NMBC_str(1,10) Node949_NMBC_str(1,11) Node949_NMBC_str(1,12);...
    Node949_NMBC_str(1,13) Node949_NMBC_str(1,14) Node949_NMBC_str(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_949_AD_NMBC_str_trans = R'*Result_S_949_AD_NMBC_str*R;


Result_S_951_AD_NMBC_str = [Node951_NMBC_str(1,7) Node951_NMBC_str(1,8) Node951_NMBC_str(1,9);Node951_NMBC_str(1,10) Node951_NMBC_str(1,11) Node951_NMBC_str(1,12);...
    Node951_NMBC_str(1,13) Node951_NMBC_str(1,14) Node951_NMBC_str(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_951_AD_NMBC_str_trans = R'*Result_S_951_AD_NMBC_str*R;


Result_S_953_AD_NMBC_str = [Node953_NMBC_str(1,7) Node953_NMBC_str(1,8) Node953_NMBC_str(1,9);Node953_NMBC_str(1,10) Node953_NMBC_str(1,11) Node953_NMBC_str(1,12);...
    Node953_NMBC_str(1,13) Node953_NMBC_str(1,14) Node953_NMBC_str(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_953_AD_NMBC_str_trans = R'*Result_S_953_AD_NMBC_str*R;


Result_S_955_AD_NMBC_str = [Node955_NMBC_str(1,7) Node955_NMBC_str(1,8) Node955_NMBC_str(1,9);Node955_NMBC_str(1,10) Node955_NMBC_str(1,11) Node955_NMBC_str(1,12);...
    Node955_NMBC_str(1,13) Node955_NMBC_str(1,14) Node955_NMBC_str(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_955_AD_NMBC_str_trans = R'*Result_S_955_AD_NMBC_str*R;


Result_S_957_AD_NMBC_str = [Node957_NMBC_str(1,7) Node957_NMBC_str(1,8) Node957_NMBC_str(1,9);Node957_NMBC_str(1,10) Node957_NMBC_str(1,11) Node957_NMBC_str(1,12);...
    Node957_NMBC_str(1,13) Node957_NMBC_str(1,14) Node957_NMBC_str(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_957_AD_NMBC_str_trans = R'*Result_S_957_AD_NMBC_str*R;


Result_S_959_AD_NMBC_str = [Node959_NMBC_str(1,7) Node959_NMBC_str(1,8) Node959_NMBC_str(1,9);Node959_NMBC_str(1,10) Node959_NMBC_str(1,11) Node959_NMBC_str(1,12);...
    Node959_NMBC_str(1,13) Node959_NMBC_str(1,14) Node959_NMBC_str(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_959_AD_NMBC_str_trans = R'*Result_S_959_AD_NMBC_str*R;


Result_S_961_AD_NMBC_str = [Node961_NMBC_str(1,7) Node961_NMBC_str(1,8) Node961_NMBC_str(1,9);Node961_NMBC_str(1,10) Node961_NMBC_str(1,11) Node961_NMBC_str(1,12);...
    Node961_NMBC_str(1,13) Node961_NMBC_str(1,14) Node961_NMBC_str(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_961_AD_NMBC_str_trans = R'*Result_S_961_AD_NMBC_str*R;

Result_S_963_AD_NMBC_str = [Node963_NMBC_str(1,7) Node963_NMBC_str(1,8) Node963_NMBC_str(1,9);Node963_NMBC_str(1,10) Node963_NMBC_str(1,11) Node963_NMBC_str(1,12);...
    Node963_NMBC_str(1,13) Node963_NMBC_str(1,14) Node963_NMBC_str(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_963_AD_NMBC_str_trans = R'*Result_S_963_AD_NMBC_str*R;


Result_S_965_AD_NMBC_str = [Node965_NMBC_str(1,7) Node965_NMBC_str(1,8) Node965_NMBC_str(1,9);Node965_NMBC_str(1,10) Node965_NMBC_str(1,11) Node965_NMBC_str(1,12);...
    Node965_NMBC_str(1,13) Node965_NMBC_str(1,14) Node965_NMBC_str(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_965_AD_NMBC_str_trans = R'*Result_S_965_AD_NMBC_str*R;

Result_S_967_AD_NMBC_str = [Node967_NMBC_str(1,7) Node967_NMBC_str(1,8) Node967_NMBC_str(1,9);Node967_NMBC_str(1,10) Node967_NMBC_str(1,11) Node967_NMBC_str(1,12);...
    Node967_NMBC_str(1,13) Node967_NMBC_str(1,14) Node967_NMBC_str(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_967_AD_NMBC_str_trans = R'*Result_S_967_AD_NMBC_str*R;

Result_S_969_AD_NMBC_str = [Node969_NMBC_str(1,7) Node969_NMBC_str(1,8) Node969_NMBC_str(1,9);Node969_NMBC_str(1,10) Node969_NMBC_str(1,11) Node969_NMBC_str(1,12);...
    Node969_NMBC_str(1,13) Node969_NMBC_str(1,14) Node969_NMBC_str(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_969_AD_NMBC_str_trans = R'*Result_S_969_AD_NMBC_str*R;

Result_S_971_AD_NMBC_str = [Node971_NMBC_str(1,7) Node971_NMBC_str(1,8) Node971_NMBC_str(1,9);Node971_NMBC_str(1,10) Node971_NMBC_str(1,11) Node971_NMBC_str(1,12);...
    Node971_NMBC_str(1,13) Node971_NMBC_str(1,14) Node971_NMBC_str(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_971_AD_NMBC_str_trans = R'*Result_S_971_AD_NMBC_str*R;

Result_S_973_AD_NMBC_str = [Node973_NMBC_str(1,7) Node973_NMBC_str(1,8) Node973_NMBC_str(1,9);Node973_NMBC_str(1,10) Node973_NMBC_str(1,11) Node973_NMBC_str(1,12);...
    Node973_NMBC_str(1,13) Node973_NMBC_str(1,14) Node973_NMBC_str(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_973_AD_NMBC_str_trans = R'*Result_S_973_AD_NMBC_str*R;

Result_S_975_AD_NMBC_str = [Node975_NMBC_str(1,7) Node975_NMBC_str(1,8) Node975_NMBC_str(1,9);Node975_NMBC_str(1,10) Node975_NMBC_str(1,11) Node975_NMBC_str(1,12);...
    Node975_NMBC_str(1,13) Node975_NMBC_str(1,14) Node975_NMBC_str(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_975_AD_NMBC_str_trans = R'*Result_S_975_AD_NMBC_str*R;

Result_S_977_AD_NMBC_str = [Node977_NMBC_str(1,7) Node977_NMBC_str(1,8) Node977_NMBC_str(1,9);Node977_NMBC_str(1,10) Node977_NMBC_str(1,11) Node977_NMBC_str(1,12);...
    Node977_NMBC_str(1,13) Node977_NMBC_str(1,14) Node977_NMBC_str(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_977_AD_NMBC_str_trans = R'*Result_S_977_AD_NMBC_str*R;

Result_S_979_AD_NMBC_str = [Node979_NMBC_str(1,7) Node979_NMBC_str(1,8) Node979_NMBC_str(1,9);Node979_NMBC_str(1,10) Node979_NMBC_str(1,11) Node979_NMBC_str(1,12);...
    Node979_NMBC_str(1,13) Node979_NMBC_str(1,14) Node979_NMBC_str(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_979_AD_NMBC_str_trans = R'*Result_S_979_AD_NMBC_str*R;

Result_S_981_AD_NMBC_str = [Node981_NMBC_str(1,7) Node981_NMBC_str(1,8) Node981_NMBC_str(1,9);Node981_NMBC_str(1,10) Node981_NMBC_str(1,11) Node981_NMBC_str(1,12);...
    Node981_NMBC_str(1,13) Node981_NMBC_str(1,14) Node981_NMBC_str(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_981_AD_NMBC_str_trans = R'*Result_S_981_AD_NMBC_str*R;

Result_S_983_AD_NMBC_str = [Node983_NMBC_str(1,7) Node983_NMBC_str(1,8) Node983_NMBC_str(1,9);Node983_NMBC_str(1,10) Node983_NMBC_str(1,11) Node983_NMBC_str(1,12);...
    Node983_NMBC_str(1,13) Node983_NMBC_str(1,14) Node983_NMBC_str(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_983_AD_NMBC_str_trans = R'*Result_S_983_AD_NMBC_str*R;

Result_S_985_AD_NMBC_str = [Node985_NMBC_str(1,7) Node985_NMBC_str(1,8) Node985_NMBC_str(1,9);Node985_NMBC_str(1,10) Node985_NMBC_str(1,11) Node985_NMBC_str(1,12);...
    Node985_NMBC_str(1,13) Node985_NMBC_str(1,14) Node985_NMBC_str(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_985_AD_NMBC_str_trans = R'*Result_S_985_AD_NMBC_str*R;

Result_S_987_AD_NMBC_str = [Node987_NMBC_str(1,7) Node987_NMBC_str(1,8) Node987_NMBC_str(1,9);Node987_NMBC_str(1,10) Node987_NMBC_str(1,11) Node987_NMBC_str(1,12);...
    Node987_NMBC_str(1,13) Node987_NMBC_str(1,14) Node987_NMBC_str(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_987_AD_NMBC_str_trans = R'*Result_S_987_AD_NMBC_str*R;

Result_S_989_AD_NMBC_str = [Node989_NMBC_str(1,7) Node989_NMBC_str(1,8) Node989_NMBC_str(1,9);Node989_NMBC_str(1,10) Node989_NMBC_str(1,11) Node989_NMBC_str(1,12);...
    Node989_NMBC_str(1,13) Node989_NMBC_str(1,14) Node989_NMBC_str(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_989_AD_NMBC_str_trans = R'*Result_S_989_AD_NMBC_str*R;

Result_S_991_AD_NMBC_str = [Node991_NMBC_str(1,7) Node991_NMBC_str(1,8) Node991_NMBC_str(1,9);Node991_NMBC_str(1,10) Node991_NMBC_str(1,11) Node991_NMBC_str(1,12);...
    Node991_NMBC_str(1,13) Node991_NMBC_str(1,14) Node991_NMBC_str(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_991_AD_NMBC_str_trans = R'*Result_S_991_AD_NMBC_str*R;

Result_S_993_AD_NMBC_str = [Node993_NMBC_str(1,7) Node993_NMBC_str(1,8) Node993_NMBC_str(1,9);Node993_NMBC_str(1,10) Node993_NMBC_str(1,11) Node993_NMBC_str(1,12);...
    Node993_NMBC_str(1,13) Node993_NMBC_str(1,14) Node993_NMBC_str(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_993_AD_NMBC_str_trans = R'*Result_S_993_AD_NMBC_str*R;

Result_S_995_AD_NMBC_str = [Node995_NMBC_str(1,7) Node995_NMBC_str(1,8) Node995_NMBC_str(1,9);Node995_NMBC_str(1,10) Node995_NMBC_str(1,11) Node995_NMBC_str(1,12);...
    Node995_NMBC_str(1,13) Node995_NMBC_str(1,14) Node995_NMBC_str(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_995_AD_NMBC_str_trans = R'*Result_S_995_AD_NMBC_str*R;

Result_S_997_AD_NMBC_str = [Node997_NMBC_str(1,7) Node997_NMBC_str(1,8) Node997_NMBC_str(1,9);Node997_NMBC_str(1,10) Node997_NMBC_str(1,11) Node997_NMBC_str(1,12);...
    Node997_NMBC_str(1,13) Node997_NMBC_str(1,14) Node997_NMBC_str(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_997_AD_NMBC_str_trans = R'*Result_S_997_AD_NMBC_str*R;

Result_S_999_AD_NMBC_str = [Node999_NMBC_str(1,7) Node999_NMBC_str(1,8) Node999_NMBC_str(1,9);Node999_NMBC_str(1,10) Node999_NMBC_str(1,11) Node999_NMBC_str(1,12);...
    Node999_NMBC_str(1,13) Node999_NMBC_str(1,14) Node999_NMBC_str(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_999_AD_NMBC_str_trans = R'*Result_S_999_AD_NMBC_str*R;

Result_S_1001_AD_NMBC_str = [Node1001_NMBC_str(1,7) Node1001_NMBC_str(1,8) Node1001_NMBC_str(1,9);Node1001_NMBC_str(1,10) Node1001_NMBC_str(1,11) Node1001_NMBC_str(1,12);...
    Node1001_NMBC_str(1,13) Node1001_NMBC_str(1,14) Node1001_NMBC_str(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_1001_AD_NMBC_str_trans = R'*Result_S_1001_AD_NMBC_str*R;

Result_S_1003_AD_NMBC_str = [Node1003_NMBC_str(1,7) Node1003_NMBC_str(1,8) Node1003_NMBC_str(1,9);Node1003_NMBC_str(1,10) Node1003_NMBC_str(1,11) Node1003_NMBC_str(1,12);...
    Node1003_NMBC_str(1,13) Node1003_NMBC_str(1,14) Node1003_NMBC_str(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_1003_AD_NMBC_str_trans = R'*Result_S_1003_AD_NMBC_str*R;

Result_S_1005_AD_NMBC_str = [Node1005_NMBC_str(1,7) Node1005_NMBC_str(1,8) Node1005_NMBC_str(1,9);Node1005_NMBC_str(1,10) Node1005_NMBC_str(1,11) Node1005_NMBC_str(1,12);...
    Node1005_NMBC_str(1,13) Node1005_NMBC_str(1,14) Node1005_NMBC_str(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_1005_AD_NMBC_str_trans = R'*Result_S_1005_AD_NMBC_str*R;

Result_S_1007_AD_NMBC_str = [Node1007_NMBC_str(1,7) Node1007_NMBC_str(1,8) Node1007_NMBC_str(1,9);Node1007_NMBC_str(1,10) Node1007_NMBC_str(1,11) Node1007_NMBC_str(1,12);...
    Node1007_NMBC_str(1,13) Node1007_NMBC_str(1,14) Node1007_NMBC_str(1,15)];
R = [sqrt(2)/2 -sqrt(2)/2 0;sqrt(2)/2 sqrt(2)/2 0;0 0 1];
Result_S_1007_AD_NMBC_str_trans = R'*Result_S_1007_AD_NMBC_str*R;


Result_S_r_teta_NMBC_str = [r_938/0.05 Result_S_938_AD_NMBC_str_trans(1,2)/P_applied;r_939/0.05  Result_S_939_AD_NMBC_str_trans(1,2)/P_applied;...
    r_941/0.05  Result_S_941_AD_NMBC_str_trans(1,2)/P_applied;r_943/0.05  Result_S_943_AD_NMBC_str_trans(1,2)/P_applied;r_945/0.05  Result_S_945_AD_NMBC_str_trans(1,2)/P_applied;...
    r_947/0.05  Result_S_947_AD_NMBC_str_trans(1,2)/P_applied;r_949/0.05  Result_S_949_AD_NMBC_str_trans(1,2)/P_applied;r_951/0.05  Result_S_951_AD_NMBC_str_trans(1,2)/P_applied;...
    r_953/0.05  Result_S_953_AD_NMBC_str_trans(1,2)/P_applied;r_955/0.05  Result_S_955_AD_NMBC_str_trans(1,2)/P_applied;r_957/0.05  Result_S_957_AD_NMBC_str_trans(1,2)/P_applied;...
    r_959/0.05  Result_S_959_AD_NMBC_str_trans(1,2)/P_applied;r_961/0.05  Result_S_961_AD_NMBC_str_trans(1,2)/P_applied;r_963/0.05  Result_S_963_AD_NMBC_str_trans(1,2)/P_applied;...
    r_965/0.05  Result_S_965_AD_NMBC_str_trans(1,2)/P_applied;r_967/0.05  Result_S_967_AD_NMBC_str_trans(1,2)/P_applied;r_969/0.05  Result_S_969_AD_NMBC_str_trans(1,2)/P_applied;...
    r_971/0.05  Result_S_971_AD_NMBC_str_trans(1,2)/P_applied;r_973/0.05  Result_S_973_AD_NMBC_str_trans(1,2)/P_applied;r_975/0.05  Result_S_975_AD_NMBC_str_trans(1,2)/P_applied;...
    r_977/0.05  Result_S_977_AD_NMBC_str_trans(1,2)/P_applied;r_979/0.05  Result_S_979_AD_NMBC_str_trans(1,2)/P_applied;r_981/0.05  Result_S_981_AD_NMBC_str_trans(1,2)/P_applied;...
    r_983/0.05  Result_S_983_AD_NMBC_str_trans(1,2)/P_applied;r_985/0.05  Result_S_985_AD_NMBC_str_trans(1,2)/P_applied;r_987/0.05  Result_S_987_AD_NMBC_str_trans(1,2)/P_applied;...
    r_989/0.05  Result_S_989_AD_NMBC_str_trans(1,2)/P_applied;r_991/0.05  Result_S_991_AD_NMBC_str_trans(1,2)/P_applied;r_993/0.05  Result_S_993_AD_NMBC_str_trans(1,2)/P_applied;...
    r_995/0.05  Result_S_995_AD_NMBC_str_trans(1,2)/P_applied;r_997/0.05  Result_S_997_AD_NMBC_str_trans(1,2)/P_applied;r_999/0.05  Result_S_999_AD_NMBC_str_trans(1,2)/P_applied;...
    r_1001/0.05  Result_S_1001_AD_NMBC_str_trans(1,2)/P_applied;r_1003/0.05  Result_S_1003_AD_NMBC_str_trans(1,2)/P_applied;r_1005/0.05  Result_S_1005_AD_NMBC_str_trans(1,2)/P_applied;...
    r_1007/0.05  Result_S_1007_AD_NMBC_str_trans(1,2)/P_applied];


Result_S_teta_r_MBC_com = [r_938/0.05 Result_S_938_AD_MBC_com_trans(2,1)/P_applied;r_939/0.05  Result_S_939_AD_MBC_com_trans(2,1)/P_applied;...
    r_941/0.05  Result_S_941_AD_MBC_com_trans(2,1)/P_applied;r_943/0.05  Result_S_943_AD_MBC_com_trans(2,1)/P_applied;r_945/0.05  Result_S_945_AD_MBC_com_trans(2,1)/P_applied;...
    r_947/0.05  Result_S_947_AD_MBC_com_trans(2,1)/P_applied;r_949/0.05  Result_S_949_AD_MBC_com_trans(2,1)/P_applied;r_951/0.05  Result_S_951_AD_MBC_com_trans(2,1)/P_applied;...
    r_953/0.05  Result_S_953_AD_MBC_com_trans(2,1)/P_applied;r_955/0.05  Result_S_955_AD_MBC_com_trans(2,1)/P_applied;r_957/0.05  Result_S_957_AD_MBC_com_trans(2,1)/P_applied;...
    r_959/0.05  Result_S_959_AD_MBC_com_trans(2,1)/P_applied;r_961/0.05  Result_S_961_AD_MBC_com_trans(2,1)/P_applied;r_963/0.05  Result_S_963_AD_MBC_com_trans(2,1)/P_applied;...
    r_965/0.05  Result_S_965_AD_MBC_com_trans(2,1)/P_applied;r_967/0.05  Result_S_967_AD_MBC_com_trans(2,1)/P_applied;r_969/0.05  Result_S_969_AD_MBC_com_trans(2,1)/P_applied;...
    r_971/0.05  Result_S_971_AD_MBC_com_trans(2,1)/P_applied;r_973/0.05  Result_S_973_AD_MBC_com_trans(2,1)/P_applied;r_975/0.05  Result_S_975_AD_MBC_com_trans(2,1)/P_applied;...
    r_977/0.05  Result_S_977_AD_MBC_com_trans(2,1)/P_applied;r_979/0.05  Result_S_979_AD_MBC_com_trans(2,1)/P_applied;r_981/0.05  Result_S_981_AD_MBC_com_trans(2,1)/P_applied;...
    r_983/0.05  Result_S_983_AD_MBC_com_trans(2,1)/P_applied;r_985/0.05  Result_S_985_AD_MBC_com_trans(2,1)/P_applied;r_987/0.05  Result_S_987_AD_MBC_com_trans(2,1)/P_applied;...
    r_989/0.05  Result_S_989_AD_MBC_com_trans(2,1)/P_applied;r_991/0.05  Result_S_991_AD_MBC_com_trans(2,1)/P_applied;r_993/0.05  Result_S_993_AD_MBC_com_trans(2,1)/P_applied;...
    r_995/0.05  Result_S_995_AD_MBC_com_trans(2,1)/P_applied;r_997/0.05  Result_S_997_AD_MBC_com_trans(2,1)/P_applied;r_999/0.05  Result_S_999_AD_MBC_com_trans(2,1)/P_applied;...
    r_1001/0.05  Result_S_1001_AD_MBC_com_trans(2,1)/P_applied;r_1003/0.05  Result_S_1003_AD_MBC_com_trans(2,1)/P_applied;r_1005/0.05  Result_S_1005_AD_MBC_com_trans(2,1)/P_applied;...
    r_1007/0.05  Result_S_1007_AD_MBC_com_trans(2,1)/P_applied];

Result_S_teta_r_MBC_rot = [r_938/0.05 Result_S_938_AD_MBC_rot_trans(2,1)/P_applied;r_939/0.05  Result_S_939_AD_MBC_rot_trans(2,1)/P_applied;...
    r_941/0.05  Result_S_941_AD_MBC_rot_trans(2,1)/P_applied;r_943/0.05  Result_S_943_AD_MBC_rot_trans(2,1)/P_applied;r_945/0.05  Result_S_945_AD_MBC_rot_trans(2,1)/P_applied;...
    r_947/0.05  Result_S_947_AD_MBC_rot_trans(2,1)/P_applied;r_949/0.05  Result_S_949_AD_MBC_rot_trans(2,1)/P_applied;r_951/0.05  Result_S_951_AD_MBC_rot_trans(2,1)/P_applied;...
    r_953/0.05  Result_S_953_AD_MBC_rot_trans(2,1)/P_applied;r_955/0.05  Result_S_955_AD_MBC_rot_trans(2,1)/P_applied;r_957/0.05  Result_S_957_AD_MBC_rot_trans(2,1)/P_applied;...
    r_959/0.05  Result_S_959_AD_MBC_rot_trans(2,1)/P_applied;r_961/0.05  Result_S_961_AD_MBC_rot_trans(2,1)/P_applied;r_963/0.05  Result_S_963_AD_MBC_rot_trans(2,1)/P_applied;...
    r_965/0.05  Result_S_965_AD_MBC_rot_trans(2,1)/P_applied;r_967/0.05  Result_S_967_AD_MBC_rot_trans(2,1)/P_applied;r_969/0.05  Result_S_969_AD_MBC_rot_trans(2,1)/P_applied;...
    r_971/0.05  Result_S_971_AD_MBC_rot_trans(2,1)/P_applied;r_973/0.05  Result_S_973_AD_MBC_rot_trans(2,1)/P_applied;r_975/0.05  Result_S_975_AD_MBC_rot_trans(2,1)/P_applied;...
    r_977/0.05  Result_S_977_AD_MBC_rot_trans(2,1)/P_applied;r_979/0.05  Result_S_979_AD_MBC_rot_trans(2,1)/P_applied;r_981/0.05  Result_S_981_AD_MBC_rot_trans(2,1)/P_applied;...
    r_983/0.05  Result_S_983_AD_MBC_rot_trans(2,1)/P_applied;r_985/0.05  Result_S_985_AD_MBC_rot_trans(2,1)/P_applied;r_987/0.05  Result_S_987_AD_MBC_rot_trans(2,1)/P_applied;...
    r_989/0.05  Result_S_989_AD_MBC_rot_trans(2,1)/P_applied;r_991/0.05  Result_S_991_AD_MBC_rot_trans(2,1)/P_applied;r_993/0.05  Result_S_993_AD_MBC_rot_trans(2,1)/P_applied;...
    r_995/0.05  Result_S_995_AD_MBC_rot_trans(2,1)/P_applied;r_997/0.05  Result_S_997_AD_MBC_rot_trans(2,1)/P_applied;r_999/0.05  Result_S_999_AD_MBC_rot_trans(2,1)/P_applied;...
    r_1001/0.05  Result_S_1001_AD_MBC_rot_trans(2,1)/P_applied;r_1003/0.05  Result_S_1003_AD_MBC_rot_trans(2,1)/P_applied;r_1005/0.05  Result_S_1005_AD_MBC_rot_trans(2,1)/P_applied;...
    r_1007/0.05  Result_S_1007_AD_MBC_rot_trans(2,1)/P_applied];

Result_S_teta_r_MBC_str = [r_938/0.05 Result_S_938_AD_MBC_str_trans(2,1)/P_applied;r_939/0.05  Result_S_939_AD_MBC_str_trans(2,1)/P_applied;...
    r_941/0.05  Result_S_941_AD_MBC_str_trans(2,1)/P_applied;r_943/0.05  Result_S_943_AD_MBC_str_trans(2,1)/P_applied;r_945/0.05  Result_S_945_AD_MBC_str_trans(2,1)/P_applied;...
    r_947/0.05  Result_S_947_AD_MBC_str_trans(2,1)/P_applied;r_949/0.05  Result_S_949_AD_MBC_str_trans(2,1)/P_applied;r_951/0.05  Result_S_951_AD_MBC_str_trans(2,1)/P_applied;...
    r_953/0.05  Result_S_953_AD_MBC_str_trans(2,1)/P_applied;r_955/0.05  Result_S_955_AD_MBC_str_trans(2,1)/P_applied;r_957/0.05  Result_S_957_AD_MBC_str_trans(2,1)/P_applied;...
    r_959/0.05  Result_S_959_AD_MBC_str_trans(2,1)/P_applied;r_961/0.05  Result_S_961_AD_MBC_str_trans(2,1)/P_applied;r_963/0.05  Result_S_963_AD_MBC_str_trans(2,1)/P_applied;...
    r_965/0.05  Result_S_965_AD_MBC_str_trans(2,1)/P_applied;r_967/0.05  Result_S_967_AD_MBC_str_trans(2,1)/P_applied;r_969/0.05  Result_S_969_AD_MBC_str_trans(2,1)/P_applied;...
    r_971/0.05  Result_S_971_AD_MBC_str_trans(2,1)/P_applied;r_973/0.05  Result_S_973_AD_MBC_str_trans(2,1)/P_applied;r_975/0.05  Result_S_975_AD_MBC_str_trans(2,1)/P_applied;...
    r_977/0.05  Result_S_977_AD_MBC_str_trans(2,1)/P_applied;r_979/0.05  Result_S_979_AD_MBC_str_trans(2,1)/P_applied;r_981/0.05  Result_S_981_AD_MBC_str_trans(2,1)/P_applied;...
    r_983/0.05  Result_S_983_AD_MBC_str_trans(2,1)/P_applied;r_985/0.05  Result_S_985_AD_MBC_str_trans(2,1)/P_applied;r_987/0.05  Result_S_987_AD_MBC_str_trans(2,1)/P_applied;...
    r_989/0.05  Result_S_989_AD_MBC_str_trans(2,1)/P_applied;r_991/0.05  Result_S_991_AD_MBC_str_trans(2,1)/P_applied;r_993/0.05  Result_S_993_AD_MBC_str_trans(2,1)/P_applied;...
    r_995/0.05  Result_S_995_AD_MBC_str_trans(2,1)/P_applied;r_997/0.05  Result_S_997_AD_MBC_str_trans(2,1)/P_applied;r_999/0.05  Result_S_999_AD_MBC_str_trans(2,1)/P_applied;...
    r_1001/0.05  Result_S_1001_AD_MBC_str_trans(2,1)/P_applied;r_1003/0.05  Result_S_1003_AD_MBC_str_trans(2,1)/P_applied;r_1005/0.05  Result_S_1005_AD_MBC_str_trans(2,1)/P_applied;...
    r_1007/0.05  Result_S_1007_AD_MBC_str_trans(2,1)/P_applied];

Result_S_teta_r_NMBC_com = [r_938/0.05 Result_S_938_AD_NMBC_com_trans(2,1)/P_applied;r_939/0.05  Result_S_939_AD_NMBC_com_trans(2,1)/P_applied;...
    r_941/0.05  Result_S_941_AD_NMBC_com_trans(2,1)/P_applied;r_943/0.05  Result_S_943_AD_NMBC_com_trans(2,1)/P_applied;r_945/0.05  Result_S_945_AD_NMBC_com_trans(2,1)/P_applied;...
    r_947/0.05  Result_S_947_AD_NMBC_com_trans(2,1)/P_applied;r_949/0.05  Result_S_949_AD_NMBC_com_trans(2,1)/P_applied;r_951/0.05  Result_S_951_AD_NMBC_com_trans(2,1)/P_applied;...
    r_953/0.05  Result_S_953_AD_NMBC_com_trans(2,1)/P_applied;r_955/0.05  Result_S_955_AD_NMBC_com_trans(2,1)/P_applied;r_957/0.05  Result_S_957_AD_NMBC_com_trans(2,1)/P_applied;...
    r_959/0.05  Result_S_959_AD_NMBC_com_trans(2,1)/P_applied;r_961/0.05  Result_S_961_AD_NMBC_com_trans(2,1)/P_applied;r_963/0.05  Result_S_963_AD_NMBC_com_trans(2,1)/P_applied;...
    r_965/0.05  Result_S_965_AD_NMBC_com_trans(2,1)/P_applied;r_967/0.05  Result_S_967_AD_NMBC_com_trans(2,1)/P_applied;r_969/0.05  Result_S_969_AD_NMBC_com_trans(2,1)/P_applied;...
    r_971/0.05  Result_S_971_AD_NMBC_com_trans(2,1)/P_applied;r_973/0.05  Result_S_973_AD_NMBC_com_trans(2,1)/P_applied;r_975/0.05  Result_S_975_AD_NMBC_com_trans(2,1)/P_applied;...
    r_977/0.05  Result_S_977_AD_NMBC_com_trans(2,1)/P_applied;r_979/0.05  Result_S_979_AD_NMBC_com_trans(2,1)/P_applied;r_981/0.05  Result_S_981_AD_NMBC_com_trans(2,1)/P_applied;...
    r_983/0.05  Result_S_983_AD_NMBC_com_trans(2,1)/P_applied;r_985/0.05  Result_S_985_AD_NMBC_com_trans(2,1)/P_applied;r_987/0.05  Result_S_987_AD_NMBC_com_trans(2,1)/P_applied;...
    r_989/0.05  Result_S_989_AD_NMBC_com_trans(2,1)/P_applied;r_991/0.05  Result_S_991_AD_NMBC_com_trans(2,1)/P_applied;r_993/0.05  Result_S_993_AD_NMBC_com_trans(2,1)/P_applied;...
    r_995/0.05  Result_S_995_AD_NMBC_com_trans(2,1)/P_applied;r_997/0.05  Result_S_997_AD_NMBC_com_trans(2,1)/P_applied;r_999/0.05  Result_S_999_AD_NMBC_com_trans(2,1)/P_applied;...
    r_1001/0.05  Result_S_1001_AD_NMBC_com_trans(2,1)/P_applied;r_1003/0.05  Result_S_1003_AD_NMBC_com_trans(2,1)/P_applied;r_1005/0.05  Result_S_1005_AD_NMBC_com_trans(2,1)/P_applied;...
    r_1007/0.05  Result_S_1007_AD_NMBC_com_trans(2,1)/P_applied];

Result_S_teta_r_NMBC_rot = [r_938/0.05 Result_S_938_AD_NMBC_rot_trans(2,1)/P_applied;r_939/0.05  Result_S_939_AD_NMBC_rot_trans(2,1)/P_applied;...
    r_941/0.05  Result_S_941_AD_NMBC_rot_trans(2,1)/P_applied;r_943/0.05  Result_S_943_AD_NMBC_rot_trans(2,1)/P_applied;r_945/0.05  Result_S_945_AD_NMBC_rot_trans(2,1)/P_applied;...
    r_947/0.05  Result_S_947_AD_NMBC_rot_trans(2,1)/P_applied;r_949/0.05  Result_S_949_AD_NMBC_rot_trans(2,1)/P_applied;r_951/0.05  Result_S_951_AD_NMBC_rot_trans(2,1)/P_applied;...
    r_953/0.05  Result_S_953_AD_NMBC_rot_trans(2,1)/P_applied;r_955/0.05  Result_S_955_AD_NMBC_rot_trans(2,1)/P_applied;r_957/0.05  Result_S_957_AD_NMBC_rot_trans(2,1)/P_applied;...
    r_959/0.05  Result_S_959_AD_NMBC_rot_trans(2,1)/P_applied;r_961/0.05  Result_S_961_AD_NMBC_rot_trans(2,1)/P_applied;r_963/0.05  Result_S_963_AD_NMBC_rot_trans(2,1)/P_applied;...
    r_965/0.05  Result_S_965_AD_NMBC_rot_trans(2,1)/P_applied;r_967/0.05  Result_S_967_AD_NMBC_rot_trans(2,1)/P_applied;r_969/0.05  Result_S_969_AD_NMBC_rot_trans(2,1)/P_applied;...
    r_971/0.05  Result_S_971_AD_NMBC_rot_trans(2,1)/P_applied;r_973/0.05  Result_S_973_AD_NMBC_rot_trans(2,1)/P_applied;r_975/0.05  Result_S_975_AD_NMBC_rot_trans(2,1)/P_applied;...
    r_977/0.05  Result_S_977_AD_NMBC_rot_trans(2,1)/P_applied;r_979/0.05  Result_S_979_AD_NMBC_rot_trans(2,1)/P_applied;r_981/0.05  Result_S_981_AD_NMBC_rot_trans(2,1)/P_applied;...
    r_983/0.05  Result_S_983_AD_NMBC_rot_trans(2,1)/P_applied;r_985/0.05  Result_S_985_AD_NMBC_rot_trans(2,1)/P_applied;r_987/0.05  Result_S_987_AD_NMBC_rot_trans(2,1)/P_applied;...
    r_989/0.05  Result_S_989_AD_NMBC_rot_trans(2,1)/P_applied;r_991/0.05  Result_S_991_AD_NMBC_rot_trans(2,1)/P_applied;r_993/0.05  Result_S_993_AD_NMBC_rot_trans(2,1)/P_applied;...
    r_995/0.05  Result_S_995_AD_NMBC_rot_trans(2,1)/P_applied;r_997/0.05  Result_S_997_AD_NMBC_rot_trans(2,1)/P_applied;r_999/0.05  Result_S_999_AD_NMBC_rot_trans(2,1)/P_applied;...
    r_1001/0.05  Result_S_1001_AD_NMBC_rot_trans(2,1)/P_applied;r_1003/0.05  Result_S_1003_AD_NMBC_rot_trans(2,1)/P_applied;r_1005/0.05  Result_S_1005_AD_NMBC_rot_trans(2,1)/P_applied;...
    r_1007/0.05  Result_S_1007_AD_NMBC_rot_trans(2,1)/P_applied];

Result_S_teta_r_NMBC_str = [r_938/0.05 Result_S_938_AD_NMBC_str_trans(2,1)/P_applied;r_939/0.05  Result_S_939_AD_NMBC_str_trans(2,1)/P_applied;...
    r_941/0.05  Result_S_941_AD_NMBC_str_trans(2,1)/P_applied;r_943/0.05  Result_S_943_AD_NMBC_str_trans(2,1)/P_applied;r_945/0.05  Result_S_945_AD_NMBC_str_trans(2,1)/P_applied;...
    r_947/0.05  Result_S_947_AD_NMBC_str_trans(2,1)/P_applied;r_949/0.05  Result_S_949_AD_NMBC_str_trans(2,1)/P_applied;r_951/0.05  Result_S_951_AD_NMBC_str_trans(2,1)/P_applied;...
    r_953/0.05  Result_S_953_AD_NMBC_str_trans(2,1)/P_applied;r_955/0.05  Result_S_955_AD_NMBC_str_trans(2,1)/P_applied;r_957/0.05  Result_S_957_AD_NMBC_str_trans(2,1)/P_applied;...
    r_959/0.05  Result_S_959_AD_NMBC_str_trans(2,1)/P_applied;r_961/0.05  Result_S_961_AD_NMBC_str_trans(2,1)/P_applied;r_963/0.05  Result_S_963_AD_NMBC_str_trans(2,1)/P_applied;...
    r_965/0.05  Result_S_965_AD_NMBC_str_trans(2,1)/P_applied;r_967/0.05  Result_S_967_AD_NMBC_str_trans(2,1)/P_applied;r_969/0.05  Result_S_969_AD_NMBC_str_trans(2,1)/P_applied;...
    r_971/0.05  Result_S_971_AD_NMBC_str_trans(2,1)/P_applied;r_973/0.05  Result_S_973_AD_NMBC_str_trans(2,1)/P_applied;r_975/0.05  Result_S_975_AD_NMBC_str_trans(2,1)/P_applied;...
    r_977/0.05  Result_S_977_AD_NMBC_str_trans(2,1)/P_applied;r_979/0.05  Result_S_979_AD_NMBC_str_trans(2,1)/P_applied;r_981/0.05  Result_S_981_AD_NMBC_str_trans(2,1)/P_applied;...
    r_983/0.05  Result_S_983_AD_NMBC_str_trans(2,1)/P_applied;r_985/0.05  Result_S_985_AD_NMBC_str_trans(2,1)/P_applied;r_987/0.05  Result_S_987_AD_NMBC_str_trans(2,1)/P_applied;...
    r_989/0.05  Result_S_989_AD_NMBC_str_trans(2,1)/P_applied;r_991/0.05  Result_S_991_AD_NMBC_str_trans(2,1)/P_applied;r_993/0.05  Result_S_993_AD_NMBC_str_trans(2,1)/P_applied;...
    r_995/0.05  Result_S_995_AD_NMBC_str_trans(2,1)/P_applied;r_997/0.05  Result_S_997_AD_NMBC_str_trans(2,1)/P_applied;r_999/0.05  Result_S_999_AD_NMBC_str_trans(2,1)/P_applied;...
    r_1001/0.05  Result_S_1001_AD_NMBC_str_trans(2,1)/P_applied;r_1003/0.05  Result_S_1003_AD_NMBC_str_trans(2,1)/P_applied;r_1005/0.05  Result_S_1005_AD_NMBC_str_trans(2,1)/P_applied;...
    r_1007/0.05  Result_S_1007_AD_NMBC_str_trans(2,1)/P_applied];


Result_phi11_intensity_MBC_com = [Node1012_SI_MBC_com(1,3)+0.25 Node1012_SI_MBC_com(1,5)/P_applied;Node1011_SI_MBC_com(1,3)+0.25 Node1011_SI_MBC_com(1,5)/P_applied;...
    Node1018_SI_MBC_com(1,3)+0.25 Node1018_SI_MBC_com(1,5)/P_applied;Node1022_SI_MBC_com(1,3)+0.25 Node1022_SI_MBC_com(1,5)/P_applied;...
    Node1022_SI_MBC_com(1,3)+0.25 Node1022_SI_MBC_com(1,5)/P_applied;Node1026_SI_MBC_com(1,3)+0.25 Node1026_SI_MBC_com(1,5)/P_applied;...
    Node1030_SI_MBC_com(1,3)+0.25 Node1030_SI_MBC_com(1,5)/P_applied;Node1034_SI_MBC_com(1,3)+0.25 Node1034_SI_MBC_com(1,5)/P_applied;...
    Node1038_SI_MBC_com(1,3)+0.25 Node1038_SI_MBC_com(1,5)/P_applied;Node1042_SI_MBC_com(1,3)+0.25 Node1042_SI_MBC_com(1,5)/P_applied;...
    Node1046_SI_MBC_com(1,3)+0.25 Node1046_SI_MBC_com(1,5)/P_applied;Node1050_SI_MBC_com(1,3)+0.25 Node1050_SI_MBC_com(1,5)/P_applied;...
    Node1054_SI_MBC_com(1,3)+0.25 Node1054_SI_MBC_com(1,5)/P_applied;Node1058_SI_MBC_com(1,3)+0.25 Node1058_SI_MBC_com(1,5)/P_applied;...
    Node1062_SI_MBC_com(1,3)+0.25 Node1062_SI_MBC_com(1,5)/P_applied;Node1066_SI_MBC_com(1,3)+0.25 Node1066_SI_MBC_com(1,5)/P_applied;...
    Node1070_SI_MBC_com(1,3)+0.25 Node1070_SI_MBC_com(1,5)/P_applied;Node1074_SI_MBC_com(1,3)+0.25 Node1074_SI_MBC_com(1,5)/P_applied;...
    Node1078_SI_MBC_com(1,3)+0.25 Node1078_SI_MBC_com(1,5)/P_applied;Node1082_SI_MBC_com(1,3)+0.25 Node1082_SI_MBC_com(1,5)/P_applied;...
    Node1086_SI_MBC_com(1,3)+0.25 Node1086_SI_MBC_com(1,5)/P_applied;Node1090_SI_MBC_com(1,3)+0.25 Node1090_SI_MBC_com(1,5)/P_applied;...
    Node1094_SI_MBC_com(1,3)+0.25 Node1094_SI_MBC_com(1,5)/P_applied;Node1098_SI_MBC_com(1,3)+0.25 Node1098_SI_MBC_com(1,5)/P_applied;...
    Node1102_SI_MBC_com(1,3)+0.25 Node1102_SI_MBC_com(1,5)/P_applied;Node1106_SI_MBC_com(1,3)+0.25 Node1106_SI_MBC_com(1,5)/P_applied;...
    Node1110_SI_MBC_com(1,3)+0.25 Node1110_SI_MBC_com(1,5)/P_applied;Node1114_SI_MBC_com(1,3)+0.25 Node1114_SI_MBC_com(1,5)/P_applied;...
    Node1118_SI_MBC_com(1,3)+0.25 Node1118_SI_MBC_com(1,5)/P_applied;Node1122_SI_MBC_com(1,3)+0.25 Node1122_SI_MBC_com(1,5)/P_applied;...
    Node1126_SI_MBC_com(1,3)+0.25 Node1126_SI_MBC_com(1,5)/P_applied;Node1130_SI_MBC_com(1,3)+0.25 Node1130_SI_MBC_com(1,5)/P_applied;...
    Node1134_SI_MBC_com(1,3)+0.25 Node1134_SI_MBC_com(1,5)/P_applied;Node1138_SI_MBC_com(1,3)+0.25 Node1138_SI_MBC_com(1,5)/P_applied;...
    Node1142_SI_MBC_com(1,3)+0.25 Node1142_SI_MBC_com(1,5)/P_applied;Node1146_SI_MBC_com(1,3)+0.25 Node1146_SI_MBC_com(1,5)/P_applied;...
    Node1150_SI_MBC_com(1,3)+0.25 Node1150_SI_MBC_com(1,5)/P_applied];


Result_phi11_intensity_MBC_rot = [Node1012_SI_MBC_rot(1,3)+0.25 Node1012_SI_MBC_rot(1,5)/P_applied;Node1011_SI_MBC_rot(1,3)+0.25 Node1011_SI_MBC_rot(1,5)/P_applied;...
    Node1018_SI_MBC_rot(1,3)+0.25 Node1018_SI_MBC_rot(1,5)/P_applied;Node1022_SI_MBC_rot(1,3)+0.25 Node1022_SI_MBC_rot(1,5)/P_applied;...
    Node1022_SI_MBC_rot(1,3)+0.25 Node1022_SI_MBC_rot(1,5)/P_applied;Node1026_SI_MBC_rot(1,3)+0.25 Node1026_SI_MBC_rot(1,5)/P_applied;...
    Node1030_SI_MBC_rot(1,3)+0.25 Node1030_SI_MBC_rot(1,5)/P_applied;Node1034_SI_MBC_rot(1,3)+0.25 Node1034_SI_MBC_rot(1,5)/P_applied;...
    Node1038_SI_MBC_rot(1,3)+0.25 Node1038_SI_MBC_rot(1,5)/P_applied;Node1042_SI_MBC_rot(1,3)+0.25 Node1042_SI_MBC_rot(1,5)/P_applied;...
    Node1046_SI_MBC_rot(1,3)+0.25 Node1046_SI_MBC_rot(1,5)/P_applied;Node1050_SI_MBC_rot(1,3)+0.25 Node1050_SI_MBC_rot(1,5)/P_applied;...
    Node1054_SI_MBC_rot(1,3)+0.25 Node1054_SI_MBC_rot(1,5)/P_applied;Node1058_SI_MBC_rot(1,3)+0.25 Node1058_SI_MBC_rot(1,5)/P_applied;...
    Node1062_SI_MBC_rot(1,3)+0.25 Node1062_SI_MBC_rot(1,5)/P_applied;Node1066_SI_MBC_rot(1,3)+0.25 Node1066_SI_MBC_rot(1,5)/P_applied;...
    Node1070_SI_MBC_rot(1,3)+0.25 Node1070_SI_MBC_rot(1,5)/P_applied;Node1074_SI_MBC_rot(1,3)+0.25 Node1074_SI_MBC_rot(1,5)/P_applied;...
    Node1078_SI_MBC_rot(1,3)+0.25 Node1078_SI_MBC_rot(1,5)/P_applied;Node1082_SI_MBC_rot(1,3)+0.25 Node1082_SI_MBC_rot(1,5)/P_applied;...
    Node1086_SI_MBC_rot(1,3)+0.25 Node1086_SI_MBC_rot(1,5)/P_applied;Node1090_SI_MBC_rot(1,3)+0.25 Node1090_SI_MBC_rot(1,5)/P_applied;...
    Node1094_SI_MBC_rot(1,3)+0.25 Node1094_SI_MBC_rot(1,5)/P_applied;Node1098_SI_MBC_rot(1,3)+0.25 Node1098_SI_MBC_rot(1,5)/P_applied;...
    Node1102_SI_MBC_rot(1,3)+0.25 Node1102_SI_MBC_rot(1,5)/P_applied;Node1106_SI_MBC_rot(1,3)+0.25 Node1106_SI_MBC_rot(1,5)/P_applied;...
    Node1110_SI_MBC_rot(1,3)+0.25 Node1110_SI_MBC_rot(1,5)/P_applied;Node1114_SI_MBC_rot(1,3)+0.25 Node1114_SI_MBC_rot(1,5)/P_applied;...
    Node1118_SI_MBC_rot(1,3)+0.25 Node1118_SI_MBC_rot(1,5)/P_applied;Node1122_SI_MBC_rot(1,3)+0.25 Node1122_SI_MBC_rot(1,5)/P_applied;...
    Node1126_SI_MBC_rot(1,3)+0.25 Node1126_SI_MBC_rot(1,5)/P_applied;Node1130_SI_MBC_rot(1,3)+0.25 Node1130_SI_MBC_rot(1,5)/P_applied;...
    Node1134_SI_MBC_rot(1,3)+0.25 Node1134_SI_MBC_rot(1,5)/P_applied;Node1138_SI_MBC_rot(1,3)+0.25 Node1138_SI_MBC_rot(1,5)/P_applied;...
    Node1142_SI_MBC_rot(1,3)+0.25 Node1142_SI_MBC_rot(1,5)/P_applied;Node1146_SI_MBC_rot(1,3)+0.25 Node1146_SI_MBC_rot(1,5)/P_applied;...
    Node1150_SI_MBC_rot(1,3)+0.25 Node1150_SI_MBC_rot(1,5)/P_applied];


Result_phi11_intensity_MBC_str = [Node1012_SI_MBC_str(1,3)+0.25 Node1012_SI_MBC_str(1,5)/P_applied;Node1011_SI_MBC_str(1,3)+0.25 Node1011_SI_MBC_str(1,5)/P_applied;...
    Node1018_SI_MBC_str(1,3)+0.25 Node1018_SI_MBC_str(1,5)/P_applied;Node1022_SI_MBC_str(1,3)+0.25 Node1022_SI_MBC_str(1,5)/P_applied;...
    Node1022_SI_MBC_str(1,3)+0.25 Node1022_SI_MBC_str(1,5)/P_applied;Node1026_SI_MBC_str(1,3)+0.25 Node1026_SI_MBC_str(1,5)/P_applied;...
    Node1030_SI_MBC_str(1,3)+0.25 Node1030_SI_MBC_str(1,5)/P_applied;Node1034_SI_MBC_str(1,3)+0.25 Node1034_SI_MBC_str(1,5)/P_applied;...
    Node1038_SI_MBC_str(1,3)+0.25 Node1038_SI_MBC_str(1,5)/P_applied;Node1042_SI_MBC_str(1,3)+0.25 Node1042_SI_MBC_str(1,5)/P_applied;...
    Node1046_SI_MBC_str(1,3)+0.25 Node1046_SI_MBC_str(1,5)/P_applied;Node1050_SI_MBC_str(1,3)+0.25 Node1050_SI_MBC_str(1,5)/P_applied;...
    Node1054_SI_MBC_str(1,3)+0.25 Node1054_SI_MBC_str(1,5)/P_applied;Node1058_SI_MBC_str(1,3)+0.25 Node1058_SI_MBC_str(1,5)/P_applied;...
    Node1062_SI_MBC_str(1,3)+0.25 Node1062_SI_MBC_str(1,5)/P_applied;Node1066_SI_MBC_str(1,3)+0.25 Node1066_SI_MBC_str(1,5)/P_applied;...
    Node1070_SI_MBC_str(1,3)+0.25 Node1070_SI_MBC_str(1,5)/P_applied;Node1074_SI_MBC_str(1,3)+0.25 Node1074_SI_MBC_str(1,5)/P_applied;...
    Node1078_SI_MBC_str(1,3)+0.25 Node1078_SI_MBC_str(1,5)/P_applied;Node1082_SI_MBC_str(1,3)+0.25 Node1082_SI_MBC_str(1,5)/P_applied;...
    Node1086_SI_MBC_str(1,3)+0.25 Node1086_SI_MBC_str(1,5)/P_applied;Node1090_SI_MBC_str(1,3)+0.25 Node1090_SI_MBC_str(1,5)/P_applied;...
    Node1094_SI_MBC_str(1,3)+0.25 Node1094_SI_MBC_str(1,5)/P_applied;Node1098_SI_MBC_str(1,3)+0.25 Node1098_SI_MBC_str(1,5)/P_applied;...
    Node1102_SI_MBC_str(1,3)+0.25 Node1102_SI_MBC_str(1,5)/P_applied;Node1106_SI_MBC_str(1,3)+0.25 Node1106_SI_MBC_str(1,5)/P_applied;...
    Node1110_SI_MBC_str(1,3)+0.25 Node1110_SI_MBC_str(1,5)/P_applied;Node1114_SI_MBC_str(1,3)+0.25 Node1114_SI_MBC_str(1,5)/P_applied;...
    Node1118_SI_MBC_str(1,3)+0.25 Node1118_SI_MBC_str(1,5)/P_applied;Node1122_SI_MBC_str(1,3)+0.25 Node1122_SI_MBC_str(1,5)/P_applied;...
    Node1126_SI_MBC_str(1,3)+0.25 Node1126_SI_MBC_str(1,5)/P_applied;Node1130_SI_MBC_str(1,3)+0.25 Node1130_SI_MBC_str(1,5)/P_applied;...
    Node1134_SI_MBC_str(1,3)+0.25 Node1134_SI_MBC_str(1,5)/P_applied;Node1138_SI_MBC_str(1,3)+0.25 Node1138_SI_MBC_str(1,5)/P_applied;...
    Node1142_SI_MBC_str(1,3)+0.25 Node1142_SI_MBC_str(1,5)/P_applied;Node1146_SI_MBC_str(1,3)+0.25 Node1146_SI_MBC_str(1,5)/P_applied;...
    Node1150_SI_MBC_str(1,3)+0.25 Node1150_SI_MBC_str(1,5)/P_applied];


Result_phi11_intensity_NMBC_com = [Node1012_SI_NMBC_com(1,3)+0.25 Node1012_SI_NMBC_com(1,5)/P_applied;Node1011_SI_NMBC_com(1,3)+0.25 Node1011_SI_NMBC_com(1,5)/P_applied;...
    Node1018_SI_NMBC_com(1,3)+0.25 Node1018_SI_NMBC_com(1,5)/P_applied;Node1022_SI_NMBC_com(1,3)+0.25 Node1022_SI_NMBC_com(1,5)/P_applied;...
    Node1022_SI_NMBC_com(1,3)+0.25 Node1022_SI_NMBC_com(1,5)/P_applied;Node1026_SI_NMBC_com(1,3)+0.25 Node1026_SI_NMBC_com(1,5)/P_applied;...
    Node1030_SI_NMBC_com(1,3)+0.25 Node1030_SI_NMBC_com(1,5)/P_applied;Node1034_SI_NMBC_com(1,3)+0.25 Node1034_SI_NMBC_com(1,5)/P_applied;...
    Node1038_SI_NMBC_com(1,3)+0.25 Node1038_SI_NMBC_com(1,5)/P_applied;Node1042_SI_NMBC_com(1,3)+0.25 Node1042_SI_NMBC_com(1,5)/P_applied;...
    Node1046_SI_NMBC_com(1,3)+0.25 Node1046_SI_NMBC_com(1,5)/P_applied;Node1050_SI_NMBC_com(1,3)+0.25 Node1050_SI_NMBC_com(1,5)/P_applied;...
    Node1054_SI_NMBC_com(1,3)+0.25 Node1054_SI_NMBC_com(1,5)/P_applied;Node1058_SI_NMBC_com(1,3)+0.25 Node1058_SI_NMBC_com(1,5)/P_applied;...
    Node1062_SI_NMBC_com(1,3)+0.25 Node1062_SI_NMBC_com(1,5)/P_applied;Node1066_SI_NMBC_com(1,3)+0.25 Node1066_SI_NMBC_com(1,5)/P_applied;...
    Node1070_SI_NMBC_com(1,3)+0.25 Node1070_SI_NMBC_com(1,5)/P_applied;Node1074_SI_NMBC_com(1,3)+0.25 Node1074_SI_NMBC_com(1,5)/P_applied;...
    Node1078_SI_NMBC_com(1,3)+0.25 Node1078_SI_NMBC_com(1,5)/P_applied;Node1082_SI_NMBC_com(1,3)+0.25 Node1082_SI_NMBC_com(1,5)/P_applied;...
    Node1086_SI_NMBC_com(1,3)+0.25 Node1086_SI_NMBC_com(1,5)/P_applied;Node1090_SI_NMBC_com(1,3)+0.25 Node1090_SI_NMBC_com(1,5)/P_applied;...
    Node1094_SI_NMBC_com(1,3)+0.25 Node1094_SI_NMBC_com(1,5)/P_applied;Node1098_SI_NMBC_com(1,3)+0.25 Node1098_SI_NMBC_com(1,5)/P_applied;...
    Node1102_SI_NMBC_com(1,3)+0.25 Node1102_SI_NMBC_com(1,5)/P_applied;Node1106_SI_NMBC_com(1,3)+0.25 Node1106_SI_NMBC_com(1,5)/P_applied;...
    Node1110_SI_NMBC_com(1,3)+0.25 Node1110_SI_NMBC_com(1,5)/P_applied;Node1114_SI_NMBC_com(1,3)+0.25 Node1114_SI_NMBC_com(1,5)/P_applied;...
    Node1118_SI_NMBC_com(1,3)+0.25 Node1118_SI_NMBC_com(1,5)/P_applied;Node1122_SI_NMBC_com(1,3)+0.25 Node1122_SI_NMBC_com(1,5)/P_applied;...
    Node1126_SI_NMBC_com(1,3)+0.25 Node1126_SI_NMBC_com(1,5)/P_applied;Node1130_SI_NMBC_com(1,3)+0.25 Node1130_SI_NMBC_com(1,5)/P_applied;...
    Node1134_SI_NMBC_com(1,3)+0.25 Node1134_SI_NMBC_com(1,5)/P_applied;Node1138_SI_NMBC_com(1,3)+0.25 Node1138_SI_NMBC_com(1,5)/P_applied;...
    Node1142_SI_NMBC_com(1,3)+0.25 Node1142_SI_NMBC_com(1,5)/P_applied;Node1146_SI_NMBC_com(1,3)+0.25 Node1146_SI_NMBC_com(1,5)/P_applied;...
    Node1150_SI_NMBC_com(1,3)+0.25 Node1150_SI_NMBC_com(1,5)/P_applied];


Result_phi11_intensity_NMBC_rot = [Node1012_SI_NMBC_rot(1,3)+0.25 Node1012_SI_NMBC_rot(1,5)/P_applied;Node1011_SI_NMBC_rot(1,3)+0.25 Node1011_SI_NMBC_rot(1,5)/P_applied;...
    Node1018_SI_NMBC_rot(1,3)+0.25 Node1018_SI_NMBC_rot(1,5)/P_applied;Node1022_SI_NMBC_rot(1,3)+0.25 Node1022_SI_NMBC_rot(1,5)/P_applied;...
    Node1022_SI_NMBC_rot(1,3)+0.25 Node1022_SI_NMBC_rot(1,5)/P_applied;Node1026_SI_NMBC_rot(1,3)+0.25 Node1026_SI_NMBC_rot(1,5)/P_applied;...
    Node1030_SI_NMBC_rot(1,3)+0.25 Node1030_SI_NMBC_rot(1,5)/P_applied;Node1034_SI_NMBC_rot(1,3)+0.25 Node1034_SI_NMBC_rot(1,5)/P_applied;...
    Node1038_SI_NMBC_rot(1,3)+0.25 Node1038_SI_NMBC_rot(1,5)/P_applied;Node1042_SI_NMBC_rot(1,3)+0.25 Node1042_SI_NMBC_rot(1,5)/P_applied;...
    Node1046_SI_NMBC_rot(1,3)+0.25 Node1046_SI_NMBC_rot(1,5)/P_applied;Node1050_SI_NMBC_rot(1,3)+0.25 Node1050_SI_NMBC_rot(1,5)/P_applied;...
    Node1054_SI_NMBC_rot(1,3)+0.25 Node1054_SI_NMBC_rot(1,5)/P_applied;Node1058_SI_NMBC_rot(1,3)+0.25 Node1058_SI_NMBC_rot(1,5)/P_applied;...
    Node1062_SI_NMBC_rot(1,3)+0.25 Node1062_SI_NMBC_rot(1,5)/P_applied;Node1066_SI_NMBC_rot(1,3)+0.25 Node1066_SI_NMBC_rot(1,5)/P_applied;...
    Node1070_SI_NMBC_rot(1,3)+0.25 Node1070_SI_NMBC_rot(1,5)/P_applied;Node1074_SI_NMBC_rot(1,3)+0.25 Node1074_SI_NMBC_rot(1,5)/P_applied;...
    Node1078_SI_NMBC_rot(1,3)+0.25 Node1078_SI_NMBC_rot(1,5)/P_applied;Node1082_SI_NMBC_rot(1,3)+0.25 Node1082_SI_NMBC_rot(1,5)/P_applied;...
    Node1086_SI_NMBC_rot(1,3)+0.25 Node1086_SI_NMBC_rot(1,5)/P_applied;Node1090_SI_NMBC_rot(1,3)+0.25 Node1090_SI_NMBC_rot(1,5)/P_applied;...
    Node1094_SI_NMBC_rot(1,3)+0.25 Node1094_SI_NMBC_rot(1,5)/P_applied;Node1098_SI_NMBC_rot(1,3)+0.25 Node1098_SI_NMBC_rot(1,5)/P_applied;...
    Node1102_SI_NMBC_rot(1,3)+0.25 Node1102_SI_NMBC_rot(1,5)/P_applied;Node1106_SI_NMBC_rot(1,3)+0.25 Node1106_SI_NMBC_rot(1,5)/P_applied;...
    Node1110_SI_NMBC_rot(1,3)+0.25 Node1110_SI_NMBC_rot(1,5)/P_applied;Node1114_SI_NMBC_rot(1,3)+0.25 Node1114_SI_NMBC_rot(1,5)/P_applied;...
    Node1118_SI_NMBC_rot(1,3)+0.25 Node1118_SI_NMBC_rot(1,5)/P_applied;Node1122_SI_NMBC_rot(1,3)+0.25 Node1122_SI_NMBC_rot(1,5)/P_applied;...
    Node1126_SI_NMBC_rot(1,3)+0.25 Node1126_SI_NMBC_rot(1,5)/P_applied;Node1130_SI_NMBC_rot(1,3)+0.25 Node1130_SI_NMBC_rot(1,5)/P_applied;...
    Node1134_SI_NMBC_rot(1,3)+0.25 Node1134_SI_NMBC_rot(1,5)/P_applied;Node1138_SI_NMBC_rot(1,3)+0.25 Node1138_SI_NMBC_rot(1,5)/P_applied;...
    Node1142_SI_NMBC_rot(1,3)+0.25 Node1142_SI_NMBC_rot(1,5)/P_applied;Node1146_SI_NMBC_rot(1,3)+0.25 Node1146_SI_NMBC_rot(1,5)/P_applied;...
    Node1150_SI_NMBC_rot(1,3)+0.25 Node1150_SI_NMBC_rot(1,5)/P_applied];


Result_phi11_intensity_NMBC_str = [Node1012_SI_NMBC_str(1,3)+0.25 Node1012_SI_NMBC_str(1,5)/P_applied;Node1011_SI_NMBC_str(1,3)+0.25 Node1011_SI_NMBC_str(1,5)/P_applied;...
    Node1018_SI_NMBC_str(1,3)+0.25 Node1018_SI_NMBC_str(1,5)/P_applied;Node1022_SI_NMBC_str(1,3)+0.25 Node1022_SI_NMBC_str(1,5)/P_applied;...
    Node1022_SI_NMBC_str(1,3)+0.25 Node1022_SI_NMBC_str(1,5)/P_applied;Node1026_SI_NMBC_str(1,3)+0.25 Node1026_SI_NMBC_str(1,5)/P_applied;...
    Node1030_SI_NMBC_str(1,3)+0.25 Node1030_SI_NMBC_str(1,5)/P_applied;Node1034_SI_NMBC_str(1,3)+0.25 Node1034_SI_NMBC_str(1,5)/P_applied;...
    Node1038_SI_NMBC_str(1,3)+0.25 Node1038_SI_NMBC_str(1,5)/P_applied;Node1042_SI_NMBC_str(1,3)+0.25 Node1042_SI_NMBC_str(1,5)/P_applied;...
    Node1046_SI_NMBC_str(1,3)+0.25 Node1046_SI_NMBC_str(1,5)/P_applied;Node1050_SI_NMBC_str(1,3)+0.25 Node1050_SI_NMBC_str(1,5)/P_applied;...
    Node1054_SI_NMBC_str(1,3)+0.25 Node1054_SI_NMBC_str(1,5)/P_applied;Node1058_SI_NMBC_str(1,3)+0.25 Node1058_SI_NMBC_str(1,5)/P_applied;...
    Node1062_SI_NMBC_str(1,3)+0.25 Node1062_SI_NMBC_str(1,5)/P_applied;Node1066_SI_NMBC_str(1,3)+0.25 Node1066_SI_NMBC_str(1,5)/P_applied;...
    Node1070_SI_NMBC_str(1,3)+0.25 Node1070_SI_NMBC_str(1,5)/P_applied;Node1074_SI_NMBC_str(1,3)+0.25 Node1074_SI_NMBC_str(1,5)/P_applied;...
    Node1078_SI_NMBC_str(1,3)+0.25 Node1078_SI_NMBC_str(1,5)/P_applied;Node1082_SI_NMBC_str(1,3)+0.25 Node1082_SI_NMBC_str(1,5)/P_applied;...
    Node1086_SI_NMBC_str(1,3)+0.25 Node1086_SI_NMBC_str(1,5)/P_applied;Node1090_SI_NMBC_str(1,3)+0.25 Node1090_SI_NMBC_str(1,5)/P_applied;...
    Node1094_SI_NMBC_str(1,3)+0.25 Node1094_SI_NMBC_str(1,5)/P_applied;Node1098_SI_NMBC_str(1,3)+0.25 Node1098_SI_NMBC_str(1,5)/P_applied;...
    Node1102_SI_NMBC_str(1,3)+0.25 Node1102_SI_NMBC_str(1,5)/P_applied;Node1106_SI_NMBC_str(1,3)+0.25 Node1106_SI_NMBC_str(1,5)/P_applied;...
    Node1110_SI_NMBC_str(1,3)+0.25 Node1110_SI_NMBC_str(1,5)/P_applied;Node1114_SI_NMBC_str(1,3)+0.25 Node1114_SI_NMBC_str(1,5)/P_applied;...
    Node1118_SI_NMBC_str(1,3)+0.25 Node1118_SI_NMBC_str(1,5)/P_applied;Node1122_SI_NMBC_str(1,3)+0.25 Node1122_SI_NMBC_str(1,5)/P_applied;...
    Node1126_SI_NMBC_str(1,3)+0.25 Node1126_SI_NMBC_str(1,5)/P_applied;Node1130_SI_NMBC_str(1,3)+0.25 Node1130_SI_NMBC_str(1,5)/P_applied;...
    Node1134_SI_NMBC_str(1,3)+0.25 Node1134_SI_NMBC_str(1,5)/P_applied;Node1138_SI_NMBC_str(1,3)+0.25 Node1138_SI_NMBC_str(1,5)/P_applied;...
    Node1142_SI_NMBC_str(1,3)+0.25 Node1142_SI_NMBC_str(1,5)/P_applied;Node1146_SI_NMBC_str(1,3)+0.25 Node1146_SI_NMBC_str(1,5)/P_applied;...
    Node1150_SI_NMBC_str(1,3)+0.25 Node1150_SI_NMBC_str(1,5)/P_applied];



Result_phi22_intensity_MBC_com = [Node1012_SI_MBC_com(1,3)+0.25 Node1012_SI_MBC_com(1,6)/P_applied;Node1011_SI_MBC_com(1,3)+0.25 Node1011_SI_MBC_com(1,6)/P_applied;...
    Node1018_SI_MBC_com(1,3)+0.25 Node1018_SI_MBC_com(1,6)/P_applied;Node1022_SI_MBC_com(1,3)+0.25 Node1022_SI_MBC_com(1,6)/P_applied;...
    Node1022_SI_MBC_com(1,3)+0.25 Node1022_SI_MBC_com(1,6)/P_applied;Node1026_SI_MBC_com(1,3)+0.25 Node1026_SI_MBC_com(1,6)/P_applied;...
    Node1030_SI_MBC_com(1,3)+0.25 Node1030_SI_MBC_com(1,6)/P_applied;Node1034_SI_MBC_com(1,3)+0.25 Node1034_SI_MBC_com(1,6)/P_applied;...
    Node1038_SI_MBC_com(1,3)+0.25 Node1038_SI_MBC_com(1,6)/P_applied;Node1042_SI_MBC_com(1,3)+0.25 Node1042_SI_MBC_com(1,6)/P_applied;...
    Node1046_SI_MBC_com(1,3)+0.25 Node1046_SI_MBC_com(1,6)/P_applied;Node1050_SI_MBC_com(1,3)+0.25 Node1050_SI_MBC_com(1,6)/P_applied;...
    Node1054_SI_MBC_com(1,3)+0.25 Node1054_SI_MBC_com(1,6)/P_applied;Node1058_SI_MBC_com(1,3)+0.25 Node1058_SI_MBC_com(1,6)/P_applied;...
    Node1062_SI_MBC_com(1,3)+0.25 Node1062_SI_MBC_com(1,6)/P_applied;Node1066_SI_MBC_com(1,3)+0.25 Node1066_SI_MBC_com(1,6)/P_applied;...
    Node1070_SI_MBC_com(1,3)+0.25 Node1070_SI_MBC_com(1,6)/P_applied;Node1074_SI_MBC_com(1,3)+0.25 Node1074_SI_MBC_com(1,6)/P_applied;...
    Node1078_SI_MBC_com(1,3)+0.25 Node1078_SI_MBC_com(1,6)/P_applied;Node1082_SI_MBC_com(1,3)+0.25 Node1082_SI_MBC_com(1,6)/P_applied;...
    Node1086_SI_MBC_com(1,3)+0.25 Node1086_SI_MBC_com(1,6)/P_applied;Node1090_SI_MBC_com(1,3)+0.25 Node1090_SI_MBC_com(1,6)/P_applied;...
    Node1094_SI_MBC_com(1,3)+0.25 Node1094_SI_MBC_com(1,6)/P_applied;Node1098_SI_MBC_com(1,3)+0.25 Node1098_SI_MBC_com(1,6)/P_applied;...
    Node1102_SI_MBC_com(1,3)+0.25 Node1102_SI_MBC_com(1,6)/P_applied;Node1106_SI_MBC_com(1,3)+0.25 Node1106_SI_MBC_com(1,6)/P_applied;...
    Node1110_SI_MBC_com(1,3)+0.25 Node1110_SI_MBC_com(1,6)/P_applied;Node1114_SI_MBC_com(1,3)+0.25 Node1114_SI_MBC_com(1,6)/P_applied;...
    Node1118_SI_MBC_com(1,3)+0.25 Node1118_SI_MBC_com(1,6)/P_applied;Node1122_SI_MBC_com(1,3)+0.25 Node1122_SI_MBC_com(1,6)/P_applied;...
    Node1126_SI_MBC_com(1,3)+0.25 Node1126_SI_MBC_com(1,6)/P_applied;Node1130_SI_MBC_com(1,3)+0.25 Node1130_SI_MBC_com(1,6)/P_applied;...
    Node1134_SI_MBC_com(1,3)+0.25 Node1134_SI_MBC_com(1,6)/P_applied;Node1138_SI_MBC_com(1,3)+0.25 Node1138_SI_MBC_com(1,6)/P_applied;...
    Node1142_SI_MBC_com(1,3)+0.25 Node1142_SI_MBC_com(1,6)/P_applied;Node1146_SI_MBC_com(1,3)+0.25 Node1146_SI_MBC_com(1,6)/P_applied;...
    Node1150_SI_MBC_com(1,3)+0.25 Node1150_SI_MBC_com(1,6)/P_applied];


Result_phi22_intensity_MBC_rot = [Node1012_SI_MBC_rot(1,3)+0.25 Node1012_SI_MBC_rot(1,6)/P_applied;Node1011_SI_MBC_rot(1,3)+0.25 Node1011_SI_MBC_rot(1,6)/P_applied;...
    Node1018_SI_MBC_rot(1,3)+0.25 Node1018_SI_MBC_rot(1,6)/P_applied;Node1022_SI_MBC_rot(1,3)+0.25 Node1022_SI_MBC_rot(1,6)/P_applied;...
    Node1022_SI_MBC_rot(1,3)+0.25 Node1022_SI_MBC_rot(1,6)/P_applied;Node1026_SI_MBC_rot(1,3)+0.25 Node1026_SI_MBC_rot(1,6)/P_applied;...
    Node1030_SI_MBC_rot(1,3)+0.25 Node1030_SI_MBC_rot(1,6)/P_applied;Node1034_SI_MBC_rot(1,3)+0.25 Node1034_SI_MBC_rot(1,6)/P_applied;...
    Node1038_SI_MBC_rot(1,3)+0.25 Node1038_SI_MBC_rot(1,6)/P_applied;Node1042_SI_MBC_rot(1,3)+0.25 Node1042_SI_MBC_rot(1,6)/P_applied;...
    Node1046_SI_MBC_rot(1,3)+0.25 Node1046_SI_MBC_rot(1,6)/P_applied;Node1050_SI_MBC_rot(1,3)+0.25 Node1050_SI_MBC_rot(1,6)/P_applied;...
    Node1054_SI_MBC_rot(1,3)+0.25 Node1054_SI_MBC_rot(1,6)/P_applied;Node1058_SI_MBC_rot(1,3)+0.25 Node1058_SI_MBC_rot(1,6)/P_applied;...
    Node1062_SI_MBC_rot(1,3)+0.25 Node1062_SI_MBC_rot(1,6)/P_applied;Node1066_SI_MBC_rot(1,3)+0.25 Node1066_SI_MBC_rot(1,6)/P_applied;...
    Node1070_SI_MBC_rot(1,3)+0.25 Node1070_SI_MBC_rot(1,6)/P_applied;Node1074_SI_MBC_rot(1,3)+0.25 Node1074_SI_MBC_rot(1,6)/P_applied;...
    Node1078_SI_MBC_rot(1,3)+0.25 Node1078_SI_MBC_rot(1,6)/P_applied;Node1082_SI_MBC_rot(1,3)+0.25 Node1082_SI_MBC_rot(1,6)/P_applied;...
    Node1086_SI_MBC_rot(1,3)+0.25 Node1086_SI_MBC_rot(1,6)/P_applied;Node1090_SI_MBC_rot(1,3)+0.25 Node1090_SI_MBC_rot(1,6)/P_applied;...
    Node1094_SI_MBC_rot(1,3)+0.25 Node1094_SI_MBC_rot(1,6)/P_applied;Node1098_SI_MBC_rot(1,3)+0.25 Node1098_SI_MBC_rot(1,6)/P_applied;...
    Node1102_SI_MBC_rot(1,3)+0.25 Node1102_SI_MBC_rot(1,6)/P_applied;Node1106_SI_MBC_rot(1,3)+0.25 Node1106_SI_MBC_rot(1,6)/P_applied;...
    Node1110_SI_MBC_rot(1,3)+0.25 Node1110_SI_MBC_rot(1,6)/P_applied;Node1114_SI_MBC_rot(1,3)+0.25 Node1114_SI_MBC_rot(1,6)/P_applied;...
    Node1118_SI_MBC_rot(1,3)+0.25 Node1118_SI_MBC_rot(1,6)/P_applied;Node1122_SI_MBC_rot(1,3)+0.25 Node1122_SI_MBC_rot(1,6)/P_applied;...
    Node1126_SI_MBC_rot(1,3)+0.25 Node1126_SI_MBC_rot(1,6)/P_applied;Node1130_SI_MBC_rot(1,3)+0.25 Node1130_SI_MBC_rot(1,6)/P_applied;...
    Node1134_SI_MBC_rot(1,3)+0.25 Node1134_SI_MBC_rot(1,6)/P_applied;Node1138_SI_MBC_rot(1,3)+0.25 Node1138_SI_MBC_rot(1,6)/P_applied;...
    Node1142_SI_MBC_rot(1,3)+0.25 Node1142_SI_MBC_rot(1,6)/P_applied;Node1146_SI_MBC_rot(1,3)+0.25 Node1146_SI_MBC_rot(1,6)/P_applied;...
    Node1150_SI_MBC_rot(1,3)+0.25 Node1150_SI_MBC_rot(1,6)/P_applied];


Result_phi22_intensity_MBC_str = [Node1012_SI_MBC_str(1,3)+0.25 Node1012_SI_MBC_str(1,6)/P_applied;Node1011_SI_MBC_str(1,3)+0.25 Node1011_SI_MBC_str(1,6)/P_applied;...
    Node1018_SI_MBC_str(1,3)+0.25 Node1018_SI_MBC_str(1,6)/P_applied;Node1022_SI_MBC_str(1,3)+0.25 Node1022_SI_MBC_str(1,6)/P_applied;...
    Node1022_SI_MBC_str(1,3)+0.25 Node1022_SI_MBC_str(1,6)/P_applied;Node1026_SI_MBC_str(1,3)+0.25 Node1026_SI_MBC_str(1,6)/P_applied;...
    Node1030_SI_MBC_str(1,3)+0.25 Node1030_SI_MBC_str(1,6)/P_applied;Node1034_SI_MBC_str(1,3)+0.25 Node1034_SI_MBC_str(1,6)/P_applied;...
    Node1038_SI_MBC_str(1,3)+0.25 Node1038_SI_MBC_str(1,6)/P_applied;Node1042_SI_MBC_str(1,3)+0.25 Node1042_SI_MBC_str(1,6)/P_applied;...
    Node1046_SI_MBC_str(1,3)+0.25 Node1046_SI_MBC_str(1,6)/P_applied;Node1050_SI_MBC_str(1,3)+0.25 Node1050_SI_MBC_str(1,6)/P_applied;...
    Node1054_SI_MBC_str(1,3)+0.25 Node1054_SI_MBC_str(1,6)/P_applied;Node1058_SI_MBC_str(1,3)+0.25 Node1058_SI_MBC_str(1,6)/P_applied;...
    Node1062_SI_MBC_str(1,3)+0.25 Node1062_SI_MBC_str(1,6)/P_applied;Node1066_SI_MBC_str(1,3)+0.25 Node1066_SI_MBC_str(1,6)/P_applied;...
    Node1070_SI_MBC_str(1,3)+0.25 Node1070_SI_MBC_str(1,6)/P_applied;Node1074_SI_MBC_str(1,3)+0.25 Node1074_SI_MBC_str(1,6)/P_applied;...
    Node1078_SI_MBC_str(1,3)+0.25 Node1078_SI_MBC_str(1,6)/P_applied;Node1082_SI_MBC_str(1,3)+0.25 Node1082_SI_MBC_str(1,6)/P_applied;...
    Node1086_SI_MBC_str(1,3)+0.25 Node1086_SI_MBC_str(1,6)/P_applied;Node1090_SI_MBC_str(1,3)+0.25 Node1090_SI_MBC_str(1,6)/P_applied;...
    Node1094_SI_MBC_str(1,3)+0.25 Node1094_SI_MBC_str(1,6)/P_applied;Node1098_SI_MBC_str(1,3)+0.25 Node1098_SI_MBC_str(1,6)/P_applied;...
    Node1102_SI_MBC_str(1,3)+0.25 Node1102_SI_MBC_str(1,6)/P_applied;Node1106_SI_MBC_str(1,3)+0.25 Node1106_SI_MBC_str(1,6)/P_applied;...
    Node1110_SI_MBC_str(1,3)+0.25 Node1110_SI_MBC_str(1,6)/P_applied;Node1114_SI_MBC_str(1,3)+0.25 Node1114_SI_MBC_str(1,6)/P_applied;...
    Node1118_SI_MBC_str(1,3)+0.25 Node1118_SI_MBC_str(1,6)/P_applied;Node1122_SI_MBC_str(1,3)+0.25 Node1122_SI_MBC_str(1,6)/P_applied;...
    Node1126_SI_MBC_str(1,3)+0.25 Node1126_SI_MBC_str(1,6)/P_applied;Node1130_SI_MBC_str(1,3)+0.25 Node1130_SI_MBC_str(1,6)/P_applied;...
    Node1134_SI_MBC_str(1,3)+0.25 Node1134_SI_MBC_str(1,6)/P_applied;Node1138_SI_MBC_str(1,3)+0.25 Node1138_SI_MBC_str(1,6)/P_applied;...
    Node1142_SI_MBC_str(1,3)+0.25 Node1142_SI_MBC_str(1,6)/P_applied;Node1146_SI_MBC_str(1,3)+0.25 Node1146_SI_MBC_str(1,6)/P_applied;...
    Node1150_SI_MBC_str(1,3)+0.25 Node1150_SI_MBC_str(1,6)/P_applied];


Result_phi22_intensity_NMBC_com = [Node1012_SI_NMBC_com(1,3)+0.25 Node1012_SI_NMBC_com(1,6)/P_applied;Node1011_SI_NMBC_com(1,3)+0.25 Node1011_SI_NMBC_com(1,6)/P_applied;...
    Node1018_SI_NMBC_com(1,3)+0.25 Node1018_SI_NMBC_com(1,6)/P_applied;Node1022_SI_NMBC_com(1,3)+0.25 Node1022_SI_NMBC_com(1,6)/P_applied;...
    Node1022_SI_NMBC_com(1,3)+0.25 Node1022_SI_NMBC_com(1,6)/P_applied;Node1026_SI_NMBC_com(1,3)+0.25 Node1026_SI_NMBC_com(1,6)/P_applied;...
    Node1030_SI_NMBC_com(1,3)+0.25 Node1030_SI_NMBC_com(1,6)/P_applied;Node1034_SI_NMBC_com(1,3)+0.25 Node1034_SI_NMBC_com(1,6)/P_applied;...
    Node1038_SI_NMBC_com(1,3)+0.25 Node1038_SI_NMBC_com(1,6)/P_applied;Node1042_SI_NMBC_com(1,3)+0.25 Node1042_SI_NMBC_com(1,6)/P_applied;...
    Node1046_SI_NMBC_com(1,3)+0.25 Node1046_SI_NMBC_com(1,6)/P_applied;Node1050_SI_NMBC_com(1,3)+0.25 Node1050_SI_NMBC_com(1,6)/P_applied;...
    Node1054_SI_NMBC_com(1,3)+0.25 Node1054_SI_NMBC_com(1,6)/P_applied;Node1058_SI_NMBC_com(1,3)+0.25 Node1058_SI_NMBC_com(1,6)/P_applied;...
    Node1062_SI_NMBC_com(1,3)+0.25 Node1062_SI_NMBC_com(1,6)/P_applied;Node1066_SI_NMBC_com(1,3)+0.25 Node1066_SI_NMBC_com(1,6)/P_applied;...
    Node1070_SI_NMBC_com(1,3)+0.25 Node1070_SI_NMBC_com(1,6)/P_applied;Node1074_SI_NMBC_com(1,3)+0.25 Node1074_SI_NMBC_com(1,6)/P_applied;...
    Node1078_SI_NMBC_com(1,3)+0.25 Node1078_SI_NMBC_com(1,6)/P_applied;Node1082_SI_NMBC_com(1,3)+0.25 Node1082_SI_NMBC_com(1,6)/P_applied;...
    Node1086_SI_NMBC_com(1,3)+0.25 Node1086_SI_NMBC_com(1,6)/P_applied;Node1090_SI_NMBC_com(1,3)+0.25 Node1090_SI_NMBC_com(1,6)/P_applied;...
    Node1094_SI_NMBC_com(1,3)+0.25 Node1094_SI_NMBC_com(1,6)/P_applied;Node1098_SI_NMBC_com(1,3)+0.25 Node1098_SI_NMBC_com(1,6)/P_applied;...
    Node1102_SI_NMBC_com(1,3)+0.25 Node1102_SI_NMBC_com(1,6)/P_applied;Node1106_SI_NMBC_com(1,3)+0.25 Node1106_SI_NMBC_com(1,6)/P_applied;...
    Node1110_SI_NMBC_com(1,3)+0.25 Node1110_SI_NMBC_com(1,6)/P_applied;Node1114_SI_NMBC_com(1,3)+0.25 Node1114_SI_NMBC_com(1,6)/P_applied;...
    Node1118_SI_NMBC_com(1,3)+0.25 Node1118_SI_NMBC_com(1,6)/P_applied;Node1122_SI_NMBC_com(1,3)+0.25 Node1122_SI_NMBC_com(1,6)/P_applied;...
    Node1126_SI_NMBC_com(1,3)+0.25 Node1126_SI_NMBC_com(1,6)/P_applied;Node1130_SI_NMBC_com(1,3)+0.25 Node1130_SI_NMBC_com(1,6)/P_applied;...
    Node1134_SI_NMBC_com(1,3)+0.25 Node1134_SI_NMBC_com(1,6)/P_applied;Node1138_SI_NMBC_com(1,3)+0.25 Node1138_SI_NMBC_com(1,6)/P_applied;...
    Node1142_SI_NMBC_com(1,3)+0.25 Node1142_SI_NMBC_com(1,6)/P_applied;Node1146_SI_NMBC_com(1,3)+0.25 Node1146_SI_NMBC_com(1,6)/P_applied;...
    Node1150_SI_NMBC_com(1,3)+0.25 Node1150_SI_NMBC_com(1,6)/P_applied];


Result_phi22_intensity_NMBC_rot = [Node1012_SI_NMBC_rot(1,3)+0.25 Node1012_SI_NMBC_rot(1,6)/P_applied;Node1011_SI_NMBC_rot(1,3)+0.25 Node1011_SI_NMBC_rot(1,6)/P_applied;...
    Node1018_SI_NMBC_rot(1,3)+0.25 Node1018_SI_NMBC_rot(1,6)/P_applied;Node1022_SI_NMBC_rot(1,3)+0.25 Node1022_SI_NMBC_rot(1,6)/P_applied;...
    Node1022_SI_NMBC_rot(1,3)+0.25 Node1022_SI_NMBC_rot(1,6)/P_applied;Node1026_SI_NMBC_rot(1,3)+0.25 Node1026_SI_NMBC_rot(1,6)/P_applied;...
    Node1030_SI_NMBC_rot(1,3)+0.25 Node1030_SI_NMBC_rot(1,6)/P_applied;Node1034_SI_NMBC_rot(1,3)+0.25 Node1034_SI_NMBC_rot(1,6)/P_applied;...
    Node1038_SI_NMBC_rot(1,3)+0.25 Node1038_SI_NMBC_rot(1,6)/P_applied;Node1042_SI_NMBC_rot(1,3)+0.25 Node1042_SI_NMBC_rot(1,6)/P_applied;...
    Node1046_SI_NMBC_rot(1,3)+0.25 Node1046_SI_NMBC_rot(1,6)/P_applied;Node1050_SI_NMBC_rot(1,3)+0.25 Node1050_SI_NMBC_rot(1,6)/P_applied;...
    Node1054_SI_NMBC_rot(1,3)+0.25 Node1054_SI_NMBC_rot(1,6)/P_applied;Node1058_SI_NMBC_rot(1,3)+0.25 Node1058_SI_NMBC_rot(1,6)/P_applied;...
    Node1062_SI_NMBC_rot(1,3)+0.25 Node1062_SI_NMBC_rot(1,6)/P_applied;Node1066_SI_NMBC_rot(1,3)+0.25 Node1066_SI_NMBC_rot(1,6)/P_applied;...
    Node1070_SI_NMBC_rot(1,3)+0.25 Node1070_SI_NMBC_rot(1,6)/P_applied;Node1074_SI_NMBC_rot(1,3)+0.25 Node1074_SI_NMBC_rot(1,6)/P_applied;...
    Node1078_SI_NMBC_rot(1,3)+0.25 Node1078_SI_NMBC_rot(1,6)/P_applied;Node1082_SI_NMBC_rot(1,3)+0.25 Node1082_SI_NMBC_rot(1,6)/P_applied;...
    Node1086_SI_NMBC_rot(1,3)+0.25 Node1086_SI_NMBC_rot(1,6)/P_applied;Node1090_SI_NMBC_rot(1,3)+0.25 Node1090_SI_NMBC_rot(1,6)/P_applied;...
    Node1094_SI_NMBC_rot(1,3)+0.25 Node1094_SI_NMBC_rot(1,6)/P_applied;Node1098_SI_NMBC_rot(1,3)+0.25 Node1098_SI_NMBC_rot(1,6)/P_applied;...
    Node1102_SI_NMBC_rot(1,3)+0.25 Node1102_SI_NMBC_rot(1,6)/P_applied;Node1106_SI_NMBC_rot(1,3)+0.25 Node1106_SI_NMBC_rot(1,6)/P_applied;...
    Node1110_SI_NMBC_rot(1,3)+0.25 Node1110_SI_NMBC_rot(1,6)/P_applied;Node1114_SI_NMBC_rot(1,3)+0.25 Node1114_SI_NMBC_rot(1,6)/P_applied;...
    Node1118_SI_NMBC_rot(1,3)+0.25 Node1118_SI_NMBC_rot(1,6)/P_applied;Node1122_SI_NMBC_rot(1,3)+0.25 Node1122_SI_NMBC_rot(1,6)/P_applied;...
    Node1126_SI_NMBC_rot(1,3)+0.25 Node1126_SI_NMBC_rot(1,6)/P_applied;Node1130_SI_NMBC_rot(1,3)+0.25 Node1130_SI_NMBC_rot(1,6)/P_applied;...
    Node1134_SI_NMBC_rot(1,3)+0.25 Node1134_SI_NMBC_rot(1,6)/P_applied;Node1138_SI_NMBC_rot(1,3)+0.25 Node1138_SI_NMBC_rot(1,6)/P_applied;...
    Node1142_SI_NMBC_rot(1,3)+0.25 Node1142_SI_NMBC_rot(1,6)/P_applied;Node1146_SI_NMBC_rot(1,3)+0.25 Node1146_SI_NMBC_rot(1,6)/P_applied;...
    Node1150_SI_NMBC_rot(1,3)+0.25 Node1150_SI_NMBC_rot(1,6)/P_applied];


Result_phi22_intensity_NMBC_str = [Node1012_SI_NMBC_str(1,3)+0.25 Node1012_SI_NMBC_str(1,6)/P_applied;Node1011_SI_NMBC_str(1,3)+0.25 Node1011_SI_NMBC_str(1,6)/P_applied;...
    Node1018_SI_NMBC_str(1,3)+0.25 Node1018_SI_NMBC_str(1,6)/P_applied;Node1022_SI_NMBC_str(1,3)+0.25 Node1022_SI_NMBC_str(1,6)/P_applied;...
    Node1022_SI_NMBC_str(1,3)+0.25 Node1022_SI_NMBC_str(1,6)/P_applied;Node1026_SI_NMBC_str(1,3)+0.25 Node1026_SI_NMBC_str(1,6)/P_applied;...
    Node1030_SI_NMBC_str(1,3)+0.25 Node1030_SI_NMBC_str(1,6)/P_applied;Node1034_SI_NMBC_str(1,3)+0.25 Node1034_SI_NMBC_str(1,6)/P_applied;...
    Node1038_SI_NMBC_str(1,3)+0.25 Node1038_SI_NMBC_str(1,6)/P_applied;Node1042_SI_NMBC_str(1,3)+0.25 Node1042_SI_NMBC_str(1,6)/P_applied;...
    Node1046_SI_NMBC_str(1,3)+0.25 Node1046_SI_NMBC_str(1,6)/P_applied;Node1050_SI_NMBC_str(1,3)+0.25 Node1050_SI_NMBC_str(1,6)/P_applied;...
    Node1054_SI_NMBC_str(1,3)+0.25 Node1054_SI_NMBC_str(1,6)/P_applied;Node1058_SI_NMBC_str(1,3)+0.25 Node1058_SI_NMBC_str(1,6)/P_applied;...
    Node1062_SI_NMBC_str(1,3)+0.25 Node1062_SI_NMBC_str(1,6)/P_applied;Node1066_SI_NMBC_str(1,3)+0.25 Node1066_SI_NMBC_str(1,6)/P_applied;...
    Node1070_SI_NMBC_str(1,3)+0.25 Node1070_SI_NMBC_str(1,6)/P_applied;Node1074_SI_NMBC_str(1,3)+0.25 Node1074_SI_NMBC_str(1,6)/P_applied;...
    Node1078_SI_NMBC_str(1,3)+0.25 Node1078_SI_NMBC_str(1,6)/P_applied;Node1082_SI_NMBC_str(1,3)+0.25 Node1082_SI_NMBC_str(1,6)/P_applied;...
    Node1086_SI_NMBC_str(1,3)+0.25 Node1086_SI_NMBC_str(1,6)/P_applied;Node1090_SI_NMBC_str(1,3)+0.25 Node1090_SI_NMBC_str(1,6)/P_applied;...
    Node1094_SI_NMBC_str(1,3)+0.25 Node1094_SI_NMBC_str(1,6)/P_applied;Node1098_SI_NMBC_str(1,3)+0.25 Node1098_SI_NMBC_str(1,6)/P_applied;...
    Node1102_SI_NMBC_str(1,3)+0.25 Node1102_SI_NMBC_str(1,6)/P_applied;Node1106_SI_NMBC_str(1,3)+0.25 Node1106_SI_NMBC_str(1,6)/P_applied;...
    Node1110_SI_NMBC_str(1,3)+0.25 Node1110_SI_NMBC_str(1,6)/P_applied;Node1114_SI_NMBC_str(1,3)+0.25 Node1114_SI_NMBC_str(1,6)/P_applied;...
    Node1118_SI_NMBC_str(1,3)+0.25 Node1118_SI_NMBC_str(1,6)/P_applied;Node1122_SI_NMBC_str(1,3)+0.25 Node1122_SI_NMBC_str(1,6)/P_applied;...
    Node1126_SI_NMBC_str(1,3)+0.25 Node1126_SI_NMBC_str(1,6)/P_applied;Node1130_SI_NMBC_str(1,3)+0.25 Node1130_SI_NMBC_str(1,6)/P_applied;...
    Node1134_SI_NMBC_str(1,3)+0.25 Node1134_SI_NMBC_str(1,6)/P_applied;Node1138_SI_NMBC_str(1,3)+0.25 Node1138_SI_NMBC_str(1,6)/P_applied;...
    Node1142_SI_NMBC_str(1,3)+0.25 Node1142_SI_NMBC_str(1,6)/P_applied;Node1146_SI_NMBC_str(1,3)+0.25 Node1146_SI_NMBC_str(1,6)/P_applied;...
    Node1150_SI_NMBC_str(1,3)+0.25 Node1150_SI_NMBC_str(1,6)/P_applied];



Result_phi_21_MBC_com = [r_938/0.05 Node938_MBC_com(1,5);r_939/0.05  Node939_MBC_com(1,5);...
    r_941/0.05  Node941_MBC_com(1,5);r_943/0.05  Node943_MBC_com(1,5);r_945/0.05  Node945_MBC_com(1,5);...
    r_947/0.05  Node947_MBC_com(1,5);r_949/0.05  Node949_MBC_com(1,5);r_951/0.05  Node951_MBC_com(1,5);...
    r_953/0.05  Node953_MBC_com(1,5);r_955/0.05  Node955_MBC_com(1,5);r_957/0.05  Node957_MBC_com(1,5);...
    r_959/0.05  Node959_MBC_com(1,5);r_961/0.05  Node961_MBC_com(1,5);r_963/0.05  Node963_MBC_com(1,5);...
    r_965/0.05  Node965_MBC_com(1,5);r_967/0.05  Node967_MBC_com(1,5);r_969/0.05  Node969_MBC_com(1,5);...
    r_971/0.05  Node971_MBC_com(1,5);r_973/0.05  Node973_MBC_com(1,5);r_975/0.05  Node975_MBC_com(1,5);...
    r_977/0.05  Node977_MBC_com(1,5);r_979/0.05  Node979_MBC_com(1,5);r_981/0.05  Node981_MBC_com(1,5);...
    r_983/0.05  Node983_MBC_com(1,5);r_985/0.05  Node985_MBC_com(1,5);r_987/0.05  Node987_MBC_com(1,5);...
    r_989/0.05  Node989_MBC_com(1,5);r_991/0.05  Node991_MBC_com(1,5);r_993/0.05  Node993_MBC_com(1,5);...
    r_995/0.05  Node995_MBC_com(1,5);r_997/0.05  Node997_MBC_com(1,5);r_999/0.05  Node999_MBC_com(1,5);...
    r_1001/0.05  Node1001_MBC_com(1,5);r_1003/0.05  Node1003_MBC_com(1,5);r_1005/0.05  Node1005_MBC_com(1,5);...
    r_1007/0.05  Node1007_MBC_com(1,5)];

Result_phi_21_MBC_str = [r_938/0.05 Node938_MBC_str(1,5);r_939/0.05  Node939_MBC_str(1,5);...
    r_941/0.05  Node941_MBC_str(1,5);r_943/0.05  Node943_MBC_str(1,5);r_945/0.05  Node945_MBC_str(1,5);...
    r_947/0.05  Node947_MBC_str(1,5);r_949/0.05  Node949_MBC_str(1,5);r_951/0.05  Node951_MBC_str(1,5);...
    r_953/0.05  Node953_MBC_str(1,5);r_955/0.05  Node955_MBC_str(1,5);r_957/0.05  Node957_MBC_str(1,5);...
    r_959/0.05  Node959_MBC_str(1,5);r_961/0.05  Node961_MBC_str(1,5);r_963/0.05  Node963_MBC_str(1,5);...
    r_965/0.05  Node965_MBC_str(1,5);r_967/0.05  Node967_MBC_str(1,5);r_969/0.05  Node969_MBC_str(1,5);...
    r_971/0.05  Node971_MBC_str(1,5);r_973/0.05  Node973_MBC_str(1,5);r_975/0.05  Node975_MBC_str(1,5);...
    r_977/0.05  Node977_MBC_str(1,5);r_979/0.05  Node979_MBC_str(1,5);r_981/0.05  Node981_MBC_str(1,5);...
    r_983/0.05  Node983_MBC_str(1,5);r_985/0.05  Node985_MBC_str(1,5);r_987/0.05  Node987_MBC_str(1,5);...
    r_989/0.05  Node989_MBC_str(1,5);r_991/0.05  Node991_MBC_str(1,5);r_993/0.05  Node993_MBC_str(1,5);...
    r_995/0.05  Node995_MBC_str(1,5);r_997/0.05  Node997_MBC_str(1,5);r_999/0.05  Node999_MBC_str(1,5);...
    r_1001/0.05  Node1001_MBC_str(1,5);r_1003/0.05  Node1003_MBC_str(1,5);r_1005/0.05  Node1005_MBC_str(1,5);...
    r_1007/0.05  Node1007_MBC_str(1,5)];

Result_phi_21_MBC_rot = [r_938/0.05 Node938_MBC_rot(1,5);r_939/0.05  Node939_MBC_rot(1,5);...
    r_941/0.05  Node941_MBC_rot(1,5);r_943/0.05  Node943_MBC_rot(1,5);r_945/0.05  Node945_MBC_rot(1,5);...
    r_947/0.05  Node947_MBC_rot(1,5);r_949/0.05  Node949_MBC_rot(1,5);r_951/0.05  Node951_MBC_rot(1,5);...
    r_953/0.05  Node953_MBC_rot(1,5);r_955/0.05  Node955_MBC_rot(1,5);r_957/0.05  Node957_MBC_rot(1,5);...
    r_959/0.05  Node959_MBC_rot(1,5);r_961/0.05  Node961_MBC_rot(1,5);r_963/0.05  Node963_MBC_rot(1,5);...
    r_965/0.05  Node965_MBC_rot(1,5);r_967/0.05  Node967_MBC_rot(1,5);r_969/0.05  Node969_MBC_rot(1,5);...
    r_971/0.05  Node971_MBC_rot(1,5);r_973/0.05  Node973_MBC_rot(1,5);r_975/0.05  Node975_MBC_rot(1,5);...
    r_977/0.05  Node977_MBC_rot(1,5);r_979/0.05  Node979_MBC_rot(1,5);r_981/0.05  Node981_MBC_rot(1,5);...
    r_983/0.05  Node983_MBC_rot(1,5);r_985/0.05  Node985_MBC_rot(1,5);r_987/0.05  Node987_MBC_rot(1,5);...
    r_989/0.05  Node989_MBC_rot(1,5);r_991/0.05  Node991_MBC_rot(1,5);r_993/0.05  Node993_MBC_rot(1,5);...
    r_995/0.05  Node995_MBC_rot(1,5);r_997/0.05  Node997_MBC_rot(1,5);r_999/0.05  Node999_MBC_rot(1,5);...
    r_1001/0.05  Node1001_MBC_rot(1,5);r_1003/0.05  Node1003_MBC_rot(1,5);r_1005/0.05  Node1005_MBC_rot(1,5);...
    r_1007/0.05  Node1007_MBC_rot(1,5)];


Result_phi_21_NMBC_com = [r_938/0.05 Node938_NMBC_com(1,5);r_939/0.05  Node939_NMBC_com(1,5);...
    r_941/0.05  Node941_NMBC_com(1,5);r_943/0.05  Node943_NMBC_com(1,5);r_945/0.05  Node945_NMBC_com(1,5);...
    r_947/0.05  Node947_NMBC_com(1,5);r_949/0.05  Node949_NMBC_com(1,5);r_951/0.05  Node951_NMBC_com(1,5);...
    r_953/0.05  Node953_NMBC_com(1,5);r_955/0.05  Node955_NMBC_com(1,5);r_957/0.05  Node957_NMBC_com(1,5);...
    r_959/0.05  Node959_NMBC_com(1,5);r_961/0.05  Node961_NMBC_com(1,5);r_963/0.05  Node963_NMBC_com(1,5);...
    r_965/0.05  Node965_NMBC_com(1,5);r_967/0.05  Node967_NMBC_com(1,5);r_969/0.05  Node969_NMBC_com(1,5);...
    r_971/0.05  Node971_NMBC_com(1,5);r_973/0.05  Node973_NMBC_com(1,5);r_975/0.05  Node975_NMBC_com(1,5);...
    r_977/0.05  Node977_NMBC_com(1,5);r_979/0.05  Node979_NMBC_com(1,5);r_981/0.05  Node981_NMBC_com(1,5);...
    r_983/0.05  Node983_NMBC_com(1,5);r_985/0.05  Node985_NMBC_com(1,5);r_987/0.05  Node987_NMBC_com(1,5);...
    r_989/0.05  Node989_NMBC_com(1,5);r_991/0.05  Node991_NMBC_com(1,5);r_993/0.05  Node993_NMBC_com(1,5);...
    r_995/0.05  Node995_NMBC_com(1,5);r_997/0.05  Node997_NMBC_com(1,5);r_999/0.05  Node999_NMBC_com(1,5);...
    r_1001/0.05  Node1001_NMBC_com(1,5);r_1003/0.05  Node1003_NMBC_com(1,5);r_1005/0.05  Node1005_NMBC_com(1,5);...
    r_1007/0.05  Node1007_NMBC_com(1,5)];

Result_phi_21_NMBC_str = [r_938/0.05 Node938_NMBC_str(1,5);r_939/0.05  Node939_NMBC_str(1,5);...
    r_941/0.05  Node941_NMBC_str(1,5);r_943/0.05  Node943_NMBC_str(1,5);r_945/0.05  Node945_NMBC_str(1,5);...
    r_947/0.05  Node947_NMBC_str(1,5);r_949/0.05  Node949_NMBC_str(1,5);r_951/0.05  Node951_NMBC_str(1,5);...
    r_953/0.05  Node953_NMBC_str(1,5);r_955/0.05  Node955_NMBC_str(1,5);r_957/0.05  Node957_NMBC_str(1,5);...
    r_959/0.05  Node959_NMBC_str(1,5);r_961/0.05  Node961_NMBC_str(1,5);r_963/0.05  Node963_NMBC_str(1,5);...
    r_965/0.05  Node965_NMBC_str(1,5);r_967/0.05  Node967_NMBC_str(1,5);r_969/0.05  Node969_NMBC_str(1,5);...
    r_971/0.05  Node971_NMBC_str(1,5);r_973/0.05  Node973_NMBC_str(1,5);r_975/0.05  Node975_NMBC_str(1,5);...
    r_977/0.05  Node977_NMBC_str(1,5);r_979/0.05  Node979_NMBC_str(1,5);r_981/0.05  Node981_NMBC_str(1,5);...
    r_983/0.05  Node983_NMBC_str(1,5);r_985/0.05  Node985_NMBC_str(1,5);r_987/0.05  Node987_NMBC_str(1,5);...
    r_989/0.05  Node989_NMBC_str(1,5);r_991/0.05  Node991_NMBC_str(1,5);r_993/0.05  Node993_NMBC_str(1,5);...
    r_995/0.05  Node995_NMBC_str(1,5);r_997/0.05  Node997_NMBC_str(1,5);r_999/0.05  Node999_NMBC_str(1,5);...
    r_1001/0.05  Node1001_NMBC_str(1,5);r_1003/0.05  Node1003_NMBC_str(1,5);r_1005/0.05  Node1005_NMBC_str(1,5);...
    r_1007/0.05  Node1007_NMBC_str(1,5)];

Result_phi_21_NMBC_rot = [r_938/0.05 Node938_NMBC_rot(1,5);r_939/0.05  Node939_NMBC_rot(1,5);...
    r_941/0.05  Node941_NMBC_rot(1,5);r_943/0.05  Node943_NMBC_rot(1,5);r_945/0.05  Node945_NMBC_rot(1,5);...
    r_947/0.05  Node947_NMBC_rot(1,5);r_949/0.05  Node949_NMBC_rot(1,5);r_951/0.05  Node951_NMBC_rot(1,5);...
    r_953/0.05  Node953_NMBC_rot(1,5);r_955/0.05  Node955_NMBC_rot(1,5);r_957/0.05  Node957_NMBC_rot(1,5);...
    r_959/0.05  Node959_NMBC_rot(1,5);r_961/0.05  Node961_NMBC_rot(1,5);r_963/0.05  Node963_NMBC_rot(1,5);...
    r_965/0.05  Node965_NMBC_rot(1,5);r_967/0.05  Node967_NMBC_rot(1,5);r_969/0.05  Node969_NMBC_rot(1,5);...
    r_971/0.05  Node971_NMBC_rot(1,5);r_973/0.05  Node973_NMBC_rot(1,5);r_975/0.05  Node975_NMBC_rot(1,5);...
    r_977/0.05  Node977_NMBC_rot(1,5);r_979/0.05  Node979_NMBC_rot(1,5);r_981/0.05  Node981_NMBC_rot(1,5);...
    r_983/0.05  Node983_NMBC_rot(1,5);r_985/0.05  Node985_NMBC_rot(1,5);r_987/0.05  Node987_NMBC_rot(1,5);...
    r_989/0.05  Node989_NMBC_rot(1,5);r_991/0.05  Node991_NMBC_rot(1,5);r_993/0.05  Node993_NMBC_rot(1,5);...
    r_995/0.05  Node995_NMBC_rot(1,5);r_997/0.05  Node997_NMBC_rot(1,5);r_999/0.05  Node999_NMBC_rot(1,5);...
    r_1001/0.05  Node1001_NMBC_rot(1,5);r_1003/0.05  Node1003_NMBC_rot(1,5);r_1005/0.05  Node1005_NMBC_rot(1,5);...
    r_1007/0.05  Node1007_NMBC_rot(1,5)];

Result_phi_12_NMBC_com = [r_938/0.05 Node938_NMBC_com(1,6);r_939/0.05  Node939_NMBC_com(1,6);...
    r_941/0.05  Node941_NMBC_com(1,6);r_943/0.05  Node943_NMBC_com(1,6);r_945/0.05  Node945_NMBC_com(1,6);...
    r_947/0.05  Node947_NMBC_com(1,6);r_949/0.05  Node949_NMBC_com(1,6);r_951/0.05  Node951_NMBC_com(1,6);...
    r_953/0.05  Node953_NMBC_com(1,6);r_955/0.05  Node955_NMBC_com(1,6);r_957/0.05  Node957_NMBC_com(1,6);...
    r_959/0.05  Node959_NMBC_com(1,6);r_961/0.05  Node961_NMBC_com(1,6);r_963/0.05  Node963_NMBC_com(1,6);...
    r_965/0.05  Node965_NMBC_com(1,6);r_967/0.05  Node967_NMBC_com(1,6);r_969/0.05  Node969_NMBC_com(1,6);...
    r_971/0.05  Node971_NMBC_com(1,6);r_973/0.05  Node973_NMBC_com(1,6);r_975/0.05  Node975_NMBC_com(1,6);...
    r_977/0.05  Node977_NMBC_com(1,6);r_979/0.05  Node979_NMBC_com(1,6);r_981/0.05  Node981_NMBC_com(1,6);...
    r_983/0.05  Node983_NMBC_com(1,6);r_985/0.05  Node985_NMBC_com(1,6);r_987/0.05  Node987_NMBC_com(1,6);...
    r_989/0.05  Node989_NMBC_com(1,6);r_991/0.05  Node991_NMBC_com(1,6);r_993/0.05  Node993_NMBC_com(1,6);...
    r_995/0.05  Node995_NMBC_com(1,6);r_997/0.05  Node997_NMBC_com(1,6);r_999/0.05  Node999_NMBC_com(1,6);...
    r_1001/0.05  Node1001_NMBC_com(1,6);r_1003/0.05  Node1003_NMBC_com(1,6);r_1005/0.05  Node1005_NMBC_com(1,6);...
    r_1007/0.05  Node1007_NMBC_com(1,6)];

Result_phi_12_NMBC_str = [r_938/0.05 Node938_NMBC_str(1,6);r_939/0.05  Node939_NMBC_str(1,6);...
    r_941/0.05  Node941_NMBC_str(1,6);r_943/0.05  Node943_NMBC_str(1,6);r_945/0.05  Node945_NMBC_str(1,6);...
    r_947/0.05  Node947_NMBC_str(1,6);r_949/0.05  Node949_NMBC_str(1,6);r_951/0.05  Node951_NMBC_str(1,6);...
    r_953/0.05  Node953_NMBC_str(1,6);r_955/0.05  Node955_NMBC_str(1,6);r_957/0.05  Node957_NMBC_str(1,6);...
    r_959/0.05  Node959_NMBC_str(1,6);r_961/0.05  Node961_NMBC_str(1,6);r_963/0.05  Node963_NMBC_str(1,6);...
    r_965/0.05  Node965_NMBC_str(1,6);r_967/0.05  Node967_NMBC_str(1,6);r_969/0.05  Node969_NMBC_str(1,6);...
    r_971/0.05  Node971_NMBC_str(1,6);r_973/0.05  Node973_NMBC_str(1,6);r_975/0.05  Node975_NMBC_str(1,6);...
    r_977/0.05  Node977_NMBC_str(1,6);r_979/0.05  Node979_NMBC_str(1,6);r_981/0.05  Node981_NMBC_str(1,6);...
    r_983/0.05  Node983_NMBC_str(1,6);r_985/0.05  Node985_NMBC_str(1,6);r_987/0.05  Node987_NMBC_str(1,6);...
    r_989/0.05  Node989_NMBC_str(1,6);r_991/0.05  Node991_NMBC_str(1,6);r_993/0.05  Node993_NMBC_str(1,6);...
    r_995/0.05  Node995_NMBC_str(1,6);r_997/0.05  Node997_NMBC_str(1,6);r_999/0.05  Node999_NMBC_str(1,6);...
    r_1001/0.05  Node1001_NMBC_str(1,6);r_1003/0.05  Node1003_NMBC_str(1,6);r_1005/0.05  Node1005_NMBC_str(1,6);...
    r_1007/0.05  Node1007_NMBC_str(1,6)];

Result_phi_12_NMBC_rot = [r_938/0.05 Node938_NMBC_rot(1,6);r_939/0.05  Node939_NMBC_rot(1,6);...
    r_941/0.05  Node941_NMBC_rot(1,6);r_943/0.05  Node943_NMBC_rot(1,6);r_945/0.05  Node945_NMBC_rot(1,6);...
    r_947/0.05  Node947_NMBC_rot(1,6);r_949/0.05  Node949_NMBC_rot(1,6);r_951/0.05  Node951_NMBC_rot(1,6);...
    r_953/0.05  Node953_NMBC_rot(1,6);r_955/0.05  Node955_NMBC_rot(1,6);r_957/0.05  Node957_NMBC_rot(1,6);...
    r_959/0.05  Node959_NMBC_rot(1,6);r_961/0.05  Node961_NMBC_rot(1,6);r_963/0.05  Node963_NMBC_rot(1,6);...
    r_965/0.05  Node965_NMBC_rot(1,6);r_967/0.05  Node967_NMBC_rot(1,6);r_969/0.05  Node969_NMBC_rot(1,6);...
    r_971/0.05  Node971_NMBC_rot(1,6);r_973/0.05  Node973_NMBC_rot(1,6);r_975/0.05  Node975_NMBC_rot(1,6);...
    r_977/0.05  Node977_NMBC_rot(1,6);r_979/0.05  Node979_NMBC_rot(1,6);r_981/0.05  Node981_NMBC_rot(1,6);...
    r_983/0.05  Node983_NMBC_rot(1,6);r_985/0.05  Node985_NMBC_rot(1,6);r_987/0.05  Node987_NMBC_rot(1,6);...
    r_989/0.05  Node989_NMBC_rot(1,6);r_991/0.05  Node991_NMBC_rot(1,6);r_993/0.05  Node993_NMBC_rot(1,6);...
    r_995/0.05  Node995_NMBC_rot(1,6);r_997/0.05  Node997_NMBC_rot(1,6);r_999/0.05  Node999_NMBC_rot(1,6);...
    r_1001/0.05  Node1001_NMBC_rot(1,6);r_1003/0.05  Node1003_NMBC_rot(1,6);r_1005/0.05  Node1005_NMBC_rot(1,6);...
    r_1007/0.05  Node1007_NMBC_rot(1,6)];

Result_phi_12_MBC_com = [r_938/0.05 Node938_MBC_com(1,6);r_939/0.05  Node939_MBC_com(1,6);...
    r_941/0.05  Node941_MBC_com(1,6);r_943/0.05  Node943_MBC_com(1,6);r_945/0.05  Node945_MBC_com(1,6);...
    r_947/0.05  Node947_MBC_com(1,6);r_949/0.05  Node949_MBC_com(1,6);r_951/0.05  Node951_MBC_com(1,6);...
    r_953/0.05  Node953_MBC_com(1,6);r_955/0.05  Node955_MBC_com(1,6);r_957/0.05  Node957_MBC_com(1,6);...
    r_959/0.05  Node959_MBC_com(1,6);r_961/0.05  Node961_MBC_com(1,6);r_963/0.05  Node963_MBC_com(1,6);...
    r_965/0.05  Node965_MBC_com(1,6);r_967/0.05  Node967_MBC_com(1,6);r_969/0.05  Node969_MBC_com(1,6);...
    r_971/0.05  Node971_MBC_com(1,6);r_973/0.05  Node973_MBC_com(1,6);r_975/0.05  Node975_MBC_com(1,6);...
    r_977/0.05  Node977_MBC_com(1,6);r_979/0.05  Node979_MBC_com(1,6);r_981/0.05  Node981_MBC_com(1,6);...
    r_983/0.05  Node983_MBC_com(1,6);r_985/0.05  Node985_MBC_com(1,6);r_987/0.05  Node987_MBC_com(1,6);...
    r_989/0.05  Node989_MBC_com(1,6);r_991/0.05  Node991_MBC_com(1,6);r_993/0.05  Node993_MBC_com(1,6);...
    r_995/0.05  Node995_MBC_com(1,6);r_997/0.05  Node997_MBC_com(1,6);r_999/0.05  Node999_MBC_com(1,6);...
    r_1001/0.05  Node1001_MBC_com(1,6);r_1003/0.05  Node1003_MBC_com(1,6);r_1005/0.05  Node1005_MBC_com(1,6);...
    r_1007/0.05  Node1007_MBC_com(1,6)];

Result_phi_12_MBC_str = [r_938/0.05 Node938_MBC_str(1,6);r_939/0.05  Node939_MBC_str(1,6);...
    r_941/0.05  Node941_MBC_str(1,6);r_943/0.05  Node943_MBC_str(1,6);r_945/0.05  Node945_MBC_str(1,6);...
    r_947/0.05  Node947_MBC_str(1,6);r_949/0.05  Node949_MBC_str(1,6);r_951/0.05  Node951_MBC_str(1,6);...
    r_953/0.05  Node953_MBC_str(1,6);r_955/0.05  Node955_MBC_str(1,6);r_957/0.05  Node957_MBC_str(1,6);...
    r_959/0.05  Node959_MBC_str(1,6);r_961/0.05  Node961_MBC_str(1,6);r_963/0.05  Node963_MBC_str(1,6);...
    r_965/0.05  Node965_MBC_str(1,6);r_967/0.05  Node967_MBC_str(1,6);r_969/0.05  Node969_MBC_str(1,6);...
    r_971/0.05  Node971_MBC_str(1,6);r_973/0.05  Node973_MBC_str(1,6);r_975/0.05  Node975_MBC_str(1,6);...
    r_977/0.05  Node977_MBC_str(1,6);r_979/0.05  Node979_MBC_str(1,6);r_981/0.05  Node981_MBC_str(1,6);...
    r_983/0.05  Node983_MBC_str(1,6);r_985/0.05  Node985_MBC_str(1,6);r_987/0.05  Node987_MBC_str(1,6);...
    r_989/0.05  Node989_MBC_str(1,6);r_991/0.05  Node991_MBC_str(1,6);r_993/0.05  Node993_MBC_str(1,6);...
    r_995/0.05  Node995_MBC_str(1,6);r_997/0.05  Node997_MBC_str(1,6);r_999/0.05  Node999_MBC_str(1,6);...
    r_1001/0.05  Node1001_MBC_str(1,6);r_1003/0.05  Node1003_MBC_str(1,6);r_1005/0.05  Node1005_MBC_str(1,6);...
    r_1007/0.05  Node1007_MBC_str(1,6)];

Result_phi_12_MBC_rot = [r_938/0.05 Node938_MBC_rot(1,6);r_939/0.05  Node939_MBC_rot(1,6);...
    r_941/0.05  Node941_MBC_rot(1,6);r_943/0.05  Node943_MBC_rot(1,6);r_945/0.05  Node945_MBC_rot(1,6);...
    r_947/0.05  Node947_MBC_rot(1,6);r_949/0.05  Node949_MBC_rot(1,6);r_951/0.05  Node951_MBC_rot(1,6);...
    r_953/0.05  Node953_MBC_rot(1,6);r_955/0.05  Node955_MBC_rot(1,6);r_957/0.05  Node957_MBC_rot(1,6);...
    r_959/0.05  Node959_MBC_rot(1,6);r_961/0.05  Node961_MBC_rot(1,6);r_963/0.05  Node963_MBC_rot(1,6);...
    r_965/0.05  Node965_MBC_rot(1,6);r_967/0.05  Node967_MBC_rot(1,6);r_969/0.05  Node969_MBC_rot(1,6);...
    r_971/0.05  Node971_MBC_rot(1,6);r_973/0.05  Node973_MBC_rot(1,6);r_975/0.05  Node975_MBC_rot(1,6);...
    r_977/0.05  Node977_MBC_rot(1,6);r_979/0.05  Node979_MBC_rot(1,6);r_981/0.05  Node981_MBC_rot(1,6);...
    r_983/0.05  Node983_MBC_rot(1,6);r_985/0.05  Node985_MBC_rot(1,6);r_987/0.05  Node987_MBC_rot(1,6);...
    r_989/0.05  Node989_MBC_rot(1,6);r_991/0.05  Node991_MBC_rot(1,6);r_993/0.05  Node993_MBC_rot(1,6);...
    r_995/0.05  Node995_MBC_rot(1,6);r_997/0.05  Node997_MBC_rot(1,6);r_999/0.05  Node999_MBC_rot(1,6);...
    r_1001/0.05  Node1001_MBC_rot(1,6);r_1003/0.05  Node1003_MBC_rot(1,6);r_1005/0.05  Node1005_MBC_rot(1,6);...
    r_1007/0.05  Node1007_MBC_rot(1,6)];


Micropolar=load('Plate_with_hole.txt');

Tteta_teta = load('Tteta_teta.txt');
Trteta_r = load('Trteta_r.txt');
Ttetar_r = load('Ttetar_r.txt');
Trtheta_r_classic = load('Trtheta_r_classic.txt');







% box on
% Figure1 = figure(1);
% hold on
% plot(Result_S11_intensity_NMBC_str(:,1)/0.05,Result_S11_intensity_NMBC_str(:,2),'--ob',Result_S11_intensity_NMBC_com(:,1)/0.05,Result_S11_intensity_NMBC_com(:,2),'--og',...
%     Result_S11_intensity_NMBC_rot(:,1)/0.05,Result_S11_intensity_NMBC_rot(:,2),'--or',Result_S11_intensity_cls(:,1)/0.05,Result_S11_intensity_cls(:,2),'--om',Micropolar(:,1),Micropolar(:,2),'--ok','LineWidth',2 )
% set(Figure1,'defaulttextinterpreter','latex');
% ylabel('$t_{\theta\theta}/p$','fontsize',16)
% xlabel('$r/a$','fontsize',16)
% lgnd=legend({'$\Phi_{11}$, $\Phi_{22}$ Micromorphic','$\Phi_{11}$, $\Phi_{22}$, $\Phi_{12}$, $\Phi_{21}$ Micromorphic'...
%     ,'$\Phi_{21}$, $\Phi_{12}$ Micromorphic','Classical Continuum','Micropolar'},'Interpreter','Latex');
% set(gca,'FontName','mwa_cmr10','FontSize',14)
% set(lgnd, 'XColor', 'w', 'YColor', 'w', 'Color', 'none');
% 
% 
% box on
% Figure2 = figure(2);
% hold on
% plot(Result_S11_intensity_MBC_str(:,1)/0.05,Result_S11_intensity_MBC_str(:,2),'--ob',Result_S11_intensity_MBC_com(:,1)/0.05,Result_S11_intensity_MBC_com(:,2),'--og',...
%     Result_S11_intensity_MBC_rot(:,1)/0.05,Result_S11_intensity_MBC_rot(:,2),'--or',Result_S11_intensity_NMBC_str(:,1)/0.05,Result_S11_intensity_NMBC_str(:,2),'--ok',Result_S11_intensity_NMBC_com(:,1)/0.05,Result_S11_intensity_NMBC_com(:,2),'--om',...
%     Result_S11_intensity_NMBC_rot(:,1)/0.05,Result_S11_intensity_NMBC_rot(:,2),'--oc','LineWidth',2 )
% set(Figure2,'defaulttextinterpreter','latex');
% ylabel('$t_{\theta\theta}/p$','fontsize',16)
% xlabel('$r/a$','fontsize',16)
% lgnd=legend({'$\Phi_{11}$, $\Phi_{22}$ Micromorphic','$\Phi_{11}$, $\Phi_{22}$, $\Phi_{12}$, $\Phi_{21}$ Micromorphic'...
%     ,'$\Phi_{21}$, $\Phi_{12}$ Micromorphic','$\Phi_{11}$, $\Phi_{22}$ No Micromorphic BC','$\Phi_{11}$, $\Phi_{22}$, $\Phi_{12}$, $\Phi_{21}$ No Micromorphic BC'...
%     ,'$\Phi_{21}$, $\Phi_{12}$ No Micromorphic BC'},'Interpreter','Latex');
% set(gca,'FontName','mwa_cmr10','FontSize',14)
% set(lgnd, 'XColor', 'w', 'YColor', 'w', 'Color', 'none');
% 
% 
% 
% box on
% Figure3 = figure(3);
% hold on
% plot(Result_S_teta_teta_NMBC_str(:,1),Result_S_teta_teta_NMBC_str(:,2),'--ob',Result_S_teta_teta_NMBC_com(:,1),Result_S_teta_teta_NMBC_com(:,2),'--og',...
%     Result_S_teta_teta_NMBC_rot(:,1),Result_S_teta_teta_NMBC_rot(:,2),'--or',Result_S_teta_teta_cls(:,1),Result_S_teta_teta_cls(:,2),'--om',Tteta_teta(:,1),Tteta_teta(:,2),'--ok','LineWidth',2 )
% set(Figure3,'defaulttextinterpreter','latex');
% ylabel('$t_{\theta\theta}/p$','fontsize',16)
% xlabel('$\theta$','fontsize',16)
% lgnd=legend({'$\Phi_{11}$, $\Phi_{22}$ Micromorphic','$\Phi_{11}$, $\Phi_{22}$, $\Phi_{12}$, $\Phi_{21}$ Micromorphic'...
%     ,'$\Phi_{21}$, $\Phi_{12}$ Micromorphic','Classical Continuum','Micropolar'},'Interpreter','Latex');
% set(gca,'FontName','mwa_cmr10','FontSize',14)
% set(lgnd, 'XColor', 'w', 'YColor', 'w', 'Color', 'none');
% 
% 
% box on
% Figure4 = figure(4);
% hold on
% plot(Result_S_teta_teta_NMBC_str(:,1),Result_S_teta_teta_NMBC_str(:,2),'--ob',Result_S_teta_teta_NMBC_com(:,1),Result_S_teta_teta_NMBC_com(:,2),'--og',...
%     Result_S_teta_teta_NMBC_rot(:,1),Result_S_teta_teta_NMBC_rot(:,2),'--or',Result_S_teta_teta_MBC_str(:,1),Result_S_teta_teta_MBC_str(:,2),'--ok',Result_S_teta_teta_MBC_com(:,1),Result_S_teta_teta_MBC_com(:,2),'--om',...
%     Result_S_teta_teta_MBC_rot(:,1),Result_S_teta_teta_MBC_rot(:,2),'--oc','LineWidth',2 )
% set(Figure4,'defaulttextinterpreter','latex');
% ylabel('$t_{\theta\theta}/p$','fontsize',16)
% xlabel('$\theta$','fontsize',16)
% lgnd=legend({'$\Phi_{11}$, $\Phi_{22}$ Micromorphic','$\Phi_{11}$, $\Phi_{22}$, $\Phi_{12}$, $\Phi_{21}$ Micromorphic'...
%     ,'$\Phi_{21}$, $\Phi_{12}$ Micromorphic','$\Phi_{11}$, $\Phi_{22}$ No Micromorphic BC','$\Phi_{11}$, $\Phi_{22}$, $\Phi_{12}$, $\Phi_{21}$ No Micromorphic BC'...
%     ,'$\Phi_{21}$, $\Phi_{12}$ No Micromorphic BC'},'Interpreter','Latex');
% set(gca,'FontName','mwa_cmr10','FontSize',14)
% set(lgnd, 'XColor', 'w', 'YColor', 'w', 'Color', 'none');
% 
% 
% box on
% Figure5 = figure(5);
% hold on
% plot(Result_S_r_teta_cls(1:23,1),Result_S_r_teta_cls(1:23,2),'--ob',Result_S_r_teta_NMBC_com(1:23,1),Result_S_r_teta_NMBC_com(1:23,2),'--og',...
%     Result_S_r_teta_NMBC_rot(1:23,1),Result_S_r_teta_NMBC_rot(1:23,2),'--or',Result_S_r_teta_NMBC_str(1:23,1),Result_S_r_teta_NMBC_str(1:23,2),'--om',Trtheta_r_classic(:,1),Trtheta_r_classic(:,2),'--ok',...
%     Ttetar_r(:,1),Ttetar_r(:,2),'--oc','LineWidth',2 )
% set(Figure5,'defaulttextinterpreter','latex');
% ylabel('$t_{r\theta}/p$','fontsize',16)
% xlabel('$r/a$','fontsize',16)
% lgnd=legend({'Classical Continuum','$\Phi_{11}$, $\Phi_{22}$, $\Phi_{12}$, $\Phi_{21}$ Micromorphic'...
%     ,'$\Phi_{21}$, $\Phi_{12}$ Micromorphic','$\Phi_{11}$, $\Phi_{22}$ Micromorphic','Analytical Solution','Micropolar'},'Interpreter','Latex');
% set(gca,'FontName','mwa_cmr10','FontSize',14)
% set(lgnd, 'XColor', 'w', 'YColor', 'w', 'Color', 'none');
% 
% 
% box on
% Figure6 = figure(6);
% hold on
% plot(Result_S_r_teta_MBC_str(1:23,1),Result_S_r_teta_MBC_str(1:23,2),'--ob',Result_S_r_teta_MBC_com(1:23,1),Result_S_r_teta_MBC_com(1:23,2),'--og',...
%     Result_S_r_teta_MBC_rot(1:23,1),Result_S_r_teta_MBC_rot(1:23,2),'--or',Result_S_r_teta_NMBC_str(1:23,1),Result_S_r_teta_NMBC_str(1:23,2),'--ok',Result_S_r_teta_NMBC_com(1:23,1),Result_S_r_teta_NMBC_com(1:23,2),'--om',...
%     Result_S_r_teta_NMBC_rot(1:23,1),Result_S_r_teta_NMBC_rot(1:23,2),'--oc','LineWidth',2 )
% set(Figure6,'defaulttextinterpreter','latex');
% ylabel('$t_{r\theta}/p$','fontsize',16)
% xlabel('$r/a$','fontsize',16)
% lgnd=legend({'$\Phi_{11}$, $\Phi_{22}$ Micromorphic','$\Phi_{11}$, $\Phi_{22}$, $\Phi_{12}$, $\Phi_{21}$ Micromorphic'...
%     ,'$\Phi_{21}$, $\Phi_{12}$ Micromorphic','$\Phi_{11}$, $\Phi_{22}$ No Micromorphic BC','$\Phi_{11}$, $\Phi_{22}$, $\Phi_{12}$, $\Phi_{21}$ No Micromorphic BC'...
%     ,'$\Phi_{21}$, $\Phi_{12}$ No Micromorphic BC'},'Interpreter','Latex');
% set(gca,'FontName','mwa_cmr10','FontSize',14)
% set(lgnd, 'XColor', 'w', 'YColor', 'w', 'Color', 'none');
% 
% 
% box on
% Figure7 = figure(7);
% hold on
% plot(Result_S_teta_r_NMBC_com(1:23,1),Result_S_teta_r_NMBC_com(1:23,2),'--og',...
%     Result_S_teta_r_NMBC_rot(1:23,1),Result_S_teta_r_NMBC_rot(1:23,2),'--or',Result_S_teta_r_NMBC_str(1:23,1),Result_S_teta_r_NMBC_str(1:23,2),'--om',...
%     Trteta_r(:,1),Trteta_r(:,2),'--oc','LineWidth',2 )
% set(Figure7,'defaulttextinterpreter','latex');
% ylabel('$t_{\theta r}/p$','fontsize',16)
% xlabel('$r/a$','fontsize',16)
% lgnd=legend({'$\Phi_{11}$, $\Phi_{22}$, $\Phi_{12}$, $\Phi_{21}$ Micromorphic'...
%     ,'$\Phi_{21}$, $\Phi_{12}$ Micromorphic','$\Phi_{11}$, $\Phi_{22}$ Micromorphic','Micropolar'},'Interpreter','Latex');
% set(gca,'FontName','mwa_cmr10','FontSize',14)
% set(lgnd, 'XColor', 'w', 'YColor', 'w', 'Color', 'none');
% 
% 
% box on
% Figure8 = figure(8);
% hold on
% plot(Result_S_teta_r_MBC_str(1:23,1),Result_S_teta_r_MBC_str(1:23,2),'--ob',Result_S_teta_r_MBC_com(1:23,1),Result_S_teta_r_MBC_com(1:23,2),'--og',...
%     Result_S_teta_r_MBC_rot(1:23,1),Result_S_teta_r_MBC_rot(1:23,2),'--or',Result_S_teta_r_NMBC_str(1:23,1),Result_S_teta_r_NMBC_str(1:23,2),'--ok',Result_S_teta_r_NMBC_com(1:23,1),Result_S_teta_r_NMBC_com(1:23,2),'--om',...
%     Result_S_teta_r_NMBC_rot(1:23,1),Result_S_teta_r_NMBC_rot(1:23,2),'--oc','LineWidth',2 )
% set(Figure8,'defaulttextinterpreter','latex');
% ylabel('$t_{\theta r}/p$','fontsize',16)
% xlabel('$r/a$','fontsize',16)
% lgnd=legend({'$\Phi_{11}$, $\Phi_{22}$ Micromorphic','$\Phi_{11}$, $\Phi_{22}$, $\Phi_{12}$, $\Phi_{21}$ Micromorphic'...
%     ,'$\Phi_{21}$, $\Phi_{12}$ Micromorphic','$\Phi_{11}$, $\Phi_{22}$ No Micromorphic BC','$\Phi_{11}$, $\Phi_{22}$, $\Phi_{12}$, $\Phi_{21}$ No Micromorphic BC'...
%     ,'$\Phi_{21}$, $\Phi_{12}$ No Micromorphic BC'},'Interpreter','Latex');
% set(gca,'FontName','mwa_cmr10','FontSize',14)
% set(lgnd, 'XColor', 'w', 'YColor', 'w', 'Color', 'none');
% 
% 
% 
% box on
% Figure9 = figure(9);
% hold on
% plot(Result_phi11_intensity_MBC_str(:,1)/0.05,Result_phi11_intensity_MBC_str(:,2),'--ob',Result_phi11_intensity_MBC_com(:,1)/0.05,Result_phi11_intensity_MBC_com(:,2),'--og',...
%     Result_phi11_intensity_NMBC_str(:,1)/0.05,Result_phi11_intensity_NMBC_str(:,2),'--ok',...
%     Result_phi11_intensity_NMBC_com(:,1)/0.05,Result_phi11_intensity_NMBC_com(:,2),'--om','LineWidth',2 )
% set(Figure9,'defaulttextinterpreter','latex');
% ylabel('$\phi_{11}$','fontsize',16)
% xlabel('$r/a$','fontsize',16)
% lgnd=legend({'$\Phi_{11}$, $\Phi_{22}$ Micromorphic','$\Phi_{11}$, $\Phi_{22}$, $\Phi_{12}$, $\Phi_{21}$ Micromorphic'...
%     ,'$\Phi_{11}$, $\Phi_{22}$ No Micromorphic BC','$\Phi_{11}$, $\Phi_{22}$, $\Phi_{12}$, $\Phi_{21}$ No Micromorphic BC'},'Interpreter','Latex');
% set(gca,'FontName','mwa_cmr10','FontSize',14)
% set(lgnd, 'XColor', 'w', 'YColor', 'w', 'Color', 'none');
% 
% 
% box on
% Figure10 = figure(10);
% hold on
% plot(Result_phi22_intensity_MBC_str(:,1)/0.05,Result_phi22_intensity_MBC_str(:,2),'--ob',Result_phi22_intensity_MBC_com(:,1)/0.05,Result_phi22_intensity_MBC_com(:,2),'--og',...
%     Result_phi22_intensity_NMBC_str(:,1)/0.05,Result_phi22_intensity_NMBC_str(:,2),'--ok',...
%     Result_phi22_intensity_NMBC_com(:,1)/0.05,Result_phi22_intensity_NMBC_com(:,2),'--om','LineWidth',2 )
% set(Figure10,'defaulttextinterpreter','latex');
% ylabel('$\phi_{22}$','fontsize',16)
% xlabel('$r/a$','fontsize',16)
% lgnd=legend({'$\Phi_{11}$, $\Phi_{22}$ Micromorphic','$\Phi_{11}$, $\Phi_{22}$, $\Phi_{12}$, $\Phi_{21}$ Micromorphic'...
%     ,'$\Phi_{11}$, $\Phi_{22}$ No Micromorphic BC','$\Phi_{11}$, $\Phi_{22}$, $\Phi_{12}$, $\Phi_{21}$ No Micromorphic BC'},'Interpreter','Latex');
% set(gca,'FontName','mwa_cmr10','FontSize',14)
% set(lgnd, 'XColor', 'w', 'YColor', 'w', 'Color', 'none');
% 
% 
% 
% box on
% Figure11 = figure(11);
% hold on
% plot(Result_phi_12_MBC_rot(:,1),Result_phi_12_MBC_rot(:,2),'--ob',Result_phi_12_MBC_com(:,1),Result_phi_12_MBC_com(:,2),'--og',...
%     Result_phi_12_NMBC_rot(:,1),Result_phi_12_NMBC_rot(:,2),'--ok',...
%     Result_phi_12_NMBC_com(:,1),Result_phi_12_NMBC_com(:,2),'--om','LineWidth',2 )
% set(Figure11,'defaulttextinterpreter','latex');
% ylabel('$\phi_{12}$','fontsize',16)
% xlabel('$r/a$','fontsize',16)
% lgnd=legend({'$\Phi_{21}$, $\Phi_{12}$ Micromorphic','$\Phi_{11}$, $\Phi_{22}$, $\Phi_{12}$, $\Phi_{21}$ Micromorphic'...
%     ,'$\Phi_{21}$, $\Phi_{12}$ No Micromorphic BC','$\Phi_{11}$, $\Phi_{22}$, $\Phi_{12}$, $\Phi_{21}$ No Micromorphic BC'},'Interpreter','Latex');
% set(gca,'FontName','mwa_cmr10','FontSize',14)
% set(lgnd, 'XColor', 'w', 'YColor', 'w', 'Color', 'none');
% 
% 
% box on
% Figure12 = figure(12);
% hold on
% plot(Result_phi_21_MBC_rot(:,1),Result_phi_21_MBC_rot(:,2),'--ob',Result_phi_21_MBC_com(:,1),Result_phi_21_MBC_com(:,2),'--og',...
%     Result_phi_21_NMBC_rot(:,1),Result_phi_21_NMBC_rot(:,2),'--ok',...
%     Result_phi_21_NMBC_com(:,1),Result_phi_21_NMBC_com(:,2),'--om','LineWidth',2 )
% set(Figure12,'defaulttextinterpreter','latex');
% ylabel('$\phi_{21}$','fontsize',16)
% xlabel('$r/a$','fontsize',16)
% lgnd=legend({'$\Phi_{21}$, $\Phi_{12}$ Micromorphic','$\Phi_{11}$, $\Phi_{22}$, $\Phi_{12}$, $\Phi_{21}$ Micromorphic'...
%     ,'$\Phi_{21}$, $\Phi_{12}$ No Micromorphic BC','$\Phi_{11}$, $\Phi_{22}$, $\Phi_{12}$, $\Phi_{21}$ No Micromorphic BC'},'Interpreter','Latex');
% set(gca,'FontName','mwa_cmr10','FontSize',14)
% set(lgnd, 'XColor', 'w', 'YColor', 'w', 'Color', 'none');
% 




box on
Figure1 = figure(1);
hold on
plot(Result_S11_intensity_cls(:,1)/0.05,Result_S11_intensity_cls(:,2),'--oc',Micropolar(:,1),Micropolar(:,2),'--ok',Result_S11_intensity_NMBC_rot(:,1)/0.05,Result_S11_intensity_NMBC_rot(:,2),'--or',...
    Result_S11_intensity_NMBC_com(:,1)/0.05,Result_S11_intensity_NMBC_com(:,2),'--og',Result_S11_intensity_NMBC_str(:,1)/0.05,Result_S11_intensity_NMBC_str(:,2),'--ob','LineWidth',2 )
set(Figure1,'defaulttextinterpreter','latex');
ylabel('$\sigma_{\theta\theta}/p$','fontsize',16)
xlabel('$r/a$','fontsize',16)
lgnd=legend({'Classical Continuum','Micropolar Bauer et al. CMAME-2010.','Case B-III Micromorphic','Case B-I Micromorphic'...
    ,'Case B-II Micromorphic'},'Interpreter','Latex');
set(gca,'FontName','Times New Roman','FontSize',14)
set(lgnd, 'XColor', 'w', 'YColor', 'w', 'Color', 'none');


box on
Figure20 = figure(20);
hold on
plot(Result_S11_intensity_cls(1:18,1)/0.05,Result_S11_intensity_cls(1:18,2),'--oc',Micropolar(1:9,1),Micropolar(1:9,2),'--ok',Result_S11_intensity_NMBC_rot(1:18,1)/0.05,Result_S11_intensity_NMBC_rot(1:18,2),'--or',...
    Result_S11_intensity_NMBC_com(1:18,1)/0.05,Result_S11_intensity_NMBC_com(1:18,2),'--og',Result_S11_intensity_NMBC_str(1:18,1)/0.05,Result_S11_intensity_NMBC_str(1:18,2),'--ob','LineWidth',2 )
set(Figure20,'defaulttextinterpreter','latex');
ylabel('$\sigma_{\theta\theta}/p$','fontsize',16)
xlabel('$r/a$','fontsize',16)
lgnd=legend({'Classical Continuum','Micropolar Bauer et al. CMAME-2010.','Case B-III Micromorphic','Case B-I Micromorphic'...
    ,'Case B-II Micromorphic'},'Interpreter','Latex');
set(gca,'FontName','Times New Roman','FontSize',14)
set(lgnd, 'XColor', 'w', 'YColor', 'w', 'Color', 'none');


box on
Figure2 = figure(2);
hold on
CC = plot(Result_S11_intensity_MBC_str(:,1)/0.05,Result_S11_intensity_MBC_str(:,2),'--ob',Result_S11_intensity_MBC_com(:,1)/0.05,Result_S11_intensity_MBC_com(:,2),'--og',...
    Result_S11_intensity_MBC_rot(:,1)/0.05,Result_S11_intensity_MBC_rot(:,2),'--or',Result_S11_intensity_NMBC_str(:,1)/0.05,Result_S11_intensity_NMBC_str(:,2),'--ok',Result_S11_intensity_NMBC_com(:,1)/0.05,Result_S11_intensity_NMBC_com(:,2),...
    Result_S11_intensity_NMBC_rot(:,1)/0.05,Result_S11_intensity_NMBC_rot(:,2),'--oc','LineWidth',2 );
set(Figure2,'defaulttextinterpreter','latex');
ylabel('$\sigma_{\theta\theta}/p$','fontsize',16)
xlabel('$r/a$','fontsize',16)
lgnd=legend({'Case A-II Micromorphic','Case A-I Micromorphic'...
    ,'Case A-III Micromorphic','Case B-II Micromorphic','Case B-I Micromorphic'...
    ,'Case B-III Micromorphic'},'Interpreter','Latex');
set(gca,'FontName','Times New Roman','FontSize',14)
set(lgnd, 'XColor', 'w', 'YColor', 'w', 'Color', 'none');
set(CC(5),'Color', [1, 0.5, 0], 'LineStyle', '--','Marker','o')


box on
Figure21 = figure(21);
hold on
CC = plot(Result_S11_intensity_MBC_str(1:18,1)/0.05,Result_S11_intensity_MBC_str(1:18,2),'--ob',Result_S11_intensity_MBC_com(1:18,1)/0.05,Result_S11_intensity_MBC_com(1:18,2),'--og',...
    Result_S11_intensity_MBC_rot(1:18,1)/0.05,Result_S11_intensity_MBC_rot(1:18,2),'--or',Result_S11_intensity_NMBC_str(1:18,1)/0.05,Result_S11_intensity_NMBC_str(1:18,2),'--ok',Result_S11_intensity_NMBC_com(1:18,1)/0.05,Result_S11_intensity_NMBC_com(1:18,2),...
    Result_S11_intensity_NMBC_rot(1:18,1)/0.05,Result_S11_intensity_NMBC_rot(1:18,2),'--oc','LineWidth',2 );
set(Figure21,'defaulttextinterpreter','latex');
ylabel('$\sigma_{\theta\theta}/p$','fontsize',16)
xlabel('$r/a$','fontsize',16)
lgnd=legend({'Case A-II Micromorphic','Case A-I Micromorphic'...
    ,'Case A-III Micromorphic','Case B-II Micromorphic','Case B-I Micromorphic'...
    ,'Case B-III Micromorphic'},'Interpreter','Latex');
set(gca,'FontName','Times New Roman','FontSize',14)
set(lgnd, 'XColor', 'w', 'YColor', 'w', 'Color', 'none');
set(CC(5),'Color', [1, 0.5, 0], 'LineStyle', '--','Marker','o')


box on
Figure3 = figure(3);
hold on
CC = plot(Result_S_teta_teta_cls(:,1),Result_S_teta_teta_cls(:,2),Tteta_teta(:,1),Tteta_teta(:,2),'--ok'...
    ,Result_S_teta_teta_NMBC_rot(:,1),Result_S_teta_teta_NMBC_rot(:,2),'--or',Result_S_teta_teta_NMBC_com(:,1),Result_S_teta_teta_NMBC_com(:,2),'--og'...
    ,Result_S_teta_teta_NMBC_str(:,1),Result_S_teta_teta_NMBC_str(:,2),'--ob','LineWidth',2 );
set(Figure3,'defaulttextinterpreter','latex');
ylabel('$\sigma_{\theta\theta}/p$','fontsize',16)
xlabel('$\theta$','fontsize',16)
lgnd=legend({'Classical Continuum','Micropolar Bauer et al. CMAME-2010','Case B-III Micromorphic','Case B-I Micromorphic'...
    ,'Case B-II Micromorphic'},'Interpreter','Latex');
set(gca,'FontName','Times New Roman','FontSize',14)
set(lgnd, 'XColor', 'w', 'YColor', 'w', 'Color', 'none');
set(CC(1),'Color', [1, 0.5, 0], 'LineStyle', '--','Marker','o')


box on
Figure4 = figure(4);
hold on
CC = plot(Result_S_teta_teta_MBC_str(:,1),Result_S_teta_teta_MBC_str(:,2),'--ok',Result_S_teta_teta_MBC_com(:,1),Result_S_teta_teta_MBC_com(:,2),Result_S_teta_teta_NMBC_str(:,1),Result_S_teta_teta_NMBC_str(:,2),'--ob',Result_S_teta_teta_NMBC_com(:,1),Result_S_teta_teta_NMBC_com(:,2),'--og',...
    Result_S_teta_teta_NMBC_rot(:,1),Result_S_teta_teta_NMBC_rot(:,2),'--or',...
    Result_S_teta_teta_MBC_rot(:,1),Result_S_teta_teta_MBC_rot(:,2),'--oc','LineWidth',2 );
set(Figure4,'defaulttextinterpreter','latex');
ylabel('$\sigma_{\theta\theta}/p$','fontsize',16)
xlabel('$\theta$','fontsize',16)
lgnd=legend({'Case A-II Micromorphic','Case A-I Micromorphic','Case B-II Micromorphic','Case B-I Micromorphic'...
    ,'Case B-III Micromorphic'...
    ,'Case A-III Micromorphic'},'Interpreter','Latex');
set(gca,'FontName','Times New Roman','FontSize',14)
set(lgnd, 'XColor', 'w', 'YColor', 'w', 'Color', 'none');
set(CC(2),'Color', [1, 0.5, 0], 'LineStyle', '--','Marker','o')

box on
Figure5 = figure(5);
hold on
CC = plot(Result_S_r_teta_cls(1:23,1),Result_S_r_teta_cls(1:23,2),'--ob',Result_S_r_teta_NMBC_com(1:23,1),Result_S_r_teta_NMBC_com(1:23,2),'--og',...
    Result_S_r_teta_NMBC_rot(1:23,1),Result_S_r_teta_NMBC_rot(1:23,2),'--or',Result_S_r_teta_NMBC_str(1:23,1),Result_S_r_teta_NMBC_str(1:23,2),Trtheta_r_classic(:,1),Trtheta_r_classic(:,2),'k',...
    Ttetar_r(:,1),Ttetar_r(:,2),'--oc','LineWidth',2 );
set(Figure5,'defaulttextinterpreter','latex');
ylabel('$\sigma_{r\theta}/p$','fontsize',16)
xlabel('$r/a$','fontsize',16)
lgnd=legend({'Classical Continuum','Case B-I Micromorphic'...
    ,'Case B-III Micromorphic','Case B-II Micromorphic','Analytical Solution','Micropolar Bauer et al. CMAME-2010'},'Interpreter','Latex');
set(gca,'FontName','Times New Roman','FontSize',14)
set(lgnd, 'XColor', 'w', 'YColor', 'w', 'Color', 'none');
set(CC(4),'Color', [1, 0.5, 0], 'LineStyle', '--','Marker','o')


box on
Figure6 = figure(6);
hold on
CC = plot(Result_S_r_teta_MBC_str(1:23,1),Result_S_r_teta_MBC_str(1:23,2),'--ob',Result_S_r_teta_MBC_com(1:23,1),Result_S_r_teta_MBC_com(1:23,2),'--og',...
    Result_S_r_teta_MBC_rot(1:23,1),Result_S_r_teta_MBC_rot(1:23,2),'--or',Result_S_r_teta_NMBC_str(1:23,1),Result_S_r_teta_NMBC_str(1:23,2),'--ok',Result_S_r_teta_NMBC_com(1:23,1),Result_S_r_teta_NMBC_com(1:23,2),...
    Result_S_r_teta_NMBC_rot(1:23,1),Result_S_r_teta_NMBC_rot(1:23,2),'--oc','LineWidth',2 );
set(Figure6,'defaulttextinterpreter','latex');
ylabel('$\sigma_{r\theta}/p$','fontsize',16)
xlabel('$r/a$','fontsize',16)
lgnd=legend({'Case A-II Micromorphic','Case A-I Micromorphic'...
    ,'Case A-III Micromorphic','Case B-II Micromorphic','Case B-I Micromorphic'...
    ,'Case B-III Micromorphic'},'Interpreter','Latex');
set(gca,'FontName','Times New Roman','FontSize',14)
set(lgnd, 'XColor', 'w', 'YColor', 'w', 'Color', 'none');
set(CC(5),'Color', [1, 0.5, 0], 'LineStyle', '--','Marker','o')

box on
Figure7 = figure(7);
hold on
plot(Result_S_teta_r_NMBC_com(1:23,1),Result_S_teta_r_NMBC_com(1:23,2),'--og',...
    Result_S_teta_r_NMBC_rot(1:23,1),Result_S_teta_r_NMBC_rot(1:23,2),'--or',Result_S_teta_r_NMBC_str(1:23,1),Result_S_teta_r_NMBC_str(1:23,2),'--ok',...
    Trteta_r(:,1),Trteta_r(:,2),'--oc','LineWidth',2 )
set(Figure7,'defaulttextinterpreter','latex');
ylabel('$\sigma_{\theta r}/p$','fontsize',16)
xlabel('$r/a$','fontsize',16)
lgnd=legend({'Case B-I Micromorphic'...
    ,'Case B-III Micromorphic','Case B-II Micromorphic','Micropolar Bauer et al. CMAME-2010'},'Interpreter','Latex');
set(gca,'FontName','Times New Roman','FontSize',14)
set(lgnd, 'XColor', 'w', 'YColor', 'w', 'Color', 'none');


box on
Figure8 = figure(8);
hold on
CC = plot(Result_S_teta_r_MBC_str(1:23,1),Result_S_teta_r_MBC_str(1:23,2),'--ob',Result_S_teta_r_MBC_com(1:23,1),Result_S_teta_r_MBC_com(1:23,2),'--og',...
    Result_S_teta_r_MBC_rot(1:23,1),Result_S_teta_r_MBC_rot(1:23,2),'--or',Result_S_teta_r_NMBC_str(1:23,1),Result_S_teta_r_NMBC_str(1:23,2),'--ok',Result_S_teta_r_NMBC_com(1:23,1),Result_S_teta_r_NMBC_com(1:23,2),...
    Result_S_teta_r_NMBC_rot(1:23,1),Result_S_teta_r_NMBC_rot(1:23,2),'--oc','LineWidth',2 );
set(Figure8,'defaulttextinterpreter','latex');
ylabel('$\sigma_{\theta r}/p$','fontsize',16)
xlabel('$r/a$','fontsize',16)
lgnd=legend({'Case A-II Micromorphic','Case A-I Micromorphic'...
    ,'Case A-III Micromorphic','Case B-II Micromorphic','Case B-I Micromorphic'...
    ,'Case B-III Micromorphic'},'Interpreter','Latex');
set(gca,'FontName','Times New Roman','FontSize',14)
set(lgnd, 'XColor', 'w', 'YColor', 'w', 'Color', 'none');
set(CC(5),'Color', [1, 0.5, 0], 'LineStyle', '--','Marker','o')


box on
Figure9 = figure(9);
hold on
plot(Result_phi11_intensity_MBC_str(:,1)/0.05,Result_phi11_intensity_MBC_str(:,2),'--ob',Result_phi11_intensity_MBC_com(:,1)/0.05,Result_phi11_intensity_MBC_com(:,2),'--og',...
    Result_phi11_intensity_NMBC_str(:,1)/0.05,Result_phi11_intensity_NMBC_str(:,2),'--ok',...
    Result_phi11_intensity_NMBC_com(:,1)/0.05,Result_phi11_intensity_NMBC_com(:,2),'--om','LineWidth',2 )
set(Figure9,'defaulttextinterpreter','latex');
ylabel('$\phi_{11}$','fontsize',16)
xlabel('$r/a$','fontsize',16)
lgnd=legend({'Case A-II Micromorphic','Case A-I Micromorphic'...
    ,'Case B-II Micromorphic','Case B-I Micromorphic'},'Interpreter','Latex');
set(gca,'FontName','Times New Roman','FontSize',14)
set(lgnd, 'XColor', 'w', 'YColor', 'w', 'Color', 'none');


box on
Figure10 = figure(10);
hold on
plot(Result_phi22_intensity_MBC_str(:,1)/0.05,Result_phi22_intensity_MBC_str(:,2),'--ob',Result_phi22_intensity_MBC_com(:,1)/0.05,Result_phi22_intensity_MBC_com(:,2),'--og',...
    Result_phi22_intensity_NMBC_str(:,1)/0.05,Result_phi22_intensity_NMBC_str(:,2),'--ok',...
    Result_phi22_intensity_NMBC_com(:,1)/0.05,Result_phi22_intensity_NMBC_com(:,2),'--om','LineWidth',2 )
set(Figure10,'defaulttextinterpreter','latex');
ylabel('$\phi_{22}$','fontsize',16)
xlabel('$r/a$','fontsize',16)
lgnd=legend({'Case A-II Micromorphic','Case A-I Micromorphic'...
    ,'Case B-II Micromorphic','Case B-I Micromorphic'},'Interpreter','Latex');
set(gca,'FontName','Times New Roman','FontSize',14)
set(lgnd, 'XColor', 'w', 'YColor', 'w', 'Color', 'none');



box on
Figure11 = figure(11);
box on
hold on
plot(Result_phi_12_MBC_rot(:,1),Result_phi_12_MBC_rot(:,2),'--ob',Result_phi_12_MBC_com(:,1),Result_phi_12_MBC_com(:,2),'--og',...
    Result_phi_12_NMBC_rot(:,1),Result_phi_12_NMBC_rot(:,2),'--ok',...
    Result_phi_12_NMBC_com(:,1),Result_phi_12_NMBC_com(:,2),'--om','LineWidth',2 )
set(Figure11,'defaulttextinterpreter','latex');
ylabel('$\phi_{12}$','fontsize',16)
xlabel('$r/a$','fontsize',16)
lgnd=legend({'Case A-III Micromorphic','Case A-I Micromorphic'...
    ,'Case B-III Micromorphic','Case B-I Micromorphic'},'Interpreter','Latex');
set(gca,'FontName','Times New Roman','FontSize',14)
set(lgnd, 'XColor', 'w', 'YColor', 'w', 'Color', 'none');



Figure12 = figure(12);
box on
hold on
plot(Result_phi_21_MBC_rot(:,1),Result_phi_21_MBC_rot(:,2),'--ob',Result_phi_21_MBC_com(:,1),Result_phi_21_MBC_com(:,2),'--og',...
    Result_phi_21_NMBC_rot(:,1),Result_phi_21_NMBC_rot(:,2),'--ok',...
    Result_phi_21_NMBC_com(:,1),Result_phi_21_NMBC_com(:,2),'--om','LineWidth',2 )
set(Figure12,'defaulttextinterpreter','latex');
ylabel('$\phi_{21}$','fontsize',16)
xlabel('$r/a$','fontsize',16)
lgnd=legend({'Case A-III Micromorphic','Case A-I Micromorphic'...
    ,'Case B-III Micromorphic','Case B-I Micromorphic'},'Interpreter','Latex');
set(gca,'FontName','Times New Roman','FontSize',14)
set(lgnd, 'XColor', 'w', 'YColor', 'w', 'Color', 'none');


