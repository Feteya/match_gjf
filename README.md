# match_gjf
match two structures saved as gjf format

编译 ifort match_gjf.f90 -o match

使用 match 1.gjf 2.gjf match_in
其中match为可执行文件，1.gjf与2.gjf为保存了连接度信息的gjf文件，match_in规定匹配原子的起点

例子：两个甲基苯结构的匹配

分子1的结构

%chk=1.chk
# hf/3-21g geom=connectivity

Title Card Required

0 1
 C                  0.33188248    1.65397182    0.00000000
 C                  1.72704248    1.65397182    0.00000000
 C                  2.42458048    2.86172282    0.00000000
 C                  1.72692648    4.07023182   -0.00119900
 C                  0.33210148    4.07015382   -0.00167800
 C                 -0.36549952    2.86194782   -0.00068200
 H                 -0.21787652    0.70165482    0.00045000
 H                  2.27655048    0.70145882    0.00131500
 H                  3.52426048    2.86180282    0.00063400
 H                 -0.21802052    5.02243482   -0.00263100
 H                 -1.46510352    2.86213082   -0.00086200
 C                  2.49743055    5.40361979   -0.00128162
 H                  2.69640187    5.70063858    1.00722730
 H                  3.42185001    5.28134019   -0.52606382
 H                  1.90938973    6.15532511   -0.48506578

 1 2 1.5 6 1.5 7 1.0
 2 3 1.5 8 1.0
 3 4 1.5 9 1.0
 4 5 1.5 12 1.0
 5 6 1.5 10 1.0
 6 11 1.0
 7
 8
 9
 10
 11
 12 13 1.0 14 1.0 15 1.0
 13
 14
 15


分子2的结构

%chk=1.chk
# hf/3-21g geom=connectivity

Title Card Required

0 1
 H                 -0.21787652    0.70165482    0.00045000
 H                  2.27655048    0.70145882    0.00131500
 H                  3.52426048    2.86180282    0.00063400
 H                 -0.21802052    5.02243482   -0.00263100
 H                 -1.46510352    2.86213082   -0.00086200
 C                  2.49743055    5.40361979   -0.00128162
 H                  2.69640187    5.70063858    1.00722730
 H                  3.42185001    5.28134019   -0.52606382
 H                  1.90938973    6.15532511   -0.48506578
 C                  0.33188248    1.65397182    0.00000000
 C                  1.72704248    1.65397182    0.00000000
 C                  2.42458048    2.86172282    0.00000000
 C                  1.72692648    4.07023182   -0.00119900
 C                  0.33210148    4.07015382   -0.00167800
 C                 -0.36549952    2.86194782   -0.00068200

 1 10 1.0
 2 11 1.0
 3 12 1.0
 4 14 1.0
 5 15 1.0
 6 7 1.0 8 1.0 9 1.0 13 1.0
 7
 8
 9
 10 11 1.5 15 1.5
 11 12 1.5
 12 13 1.5
 13 14 1.5
 14 15 1.5
 15


match_in 输入
1 10

运行
match 1.gjf 2.gjf match_in
