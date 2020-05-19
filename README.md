# match_gjf
match two structures saved as gjf format

编译 ifort match_gjf.f90 -o match

使用 match 1.gjf 2.gjf match_in
其中match为可执行文件，1.gjf与2.gjf为保存了连接度信息的gjf文件，match_in规定匹配原子的起点

例子：两个甲基苯结构的匹配
