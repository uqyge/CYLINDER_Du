%% 每次调用前清空全部变量
clear;clear global;clc;
%% 参数说明(model)
% .E: 弹性模量
% .mu: 泊松比
% .t: 膜厚度
% .presload: 预压力载荷
% .L: 管长度
% .R: 管半径
% .pw: 横向载荷加载宽度
% .elesize: 单元尺寸
% .FAXIS: 轴向载荷
% .MAXDISP: 最大失效挠度
% .LPRES: 横向压力
% .num: 模型文件编号
%% 模型参数定义
model.E = 4.9E9;
model.mu = 0.34;
model.t = 25E-6;
model.presload = 1E4;
model.L = 0.6;
model.R = 0.03;
model.pw = 0.01;
model.elesize = 0.005;
model.FAXIS = 20;
model.num = 1;
%% 路径定义
path.ANSYSPATH = 'D:\programs\ANSYS Inc\v150\ansys\bin\winx64\ANSYS150.exe';
%% 调用函数

DOE.num = [1, 2, 3];
DOE.t = [25E-6, 25E-6, 25E-6];
DOE.R = [0.03, 0.04, 0.05]
DOE.presload = [5E4, 2E5]
csv_in = csvread('./geo.csv')

for i = 1:3
    model.num = i;
    model.t = csv_in(i, 1);
    model.R = csv_in(i, 2);
    createcdb(model, path)
end
