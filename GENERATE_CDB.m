%% ÿ�ε���ǰ���ȫ������
clear;clear global;clc;
%% ����˵��(model)
% .E: ����ģ��
% .mu: ���ɱ�
% .t: Ĥ���
% .presload: Ԥѹ���غ�
% .L: �ܳ���
% .R: �ܰ뾶
% .pw: �����غɼ��ؿ��
% .elesize: ��Ԫ�ߴ�
% .FAXIS: �����غ�
% .MAXDISP: ���ʧЧ�Ӷ�
% .LPRES: ����ѹ��
% .num: ģ���ļ����
%% ģ�Ͳ�������
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
%% ·������
path.ANSYSPATH = 'D:\programs\ANSYS Inc\v150\ansys\bin\winx64\ANSYS150.exe';
%% ���ú���
createcdb(model,path)

