%% ÿ�ε���ǰ���ȫ������
clear;clear global;clc;
%% ����˵��(model)
% .E: ����ģ��
% .mu: ���ɱ�
% .t: Ĥ���
% .presload: Ԥѹ���غ�
% .L: �ܳ���
% .R: �ܰ뾶
% .pw: �����غɼ��ؿ���
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
% path.ANSYSPATH = 'D:\programs\ANSYS Inc\v150\ansys\bin\winx64\ANSYS150.exe';
path.ANSYSPATH = 'C:\Program Files\ANSYS Inc\v180\ansys\bin\winx64\ANSYS180.exe'
%% ���ú���

DOE.num = [1, 2, 3];
DOE.t = [25E-6, 25E-6, 25E-6];
DOE.R = [0.03, 0.04, 0.05]
DOE.presload = [5E4, 2E5]
csv_in = csvread('./geo.csv')

for i = 1:size(csv_in, 1)
    model.num = i;
    model.t = csv_in(i, 1);
    model.R = csv_in(i, 2);
    createcdb(model, path)
end
