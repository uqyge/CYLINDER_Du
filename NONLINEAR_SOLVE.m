%% �������
clear;clear global;clc;
%% �����ļ���
num = 1;
global FILE
FILE.FIRST = ['CYLINDER_WANG_SS',num2str(num),'.cdb'];
FILE.SECOND = ['LATERAL',num2str(num),'.cdb'];
%% ���ú����غ�(��ֵ��
global LPRES
LPRES = 3000;
%% ���ü���
main;
%% ��ȡ���
global RESULT
% RESULT.WRINKLELOAD ��һ�����飬������СֵΪ��һ�γ���������غ�
% RESULT.FAILLOAD ��ʧЧ�غ�
% RESULT.FAILDISP ��ʧЧλ��
