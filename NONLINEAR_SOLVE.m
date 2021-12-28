%% �������?
clear;clear global;clc;
%% �����ļ���
geo = csvread('./geo.csv');

for num = 1:size(geo, 1)
    global FILE
    FILE.FIRST = ['CYLINDER_WANG_SS', num2str(num), '.cdb'];
    FILE.SECOND = ['LATERAL', num2str(num), '.cdb'];
    %% ���ú����غ�(��ֵ��
    global LPRES
    LPRES = 3000;
    %% ���ü���
    main;
    %% ��ȡ���?
    global RESULT
    % RESULT.WRINKLELOAD ��һ�����飬������СֵΪ��һ�γ���������غ�?
    % RESULT.FAILLOAD ��ʧЧ�غ�
    % RESULT.FAILDISP ��ʧЧλ��
    doe_res(num).res = RESULT;
    doe_res(num).u = RESU;
    doe_res(num).load = RESLOAD;
    save('doe_res', 'doe_res')
end

% save('doe_res', 'doe_res')
% uisave({'doe_res'}, 'doe_res')
