%% �������
clear;clear global;clc;
%% �����ļ���
for num = 1:1
    global FILE
    FILE.FIRST = ['CYLINDER_WANG_SS', num2str(num), '.cdb'];
    FILE.SECOND = ['LATERAL', num2str(num), '.cdb'];
    %% ���ú����غ�(��ֵ��
    global LPRES
    LPRES = 1000;
    %% ���ü���
    main;
    %% ��ȡ���
    global RESULT
    % RESULT.WRINKLELOAD ��һ�����飬������СֵΪ��һ�γ���������غ�
    % RESULT.FAILLOAD ��ʧЧ�غ�
    % RESULT.FAILDISP ��ʧЧλ��
    doe_res(num).res = RESULT;
end

% uisave({'doe_res'}, 'doe_res')
