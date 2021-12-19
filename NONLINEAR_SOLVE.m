%% 清空数据
clear;clear global;clc;
%% 设置文件名
for num = 1:1
    global FILE
    FILE.FIRST = ['CYLINDER_WANG_SS', num2str(num), '.cdb'];
    FILE.SECOND = ['LATERAL', num2str(num), '.cdb'];
    %% 设置横向载荷(估值）
    global LPRES
    LPRES = 1000;
    %% 调用计算
    main;
    %% 获取结果
    global RESULT
    % RESULT.WRINKLELOAD 是一个数组，其中最小值为第一次出现褶皱的载荷
    % RESULT.FAILLOAD 是失效载荷
    % RESULT.FAILDISP 是失效位移
    doe_res(num).res = RESULT;
end

% uisave({'doe_res'}, 'doe_res')
