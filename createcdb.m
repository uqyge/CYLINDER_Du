%% 驱动脚本
function [] = createcdb(model,path)
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
%% 路径说明
% .ANSYSPATH: ansys软件路径
%% 生成APDL文件
fid = fopen('MODEL_AUTO.mac','w+');
fprintf(fid,'finish\n/clear\n/prep7\net,1,shell63\nKEYOPT,1,1,1\nKEYOPT,1,2,0\nKEYOPT,1,3,1\nKEYOPT,1,5,0\nKEYOPT,1,6,0\nKEYOPT,1,7,0\nKEYOPT,1,8,0\nKEYOPT,1,9,0\nKEYOPT,1,11,0\n\n');
fprintf(fid,'E = %8.3e\n',model.E);
fprintf(fid,'mu = %8.3e\n',model.mu);
fprintf(fid,'t = %8.3e\n',model.t);
fprintf(fid,'presload = %8.3e\n',model.presload);
fprintf(fid,'L = %8.3e\n',model.L);
fprintf(fid,'R = %8.3e\n',model.R);
fprintf(fid,'pw = %8.3e\n',model.pw);
fprintf(fid,'elesize = %8.3e\n',model.elesize);
fprintf(fid,'FAXIS = %8.3e\n',model.FAXIS);
fprintf(fid,'MPTEMP,,,,,,,,\nMPTEMP,1,0\nMPDATA,EX,1,,E\nMPDATA,NUXY,1,,mu\nMPDATA,PRXY,1,,mu\n\n');
fprintf(fid,'R,1,t,t,t,t, , ,\nRMORE, , , ,\nRMORE\nRMORE,,\n\n');
fprintf(fid,'K, ,R,0,0,\nK, ,R,L/2-pw/2,0,\nK, ,R,L/2+pw/2,0,\nK, ,R,L,0,\nLSTR,1,2\nLSTR,2,3\nLSTR,3,4\n');
fprintf(fid,'K,,0,0,0\nK,,0,0.1,0\n\n');
fprintf(fid,'FLST,2,3,4,ORDE,2\nFITEM,2,1\nFITEM,2,-3\nFLST,2,3,4,ORDE,2\nFITEM,2,1\nFITEM,2,-3\nFLST,8,2,3\nFITEM,8,5\nFITEM,8,6\nAROTAT,P51X, , , , , ,P51X, ,360,4,\n\n');
fprintf(fid,"MSHAPE,1,2D\nMSHKEY,1\nESIZE,elesize,0,\nFLST,5,20,5,ORDE,4\nFITEM,5,1\nFITEM,5,-12\nFITEM,5,21\nFITEM,5,-28\nCM,_Y,AREA\nASEL, , , ,P51X\nCM,_Y1,AREA\nCHKMSH,\'AREA\'\nCMSEL,S,_Y\nAMESH,_Y1\nCMDELE,_Y\nCMDELE,_Y1\nCMDELE,_Y2\n\n");
fprintf(fid,'NUMMRG,ALL, , , ,LOW\nNUMCMP,ALL\n\n');
fprintf(fid,'FINISH\n/SOL\nANTYPE,0\nCNVTOL,F, ,1E-6,2, ,\nANTYPE,0\nNLGEOM,1\nNSUBST,10,0,0\nOUTRES,ERASE\nOUTRES,ALL,ALL\nAUTOTS,0\nLNSRCH,0\nNEQIT,100000\nPRED,0\nTIME,1\n\n');
fprintf(fid,'NSEL,S,LOC,Y,-1E-6,1E-6\nNSEL,R,LOC,Z,-R-1E-6,1E-6\nD,ALL,ALL\nALLSEL,ALL\n\n');
fprintf(fid,'NSEL,S,LOC,Y,L-1E-6,L+1E-6\nNSEL,R,LOC,Z,-R-1E-6,1E-6\nD,ALL, , , , , ,UX,UZ, , , ,\nALLSEL,ALL\n\n');
fprintf(fid,'ASEL,S,LOC,Y,1E-6,L-1E-6\nSFA,ALL,1,PRES,presload\nALLSEL,ALL\n\n');
fprintf(fid,'NSEL,S,LOC,Y,-1E-6,1E-6\n*get,nn,node,,count\nALLSEL,ALL\nfnn = FAXIS/nn\n\n');
fprintf(fid,'NSEL,S,LOC,Y,-1E-6,1E-6\nF,ALL,FY,-fnn\nALLSEL,ALL\n\n');
fprintf(fid,'NSEL,S,LOC,Y,L-1E-6,L+1E-6\nF,ALL,FY,fnn\nALLSEL,ALL\n\n');
fprintf(fid,"FINISH\n/PREP7\nCDWRITE,DB,'CYLINDER_WANG_SS%d','cdb',,'',''\n\n",model.num);
fprintf(fid,'FLST,5,2,5,ORDE,2\nFITEM,5,8\nFITEM,5,11\nASEL,S, , ,P51X\nESLA,S\n\n');
fprintf(fid,"CDWRITE,DB,'LATERAL%d','cdb',,'',''\n\n",model.num);
fclose(fid);
%% 调用ANSYS批处理
ws = pwd;
command = ['"',path.ANSYSPATH,'"',' -dir ','"',ws,'"',' -p ','ansys',' -j ','"','file','"',' -s read -l en-us -b ','-i ','"',ws,'\','MODEL_AUTO.mac','"',' -o ', '"',ws,'\file.out','"'];
system(command);
end

