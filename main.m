global FILE
global RESULT
RESULT.WRINKLELOAD = [];

format short e
importcdb_shell(FILE.FIRST);
[s, DI, e1, e2] = initialize_trishell();
calc_quadrature();
calc_localcoord_initial();
calc_localcoord_current();
[Be, te, Ae] = calc_local_stiffness([1E-10; 1E-10; 0]); % Change INITIAL STRAIN
plotvtk_trishell(0, 1);

%% Newton-Loop
global PARA
global FEXT
global FINT
global K
global U
global NODE
global ELEM
global FEXTL
global FG
global MAXDISP

PARA.EPS = 1e-6;

%  FEXTL = FEXTL ;

RESITER = zeros(0);
RESU = zeros(PARA.SUBSTEPS, 1);
RESLOAD = zeros(PARA.SUBSTEPS, 1);
RESDOF = (3143 - 1) * 3 + 3;

%% Modify Additional Tip force
FEXT = FEXT * 0.84;

%% Modify Gravity
FG = FEXTL * 0.0;
A = zeros(2);
%% First Load Step - Inflation
tic;
SUBSTEPS_LATERAL = PARA.SUBSTEPS;
PARA.SUBSTEPS = 1;

for incr = 1:PARA.SUBSTEPS
    display(incr);
    lambda = incr / PARA.SUBSTEPS;
    FEXTi = FEXT * lambda + FG;
    iter = 1;

    while (1)
        e2 = e1;
        [lambdas, strain_matrix] = calc_stiff_intforce_Du(e1, e2, s, DI, Ae, Be, te);
        %          [ lambdas, strain_matrix ] = calc_stiff_intforce(e1,e2,s,DI,Ae,Be,te);
        [Fpres] = calcPressureForce(lambda);
        calcPressureStiff(lambda);
        Fubl = FEXTi + Fpres - FINT;
        dU = K \ (Fubl);
        U_TEMP = U;
        ss = 1.0;

        U = U_TEMP + ss * dU;
        update_config();

        if (iter ~= 1)
            iter;
            eps = norm(Fubl) / norm(FEXT + Fpres + FG);
            RESITER(end + 1) = eps;
            semilogy(RESITER);
            %             drawnow();
            display(eps);

            if (eps < PARA.EPS)
                %                 postdisposal(NODE(1:3,:)',U,ELEM(1:3,:)',PARA.NELEM,PARA.NNODE,s,lambdas,DI,strain_matrix,incr);
                %                 plotvtk_trishell(incr,1);
                RESU(incr) = U(RESDOF);
                RESLOAD(incr) = lambda;
                break;
            end

        end

        calc_localcoord_current();
        iter = iter + 1;

        % if (iter > 20)
        %     break
        % end

    end

end

a = toc;
A(1) = a;
tic;
PARA.SUBSTEPS = SUBSTEPS_LATERAL * 3;
%% Second Load Step - Tip Load
incr = 1;
endcalculation = 0;

% if (iter > 20)
%     endcalculation = 1
% end

while (~endcalculation)
    display(incr);
    lambda = incr / PARA.SUBSTEPS;
    FEXTi = FEXT + FG + FEXTL * lambda;
    iter = 1;

    while (1)
        % 每次迭代前清空RESULT中的褶皱判据
        RESULT.ISWRINKLE = 0;
        e2 = 0.000001; % Material Modification
        %         e2=e1;
        [lambdas, strain_matrix] = calc_stiff_intforce_Du(e1, e2, s, DI, Ae, Be, te);
        %          [ lambdas, strain_matrix ] = calc_stiff_intforce(e1,e2,s,DI,Ae,Be,te);
        [Fpres] = calcPressureForce(1);
        calcPressureStiff(1);
        Fubl = FEXTi + Fpres - FINT;

        dU = K \ (Fubl);
        U = U + dU;
        update_config();

        if (iter ~= 1)
            iter;
            eps = norm(Fubl) / norm(FEXT + Fpres + FEXTL + FG);
            RESITER(end + 1) = eps;
            semilogy(RESITER);
            %             drawnow();
            display(eps);

            if (eps < PARA.EPS)
                %                 postdisposal(NODE(1:3,:)',U,ELEM(1:3,:)',PARA.NELEM,PARA.NNODE,s,lambdas,DI,strain_matrix,incr);
                %                 plotvtk_trishell(incr+1,1);
                RESU(incr + 1) = U(RESDOF);
                RESLOAD(incr + 1) = lambda * sum(FEXTL);

                if (iter >= 3)
                    RESULT.WRINKLELOAD(end + 1) = abs(lambda * sum(FEXTL));
                end

                % if (RESULT.ISWRINKLE == 1)
                %     RESULT.WRINKLELOAD(end + 1) = abs(lambda * sum(FEXTL));
                % end

                break;
            elseif (eps > 1E4 || iter > 20) % 终止条件，不收敛，退出计算
                RESULT.FAILLOAD = abs(sum(RESLOAD(incr)));
                RESULT.FAILDISP = max(abs(RESU(incr)));
                endcalculation = 1;
                break;
            end

        end

        calc_localcoord_current();
        iter = iter + 1;
    end

    incr = incr + 1;
end

b = toc;
A(2) = b;

% %% Third Load Step - Recover Force
% recovlist = [33,34,35,36,37,38,198,199,200,201,202,203,264,265,266,267,268,269,330,331,332,333,334,335,396,397,398,399,400,401,1342,1343,1344,1345,1346,1347,1408,1409,1410,1411,1412,1413,1474,1475,1476,1477,1478,1479,1540,1541,1542,1543,1544,1545];
% noderecf = -30/size(recovlist,2);
% FREC = zeros(PARA.NNODE*3,1);
% for i=1:size(recovlist,2)
%     FREC(3*recovlist(i))=noderecf;
% end
%
% for incr=1:PARA.SUBSTEPS
%     display(incr);
%     lambda = incr/PARA.SUBSTEPS;
%     FEXTi=FEXT+FG+FEXTL+FREC*lambda;
%     iter=1;
%     while(1)
%         e2=1E-4; % Material Modification
%         [ lambdas, strain_matrix ] = calc_stiff_intforce_Dong(e1,e2,s,DI,Ae,Be,te);
%         [ Fpres ] = calcPressureForce( 1 );
%         calcPressureStiff( 1 );
%         Fubl=FEXTi+Fpres-FINT;
%
%         dU = K\(Fubl);
%         U = U + dU;
%         update_config();
%
%         if(iter~=1)
%             iter;
%             eps=norm(Fubl)/norm(FEXT+Fpres+FEXTL+FG+FREC);
%             RESITER(end+1)=eps;
%             semilogy(RESITER);
%             drawnow();
%             display(eps);
%
%             if(eps<PARA.EPS)
% %                 postdisposal(NODE(1:3,:)',U,ELEM,PARA.NELEM,PARA.NNODE,s,lambdas,DI,strain_matrix);
%                 postdisposal3(NODE(1:3,:)',U,ELEM(1:3,:)',PARA.NELEM,PARA.NNODE,s,lambdas,DI,strain_matrix,incr);
%                 plotvtk_trishell(incr+PARA.SUBSTEPS*2,1);
%                 RESU(incr+PARA.SUBSTEPS*2)=U(RESDOF);
%                 RESLOAD(incr+PARA.SUBSTEPS*2)=lambda+2;
%                 break;
%             end
%         end
%
%         calc_localcoord_current();
%         iter=iter+1;
%     end
% end
