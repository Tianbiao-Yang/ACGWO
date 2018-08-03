%_____________________________________________________________________%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% An improved chaotic gray wolf optimization algorithm (ACGWO)                                         
% 作者和程序员:Tianbiao Yang           
% 主要论文：Bioactive assay and hyphenated chromatography detection for
%		   complex supercritical CO2 extract from Chaihu Shugan San 
% 		   using an experimental design approach                         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 你需要的初始参数：
%__________________________________________
% fobj = @YourCostFunction
% dim=变量的数量
% Max_iteration = 最大迭代数
% SearchAgents_no = 狼群数量
% lb=[lb1,lb2,...,lbn] 参数取值下界
% ub=[ub1,ub2,...,ubn] 参数取值上界
% 运行GWO: 
%[Best_score,Best_pos,GWO_cg_curve]
% 				=ACGWO(SearchAgents_no,Max_iteration,lb,ub,dim,fobj)
%_____________________________________________________________________%
clear all 
clc
SearchAgents_no=20; % 狼群数量30
Function_name='F19'; %测试函数的名称,可以用论文中的从F1到F23(表1、2、3)
Max_iteration=100; % 最大迭代数为500
% 加载选择的基准函数的详细信息 
[lb,ub,dim,fobj]=Get_Functions_details(Function_name);
%基本灰狼优化GWO
[Best_score,Best_pos,GWO_cg_curve,K]=GWO(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);
x1=Best_pos(1,1)*20/2+50;
x2=Best_pos(1,2)*100+250;
x3=Best_pos(1,3)*0.225+0.375;
X=[x1,x2,x3];
YMAX=max(GWO_cg_curve);
YMIN=min(GWO_cg_curve);
YMAIN=YMAX-YMIN;
%改善收敛因子的灰狼优化AGWO
[ABest_score,ABest_pos,AGWO_cg_curve,K]=AGWO(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);
ax1=ABest_pos(1,1)*20/2+50;
ax2=ABest_pos(1,2)*100+250;
ax3=ABest_pos(1,3)*0.225+0.375;
X=[ax1,ax2,ax3];
AYMAX=max(AGWO_cg_curve);
AYMIN=min(AGWO_cg_curve);
AYMAIN=AYMAX-AYMIN;
%混沌灰狼优化CGWO
[CBest_score,CBest_pos,CGWO_cg_curve,K]=CGWO(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);
cx1=Best_pos(1,1)*20/2+50;
cx2=Best_pos(1,2)*100+250;
cx3=Best_pos(1,3)*0.225+0.375;
X=[cx1,cx2,cx3];
CYMAX=max(CGWO_cg_curve);
CYMIN=min(CGWO_cg_curve);
CYMAIN=CYMAX-CYMIN;
%最终改进后的灰狼优化ACGWO
[ACBest_score,ACBest_pos,ACGWO_cg_curve,K]=ACGWO(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);
acx1=Best_pos(1,1)*20/2+50;
acx2=Best_pos(1,2)*100+250;
acx3=Best_pos(1,3)*0.225+0.375;
ACX=[acx1,acx2,acx3];
ACYMAX=max(ACGWO_cg_curve);
ACYMIN=min(ACGWO_cg_curve);
ACYMAIN=ACYMAX-ACYMIN;
[gBestScore,PSO_cg_curve]=PSO(SearchAgents_no,Max_iteration,lb,ub,dim,fobj); % run PSO to compare to results
%绘制搜索空间 
figure('Position',[500 500 500 290])
subplot(1,1,1);
semilogy(GWO_cg_curve,'Color','r')
hold on
semilogy(CGWO_cg_curve,'Color','b')
hold on
semilogy(AGWO_cg_curve,'Color','k')
hold on
semilogy(ACGWO_cg_curve,'Color','g')
hold on
semilogy(PSO_cg_curve,'Color','m')
%title('Objective space')
xlabel('迭代次数');
ylabel('目标函数值');
axis tight
grid on
box on
legend('GWO','CGWO','AGWO','ACGWO','PSO')
display(['The best solution obtained by ACGWO is : ', num2str(ACX)]);
display(['The best optimal value of the objective funciton found by ACGWO is : ', num2str(ACBest_score)]);

        



