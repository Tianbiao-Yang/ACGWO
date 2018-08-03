%___________________________________________________________________%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% An improved chaotic gray wolf optimization algorithm (ACGWO)                                         
% 作者和程序员:Tianbiao Yang           
% 主要论文：Bioactive assay and hyphenated chromatography detection for
%          complex supercritical CO2 extract from Chaihu Shugan San 
%          using an experimental design approach                         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Alpha_score,Alpha_pos,Convergence_curve,K]=AGWO(SearchAgents_no,Max_iter,lb,ub,dim,fobj)

% 初始化α,β,Delta_pos
Alpha_pos=zeros(1,dim);% 初始化Alpha狼的位置
Alpha_score=inf; %初始化Alpha狼的目标函数值,change this to -inf for maximization problems

Beta_pos=zeros(1,dim);%初始化Beta狼的位置
Beta_score=inf; %初始化Beta狼的目标函数值,change this to -inf for maximization problems

Delta_pos=zeros(1,dim);%初始化Delta狼的位置
Delta_score=inf; %初始化Delta狼的目标函数值,change this to -inf for maximization problems

%初始化的位置搜索
Positions=Ainitialization(SearchAgents_no,dim,ub,lb);

Convergence_curve=zeros(1,Max_iter);
K=zeros(1,Max_iter);

l=0;% 循环计数器

% 主循环
while l<Max_iter                  % 对迭代次数循环
    for i=1:size(Positions,1)     % 遍历每个狼  
        
        %若搜索位置超过了搜索空间，需要重新回到搜索空间
        Flag4ub=Positions(i,:)>ub;
        Flag4lb=Positions(i,:)<lb;
        %若狼的位置在最大值和最小值之间，则位置不需要调整，若超出最大值，最回到最大值边界；
        %若超出最小值，最回答最小值边界
        Positions(i,:)=(Positions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;               
        %~表示取反 
        
        % 计算适应度函数值
        fitness=fobj(Positions(i,:));
        
        % 更新Alpha, Beta和Delta的值
        if fitness<Alpha_score       %如果目标函数值小于Alpha狼的目标函数值
            Alpha_score=fitness;     %则将Alpha狼的目标函数值更新为最优目标函数值
            Alpha_pos=Positions(i,:);%则将Alpha狼的目标函数值更新为最优目标函数值
        end
        
        if fitness>Alpha_score && fitness<Beta_score %如果目标函数值介于于Alpha狼和Beta狼的目标函数值之间
            Beta_score=fitness;                      %则将Beta狼的目标函数值更新为最优目标函数值
            Beta_pos=Positions(i,:);                 %同时更新Beta狼的位置
        end
        
        if fitness>Alpha_score && fitness>Beta_score && fitness<Delta_score  % 如果目标函数值介于于Beta狼和Delta狼的目标函数值之间
            Delta_score=fitness;                                             % 则将Delta狼的目标函数值更新为最优目标函数值
            Delta_pos=Positions(i,:);                                        % 同时更新Delta狼的位置
        end
        K=Alpha_score;
    end
    
    
    %a=2-2*((i)/Max_iter); % 对每一次迭代，计算相应的a值，a decreases linearly fron 2 to 0
    a=2-2*((1/(exp(1)-1))*(exp(i/Max_iter)-1));
    
    % 包围猎物，位置更新
    for i=1:size(Positions,1)
        for j=1:size(Positions,2)     
                       
            r1=rand(); % r1 is a random number in [0,1]
            r2=rand(); % r2 is a random number in [0,1]
            
            A1=2*a*r1-a; % 计算系数A，Equation (3.3)
            C1=2*r2;     % 计算系数C，Equation (3.4)
            
            % Alpha狼位置更新
            D_alpha=abs(C1*Alpha_pos(j)-Positions(i,j)); % Equation (3.5)-part 1
            X1=Alpha_pos(j)-A1*D_alpha; % Equation (3.6)-part 1
                       
            r1=rand();
            r2=rand();
            
            A2=2*a*r1-a; % Equation (3.3)
            C2=2*r2;     % Equation (3.4)
            
            % Beta狼位置更新
            D_beta=abs(C2*Beta_pos(j)-Positions(i,j)); % Equation (3.5)-part 2
            X2=Beta_pos(j)-A2*D_beta;                  % Equation (3.6)-part 2       
            
            r1=rand();
            r2=rand(); 
            
            A3=2*a*r1-a; % 计算系数A，Equation (3.3)
            C3=2*r2;     %计算系数C， Equation (3.4)
            
            % Delta狼位置更新
            D_delta=abs(C3*Delta_pos(j)-Positions(i,j)); % Equation (3.5)-part 3
            X3=Delta_pos(j)-A3*D_delta; % Equation (3.5)-part 3             
            
            % 位置更新
            Positions(i,j)=(X1+X2+X3)/3;% Equation (3.7)
            
        end
    end
    l=l+1;    
    Convergence_curve(l)=Alpha_score;
end



