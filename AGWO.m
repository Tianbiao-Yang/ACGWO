%___________________________________________________________________%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% An improved chaotic gray wolf optimization algorithm (ACGWO)                                         
% ���ߺͳ���Ա:Tianbiao Yang           
% ��Ҫ���ģ�Bioactive assay and hyphenated chromatography detection for
%          complex supercritical CO2 extract from Chaihu Shugan San 
%          using an experimental design approach                         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Alpha_score,Alpha_pos,Convergence_curve,K]=AGWO(SearchAgents_no,Max_iter,lb,ub,dim,fobj)

% ��ʼ����,��,Delta_pos
Alpha_pos=zeros(1,dim);% ��ʼ��Alpha�ǵ�λ��
Alpha_score=inf; %��ʼ��Alpha�ǵ�Ŀ�꺯��ֵ,change this to -inf for maximization problems

Beta_pos=zeros(1,dim);%��ʼ��Beta�ǵ�λ��
Beta_score=inf; %��ʼ��Beta�ǵ�Ŀ�꺯��ֵ,change this to -inf for maximization problems

Delta_pos=zeros(1,dim);%��ʼ��Delta�ǵ�λ��
Delta_score=inf; %��ʼ��Delta�ǵ�Ŀ�꺯��ֵ,change this to -inf for maximization problems

%��ʼ����λ������
Positions=Ainitialization(SearchAgents_no,dim,ub,lb);

Convergence_curve=zeros(1,Max_iter);
K=zeros(1,Max_iter);

l=0;% ѭ��������

% ��ѭ��
while l<Max_iter                  % �Ե�������ѭ��
    for i=1:size(Positions,1)     % ����ÿ����  
        
        %������λ�ó����������ռ䣬��Ҫ���»ص������ռ�
        Flag4ub=Positions(i,:)>ub;
        Flag4lb=Positions(i,:)<lb;
        %���ǵ�λ�������ֵ����Сֵ֮�䣬��λ�ò���Ҫ���������������ֵ����ص����ֵ�߽磻
        %��������Сֵ����ش���Сֵ�߽�
        Positions(i,:)=(Positions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;               
        %~��ʾȡ�� 
        
        % ������Ӧ�Ⱥ���ֵ
        fitness=fobj(Positions(i,:));
        
        % ����Alpha, Beta��Delta��ֵ
        if fitness<Alpha_score       %���Ŀ�꺯��ֵС��Alpha�ǵ�Ŀ�꺯��ֵ
            Alpha_score=fitness;     %��Alpha�ǵ�Ŀ�꺯��ֵ����Ϊ����Ŀ�꺯��ֵ
            Alpha_pos=Positions(i,:);%��Alpha�ǵ�Ŀ�꺯��ֵ����Ϊ����Ŀ�꺯��ֵ
        end
        
        if fitness>Alpha_score && fitness<Beta_score %���Ŀ�꺯��ֵ������Alpha�Ǻ�Beta�ǵ�Ŀ�꺯��ֵ֮��
            Beta_score=fitness;                      %��Beta�ǵ�Ŀ�꺯��ֵ����Ϊ����Ŀ�꺯��ֵ
            Beta_pos=Positions(i,:);                 %ͬʱ����Beta�ǵ�λ��
        end
        
        if fitness>Alpha_score && fitness>Beta_score && fitness<Delta_score  % ���Ŀ�꺯��ֵ������Beta�Ǻ�Delta�ǵ�Ŀ�꺯��ֵ֮��
            Delta_score=fitness;                                             % ��Delta�ǵ�Ŀ�꺯��ֵ����Ϊ����Ŀ�꺯��ֵ
            Delta_pos=Positions(i,:);                                        % ͬʱ����Delta�ǵ�λ��
        end
        K=Alpha_score;
    end
    
    
    %a=2-2*((i)/Max_iter); % ��ÿһ�ε�����������Ӧ��aֵ��a decreases linearly fron 2 to 0
    a=2-2*((1/(exp(1)-1))*(exp(i/Max_iter)-1));
    
    % ��Χ���λ�ø���
    for i=1:size(Positions,1)
        for j=1:size(Positions,2)     
                       
            r1=rand(); % r1 is a random number in [0,1]
            r2=rand(); % r2 is a random number in [0,1]
            
            A1=2*a*r1-a; % ����ϵ��A��Equation (3.3)
            C1=2*r2;     % ����ϵ��C��Equation (3.4)
            
            % Alpha��λ�ø���
            D_alpha=abs(C1*Alpha_pos(j)-Positions(i,j)); % Equation (3.5)-part 1
            X1=Alpha_pos(j)-A1*D_alpha; % Equation (3.6)-part 1
                       
            r1=rand();
            r2=rand();
            
            A2=2*a*r1-a; % Equation (3.3)
            C2=2*r2;     % Equation (3.4)
            
            % Beta��λ�ø���
            D_beta=abs(C2*Beta_pos(j)-Positions(i,j)); % Equation (3.5)-part 2
            X2=Beta_pos(j)-A2*D_beta;                  % Equation (3.6)-part 2       
            
            r1=rand();
            r2=rand(); 
            
            A3=2*a*r1-a; % ����ϵ��A��Equation (3.3)
            C3=2*r2;     %����ϵ��C�� Equation (3.4)
            
            % Delta��λ�ø���
            D_delta=abs(C3*Delta_pos(j)-Positions(i,j)); % Equation (3.5)-part 3
            X3=Delta_pos(j)-A3*D_delta; % Equation (3.5)-part 3             
            
            % λ�ø���
            Positions(i,j)=(X1+X2+X3)/3;% Equation (3.7)
            
        end
    end
    l=l+1;    
    Convergence_curve(l)=Alpha_score;
end



