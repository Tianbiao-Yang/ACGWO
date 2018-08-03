%_____________________________________________________________________%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% An improved chaotic gray wolf optimization algorithm (ACGWO)                                         
% 作者和程序员:Tianbiao Yang           
% 主要论文：Bioactive assay and hyphenated chromatography detection for
%      complex supercritical CO2 extract from Chaihu Shugan San 
%        using an experimental design approach                         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function initialize the first population of search agents
function Positions=initialization1(SearchAgents_no,dim,ub,lb)
Boundary_no= size(ub,2); % numnber of boundaries
for t=1:20 %次数
      %1生成
      cxl=rand(SearchAgents_no,dim);
      for j=1:dim
          if cxl(j)==0
            cxl(j)=0.1;
          end
          if cxl(j)==0.25
             cxl(j)=0.26;
          end
          if cxl(j)==0.5
             cxl(j)=0.51;
          end
          if cxl(j)==0.75
             cxl(j)=0.76;
          end
          if cxl(j)==1
             cxl(j)=0.9;
          end
      end
end
if Boundary_no==1
    Positions=cxl.*(ub-lb)+lb;
end
% If each variable has a different lb and ub
if Boundary_no>1
    for i=1:dim
        ub_i=ub(i);
        lb_i=lb(i);
        Positions(:,i)=cxl.*(ub_i-lb_i)+lb_i;
    end
end