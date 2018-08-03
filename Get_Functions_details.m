%_____________________________________________________________________%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% An improved chaotic gray wolf optimization algorithm (ACGWO)                                         
% 作者和程序员:Tianbiao Yang           
% 主要论文：Bioactive assay and hyphenated chromatography detection for
%		   complex supercritical CO2 extract from Chaihu Shugan San 
% 		   using an experimental design approach                         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [lb,ub,dim,fobj] = Get_Functions_details(F)


switch F      
    case 'F19'
        fobj = @F19;
        lb=-1;
        ub=1;
        dim=3;         
    end

end

% F19

function o = F19(x)
o=-(2.88-0.18.*x(1)+0.44.*x(2)-0.14.*x(3)+0.11.*x(1).*x(2)-0.11.*x(1).*x(3)-0.58.*x(2).*x(3)+0.12.*x(1).*x(1)-0.37.*x(2).*x(2)-0.83.*x(3).*x(3));
end
