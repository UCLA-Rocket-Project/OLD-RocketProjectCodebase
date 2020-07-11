clear all; close all; clc;

%% Inputs
Fs = 14;
Ms = 12;

G_g = 0.9669;
gam_g = 1.4;

Cv_reg = 0.8;
Cv_iso = 6.45;
Cv_check = 1.47;
Cv_relief = 1.65;

T_g = 529;
p_bottle = 2000;
p_set = 500;
p_relief = p_set*1.125;
p_atm = 14.7;

%% Calculations
p_reg2iso = p_set;
p_iso2check = p_set;
p_check2relief = p_set;

q_reg = ValveGasFlow(Cv_reg,p_bottle,p_reg2iso,G_g,T_g,gam_g);
q_iso = ValveGasFlow(Cv_iso,p_reg2iso,p_iso2check,G_g,T_g,gam_g);
q_check = ValveGasFlow(Cv_check,p_iso2check,p_check2relief,G_g,T_g,gam_g)*2;
q_relief = ValveGasFlow(Cv_relief,p_check2relief,p_atm,G_g,T_g,gam_g)*2;

q_reg_arr(1) = q_reg;
q_iso_arr(1) = q_iso;
q_check_arr(1) = q_check;
q_relief_arr(1) = q_relief;

i = 1;

while ((abs(q_reg - q_iso) > 1) && (abs(q_reg - q_check) > 1) && (abs(q_reg - q_relief) > 1))
    while (abs(q_reg - q_relief) > 1)
        if (q_reg > q_relief)
            p_check2relief = p_check2relief + 0.01*abs(q_reg - q_relief);
            p_iso2check = p_check2relief;
        else
            p_check2relief = p_check2relief -  0.01*abs(q_reg - q_relief);
        end
        if (p_check2relief > p_bottle)
            error('Relief Valve Cv too small')
        end
        q_relief = ValveGasFlow(Cv_relief,p_check2relief,p_atm,G_g,T_g,gam_g)*2;
    end
    q_check = ValveGasFlow(Cv_check,p_iso2check,p_check2relief,G_g,T_g,gam_g);
    while (abs(q_reg - q_check) > 1)
        if (q_reg > q_check)
            p_iso2check = p_iso2check + 0.01*abs(q_reg - q_check);
            p_reg2iso = p_iso2check;
        else
            p_iso2check = p_iso2check - 0.01*abs(q_reg - q_check);
        end
        if (p_iso2check > p_bottle)
            error('Check Valve Cv too small')
        end
        q_check = ValveGasFlow(Cv_check,p_iso2check,p_check2relief,G_g,T_g,gam_g)*2;
    end
    q_iso = ValveGasFlow(Cv_iso,p_reg2iso,p_iso2check,G_g,T_g,gam_g);
    while (abs(q_reg - q_iso) > 1)
        if (q_reg > q_iso)
            p_reg2iso = p_reg2iso + 1.01*abs(q_reg - q_iso);
        else
            p_reg2iso = p_reg2iso - 0.01*abs(q_reg - q_iso);
        end
        q_iso = ValveGasFlow(Cv_iso,p_reg2iso,p_iso2check,G_g,T_g,gam_g);
    end
    if (p_reg2iso > p_bottle)
        error('Iso Valve Cv too small')
    end
    q_reg = ValveGasFlow(Cv_reg,p_bottle,p_reg2iso,G_g,T_g,gam_g);
end

%     if i == 1
%         figure(1)
%             set(gcf,'position',[100,100,750,500])
%             set(gca,'fontsize',Fs)
%             ylabel('Std Volume Flow Rate [SCFM]')
%             xlabel('Iteration')
%             title('Convergence Plot for Pressurant Flow Rate')
%             hold on
%             plot(i,q_reg_arr(i),'.k','MarkerSize',Ms)
%             plot(i,q_iso_arr(i),'.g','MarkerSize',Ms)
%             plot(i,q_check_arr(i),'.b','MarkerSize',Ms)
%             plot(i,q_relief_arr(i),'.r','MarkerSize',Ms)
%             hold off
%             legend('Reg','Iso','Check','Relief')
%     end
% 
%     q_reg_arr(i) = q_reg;
%     q_iso_arr(i) = q_iso;
%     q_check_arr(i) = q_check;
%     q_relief_arr(i) = q_relief;
%     
% 
%     if i > 1
%         figure(1)
%             hold on
%             plot(i,q_reg_arr(i),'.b','MarkerSize',Ms)
%             plot(i,q_iso_arr(i),'.g','MarkerSize',Ms)
%             plot(i,q_check_arr(i),'b.','MarkerSize',Ms)
%             plot(i,q_relief_arr(i),'.r','MarkerSize',Ms)
%             hold off
%             drawnow
%     end
%     i = i + 1;