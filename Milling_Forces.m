% Milling_Force.m
% Code for milling tool force model presented in T. A. Dow, E. L. Miller, and K. Garrard, "Tool force and deflection compensation for small milling tools," 
% Precision Engineering, vol. 28, no. 1, pp. 31–45, Jan. 2004.
%
% Written by Kazi Sher Ahmed (kazisherahmed@gmail.com) - July 2024
% License: MIT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear; clc;

%% USER-DEFINED PARAMETERS  

    H    = 5.5e9;           % Hardness of workpiece material [Pa]
    E    = 210e9;           % Elastic modulus of workpiece [Pa]    
    mu   = 0.46;            % Friction coefﬁcient between tool rake face and workpiece
    mu_f = 0.2;             % Friction coefﬁcient between ﬂank face and workpiece
    ga   = 50*pi/180;       % Shear angle of the chip [rad]
    
    d    = 5e-6;           % Depth of cut [m]        
    f    = 7*1e-6;          % Up-feed per flute [m/flute]
    R_t  = 0.2e-3;          % Radius of tool [m]
    W_l  = 8e-6;            % Width of wear land [m]
    alp  = 0;               % Sweep angle is between the rotation axis of the tool and a line normal to the workpiece surface [rad]
    zeta = 0;               % Offset of one flute w.r.t. the other caused by runout

%% DOW ET AL. PARAMETERS
%  Note: Uncomment this section to reproduce Fig. 7 of Dow et al. paper

    % H    = 5.5e9;           % Hardness of workpiece material [Pa] [for S-7 steel as used by Dow]
    % E    = 210e9;           % Elastic modulus of workpiece [Pa]     
    % mu   = 0.46;            % Friction coefﬁcient between tool rake face and workpiece
    % mu_f = 0.2;             % Friction coefﬁcient between ﬂank face and workpiece
    % ga   = 50*pi/180;       % Shear angle of the chip [rad]
    % 
    % d    = 25e-6;           % Depth of cut [m]     
    % f    = 18.75 * 1e-6;    % Up-feed per flute [m/flute]
    % R_t  = 0.4e-3;          % Radius of tool [m]
    % W_l  = 7e-6;            % Width of wear land [m]
    % alp  = 0;               % Sweep angle is between the rotation axis of the tool and a line normal to the workpiece surface [rad]
    % zeta = 0;               % Offset of one flute w.r.t. the other caused by runout

 %% CALCULATIONS AND PLOTS

    phi   = acos(1-d/R_t);             % Contact angle between the tool and the workpiece (deﬁned by the tool radius and the depth of cut) [rad]
    theta = linspace(0,2*pi,360);   
    s1    = sin(theta);                % flute 1
    s2    = sin(theta + pi);           % flute 2 (phase-shifted by 180 degrees)

    % Gates: 1 when sinθ >= 0  else 0
    gate1    = double(s1 >= 0);
    gate2    = double(s2 >= 0);
    gate2(1) = 0;                                 % To address the unintended 1 at position 0
          
    A_c1  = d * f * s1 .* gate1 * (1 + zeta);     % Face area of chip
    A_c2  = d * f * s2 .* gate2 * (1 - zeta);     % Face area of chip
    A_f1  = R_t * W_l * phi ./ (2-s1) .* gate1;   % Flank area of tool (width times length of contact region on ﬂank face)
    A_f2  = R_t * W_l * phi ./ (2-s2) .* gate2;   % Flank area of tool (width times length of contact region on ﬂank face)
    
    Fc1   = H*A_c1/3*( (cot(ga)/sqrt(3))+1 ) + mu_f*A_f1*( 0.62*H*sqrt(43*H/E) );
    Ft1   = mu*H*A_c1/3*( (cot(ga)/sqrt(3))+1 ) + A_f1*( 0.62*H*sqrt(43*H/E) );

    Fc2   = H*A_c2/3*( (cot(ga)/sqrt(3))+1 ) + mu_f*A_f2*( 0.62*H*sqrt(43*H/E) );
    Ft2   = mu*H*A_c2/3*( (cot(ga)/sqrt(3))+1 ) + A_f2*( 0.62*H*sqrt(43*H/E) );
        
    Fx1 =  Fc1 * cos(alp) .* sin(theta) - Ft1 * sin(0.5*phi) .* cos(theta);
    Fy1 = -Fc1 * cos(alp) .* cos(theta) - Ft1 * sin(0.5*phi) .* sin(theta);
    Fz1 =  Fc1 * sin(alp) + Ft1 * cos(0.5*phi);

    Fx2 =  Fc2 * cos(alp) .* sin(theta+pi) - Ft2 * sin(0.5*phi) .* cos(theta+pi);
    Fy2 = -Fc2 * cos(alp) .* cos(theta+pi) - Ft2 * sin(0.5*phi) .* sin(theta+pi);
    Fz2 =  Fc2 * sin(alp) + Ft2 * cos(0.5*phi);

    Fx_ext = Fx1 + Fx2;
    Fy_ext = Fy1 + Fy2;
    Fz_ext = Fz1 + Fz2;

% PRINTING THE MEAN VALUES
    mean(Fx_ext)
    mean(Fy_ext)
    mean(Fz_ext)

% MILLING FORCE PLOTS
    figure;
    plot(theta,Fx_ext,'b','LineWidth',2); hold on;
    plot(theta,Fy_ext,'-.','LineWidth',2); hold on;
    plot(theta,Fz_ext,'--k','LineWidth',2); hold off;
    grid on;
    legend('F_x','F_y','F_z');
    ylabel('Force [N]');
    xlabel('Tool rotation angle [Rad]');
    % xlim([0 6.3]);
    % ylim([-0.4 1.1]);
    set(gca,'Box','off','FontSize',15,'LineWidth',3); 
    set(gcf,'Color','w');
