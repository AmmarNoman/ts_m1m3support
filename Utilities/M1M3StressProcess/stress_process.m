% INPUT before running "stress_process.m": stress, scale, coeff, elemID,
% index
% 
% Run "stress_process.m";
% 
% m file calculates algebraic sum of all plane stress values across
% combined mode shapes as defined by coeff, i.e. sumXtop sumYtop, sumXbot,
% sumYbot, sumXYtop, sumXYbot.
% 
% From these, the major/minor principle stress values are calculated for
% every element over the combined mode shapes, A_majorPrinTop,
% B_minorPrinTop, C_majorPrinBot,D_minorPrinBot.

% Next, the max principle stress values for any element in the mirror are
% determined, along with their element index number(row number).
% 
% 
% % stress: Stress values,top face fx,fy,fxy, bottom face fx,fy,fxy; for all
% % 27 % modes, for each element of the M1M3 glass. From FEA output for modes
% % based on 1 Newtom RMS.
% 
% % scale: scaling factor to convert stress from 1N RMS mode force to 1
% % micron RMS shape.
% 
% % factor: matrix of all scale factors(including multiplying by coeff)
% % used to % piecewise (.*) multiply times stress to get scaled_stress.
% 
% % scaled_stress: stress values scaled by factor to be consistent with 1
% % micron RMS mode shape and include scaling by coeff.
% 
% % stressXtop, stressYtop, stressXYtop, stressXbot, stressYbot, stressXYbot:
% % scaled stress values for each mode by stress type.
% 
% % sumXtop, sumYtop, sumXYtop, sumXbot, sumYbot, sumXYbot: algebraic sum of
% % element stress across all 27 modes.
% 
% % A_majorPrinTop, B_minorPrinTop, C_majorPrinBot, D_minorPrinBot; principle
% % major and minor stress values calculated from the plane normal stresses
% % x,y, and shear xy for each element, all modes combined.
% 
% % Amax, Bmin, Cmax, Dmin: Max tensile stress for +mode shape(Amax,Cmax),
% % Max tensile stress for -mode shape(Bmin,Dmin ie. minor principle stress
% % min).
% 
% % IAmax, IBmin, ICmax, IDmin; indicies (row address of max/min stress
% % values)
% 
% % interval: number of stress values in stress table for each element of
% % each mode shape. Always 6.
% 
% % coeff: coefficient of each mode shape used to describe the total force
% % set under investigation.


% Columns of the stress matrix represent the 6 plane stresses of each
% bending mode for a total of 27*6=162 columns. Each row of the stress
% matrix represents an element of the mirror model. Below, we set up vector
% A, which is used as a way to index from stresses for one bending mode to
% the next as we generate the segregated individual stress matricies
% stressXtop, stressXbot, stressYtop, ect.

interval=6;
A=1:interval:27*interval;

% Use the repeat elements command to generate the factor matrix. The stress
% matrix values (directly from the FEA output) are based on bending modes
% loads that were normalized to 1N RMS. The bending mode coefficients we
% use to define wavefront corrections in the control software are defined
% on the basis of loads consistent with 1 micron RMS displacement of the
% optical surface. A scale factor is used to transition between forces
% normalized to 1 N RMS and displacements normalized to 1 micron RMS. Each
% mode has a unique scale factor, hence the scale vector has 27 scale
% values, one for each bending mode. Also, we need to scale the stresses by
% the desired bending mode coefficients. The factor matrix includes both
% the scaling factors and the bending mode coeff product for all 27 modes,
% repeated row wise for all 41254 elements. 

% When the stress matrix is piecewise multiplied by the factor matrix, we
% obtain the scaled_stress matrix, wherein each stress value is multiplied
% by the appropriate coeff*scale factor to obtain stress values
% corresponding to the input bending mode coefficients. 

factor=repelem(coeff.*scale,41254,[6]);
scaled_stress=factor.*stress;

% Generate segregated 2D normal stress values for sigma X, sigma Y, and sigma XY on
% both the plate top surface and bottom surface.

stressXtop=scaled_stress(:,A);
stressYtop=scaled_stress(:,A+1);
stressXYtop=scaled_stress(:,A+2);
stressXbot=scaled_stress(:,A+3);
stressYbot=scaled_stress(:,A+4);
stressXYbot=scaled_stress(:,A+5);

% Sum together the normal stresses across all bending modes for every
% element. Each row corresponds to an element.

sumXtop=sum(stressXtop,2);  
sumYtop=sum(stressYtop,2);
sumXYtop=sum(stressXYtop,2);
sumXbot=sum(stressXbot,2);
sumYbot=sum(stressYbot,2);
sumXYbot=sum(stressXYbot,2);

% Calculate the major principal stress for each element corresponding to
% the applied bending mode coefficients. The major principal stress
% represents the maximum glass tensile stress, regardless of stress
% direction. This is calculated from the plane normal and shear stress
% values for each element, both on the plate upper surface (top) and lower
% surface (bot). They are not the same due to bending of the plates.
% Maximum values may come from either the top or bottom surface, depending
% on local element bending strains, so all must be checked. 
% The minor principal stress values represent the maximum compressive
% stress in the glass. These are not needed because the glass is far
% stronger in compression than in tension, and are therefore not used.
% 
% We sort thru each element principal stress values to find the maximum
% value (Amax, Cmax) and ID of the array row where it occurs (IAmax,
% ICmax). Find the largest major principal stress between plate top (Amax)
% or plate bot (Cmax), and assign that value to Max_Tensile_Stress.

A_majorPrinTop=0.5*((sumXtop+sumYtop)+sqrt((sumXtop-sumYtop).^2+(4*sumXYtop.^2)));
% B_minorPrinTop=0.5*((sumXtop+sumYtop)-sqrt((sumXtop-sumYtop).^2+(4*sumXYtop.^2)));
C_majorPrinBot=0.5*((sumXbot+sumYbot)+sqrt((sumXbot-sumYbot).^2+(4*sumXYbot.^2)));
% D_minorPrinBot=0.5*((sumXbot+sumYbot)-sqrt((sumXbot-sumYbot).^2+(4*sumXYbot.^2)));
[Amax,IAmax]=max(A_majorPrinTop);
% [Bmin,IBmin]=min(B_minorPrinTop);
[Cmax,ICmax]=max(C_majorPrinBot);
% [Dmin,IDmin]=min(D_minorPrinBot);
Max_Tensile_Stress=max(Amax,Cmax);

% Find the element number in the FEA corresponding to IAmax and ICmax.
% These are useful for comparison with FEA results.

EIDa=elemID(IAmax);
EIDc=elemID(ICmax);