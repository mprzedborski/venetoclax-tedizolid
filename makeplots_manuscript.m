%==========================================================================
% Plotting function. 
%==========================================================================
% This file creates the plots of the protein levels in the coupled 
% ISR-apoptosis pathway with ABT-199 and Tedizolid treatment. 

function makeplots_manuscript(t_equil,protein_equil,t_treatment,protein_treatment,protein_expt,P_0,sw,tr) 
ABT_0=P_0(1); mcl1_0=P_0(2); bcl2_0=P_0(3); bim_0=P_0(4);
bax_0=P_0(5); bak_0=P_0(6);

%==========================================================================
% Read in the equilibrium protein data: 
%==========================================================================
%mcl=x(1); bcl2=x(2); bim=x(3); bax=x(4); bak=x(5); 
%abax=x(6); abak=x(7); casp3=x(8); acasp3=x(9); 
%mcli=x(10); 
%mclbim=x(11); bcl2bim=x(12); 
%bimbax=x(13); bimbak=x(14);        
%mclabax=x(15); mclabak=x(16); bcl2abax=x(17); bcl2abak=x(18); 
%bcl2i=x(19); 
%ABT-Bcl-2 complex=x(24);

% %--------------------------------------------------------------------------
% % (1) Free Mcl-1, Bcl-2:
% %-----------------------------------
% mcl1_equil = protein_equil(:,1);
% mcl1i_equil = protein_equil(:,10);
% 
% bcl2_equil = protein_equil(:,2);
% bcl2i_equil = protein_equil(:,19);
% 
% % (2) Free Bim, Bak and Bax:
% %---------------------------------
% bim_equil = protein_equil(:,3);
% 
% bax_equil = protein_equil(:,4);
% abax_equil = protein_equil(:,6);
% 
% bak_equil = protein_equil(:,5);
% abak_equil = protein_equil(:,7);
% 
% % (3) Casp-3 and aCasp-3:
% %------------------------
% casp3_equil = protein_equil(:,8); 
% acasp3_equil = protein_equil(:,9); 
% 
% % (4) Mcl-1, Bcl-2 with Bim:
% %-----------------------------------------------
% mclbim_equil = protein_equil(:,11);
% 
% bcl2bim_equil = protein_equil(:,12);
% 
% % (5) Bak and Bax with Bim and Bid:
% %----------------------------------
% bimbax_equil = protein_equil(:,13);
% 
% bimbak_equil = protein_equil(:,14);
% 
% % (6) Bak and Bax with Bcl-2, Mcl-1:
% %-----------------------------------------------      
% mcl1abax_equil = protein_equil(:,15);
% mcl1abak_equil = protein_equil(:,16);
% 
% bcl2abax_equil = protein_equil(:,17);
% bcl2abak_equil = protein_equil(:,18);
% 
% 
% % (8) Transcription factors Myc and Chop:
% %----------------------------------------
% myc_equil = protein_treatment(1,20);      % Myc and Chop levels are constant during
% chop_equil = protein_treatment(1,21);      % equilibration simulations.

% Total untreated protein levels:
%----------------------------------
Total_protein_equil = zeros(8,max(size(protein_equil)));
max_ind = max(size(t_equil));
Calculated_protein_SS = protein_equil(max_ind,:);


%Mcl-1: 
Total_protein_equil(1) = Calculated_protein_SS(1) + Calculated_protein_SS(10)...
                       + Calculated_protein_SS(11)...
                       + Calculated_protein_SS(15) + Calculated_protein_SS(16);   
%Bcl-2:
Total_protein_equil(2) = Calculated_protein_SS(2) + Calculated_protein_SS(19)...
                       + Calculated_protein_SS(12)...
                       + Calculated_protein_SS(17) + Calculated_protein_SS(18);
%Bim:
Total_protein_equil(3) = Calculated_protein_SS(3)...
                       + Calculated_protein_SS(11)*mcl1_0/bim_0 + Calculated_protein_SS(12)*bcl2_0/bim_0...
                       + Calculated_protein_SS(13)*bax_0/bim_0 + Calculated_protein_SS(14)*bak_0/bim_0;
%Bax:
Total_protein_equil(4) = Calculated_protein_SS(4) + Calculated_protein_SS(6)...
                       + Calculated_protein_SS(13)...
                       + Calculated_protein_SS(15)*mcl1_0/bax_0 + Calculated_protein_SS(17)*bcl2_0/bax_0;
%Bak:
Total_protein_equil(5) = Calculated_protein_SS(5) + Calculated_protein_SS(7)...
                       + Calculated_protein_SS(14)...
                       + Calculated_protein_SS(16)*mcl1_0/bak_0 + Calculated_protein_SS(18)*bcl2_0/bak_0;
%aCasp-3:                
Total_protein_equil(6) = Calculated_protein_SS(9);    
   
                
%==========================================================================
% Read in the ABT-199 treatment protein data: 
%==========================================================================
%mcl=x(1); bcl2=x(2); bclx=x(3); bid=x(4); bim=x(5); bax=x(6);
%bak=x(7); abax=x(8); abak=x(9); casp3=x(10); acasp3=x(11); 
%mcli=x(12); mclbim=x(13); mclbid=x(14);
%bcl2bim=x(15); bcl2bid=x(16); bclxbim=x(17); bclxbid=x(18);
%bidbax=x(19); bimbax=x(20); bidbak=x(21); bimbak=x(22);        
%mclabax=x(23); mclabak=x(24); bcl2abax=x(25); bcl2abak=x(26); 
%bclxabax=x(27); bclxabak=x(28); bcl2i=x(29); bclxi=x(30);
%free ABT=x(31); ABT:bcl-2 complex=x(32);
%cmyc=x(33); chop=x(34);

% %--------------------------------------------------------------------------
% % (1) Free Mcl-1, Bcl-2 and Bcl-XL:
% %-----------------------------------
% mcl1_treatment = protein_treatment(:,1);
% mcl1i_treatment = protein_treatment(:,10);
% 
% bcl2_treatment = protein_treatment(:,2);
% bcl2i_treatment = protein_treatment(:,19);
% 
% % (2) Free Bim, Bak and Bax:
% %---------------------------------
% bim_treatment = protein_treatment(:,3);
% 
% bax_treatment = protein_treatment(:,4);
% abax_treatment = protein_treatment(:,6);
% 
% bak_treatment = protein_treatment(:,5);
% abak_treatment = protein_treatment(:,7);
% 
% % (3) Casp-3 and aCasp-3:
% %------------------------
% casp3_treatment = protein_treatment(:,8); 
% acasp3_treatment = protein_treatment(:,9); 
% 
% % (4) Mcl-1, Bcl-2 with Bim:
% %-----------------------------------------------
% mcl1bim_treatment = protein_treatment(:,11);
% 
% bcl2bim_treatment = protein_treatment(:,12);
% 
% 
% % (5) Bak and Bax with Bim:
% %----------------------------------
% bimbax_treatment = protein_treatment(:,13);
% 
% bimbak_treatment = protein_treatment(:,14);
% 
% % (6) Bak and Bax with Bcl-2, Mcl-1:
% %-----------------------------------------------      
% mcl1abax_treatment = protein_treatment(:,15);
% mcl1abak_treatment = protein_treatment(:,16);
% 
% bcl2abax_treatment = protein_treatment(:,17);
% bcl2abak_treatment = protein_treatment(:,18);
% 
% 
% % (8) Transcription factors Myc and Chop:
% %----------------------------------------
% myc_treatment = protein_treatment(:,20);      % Myc and Chop levels are constant during
% chop_treatment = protein_treatment(:,21);      % equilibration simulations.

% Total untreated protein levels:
%----------------------------------
Total_protein_treatment = zeros(8,max(size(protein_treatment)));


if sw == 1
   Bcl2_complex = protein_treatment(:,24);
else
   Bcl2_complex = 0;
end


%Mcl-1:
Total_protein_treatment(1,:) = protein_treatment(:,1) + protein_treatment(:,10)...
                             + protein_treatment(:,11)...
                             + protein_treatment(:,15) + protein_treatment(:,16);   

%Bcl-2:
Total_protein_treatment(2,:) = protein_treatment(:,2) + protein_treatment(:,19)...
                             + protein_treatment(:,12)...
                             + protein_treatment(:,17) + protein_treatment(:,18)...
                             + Bcl2_complex*ABT_0/bcl2_0;                 % Bcl-2 complex levels 
%Bim:
Total_protein_treatment(3,:) = protein_treatment(:,3)...
                             + protein_treatment(:,11)*mcl1_0/bim_0...
                             + protein_treatment(:,12)*bcl2_0/bim_0......
                             + protein_treatment(:,13)*bax_0/bim_0... 
                             + protein_treatment(:,14)*bak_0/bim_0;

%Bax:
Total_protein_treatment(4,:) = protein_treatment(:,4) + protein_treatment(:,6)...
                             + protein_treatment(:,13)...
                             + protein_treatment(:,15)*mcl1_0/bax_0...
                             + protein_treatment(:,17)*bcl2_0/bax_0;

%Bak:
Total_protein_treatment(5,:) = protein_treatment(:,5) + protein_treatment(:,7)...
                             + protein_treatment(:,14)...
                             + protein_treatment(:,16)*mcl1_0/bak_0... 
                             + protein_treatment(:,18)*bcl2_0/bak_0;

%aCasp-3:                
Total_protein_treatment(6,:) = protein_treatment(:,9);              

%c-Myc:
Total_protein_treatment(7,:) = protein_treatment(:,20);
   
%Chop:
Total_protein_treatment(8,:) = protein_treatment(:,21);


%Total intracellular ABT levels:
%Total_ABT = protein_treatment(:,23) + protein_treatment(:,24);


%==========================================================================
% Now create the plots to illustrate protein dynamics, and save them to the 
% working directory specified above:
%==========================================================================   
% Screensize for plots
pos_fig = [0 0 1920 1080];

%--------------------------------------------------------------------------
% Plots of simulation results for time evolution of important protein 
% concentrations:
%--------------------------------------------------------------------------
%xmax_equil = max(t_equil);
%if xmax_equil > 480
%   t_equil = t_equil/60.0;
%   timescale_equil = 'hours';
%   xmax_equil = max(t_equil);
%else 
%   timescale_equil = 'minutes';
%end


%--------------------------------------------------------------------------
% Treatment plots:
%--------------------------------------------------------------------------
t_treatment = t_treatment/60.0;
%xmax_treatment = max(t_treatment)+1;
%timescale_treatment = 'hours';

%label_str_treatment = ["Time (",timescale_treatment,")"];
%xlab_treatment = join(label_str_treatment); 

t_expt = 3.0*24;
%Find the time index that corresponds to t_expt:
[dist, ind] = min(abs(t_treatment-t_expt));

if tr == 1
   Mcl1_expt = protein_expt(1,1); 
   Bcl2_expt = protein_expt(2,1);
   Bim_expt = protein_expt(3,1);
   Bax_expt = protein_expt(4,1);
   Bak_expt = protein_expt(5,1);
   cMyc_expt = protein_expt(7,1);
   Chop_expt = protein_expt(8,1);

   err_Mcl1 = protein_expt(1,2);
   err_Bcl2 = protein_expt(2,2);
   err_Bim = protein_expt(3,2);
   err_Bax = protein_expt(4,2);
   err_Bak = protein_expt(5,2);
   err_Casp = 0.5; 
   err_cMyc = protein_expt(7,2);   
   err_Chop = protein_expt(8,2);
elseif tr == 2
   Mcl1_expt = protein_expt(1,3); 
   Bcl2_expt = protein_expt(2,3);
   Bim_expt = protein_expt(3,3);
   Bax_expt = protein_expt(4,3);
   Bak_expt = protein_expt(5,3);
   cMyc_expt = protein_expt(7,3);
   Chop_expt = protein_expt(8,3);

   err_Mcl1 = protein_expt(1,4);
   err_Bcl2 = protein_expt(2,4);
   err_Bim = protein_expt(3,4);
   err_Bax = protein_expt(4,4);
   err_Bak = protein_expt(5,4);
   err_Casp = 0.5;
   err_cMyc = protein_expt(7,4);   
   err_Chop = protein_expt(8,4);
elseif tr == 3
   Mcl1_expt = protein_expt(1,5); 
   Bcl2_expt = protein_expt(2,5);
   Bim_expt = protein_expt(3,5);
   Bax_expt = protein_expt(4,5);
   Bak_expt = protein_expt(5,5);
   cMyc_expt = protein_expt(7,5);
   Chop_expt = protein_expt(8,5);

   err_Mcl1 = protein_expt(1,6);
   err_Bcl2 = protein_expt(2,6);
   err_Bim = protein_expt(3,6);
   err_Bax = protein_expt(4,6);
   err_Bak = protein_expt(5,6);
   err_Casp = protein_expt(6,6);
   err_cMyc = protein_expt(7,6);   
   err_Chop = protein_expt(8,6);
end

if (sw==1)
   if tr == 1
      aCasp_expt = 1.5;
    else
      aCasp_expt = protein_expt(6,5);
   end       
else
   aCasp_expt = 1.25;
end   

Expt_proteins = [Mcl1_expt, Bcl2_expt, Bim_expt, Bax_expt, Bak_expt, aCasp_expt, cMyc_expt, Chop_expt];
Expt_err = [err_Mcl1, err_Bcl2, err_Bim, err_Bax, err_Bak, err_Casp, err_cMyc, err_Chop];
Sim_proteins = [Total_protein_treatment(1,ind)/Total_protein_treatment(1,1),...
                Total_protein_treatment(2,ind)/Total_protein_treatment(2,1),...
                Total_protein_treatment(3,ind)/Total_protein_treatment(3,1),...
                Total_protein_treatment(4,ind)/Total_protein_treatment(4,1),...
                Total_protein_treatment(5,ind)/Total_protein_treatment(5,1),...
                Total_protein_treatment(6,ind)/Total_protein_treatment(6,1),...
                Total_protein_treatment(7,ind)/Total_protein_treatment(7,1),...
                Total_protein_treatment(8,ind)/Total_protein_treatment(8,1)];
Protein_names = {'Mcl-1', 'Bcl-2', 'Bim', 'Bax', 'Bak', 'Casp-3*', 'c-Myc', 'Chop'};

figure('Position',pos_fig);
hold on;
errorbar(Expt_proteins, Expt_err,'bo', 'LineWidth',3, 'MarkerSize',10,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor','b');
plot(Sim_proteins, 'rs','MarkerSize',15,'MarkerEdgeColor','r','MarkerFaceColor','r');
box on;
set(gca,'fontsize',36);
set(gca,'LineWidth',2,'TickLength',[0.015 0.015]);
set(gca,'YMinorTick','on')
set(gca, 'XTick', 1:1:8, 'XTickLabel', Protein_names)
xlim([0.5 8.5]);
if tr == 3
    ylim([0 10]);
else
    ylim([0 2]);
end
ylabel('Protein level normalized to untreated control');
xlabel('Protein');
hold off;