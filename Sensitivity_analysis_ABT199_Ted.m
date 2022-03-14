% Script to run local sensitivity analysis
function Sensitivity_analysis_ABT199_Ted
close all;    
format long;
t_scale = 12.0*60.0;   %12 hours

%Coupled ISR-Bcl2 pathway parameters: 
myPars = [13.669755696602916, 485.697261747894288, 1044.141691325814236,...
          1153.185096086854855, 7.192553700040385, 82.687472599671310,...
          452.154079358190415, 0.161810685573990, 3808.779122436302714,...
          129.902673437969696, 2.320021253551643, 11.452565945560234,...
          84.164288956134570, 1.729787988695099, 3.151231128658153,...
          0.602022914394669, 3.520635994758134, 1.440200545088872,...
          20.192519980890044, 14.705005047381947, 25.637206631177463,...
          1.950200147859286, 0.735094861939094, 0.846708465943897,...
          1.726707485864304, 0.530205713378790, 1.638885619999066,...
          1.858928933236780, 3.438405173378365, 8.080797703919686,...
          10.536611750584026, 3.706640794850193, 2.721147803974464,...
          3.922800957032131, 2.377366569413438, 4.620362679300571,...
          6.817649060433075, 2.060800458146430, 10.993215713647745,...
          0.521150510527861, 2.002345966157941, 2.727210293394596,...
          2.271195899462174, 6.699538751856998, 3.484730980617899,...
          3.721176559333263, 1.025447266055947, 1.702891500541281,...
          0.140738806893586, 1.019986323110426, 0.014757110602751];

%Transcription factor parameters:
TFPars = [1.159914, 0.027393, 0.972924, 0.106661, 11.425474,...
          2.146729, 3.284086, 4.849403, 2.669784, 3.7034320,...
          1.682787, 1.568424, 6.943936, 4.620715, 0.0001039,...
          2.423171, 208.89961, 0.863131, 2.285103, 0.499996,...
          0.418970, 382.49187, 1.288237, 20.22993, 5.736299,...
          0.920794, 1.5560653, 10.74388, 8.771861, 0.634600,...
          2.107549, 113.76057];     

all_params = [myPars, TFPars];

%Names of the parameter values (I have taken the order here and translated 
%according to the names in Table S6 of the paper, which is why they do not
%seem to line up with the ga file or nomenclature in simulations):
Parameter_names = {'k_{9}', 'k_{11}',...
                   '[Mcl-1]_{0}', '[Bcl-2]_{0}', '[Bim]_{0}',...
                   '[Bax]_{0}', '[Bak]_{0}', '[Casp-3]_{0}',...
                   'k_{0,Mcl-1}', 'k_{0,Bcl-2}', 'k_{0,Bim}',...
                   'k_{0,Bax}', 'k_{0,Bak}',...
                   'K_{1}', 'K_{4}', 'K_{5}', 'K_{6}', 'K_{7}',...
                   'K_{9}', 'K_{10}', 'K_{11}', 'K_{12}',...
                   'K_{13}',...
                   '\lambda_{Bcl-2}', '\lambda_{c-Myc}', '\lambda_{0}',...
                   '\Omega_{Casp-3}',...
                   'K_{14}',...
                   'n_{1}', 'n_{4}', 'n_{9}', 'n_{10}', 'n_{5}',...
                   'n_{11}', 'n_{6}', 'n_{12}', 'n_{7}', 'n_{13}',...
                   'K_{15}', 'K_{17}',...
                   'K_{8}', 'n_{8}',...
                   'n_{17}', 'n_{14}', 'n_{15}',...
                   'K_{16}', 'n_{16}',...
                   'n_{3}', 'K_{3}', 'n_{2}', 'K_{2}',... 
                   '\delta_{T_{i,1}}', '\delta_{T_{E,1}}',... 
                   '\delta_{T_{i,2}}', '\delta_{T_{E,2}}',... 
                   'k_{T_{1,1}}',...
                   'K_{1,1,1}', 'n_{1,1,1}', 'K_{1,2,1}', 'n_{1,2,1}',...
                   'K_{1,1,4}', 'n_{1,1,4}', 'K_{1,2,4}', 'n_{1,2,4}',...
                   '\delta_{T_{1,0}}',...
                   'w_{1,1,1}', 'w_{1,2,1}',...
                   'q_{1,4}',...
                   '\tau_{1,1,1}', '\tau_{1,2,1}',...
                   '\tau_{2,2,4}', '\tau_{2,2,1}',...
                   'k_{T_{2,1}}',...
                   'K_{2,1,1}', 'n_{2,1,1}', 'K_{2,1,4}', 'n_{2,1,4}',...
                   'K_{2,2,1}', 'n_{2,2,1}', 'K_{2,2,4}', 'n_{2,2,4}',...
                   '\delta_{T_{2,0}}',...
                   'q_{2,1}'};

rng default
%==========================================================================
% (1) Read in the experimental data
%==========================================================================
fname_protein = 'Protein-April-2019.dat';         % Protein data file.

data_format = '%f';                             % data file formats (or e.g. %16f)
fid_protein = fopen(fname_protein, 'r');

m = 1;                     % loop counter variable
num_proteins = 8;            % number of proteins in the Western Blot data file.
num_protocols_protein = 6;   % x 2 because there is an average and standard deviation.

protein_expt = zeros(num_proteins,num_protocols_protein);

while true                % infinite loop until end of all files is reached
   data_protein = textscan(fid_protein, data_format, num_protocols_protein); % (num_protocols-1) data columns 
   if isempty(data_protein{1}(:)) 
       break
   end
      
   if ~isempty(data_protein{1}(:)) 
      for i=1:num_protocols_protein
        protein_expt(m,i) = data_protein{1}(i);            % Protein file goes as: ABT-199, Tedizolid, Combination.
      end
   end                 
   m = m + 1;
end
fclose(fid_protein);


%==========================================================================
% (2) Run sensitivity analysis:
%==========================================================================
%Do a single reference run:
[fval, via_orig] = obj_ABT199_Ted_sensitivity(myPars,TFPars,protein_expt,t_scale);
%disp(fval)

% Now perturb the parameters and do the sensitivity analysis:
num_params = max(size(all_params));
Sensitivity = ones(4,num_params,2);
Perturbation = 0.01;

% Monotherapies:
y_un = zeros(num_params,2);
y_abt = zeros(num_params,2);
y_ted = zeros(num_params,2);
y_at = zeros(num_params,2);
all_params_nom = all_params;

for j=1:num_params
    all_params = all_params_nom;     %Reset the parameter array
    all_params(j) = (1+Perturbation)*all_params_nom(j);   
    myPars = all_params(1:51);
    TFPars = all_params(52:end);
    [fval, via_pert] = obj_ABT199_Ted_sensitivity(myPars,TFPars,protein_expt,t_scale);
    if fval >= 1.0e8
       disp('fval too large, parameter number: ')
       disp(j)
    end
    
    for k=1:4
        Sensitivity(k,j,1) = ((via_pert(k) - via_orig(k))/via_orig(k))/Perturbation;          
    end    

    %Now do the -Perturbation:      
    all_params(j) = (1-Perturbation)*all_params_nom(j);
    myPars = all_params(1:51);
    TFPars = all_params(52:end);
    [fval, via_pert] = obj_ABT199_Ted_sensitivity(myPars,TFPars,protein_expt,t_scale);
    if fval >= 1.0e8
       disp('fval too large, parameter number: ')
       disp(j)
    end

    for k=1:4    
        Sensitivity(k,j,2) = ((via_pert(k) - via_orig(k))/via_orig(k))/Perturbation;          
    end    
    
     y_un(j,:) = [Sensitivity(1,j,1), Sensitivity(1,j,2)];
     y_abt(j,:) = [Sensitivity(2,j,1), Sensitivity(2,j,2)];
     y_ted(j,:) = [Sensitivity(3,j,1), Sensitivity(3,j,2)];
     y_at(j,:) = [Sensitivity(4,j,1), Sensitivity(4,j,2)];
end

% Make plots for each parameter sensitivity:
figure;
y_un_plus = y_un(:,1);
[~, sort_id] = sort(abs(y_un_plus),'descend');
Param_names = Parameter_names(sort_id);
y_un_new = y_un(sort_id,:);
y_un_plot = y_un_new(1:25,:);
b = bar(y_un_plot, 'Facecolor', 'flat');
%set (gca, 'XTickLabel', Parameter_names(1:end))
set(gca, 'XTick', 1:1:25, 'XTickLabel', Param_names(1:25))
set(gca,'LineWidth',2,'TickLength',[0.015 0.015]);
xtickangle(45)
title('Untreated')
for k=1:size(y_un_plot,2)
    b(k).CData = k;
end

figure;
y_abt_plus = y_abt(:,1);
[~, sort_id] = sort(abs(y_abt_plus),'descend');
Param_names = Parameter_names(sort_id);
y_abt_new = y_abt(sort_id,:);
y_abt_plot = y_abt_new(1:25,:);
b = bar(y_abt_plot, 'Facecolor', 'flat');
%set (gca, 'XTickLabel', Parameter_names(1:end))
set(gca, 'XTick', 1:1:25, 'XTickLabel', Param_names(1:25))
set(gca,'LineWidth',2,'TickLength',[0.015 0.015]);
xtickangle(45)
title('ABT-199 monotherapy')
for k=1:size(y_abt_plot,2)
    b(k).CData = k;
end

figure;
y_ted_plus = y_ted(:,1);
[~, sort_id] = sort(abs(y_ted_plus),'descend');
Param_names = Parameter_names(sort_id);
y_ted_new = y_ted(sort_id,:);
y_ted_plot = y_ted_new(1:25,:);
b = bar(y_ted_plot, 'Facecolor', 'flat');
%set (gca, 'XTickLabel', Parameter_names(1:end))
set(gca, 'XTick', 1:1:25, 'XTickLabel', Param_names(1:25))
set(gca,'LineWidth',2,'TickLength',[0.015 0.015]);
xtickangle(45)
title('Tedizolid monotherapy')
for k=1:size(y_ted_plot,2)
    b(k).CData = k;
end    

figure;
y_at_plus = y_at(:,1);
[~, sort_id] = sort(abs(y_at_plus),'descend');
Param_names = Parameter_names(sort_id);
y_at_new = y_at(sort_id,:);
y_at_plot = y_at_new(1:25,:);
b = bar(y_at_plot, 'Facecolor', 'flat');
%set (gca, 'XTickLabel', Parameter_names(1:end))
set(gca, 'XTick', 1:1:25, 'XTickLabel', Param_names(1:25))
set(gca,'LineWidth',2,'TickLength',[0.015 0.015]);
xtickangle(45)
title('Combination monotherapy')
for k=1:size(y_at_plot,2)
    b(k).CData = k;
end    

end