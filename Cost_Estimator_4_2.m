clear

%% Define Model Characteristics
n_iter = 50000; % Define the number of iterations in the model [size: n]
Neg_Eff_Max = 2500 ;% Define the upper boundary for the negative effect for the optimisation
%% Import Data
AllData=readcell('MitCon_Test_Case.xlsx');

%% Extract the Activity data from datafile
ActivityNames_m = rmmissing(string(AllData(2:end,1))); % All names of activities
m_Activities = size(ActivityNames_m,1); %Define the number of activities in the project [Size: m]

Act_PERT_Labour_m3 =rmmissing(cell2mat(AllData(2:end,2:4)));
Act_PERT_LabourCost_m3 = rmmissing(cell2mat(AllData(2:end,5:7)));
Act_PERT_Quantity_m3 = rmmissing(cell2mat(AllData(2:end,8:10)));
Act_PERT_MatRate_m3 = rmmissing(cell2mat(AllData(2:end,11:13)));
Act_PERT_EqRate_m3 = rmmissing(cell2mat(AllData(2:end,14:16)));
Act_PERT_SubRate_m3 = rmmissing(cell2mat(AllData(2:end,17:19)));

Act_TimeStep_m = [AllData{2:end,20}].';% Timestep in which an activity takes place
%% Extract Risk data from Datafile

RiskNames_o = rmmissing(string(AllData(2:end,21))); % All names of the risks
o_Risk = size(RiskNames_o,1); % Define the number of risks in the project [Size: o]

Risk_Prob_o = rmmissing([AllData{2:end,22}]); % Probabilities for all risk events

Risk_PERT_Effect_o3 = rmmissing(cell2mat((AllData(2:o_Risk+1,23:25))));

%% Extract Chances data from Datafile

Chances_Names_q = rmmissing(string(AllData(2:end,26))); % All names of the chances
q_Chance = size(Chances_Names_q,1);

Chance_Prob_q=rmmissing([AllData{2:end,27}]); % Probabilities for all chance events
Chance_PERT_Effect_o3 = rmmissing(cell2mat((AllData(2:q_Chance+1,28:30))));
%% Extract Mitigation Data from datafile

MitNames_p = rmmissing(string(AllData(2:end,31))); % All Mitigation names
p_Mitigations = size(MitNames_p,1); % Defines the number of Mitigations in the project [Size: p]

Mit_PERT_Labour_Con_p3 = cell2mat((AllData(2:p_Mitigations+1,32:34)));
Mit_Labour_Var_p = [AllData{2:(p_Mitigations+1),35}].'; %Variable element of mitigation on labour aspect

Mit_PERT_LabourCost_Con_p3 = cell2mat((AllData(2:p_Mitigations+1,36:38)));
Mit_LabourCost_Var_p = [AllData{2:(p_Mitigations+1),39}].'; % Effects of all mitigations on LabourCost

Mit_PERT_Quantity_Con_p3 = cell2mat((AllData(2:p_Mitigations+1,40:42)));
Mit_Quantity_Var_p = [AllData{2:(p_Mitigations+1),43}].'; % Effects of all mitigations on Material aspect

Mit_PERT_MatRate_Con_p3 = cell2mat((AllData(2:p_Mitigations+1,44:46)));
Mit_MatRate_Var_p = [AllData{2:(p_Mitigations+1),47}].'; % Effects of all mitigations on UnitCost

Mit_PERT_EqRate_Con_p3 = cell2mat(AllData(2:p_Mitigations+1,48:50)); %Negative effect of all mitigations
Mit_EqRate_Var_p = [AllData{2:(p_Mitigations+1),51}].'; % Timestep in which this mitigation has effect

Mit_PERT_SubRate_Con_p3 = cell2mat(AllData(2:p_Mitigations+1,52:54)); %Negative effect of all mitigations
Mit_SubRate_Var_p = [AllData{2:(p_Mitigations+1),55}].'; % Timestep in which this mitigation has effect

Mit_Neg_Effect_p = cell2mat(AllData(2:p_Mitigations+1,56)); %Negative effect of all mitigations
Mit_Neg_Effect_p_inv = Mit_Neg_Effect_p.';
Mit_Effected_Activity_p = [AllData{2:(p_Mitigations+1),57}].'; % Timestep in which this mitigation has effect

%% Drawing Values for Activities
Act_LabourDraw_nm = zeros(n_iter,m_Activities); %Pre-Allocate matrix for drawn labour values
Act_LabourCostDraw_nm= zeros(n_iter,m_Activities); %Pre-Allocate matrix for drawn labour cost values
Act_Quantity_Draw_nm= zeros(n_iter,m_Activities); %Pre-Allocate matrix for drawn Quantity values
Act_MaterialRate_Draw_nm= zeros(n_iter,m_Activities); %Pre-Allocate matrix for drawn MaterialRate values
Act_EquipmentRate_Draw_nm= zeros(n_iter,m_Activities); %Pre-Allocate matrix for drawn EquipmentRate values
Act_SubcontractRate_Draw_nm= zeros(n_iter,m_Activities); %Pre-Allocate matrix for drawn SubcontractorRate values
Act_CostDraw_nm = zeros(n_iter,m_Activities); %Pre-Allocate matrix for drawn cost item values

for act = 1:m_Activities
    if Act_PERT_Labour_m3(act,2) == Act_PERT_Labour_m3(act,1) && Act_PERT_Labour_m3(act,2)>0
        [labour] = Act_PERT_Labour_m3(act,2);
    elseif Act_PERT_Labour_m3(act,2)>0
        [labour]=RandPert (Act_PERT_Labour_m3(act,1),Act_PERT_Labour_m3(act,2),Act_PERT_Labour_m3(act,3),n_iter); %Draw value [w] from distribution based on Min/ML/Max for Labour
    else
        [labour] = 0;
    end
    Act_LabourDraw_nm(:,act) = labour;
    
    if Act_PERT_LabourCost_m3(act,2) == Act_PERT_LabourCost_m3(act,1) && Act_PERT_LabourCost_m3(act,2)>0
        [labourcost] = Act_PERT_LabourCost_m3(act,2);
    elseif Act_PERT_LabourCost_m3(act,2)>0
        Act_LabourCostDraw_nm(:,act)=RandPert (Act_PERT_LabourCost_m3(act,1),Act_PERT_LabourCost_m3(act,2),Act_PERT_LabourCost_m3(act,3),n_iter); %Draw value [x] from distribution based on Min/ML/Max for Labour Cost
    else
        [labourcost] = 0;
    end
    Act_LabourCostDraw_nm(:,act) = labourcost;
    
    if Act_PERT_Quantity_m3(act,2) == Act_PERT_Quantity_m3(act,1) && Act_PERT_Quantity_m3(act,2)>0
        [quantity] = Act_PERT_Quantity_m3(act,2);
    elseif Act_PERT_Quantity_m3(act,2)>0
        [quantity]=RandPert (Act_PERT_Quantity_m3(act,1),Act_PERT_Quantity_m3(act,2),Act_PERT_Quantity_m3(act,3),n_iter); %Draw value [y] from distribution based on Min/ML/Max for Quantity
    else
        [quantity] = 0;
    end
    Act_Quantity_Draw_nm(:,act) = quantity;
    
    if Act_PERT_MatRate_m3(act,2) == Act_PERT_MatRate_m3(act,1) && Act_PERT_MatRate_m3(act,2)>0
        [material] = Act_PERT_MatRate_m3(act,2);
    elseif Act_PERT_MatRate_m3(act,2)>0
        [material]=RandPert (Act_PERT_MatRate_m3(act,1),Act_PERT_MatRate_m3(act,2),Act_PERT_MatRate_m3(act,3),n_iter); %Draw value [z] from distribution based on Min/ML/Max for EquipmetRate
    else
        [material] = 0;
    end
    Act_MaterialRate_Draw_nm(:,act) = material;
    
    if Act_PERT_EqRate_m3(act,2) == Act_PERT_EqRate_m3(act,1) && Act_PERT_EqRate_m3(act,2)>0
        [equipment] = Act_PERT_EqRate_m3(act,2);
    elseif Act_PERT_EqRate_m3(act,2)>0
        [equipment]=RandPert (Act_PERT_EqRate_m3(act,1),Act_PERT_EqRate_m3(act,2),Act_PERT_EqRate_m3(act,3),n_iter); %Draw value [z] from distribution based on Min/ML/Max for EquipmetRate
    else
        [equipment] = 0;
    end
    Act_EquipmentRate_Draw_nm(:,act) = equipment;
    
    if Act_PERT_SubRate_m3(act,2) == Act_PERT_SubRate_m3(act,1) && Act_PERT_SubRate_m3(act,2)>0
        [subcontractor] = Act_PERT_SubRate_m3(act,2);
    elseif Act_PERT_SubRate_m3(act,2)>0
        [subcontractor]=RandPert (Act_PERT_SubRate_m3(act,1),Act_PERT_SubRate_m3(act,2),Act_PERT_SubRate_m3(act,3),n_iter); %Draw value [z] from distribution based on Min/ML/Max for SubcontractorRate
    else
        [subcontractor] = 0;
    end
    Act_SubcontractRate_Draw_nm(:,act) = subcontractor;
    
    Act_CostDraw_nm(:,act)=(Act_LabourDraw_nm(:,act).*Act_LabourCostDraw_nm(:,act)+(Act_Quantity_Draw_nm(:,act).*(Act_MaterialRate_Draw_nm(:,act)+ Act_EquipmentRate_Draw_nm(:,act)+Act_SubcontractRate_Draw_nm(:,act))));
end

Act_ProjectCostDraw=sum(Act_CostDraw_nm,2);


%% Drawing Values for Risk events
Risk_Drawn_Values_no = zeros(n_iter,o_Risk); %Pre-allocate matrix for all drawn values for risk
Risk_Fired_no = rand(n_iter,o_Risk)<Risk_Prob_o; % Generate logical matrix of random values being lower than risk probability
Risk_Effect = zeros(n_iter,o_Risk);
for risks = 1:o_Risk
    Risk_Drawn_Values_no(:,risks)=RandPert(Risk_PERT_Effect_o3(risks,1),Risk_PERT_Effect_o3(risks,2),Risk_PERT_Effect_o3(risks,3),n_iter);
    Risk_Effect(:,risks) = Risk_Drawn_Values_no(:,risks).*Risk_Fired_no(:,risks);
end

Risk_Act_Total_Effect = sum(Risk_Effect,2); 
%% Drawing Values for Chance events
Chance_Drawn_Values_no = zeros(n_iter,q_Chance); %Pre-allocate matrix for all drawn values for risk
Chance_Fired_no = rand(n_iter,q_Chance)<Chance_Prob_q; % Generate logical matrix of random values being lower than risk probability
Chance_Effect = zeros(n_iter,q_Chance);
for chance = 1:q_Chance
    Chance_Drawn_Values_no(:,chance)=RandPert(Chance_PERT_Effect_o3(chance,1),Chance_PERT_Effect_o3(chance,2),Chance_PERT_Effect_o3(chance,3),n_iter);
    Chance_Effect(:,chance) = Chance_Drawn_Values_no(:,chance).*Chance_Fired_no(:,chance);
end

%% Generate item costs with Risks
Act_Project_Cost_Risk = sum(Act_ProjectCostDraw,2)-sum(Chance_Effect,2)+sum(Risk_Effect,2);

%% P70 Element

[f,x] = ecdf(Act_Project_Cost_Risk);
A = [f,x];
P70 = interp1(f,x,0.7,'spline');
TargetCost = P70;

%% Re-adusting probabilities Risks & Chances for reality check

Real_Chance_fired = Chance_Fired_no;
Real_Chance_fired(:,1) = zeros(n_iter,1);
Real_Chance_Effect = Chance_Drawn_Values_no.*Real_Chance_fired;

Real_Risk_fired = Risk_Fired_no;
Real_Risk_fired(:,[3,12]) = ones(n_iter,2);
Real_Risk_Effect = Risk_Drawn_Values_no.*Real_Risk_fired; 

Current_Project_Cost = sum(Act_ProjectCostDraw,2)-sum(Real_Chance_Effect,2)+sum(Real_Risk_Effect,2);
%% Drawing Random Values for the mitigations

Mit_Drawn_Labour_Values_np = zeros(n_iter,p_Mitigations); %Pre-Allocate matrix for drawn mitigation labour values
Mit_Drawn_LabourCost_Values_np = zeros(n_iter,p_Mitigations); %Pre-Allocate matrix for drawn mitigation labour cost values
Mit_Drawn_Quantity_Values_np = zeros(n_iter,p_Mitigations); %Pre-Allocate matrix for drawn mitigation materials values
Mit_Drawn_MatRate_Values_np = zeros(n_iter,p_Mitigations); %Pre-Allocate matrix for drawn mitigation unitcost values
Mit_Drawn_EqRate_Values_np = zeros(n_iter,p_Mitigations); %Pre-Allocate matrix for drawn mitigation unitcost values
Mit_Drawn_SubRate_Values_np = zeros(n_iter,p_Mitigations); %Pre-Allocate matrix for drawn mitigation unitcost values

for mitigations=1:p_Mitigations
    if Mit_PERT_Labour_Con_p3(mitigations,2) > 0
        [Mitigation_Labour_Draw]=RandPert(Mit_PERT_Labour_Con_p3(mitigations,1),Mit_PERT_Labour_Con_p3(mitigations,2),Mit_PERT_Labour_Con_p3(mitigations,3),n_iter);
    else
        [Mitigation_Labour_Draw]=0;
    end
    Mit_Drawn_Labour_Values_np(:,mitigations)= Mitigation_Labour_Draw;
end
for mitigations=1:p_Mitigations
    if Mit_PERT_LabourCost_Con_p3(mitigations,2) > 0
        [Mitigation_LabourCost_Draw]=RandPert(Mit_PERT_LabourCost_Con_p3(mitigations,1),Mit_PERT_LabourCost_Con_p3(mitigations,2),Mit_PERT_LabourCost_Con_p3(mitigations,3),n_iter);
    else
        [Mitigation_LabourCost_Draw]= 0;
    end
    Mit_Drawn_LabourCost_Values_np(:,mitigations)= Mitigation_LabourCost_Draw ;
end

for mitigations=1:p_Mitigations
    if Mit_PERT_Quantity_Con_p3(mitigations,2) > 0
        [Mitigation_Quantity_Draw]=RandPert(Mit_PERT_Quantity_Con_p3(mitigations,1),Mit_PERT_Quantity_Con_p3(mitigations,2),Mit_PERT_Quantity_Con_p3(mitigations,3),n_iter);
    else
        [Mitigation_Quantity_Draw]=0;
    end
    Mit_Drawn_Quantity_Values_np(:,mitigations)= Mitigation_Quantity_Draw;
end
for mitigations=1:p_Mitigations
    if Mit_PERT_MatRate_Con_p3(mitigations,2) > 0
        [Mitigation_MatRate_Draw]=RandPert(Mit_PERT_MatRate_Con_p3(mitigations,1),Mit_PERT_MatRate_Con_p3(mitigations,2),Mit_PERT_MatRate_Con_p3(mitigations,3),n_iter);
    else
        [Mitigation_MatRate_Draw]=0;
    end
    Mit_Drawn_MatRate_Values_np(:,mitigations)=Mitigation_MatRate_Draw;
end
for mitigations=1:p_Mitigations
    if Mit_PERT_EqRate_Con_p3(mitigations,2) > 0
        [Mitigation_EqRate_Draw]=RandPert(Mit_PERT_EqRate_Con_p3(mitigations,1),Mit_PERT_EqRate_Con_p3(mitigations,2),Mit_PERT_EqRate_Con_p3(mitigations,3),n_iter);
    else
        [Mitigation_EqRate_Draw]=0;
    end
    Mit_Drawn_EqRate_Values_np(:,mitigations)=Mitigation_EqRate_Draw;
end
for mitigations=1:p_Mitigations
    if Mit_PERT_SubRate_Con_p3(mitigations,2) > 0
        [Mitigation_SubRate_Draw]=RandPert(Mit_PERT_SubRate_Con_p3(mitigations,1),Mit_PERT_SubRate_Con_p3(mitigations,2),Mit_PERT_SubRate_Con_p3(mitigations,3),n_iter);
    else
        [Mitigation_SubRate_Draw]=0;
    end
    Mit_Drawn_SubRate_Values_np(:,mitigations)=Mitigation_SubRate_Draw;
end

%% Allocate mitigation to activity it has effect on
Mit_Match_Matrix_pm = zeros(p_Mitigations,m_Activities); %Pre-Allocate matrix for matching Mitigations and activity

for mit = 1:p_Mitigations
    for act=1:m_Activities
        if Mit_Effected_Activity_p(mit) == act
            Mit_Match_Matrix_pm(mit,act)=1;
        else
            Mit_Match_Matrix_pm(mit,act)=0;
        end
    end
end  %Define which mitigation has an effect on which activity


%% Generate Effects mitigations
Mitigated_Value_Labour_pm = zeros(p_Mitigations,m_Activities);
Mitigated_Value_LabourCost_pm = zeros(p_Mitigations,m_Activities);
Mitigated_Value_Quantity_pm = zeros(p_Mitigations,m_Activities);
Mitigated_Value_MatRate_pm = zeros(p_Mitigations,m_Activities);
Mitigated_Value_EqRate_pm = zeros(p_Mitigations,m_Activities);
Mitigated_Value_SubRate_pm = zeros(p_Mitigations,m_Activities);
Mitigation_Effect_Total_pm = zeros(p_Mitigations,m_Activities);

Tussenstap_Labour = zeros(p_Mitigations,m_Activities);
Tussenstap_LabourCost = zeros(p_Mitigations,m_Activities);
Tussenstap_Quantity = zeros(p_Mitigations,m_Activities);
Tussenstap_MatRate = zeros(p_Mitigations,m_Activities);
Tussenstap_EqRate = zeros(p_Mitigations,m_Activities);
Tussenstap_SubRate = zeros(p_Mitigations,m_Activities);
Tussenstap_Overall = zeros(p_Mitigations,m_Activities);
Tussenstap_Overall_Neg_eff= zeros(p_Mitigations,m_Activities);

Mit_Effect_Factor = zeros(n_iter,p_Mitigations);
Mit_Effect_Per_Iteration_np = zeros(n_iter,p_Mitigations);
Mitigated_Effect_Per_Mitigation = zeros(n_iter,p_Mitigations);
Item_Cost_With_All_mits = Act_CostDraw_nm;
Project_Cost_After_All_Mitigations = zeros (n_iter,1);

for iter = 1:n_iter
    for mit = 1:p_Mitigations
        Tussenstap_Labour(mit,:) = Mit_Match_Matrix_pm(mit,:)*(Mit_Drawn_Labour_Values_np(iter,mit)+Mit_Labour_Var_p(mit)*Act_LabourDraw_nm(iter,Mit_Effected_Activity_p(mit)));
        Tussenstap_LabourCost(mit,:) = Mit_Match_Matrix_pm(mit,:)*(Mit_Drawn_LabourCost_Values_np(iter,mit)+Mit_LabourCost_Var_p(mit)*Act_LabourCostDraw_nm(iter,Mit_Effected_Activity_p(mit)));
        Tussenstap_Quantity(mit,:) = Mit_Match_Matrix_pm(mit,:)*(Mit_Drawn_Quantity_Values_np(iter,mit)+Mit_Quantity_Var_p(mit)*Act_Quantity_Draw_nm(iter,Mit_Effected_Activity_p(mit)));
        Tussenstap_MatRate(mit,:) = Mit_Match_Matrix_pm(mit,:)*(Mit_Drawn_MatRate_Values_np(iter,mit)+Mit_MatRate_Var_p(mit)*Act_MaterialRate_Draw_nm(iter,Mit_Effected_Activity_p(mit)));
        Tussenstap_EqRate(mit,:) = Mit_Match_Matrix_pm(mit,:)*(Mit_Drawn_EqRate_Values_np(iter,mit)+Mit_EqRate_Var_p(mit)*Act_EquipmentRate_Draw_nm(iter,Mit_Effected_Activity_p(mit)));
        Tussenstap_SubRate(mit,:) = Mit_Match_Matrix_pm(mit,:)*(Mit_Drawn_SubRate_Values_np(iter,mit)+Mit_SubRate_Var_p(mit)*Act_SubcontractRate_Draw_nm(iter,Mit_Effected_Activity_p(mit)));
        
        Tussenstap_Overall(mit,:) = Tussenstap_Labour(mit,:)+Tussenstap_LabourCost(mit,:)+ Tussenstap_Quantity(mit,:)+ Tussenstap_MatRate(mit,:)+ Tussenstap_EqRate(mit,:)+Tussenstap_SubRate(mit,:);
        Tussenstap_Overall_Neg_eff(mit,:)= Tussenstap_Overall(mit,:)/(10^3*Mit_Neg_Effect_p(mit));
        
        
        Mitigated_Value_Labour_pm(mit,:) =Act_LabourDraw_nm(iter,:) - Tussenstap_Labour(mit,:);
        Mitigated_Value_LabourCost_pm(mit,:) = Act_LabourCostDraw_nm(iter,:) - Tussenstap_LabourCost(mit,:) ;
        Mitigated_Value_Quantity_pm(mit,:) = Act_Quantity_Draw_nm(iter,:) - Tussenstap_Quantity(mit,:);
        Mitigated_Value_MatRate_pm(mit,:) = Act_MaterialRate_Draw_nm(iter,:) - Tussenstap_MatRate(mit,:);
        Mitigated_Value_EqRate_pm(mit,:) = Act_EquipmentRate_Draw_nm(iter,:) - Tussenstap_EqRate(mit,:);
        Mitigated_Value_SubRate_pm(mit,:) = Act_SubcontractRate_Draw_nm(iter,:) - Tussenstap_SubRate(mit,:);
        for act = 1:m_Activities
            Mitigation_Effect_Total_pm (mit,act) = (Mitigated_Value_Labour_pm(mit,act).* Mitigated_Value_LabourCost_pm(mit,act)+Mitigated_Value_Quantity_pm(mit,act)*(Mitigated_Value_MatRate_pm(mit,act)+Mitigated_Value_EqRate_pm(mit,act)+Mitigated_Value_SubRate_pm(mit,act))) ;
        end
        Mitigated_Effect_Per_Mitigation(iter,mit) = sum(Act_CostDraw_nm(iter,:)-Mitigation_Effect_Total_pm(mit,:),2);
        Mit_Effect_Per_Iteration_np(iter,:) = sum(Mitigation_Effect_Total_pm,2).'; % eventually this will state the new value of all four elements after
        
        Mit_Effect_Factor(iter,mit) = Mitigated_Effect_Per_Mitigation(iter,mit)/(10^3*Mit_Neg_Effect_p(mit));
    end
    Item_Cost_With_All_mits(iter,:) = Item_Cost_With_All_mits(iter,:) - sum((Mitigated_Effect_Per_Mitigation(iter,:).'.*Mit_Match_Matrix_pm),1);
    Project_Cost_After_All_Mitigations(iter)= Current_Project_Cost(iter) - sum(Mitigated_Effect_Per_Mitigation(iter,:),2);
end %The cost reduction for all mitigation measures are defined for all iterations as well as effectiveness factor
B = ones(n_iter,p_Mitigations)./Mit_Effect_Factor;

%% What happens with cost overrun
Test=zeros(n_iter,1);
Project_Cost_After_Mitigations = zeros (n_iter,1);
Strat_Optim = zeros(n_iter,p_Mitigations+1);
Optim_Cost_Red = zeros(n_iter,p_Mitigations);
Cost_Reduction_Total = zeros(n_iter,1);

F = zeros(n_iter,1);
delta = zeros(n_iter,1);
for iter = 1:n_iter
    if Current_Project_Cost(iter) > TargetCost
        F(iter) = 1;
        Neg_Eff_Dur = Mit_Neg_Effect_p.';
        [S,OCR,exitflag]=opt_mit_lin(Mitigated_Effect_Per_Mitigation,B,Current_Project_Cost,TargetCost,p_Mitigations,Neg_Eff_Dur,iter,Neg_Eff_Max);
        delta(iter)=S(p_Mitigations+1);
        Strat_Optim(iter,:)=S.';
        Optim_Cost_Red(iter,:)=S(1:p_Mitigations).'.*Mitigated_Effect_Per_Mitigation(iter,:);
        Cost_Reduction_Total(iter)=sum(Optim_Cost_Red(iter,:),2);
        Project_Cost_After_Mitigations(iter)= Current_Project_Cost(iter)-Cost_Reduction_Total(iter);
    else
        Project_Cost_After_Mitigations(iter)= Current_Project_Cost(iter);
    end
end
%% Mitigation Use
Strat_Optim = round(Strat_Optim,0);
Count_Mits = sum(Strat_Optim,1);
Count_Over_Budget = sum(F);

Delta_Average = sum(delta) / Count_Over_Budget;

Count_Used = [MitNames_p,Count_Mits(1:p_Mitigations).'];
Percent_Used = [MitNames_p,(Count_Mits(1:p_Mitigations)./Count_Over_Budget).'];

%% Calculating Effect of strategies
%In here we calculate what the actual effect is of the optimal mitigation
%strategy. As there should be an advice at the end of the model for a
%strategy, one should compare the strategy to the already existing
%strategies that are being used. As a PM will choose a certain strategy,
%this will be implemented for all runs as this will aid the probability to
%stay within budget. 

% Define Optimal strategy based on optimisation count
[strat,fval, exitflag2]=intlinprog(-Count_Mits(1:p_Mitigations),1:p_Mitigations,Mit_Neg_Effect_p_inv,Neg_Eff_Max,[],[],zeros(p_Mitigations,1),ones(p_Mitigations,1),[],optimoptions('intlinprog','MaxTime',5,'Display','off')); 
Strat_Final = strat.'; 
Strat_Final_Names = [strat,MitNames_p];
Strat_Final_Names = sortrows(Strat_Final_Names,1,'descend');
Strat_Final_Names = Strat_Final_Names(1:int8(sum(strat)),:);
Partial_Effect_Strat_Final = Mitigated_Effect_Per_Mitigation.*Strat_Final;
Total_Effect_Strat_Final = sum(Partial_Effect_Strat_Final,2); 
Project_Cost_Reducted_Final_Strat = Current_Project_Cost - Total_Effect_Strat_Final;

% Strategy testing
Oefen = Count_Mits(1:p_Mitigations).';
Count_Neg_Mits = [Oefen,Mit_Neg_Effect_p];
[Sort_Count_Mits, index1] = sortrows(Count_Neg_Mits,1,'descend'); 
Sort_Names_Opt = MitNames_p(index1,:);
mit = 1; 
neg = 0;
while neg + Sort_Count_Mits(mit,2) < Neg_Eff_Max && mit < p_Mitigations
        neg = neg + Sort_Count_Mits(mit,2);
        mit = mit + 1;
end
Test_Strat = [Sort_Names_Opt(1:mit)];
Test_Mit_Numbers = [index1(1:mit)];

Mit_FinalStrat_Reduction = sum(Mitigated_Effect_Per_Mitigation(:,Test_Mit_Numbers),2);
Project_Cost_Reduced_Test = Current_Project_Cost - Mit_FinalStrat_Reduction;


% Define Maximise reduction (given max negative effect) -> top values first
Mitigation_Effect_Average = (sum(Mitigated_Effect_Per_Mitigation)/n_iter).';
Mitigation_Av_Neg = [Mitigation_Effect_Average,Mit_Neg_Effect_p];
[Sort_Mits, index2] = sortrows(Mitigation_Av_Neg, 1 , 'descend'); 
Sorted_Names_Mit = MitNames_p(index2,:);
Old_Strat = strings(p_Mitigations);
mit = 1; 
neg = 0;
red = zeros(p_Mitigations,1);
while neg + Sort_Mits(mit,2) < Neg_Eff_Max && mit < p_Mitigations
        red(mit) = Sort_Mits(mit,1);
        neg = neg + Sort_Mits(mit,2);
        mit = mit + 1;
end
red(red==0) = [];
Strategy_Old = [red,Sorted_Names_Mit(1:size(red))];
Total_Effect_Orig_Strat = sum(str2double(Strategy_Old(:,1)));
%[Res,fval2, exitflag3]=intlinprog(-Mitigation_Effect_Average,1:p_Mitigations,[Mit_Neg_Effect_p_inv],Neg_Eff_Max,[],[],zeros(p_Mitigations,1),ones(p_Mitigations,1),[],optimoptions('intlinprog','MaxTime',5,'Display','off')); 
%Partial_Effect_Orig_Strat = Mitigated_Effect_Per_Mitigation.*Res.';
%Total_Effect_Orig_Strat = sum(Partial_Effect_Orig_Strat,2); 
Project_Cost_Reducted_Orig_Strat = Current_Project_Cost - Total_Effect_Orig_Strat;
%% Results for Changing probability for original P70 value
[g,y] = ecdf(Project_Cost_Reducted_Orig_Strat); %Define the ecdf for the Original mitigation defining strategy
[h,z] = ecdf(Project_Cost_Reducted_Final_Strat); %Define the ecdf for the Optimal mitigation defining strategy
[h2,z2] = ecdf(Project_Cost_After_All_Mitigations); %Define the ecdf for the All mitigations
[h3,z3] = ecdf(Current_Project_Cost);

g= g(2:end,:);
y= y(2:end,:);
h = h(2:end,:);
z = z(2:end,:);
h2 = h2(2:end,:);
z2 = z2(2:end,:);
h3 = h3(2:end,:);
z3 = z3(2:end,:);

Prob_All_Mits = interp1(z2,h2,TargetCost) ; 
Prob_Old_Mits = interp1(y,g,TargetCost) ;
Prob_Opt_Mits = interp1(z,h,TargetCost) ;
Prob_Cur_Sit = interp1(z3,h3,TargetCost);

%% Statistical Tests
%t-test current vs optimised
[sig_1,p1,ci1,stats1] = ttest(Current_Project_Cost,Project_Cost_Reducted_Final_Strat); 
[d_cohen_Curr_Opt] = (mean(Current_Project_Cost)-mean(Project_Cost_Reducted_Final_Strat))/std(Current_Project_Cost);
r_Curr_Opt = d_cohen_Curr_Opt/sqrt(d_cohen_Curr_Opt^2+4);

%Optimised vs Original method
[sig_2,p2,ci2,stats2] = ttest(Project_Cost_Reducted_Orig_Strat,Project_Cost_Reducted_Final_Strat); 
[d_cohen_Opt_Orig] = (mean(Project_Cost_Reducted_Orig_Strat)-mean(Project_Cost_Reducted_Final_Strat))/std(Project_Cost_Reducted_Final_Strat);
r_Opt_Orig = d_cohen_Opt_Orig/sqrt(d_cohen_Opt_Orig^2+4);
%% Functions
function [x]=RandPert(a,m,b,n_iter) % this function is for drawing random number from Pert-Beta dist
%mu=(a+4*m+b)/6; % mean value
%sd=(b-a)/2; %standard deviation
alpha=1+4*(m-a)/(b-a);% First Beta parameter computed using the three points optimistic, most likely, and pessimistic values
beta=1+4*(b-m)/(b-a);% Second Beta parameter computed using the three points optimistic, most likely, and pessimistic values
x=randraw('Beta',[a b alpha beta],n_iter); %draw a random number according to the Beta-Pert distribution
end

function [x,fval, exitflag]=opt_mit_lin(Mitigated_Effect_Per_Mitigation,B,Current_Project_Cost,TargetCost,p_Mitigations,Neg_Eff_Dur,iter,Neg_Eff_Max)
%--- Input
%--- Definition of constrains
A = [-Mitigated_Effect_Per_Mitigation(iter,:),-1;Neg_Eff_Dur,0];
b = [(-Current_Project_Cost(iter)+TargetCost);Neg_Eff_Max];
%--- Boundary constraints
lb = zeros(1,p_Mitigations+1) ; %All zeros
ub = [ones(1,p_Mitigations),20e6]; %All ones for the mitigations
%--- Objective function
Obj_Function = [1./B(iter,:),10^7];
%--- Options
intcon = 1:p_Mitigations;
options = optimoptions('intlinprog','MaxTime',60,'Display','off');
%--- Optimatization function
[x,fval, exitflag] = intlinprog(Obj_Function,intcon,A,b,[],[],lb,ub,[],options); %multi interger linear optimization function
end
