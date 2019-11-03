%load('Figure2C.mat')
%load('exchanges_to_close.mat')
ProtMolMass = [44702; 49313; 21307; 51234.9; 267009; 22294; 143859.7; 62917];
Prot = {'BMP2'; 'BMP7'; 'EPO'; 'ETANERCEPT'; 'F8'; 'IFNB1'; 'RITUXIMAB'; 'TPA'};
k = 1;
results_mu = zeros(20,8);
results_qp = zeros(20,8);

%%%
model = choBMP2;
model = changeRxnBounds(model,exchanges_to_close,0,'l');

for i=1:length(selvarasu_rxns)
    if selvarasu_bounds(i,1) < 0
        model = changeRxnBounds(model,selvarasu_rxns(i),selvarasu_bounds(i,1),'l');
        model = changeRxnBounds(model,selvarasu_rxns(i),0,'u');
    else
        model = changeRxnBounds(model,selvarasu_rxns(i),selvarasu_bounds(i,1),'u');
        model = changeRxnBounds(model,selvarasu_rxns(i),0,'l');
    end
end
[results_mu(:,k), results_qp(:,k)] = robustnessAnalysis(model,'biomass_cho_producing');
k = k + 1;
%%%
model = choBMP7;
model = changeRxnBounds(model,exchanges_to_close,0,'l');

for i=1:length(selvarasu_rxns)
    if selvarasu_bounds(i,1) < 0
        model = changeRxnBounds(model,selvarasu_rxns(i),selvarasu_bounds(i,1),'l');
        model = changeRxnBounds(model,selvarasu_rxns(i),0,'u');
    else
        model = changeRxnBounds(model,selvarasu_rxns(i),selvarasu_bounds(i,1),'u');
        model = changeRxnBounds(model,selvarasu_rxns(i),0,'l');
    end
end
[results_mu(:,k), results_qp(:,k)] = robustnessAnalysis(model,'biomass_cho_producing');
k = k + 1;%%%
model = choEPO;
model = changeRxnBounds(model,exchanges_to_close,0,'l');

for i=1:length(selvarasu_rxns)
    if selvarasu_bounds(i,1) < 0
        model = changeRxnBounds(model,selvarasu_rxns(i),selvarasu_bounds(i,1),'l');
        model = changeRxnBounds(model,selvarasu_rxns(i),0,'u');
    else
        model = changeRxnBounds(model,selvarasu_rxns(i),selvarasu_bounds(i,1),'u');
        model = changeRxnBounds(model,selvarasu_rxns(i),0,'l');
    end
end
[results_mu(:,k), results_qp(:,k)] = robustnessAnalysis(model,'biomass_cho_producing');
k = k + 1;
%%%
model = choETANERCEPT;
model = changeRxnBounds(model,exchanges_to_close,0,'l');

for i=1:length(selvarasu_rxns)
    if selvarasu_bounds(i,1) < 0
        model = changeRxnBounds(model,selvarasu_rxns(i),selvarasu_bounds(i,1),'l');
        model = changeRxnBounds(model,selvarasu_rxns(i),0,'u');
    else
        model = changeRxnBounds(model,selvarasu_rxns(i),selvarasu_bounds(i,1),'u');
        model = changeRxnBounds(model,selvarasu_rxns(i),0,'l');
    end
end
[results_mu(:,k), results_qp(:,k)] = robustnessAnalysis(model,'biomass_cho_producing');
k = k + 1;
%%%
model = choF8;
model = changeRxnBounds(model,exchanges_to_close,0,'l');

for i=1:length(selvarasu_rxns)
    if selvarasu_bounds(i,1) < 0
        model = changeRxnBounds(model,selvarasu_rxns(i),selvarasu_bounds(i,1),'l');
        model = changeRxnBounds(model,selvarasu_rxns(i),0,'u');
    else
        model = changeRxnBounds(model,selvarasu_rxns(i),selvarasu_bounds(i,1),'u');
        model = changeRxnBounds(model,selvarasu_rxns(i),0,'l');
    end
end
[results_mu(:,k), results_qp(:,k)] = robustnessAnalysis(model,'biomass_cho_producing');
k = k + 1;
%%%
model = choIFNB1;
model = changeRxnBounds(model,exchanges_to_close,0,'l');

for i=1:length(selvarasu_rxns)
    if selvarasu_bounds(i,1) < 0
        model = changeRxnBounds(model,selvarasu_rxns(i),selvarasu_bounds(i,1),'l');
        model = changeRxnBounds(model,selvarasu_rxns(i),0,'u');
    else
        model = changeRxnBounds(model,selvarasu_rxns(i),selvarasu_bounds(i,1),'u');
        model = changeRxnBounds(model,selvarasu_rxns(i),0,'l');
    end
end
[results_mu(:,k), results_qp(:,k)] = robustnessAnalysis(model,'biomass_cho_producing');
k = k + 1;
%%%
model = choRITUXIMAB;
model = changeRxnBounds(model,exchanges_to_close,0,'l');

for i=1:length(selvarasu_rxns)
    if selvarasu_bounds(i,1) < 0
        model = changeRxnBounds(model,selvarasu_rxns(i),selvarasu_bounds(i,1),'l');
        model = changeRxnBounds(model,selvarasu_rxns(i),0,'u');
    else
        model = changeRxnBounds(model,selvarasu_rxns(i),selvarasu_bounds(i,1),'u');
        model = changeRxnBounds(model,selvarasu_rxns(i),0,'l');
    end
end
[results_mu(:,k), results_qp(:,k)] = robustnessAnalysis(model,'biomass_cho_producing');
k = k + 1;
%%%
model = choTPA;
model = changeRxnBounds(model,exchanges_to_close,0,'l');

for i=1:length(selvarasu_rxns)
    if selvarasu_bounds(i,1) < 0
        model = changeRxnBounds(model,selvarasu_rxns(i),selvarasu_bounds(i,1),'l');
        model = changeRxnBounds(model,selvarasu_rxns(i),0,'u');
    else
        model = changeRxnBounds(model,selvarasu_rxns(i),selvarasu_bounds(i,1),'u');
        model = changeRxnBounds(model,selvarasu_rxns(i),0,'l');
    end
end
[results_mu(:,k), results_qp(:,k)] = robustnessAnalysis(model,'biomass_cho_producing');
k = k + 1;

% PLOT FIGURES
co = [0 0 1;
      0 0.5 0;
      1 0 0;
      0 0.75 0.75;
      0.75 0 0.75;
      0.75 0.75 0;
      0.25 0.25 0.25;
      1 0.75 0.5];
set(groot,'defaultAxesColorOrder',co)

figure
plot(results_mu,results_qp,'-o','lineWidth',1.5)
xlabel('Growth rate (1/h)')
ylabel('Qp (mmol/gDW/h)')
legend(Prot)

figure
hold on

results_pcd = zeros(20,8);
for i=1:length(Prot)
    y = results_qp(:,i) .* ProtMolMass(i) ./ 1000 .* 456 .* 24; % Assume 1 cell = 456 pg
    results_pcd(:,i) = y;
end

plot(results_mu,results_pcd,'-','lineWidth',1.5)
xlabel('Growth rate (1/h)')
ylabel('PCD (pg/cell/day)')
legend(Prot)
