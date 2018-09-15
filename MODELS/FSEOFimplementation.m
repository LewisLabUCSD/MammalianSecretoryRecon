load('FSEOFmodels.mat')
initCobraToolbox
changeCobraSolver('glpk');

%% Bidirectional FSEOF

% Create new choRITUXIMAB model
new_model = choEPO;

% Set Deadend reactions to zero flux
new_model = changeRxnBounds(new_model,removedRxns,0,'b');
new_model = changeRxnBounds(new_model,exchanges_to_close,'l',0);

% Constrain uptake rates only
new_model = changeRxnBounds(new_model,uptakeRxns,uptakeBounds,'b');

% Unbound excretion reactions
new_model = changeRxnBounds(new_model,excretionRxns,0,'l');
new_model = changeRxnBounds(new_model,excretionRxns,excretionBounds*2,'u');

enforced_obj = new_model.rxns(new_model.c==1);

% Calculate maximum theoretical IgG secretion
max_igg = optimizeCbModel(new_model);
max_igg = max_igg.f;

% Set min_IgG to an nth of maximum
n = 4;
min_igg = max_igg/n;

% Set objective to growth
new_model = changeObjective(new_model,'biomass_cho_producing');

% Enforce production of IgG to its minimum
new_model = changeRxnBounds(new_model,enforced_obj,min_igg,'b');
 
% Scan fluxes n times from igg_min to igg_max with FBA
FSEOF_matrix = zeros(length(new_model.rxns),n);

% Initialize enforced constraint
k = min_igg;

% Compute values in FSEOF_matrix
qp = min_igg:max_igg/n:max_igg;
mu = zeros(1,n);

for i=1:n
    fseof = optimizeCbModel(new_model,'max',0,1); % Allow loops
    mu(i) = fseof.f; 
    FSEOF_matrix(:,i) = fseof.x;
    k = k + (max_igg / n);
    new_model = changeRxnBounds(new_model,enforced_obj,k,'b');
end

% Process FSEOF_matrix to Identify relevant reactions coupled to production
fseof_target_rxns_ids = [];
fseof_target_fluxes = [];

for i=1:size(FSEOF_matrix,1)
    v = FSEOF_matrix(i,:);
    if abs(v(n)) > abs(v(1)) && v(n)*v(1) > 0 && sum(abs(v)) < 1000
        fseof_target_rxns_ids = [fseof_target_rxns_ids;i];
        fseof_target_fluxes = [fseof_target_fluxes;v];
    end
end

% Keep metabolic reactions only
% % k = length(model.rxns); % Number of metabolic reactions
% % fseof_target_fluxes = fseof_target_fluxes(fseof_target_rxns_ids<=k,:);
% % fseof_target_rxns_ids = fseof_target_rxns_ids(fseof_target_rxns_ids<=k);

fseof_target_rxns = new_model.rxns(fseof_target_rxns_ids);
fseof_target_rxns_pathway = new_model.subSystems(fseof_target_rxns_ids);

%% FVA analysis of identified reactions

% Enforce production of IgG to its minimum
new_model = changeRxnBounds(new_model,enforced_obj,min_igg,'b');

% Initialize fseo_fva_max and fseo_fva_min values
fseof_fva_max = zeros(length(fseof_target_rxns),n);
fseof_fva_min = zeros(length(fseof_target_rxns),n);

% Initialize enforced constraint
k = min_igg;

c = clock;
tic %Time it
for i=1:n
    [minFluxLoopy,maxFluxLoopy] = fluxVariability(new_model,100,'max',fseof_target_rxns,0,1);
    fseof_fva_min(:,i) = minFluxLoopy; 
    fseof_fva_max(:,i) = maxFluxLoopy;
    k = k + (max_igg / n);
    new_model = changeRxnBounds(new_model,enforced_obj,k,'b');
end
toc
d = clock;
fix(d) - fix(c)

