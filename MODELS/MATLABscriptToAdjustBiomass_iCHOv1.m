% You must load a new iCHO1921s model generated with the Jupyter Notebooks before running this script
aas = {'pro_L[c]','met_L[c]','gly[c]','phe_L[c]','glu_L[c]','ala_L[c]','asn_L[c]','cys_L[c]','gln_L[c]','ser_L[c]','thr_L[c]','arg_L[c]','asp_L[c]','his_L[c]','leu_L[c]','ile_L[c]','lys_L[c]','val_L[c]','trp_L[c]','tyr_L[c]'}';
aas_ids = findMetIDs(model,aas);
biomass_rxn_id = findRxnIDs(iCHO1921s,'biomass_cho_producing');
new = iCHO1921s;
for i=1:length(aas_ids)
    
    new.S(aas_ids(i),biomass_rxn_id) = full(model.S(aas_ids(i),biomass_rxn_id))*0.8;

end

% The variable "new" now contains an iCHO1921s model with an adjusted biomass protein