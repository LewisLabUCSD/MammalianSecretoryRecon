function new_model = addSecretionReactions(model, rxnFile, rxnNames, rxnGPRs, geneList)
% addSecretionReactions Expands a metabolic model to include secretory
% reactions for a given protein
%
%INPUTS
% model = COBRA model structure
% rxnFile = String of .txt file containing reactions list
% rxnNames = String of .txt file containing reactions names
% rxnGPRs = String of .txt file containing reactions GPRs
% geneList = String of .txt file containing genes to be added to model

%OUTPUT
% new_model = Expanded model with Metabolic and Secretory capabilities

%EXAMPLES
% for CHO:
% choEPO = addSecretionReactions(model, 'CHO_P01588_reactions.txt', 'CHO_P01588_names.txt', 'CHO_P01588_GPRs.txt', 'geneIDs_CHO.txt')

% for HUMAN:
% humanEPO = addSecretionReactions(recon1, 'HUMAN_P01588_reactions.txt', 'HUMAN_P01588_names.txt', 'HUMAN_P01588_GPRs.txt', 'geneIDs_HUMAN.txt') 

% NOTE: You need to run 'initCobraToolbox' first!

% Jahir M. Gutierrez 08/24/2016

warning off all

new_model = model;
name = {};
reaction = {};
gpr = {};
gene = {};
startGPRs = length(new_model.rxns);
startGenes = length(new_model.genes);

% Import reactions into a cell array
fileID = fopen(rxnNames);
i=1;
while feof(fileID)==0
    name{i} = fgetl(fileID);
    i = i+1;
end
fclose(fileID);

% Import reaction names into a cell array
fileID = fopen(rxnFile);
i=1;
while feof(fileID)==0
    reaction{i} = fgetl(fileID);
    i=i+1;
end
fclose(fileID);

% Import reaction GPRs into a cell array
fileID = fopen(rxnGPRs);
i=1;
while feof(fileID)==0
    gpr{i} = fgetl(fileID);
    i=i+1;
end
fclose(fileID);

% Import gene IDs into a cell array
fileID = fopen(geneList);
i=1;
while feof(fileID)==0
    gene{i} = fgetl(fileID);
    i=i+1;
end
fclose(fileID);

% Add protein secretion reactions to iCHO model and modify grule Matrix
for i=1:length(reaction)
    new_model = addReaction(new_model,name{i},reaction{i});
    new_model = changeGeneAssociation(new_model,name{i},gpr{i});
end

% % Add GPRs
% k=1;
% for i=(startGPRs+1):(startGPRs+length(gpr))
%     new_model.grRules{i} = gpr{k};
%     k = k+1;
% end
% 
% % Add Genes
% k=1;
% for i=(startGenes+1):(startGenes+length(gene))
%     new_model.genes{i} = gene{k};
%     k = k+1;
% end
% 
% % Modify grule Matrix
% for i=1:length(gpr)
%     new_model = changeGeneAssociation(new_model,name{i},gpr{i});
% end
obj_identifier = findRxnIDs(new_model,'SEC61C') - 1;
new_model = changeObjective(new_model,new_model.rxns(obj_identifier));

warning on all


