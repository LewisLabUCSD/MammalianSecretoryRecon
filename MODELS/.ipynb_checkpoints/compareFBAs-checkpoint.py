import cobra
import pandas as pd

import iCHO2048s_Builder

# Load CHO metabolic model
model = cobra.io.load_matlab_model('../iCHO2048s_matlab/iCHOv1.mat')

# Load CHO PSIM
psim = pd.read_csv('PSIM_CHO.tab', sep='\t')

# Get all protein IDs of secreted proteins
psim = psim[psim['SP']==1].reset_index(drop=True)

# Open output file
f = open('compareFBA_results.csv', 'w')
f.write('FBA,secretoryFBA\n')

for entryID in list(psim['Entry']):
    # Generate lists for reactions (rxns), reaction names (rxnNames), and Gene-Protein-Reactions (GPRs)
    [rxns, rxnNames, GPRs] = iCHO2048s_Builder.generateProteinSpecificRxns_A(entryID)
    # Initialize new model
    secretory_model = model.copy()
    # Add new reactions to metabolic model
    for i in range(len(rxns)):
        # Create place holder for creating reactions
        r = secretory_model.reactions[0].copy()
        r.name = rxnNames[i]
        r.id = rxnNames[i]
        r.gene_reaction_rule = GPRs[i]
        secretory_model.add_reaction(r)
        r.build_reaction_from_string(reaction_str=rxns[i],rev_arrow='<-',fwd_arrow='->',reversible_arrow='<=>')
    
    # Add regular-FBA demand reaction
    metaboliteStr = entryID + '[c]'
    r = secretory_model.reactions[0].copy()
    r.name = metaboliteStr + '_demand'
    r.id = metaboliteStr + '_demand'
    r.gene_reaction_rule = ''
    secretory_model.add_reaction(r)
    r.build_reaction_from_string(reaction_str=metaboliteStr+' -> ',rev_arrow='<-',fwd_arrow='->',reversible_arrow='<=>')
    
    # Set objective to regular-FBA demand reaction
    secretory_model.objective = metaboliteStr + '_demand'
    f.write(str(secretory_model.slim_optimize())+',')
    # Set the objective to secretion of target protein
    secretory_model.objective = entryID + '_Final_demand'
    f.write(str(secretory_model.slim_optimize())+'\n')
f.close()
