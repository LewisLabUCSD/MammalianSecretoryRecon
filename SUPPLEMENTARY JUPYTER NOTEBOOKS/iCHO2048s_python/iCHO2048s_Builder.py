"""
This script contains functions to create a list of protein-specific reactions
for the protein secretion pathway in iCHO1921s.
"""
# Read files with appropriate template iCHO1921s reactions and remove "\n"
f = open("rxnFormula_CHO.txt")
rxnFormula = f.read().splitlines()
f.close()

f = open("rxnAbbreviation_CHO.txt")
rxnAbbreviation = f.read().splitlines()
f.close()

f = open("rxnPathway_CHO.txt")
rxnPathway = f.read().splitlines()
f.close()

f = open("rxnConditions_CHO.txt")
rxnConditions = f.read().splitlines()
f.close()

f = open("rxnGPRs_CHO.txt")
rxnGPR = f.read().splitlines()
f.close()

# Load PSIM matrix
f = open('PSIM_CHO.tab','r')
PSIM=f.read().splitlines()
f.close()

#Extract the Uniprot IDs
PSIM_entries = []
for line in PSIM:
    PSIM_entries.append(line.split('\t')[0])
    
# Define Basic Functions

#Define a functions that counts of amino acids
def count_AAs(sequence):
    AAcounts = []    
    AAs = ['G', 'A', 'V', 'L', 'I', 'M', 'W', 'F', 'P', 'S', 'T', 'C', 'Y', 'N', 'Q', 'E', 'D', 'K', 'R', 'H']
    for aa in AAs:
        AAcounts.append(sequence.count(aa))
    return AAcounts

#Define a function that substitutes question marks for AA counts in a formula (string)
#Useful for translation and protein degradation pathway

def substitute_AAs_count(formula,AAcounts):
    template = "? gly[c] + ? ala_L[c] + ? val_L[c] + ? leu_L[c] + ? ile_L[c] + ? met_L[c] + ? trp_L[c] + ? phe_L[c] + ? pro_L[c] + ? ser_L[c] + ? thr_L[c] + ? cys_L[c] + ? tyr_L[c] + ? asn_L[c] + ? gln_L[c] + ? glu_L[c] + ? asp_L[c] + ? lys_L[c] + ? arg_L[c] + ? his_L[c]"
    new = ''
    index = 0
    for character in template:
        if character == "?":
            new = new + str(AAcounts[index])
            index = index + 1
        else:
            new = new + character
    newFormula = formula.replace(template,new)
    return newFormula
    
# Function that generates the translation reaction of a protein given its UniProt ID
def translate_protein(entryID):
    #Obtain protein sequence and length
    PSI_row = PSIM[PSIM_entries.index(str(entryID))] 
    PSI_row = PSI_row.split('\t')
    sequence = PSI_row[11]
    AAcounts = count_AAs(sequence)
    N = len(sequence) #Number to replace atp's, gtp's, ppi's, pi's, h2o's, amp's, and gpd's
    templateFormula = "? h2o[c] + ? atp[c] + ? gtp[c] + ? gly[c] + ? ala_L[c] + ? val_L[c] + ? leu_L[c] + ? ile_L[c] + ? met_L[c] + ? trp_L[c] + ? phe_L[c] + ? pro_L[c] + ? ser_L[c] + ? thr_L[c] + ? cys_L[c] + ? tyr_L[c] + ? asn_L[c] + ? gln_L[c] + ? glu_L[c] + ? asp_L[c] + ? lys_L[c] + ? arg_L[c] + ? his_L[c] -> ? h[c] + ? amp[c] + adp[c] + ? pi[c] + ? gdp[c] + ? ppi[c] + XXX[c]"
    translationFormula = substitute_AAs_count(templateFormula, AAcounts)
    translationFormula = translationFormula.replace("? h2o[c]",str(2*N-1)+" h2o[c]")
    translationFormula = translationFormula.replace("? h[c]",str(2*N-1)+" h[c]")
    translationFormula = translationFormula.replace("? ppi[c]",str(N)+" ppi[c]")
    translationFormula = translationFormula.replace("? pi[c]",str(2*N-1)+" pi[c]")
    translationFormula = translationFormula.replace("? gdp[c]",str(2*N-2)+" gdp[c]")
    translationFormula = translationFormula.replace("? amp[c]",str(N)+" amp[c]")
    translationFormula = translationFormula.replace("? gtp[c]",str(2*N-2)+" gtp[c]")
    translationFormula = translationFormula.replace("? atp[c]",str(N+1)+" atp[c]")
    translationFormula = translationFormula.replace("XXX[c]",str(entryID)+"[c]")
    return translationFormula

# Function that replaces the XXX template name for the UniProt ID
def insert_prot_name_in_rxnFormula(formula,entryID):
    newFormula = formula.replace("XXX",str(entryID))
    return newFormula

# Function that replaces the XXX template name for the Reaction abbreviation
def insert_prot_name_in_rxnName(rxnAbbrev,entryID):
    newAbbreviation = str(entryID) + "_" + rxnAbbrev
    return newAbbreviation

#Function for adding the reactions of a given PathwayName
def addPathway(pathwayName,listOfRxns,listOfRxnsNames):   
    newList = listOfRxns    
    newList2 = listOfRxnsNames
    for i in range(len(rxnPathway)):
        if rxnPathway[i] == pathwayName:
            newList.append(rxnFormula[i])
            newList2.append(rxnAbbreviation[i])
    return newList,newList2

#Function for adding the reactions of a given PathwayName given a condition
def addPathwayFromCondition(conditionName,listOfRxns,listOfRxnsNames):   
    newList = listOfRxns 
    newList2 = listOfRxnsNames
    for i in range(len(rxnConditions)):
        if rxnConditions[i] == conditionName:
            newList.append(rxnFormula[i])
            newList2.append(rxnAbbreviation[i])
    return newList,newList2

#Function for creating a list of GPRs given a list of Rxn names
def getGPRsFromRxnNames(listOfRxnsNames):
    GPR_list = []
    for reaction in listOfRxnsNames:
        if reaction in rxnAbbreviation:
            GPR_list.append(rxnGPR[rxnAbbreviation.index(reaction)])
        else:
            GPR_list.append('')
    return GPR_list
     
#Function that adds canonical reactions to overall list
def addCanonicalRxns(listOfRxns,listOfRxnNames,listOfGPRs):
    newListRxns = listOfRxns
    newListNames = listOfRxnNames
    newListGPRs = listOfGPRs
    #Add  canonical reactions
    [newListRxns, newListNames] = addPathway("Canonical",listOfRxns,listOfRxnNames)
    #Add canonical GPRs
    r = []
    n = []
    [r,n] = addPathway("Canonical",r,n)
    newGPRs = getGPRsFromRxnNames(n)
    for gpr in newGPRs:
        newListGPRs.append(gpr)
    return newListRxns, newListNames, newListGPRs


	
# Function that generates protein-specific reactions list given a string of a Uniprot ID (use it for CHO proteins)
def generateProteinSpecificRxns_A(entryID):
    #Obtain protein sequence and length
    PSI_row = PSIM[PSIM_entries.index(str(entryID))] 
    PSI_row = PSI_row.split('\t')
    sequence = PSI_row[11]
    L = float(PSI_row[2]) # Protein Length
    MW = float(PSI_row[3]) # Molecular weight
    PSI = []
    Kv = 0.7
    V = MW * 1.21 / 1000.0 # Protein Volume in nm^3
    clathrin_coeff = int(round(29880.01 * Kv / V)) # Number of proteins per clathrin vesicle  
    copi_coeff = int(round(143793.19 * Kv / V))
    copii_coeff = int(round(268082.35 * Kv / V))
    connector = ''
    #Prepare vectors that will store reactions and components
    protName = str(entryID)
    rxns = []
    rxnNames = []   

    #Add translation reaction
    translation_reaction = translate_protein(protName)
    rxns.append(translation_reaction)
    rxnNames.append("TRANSLATION_protein") 
    for i in [0,4,5,6,7,8,9,10]:
        PSI.append(PSI_row[i])  #This is the vector from the PSIM that corresponds to the given protein 
										  #[entry,SP,DSB,GPI,NG,OG,TMD,SubCellLoc] 
   
#------------------------------------------------------------------------------ 
    #Add translocation reactions PSI[1] = "Signal Peptide"
    
    if PSI[1] == '0': #If it doesn't have signal peptide then ignore
        GPRs = getGPRsFromRxnNames(rxnNames)
        #Change the reaction names and their formulas to include the UniProt Identifier
        for i in range(len(rxns)):
            rxns[i] = insert_prot_name_in_rxnFormula(rxns[i],protName)
        for i in range(len(rxnNames)):
            rxnNames[i] = insert_prot_name_in_rxnName(rxnNames[i],protName)
        rxns.append(protName + '[c] -> ')
        rxnNames.append(protName + '_Final_demand')
        GPRs.append('')
        return rxns,rxnNames,GPRs
        
    elif PSI[1] == '1': #Translocate protein if it has signal peptide
        if L <= 160:
            [rxns,rxnNames] = addPathway("Post-translational Translocation",rxns,rxnNames)
            if PSI[7] == '[e]' or PSI[7] == '':
                [rxns,rxnNames] = addPathway("Post-translational Translocation (Secretory protein)",rxns,rxnNames)
            if PSI[7] == "[pm]" or PSI[6] != '0':
                [rxns,rxnNames] = addPathway("Post-translational Translocation (Tail anchored membrane protein)",rxns,rxnNames)
        else:
            [rxns,rxnNames] = addPathway("Translocation",rxns,rxnNames)
        
        number_BiP = L/40 #Number of BiPs depends on protein length http://www.cshperspectives.com/content/5/5/a013201.full
        for i in range(len(rxns)):
            rxns[i] = rxns[i].replace("!",str(number_BiP))
        connector = 'XXX[r]'
        
#------------------------------------------------------------------------------
    #Add Disulphide Bond Reactions
    if PSI[2] != '0':
        number_DSB = PSI[2]
        DSBrxns = []
        DSBrxnNames = []
        DSBrxns.append(connector + ' -> XXX_preDSB[r]')
        DSBrxnNames.append('Start_DSB')
        [DSBrxns, DSBrxnNames] = addPathwayFromCondition('DSB>0',DSBrxns,DSBrxnNames)
        for i in range(len(DSBrxns)):
            if '?' in DSBrxns[i]:
                if number_DSB == '1':
                    DSBrxns[i] = DSBrxns[i].replace('? ','')
                else:
                    DSBrxns[i] = DSBrxns[i].replace('?',number_DSB)
        
        connector = 'XXX_DSB[r]'

        for reaction in DSBrxns:
                rxns.append(reaction)
        for reactionName in DSBrxnNames:
                rxnNames.append(reactionName)
#------------------------------------------------------------------------------              
    #Add GPI reactions
    if PSI[3] == '1':
        rxns.append(connector + ' -> XXX_preGPI[r]')
        rxnNames.append('Start_GPI')
        [rxns,rxnNames] = addPathwayFromCondition('GPI=1',rxns,rxnNames)
        connector = 'XXX-dgpi_cho[r]'
#------------------------------------------------------------------------------    
    #Add N-glycosylation reactions
    if PSI[4] != '0': #Has N-glycans?
        rxns.append(connector + ' -> XXX_preNG[r]')
        rxnNames.append('Start_NG')
        number_Nglycans = PSI[4] #Get number of N-Glycans
        NGlyrxns = []
        NGlyrxnNames = []
        
        [NGlyrxns,NGlyrxnNames] = addPathwayFromCondition('NG>0',NGlyrxns,NGlyrxnNames)
        [NGlyrxns,NGlyrxnNames] = addPathway('Golgi processing N',NGlyrxns,NGlyrxnNames)
        for i in range(len(NGlyrxns)): #Change the '?' for the number of N-glycans
            if NGlyrxns[i] == 'XXX-M5-unfold-UBIQP[c] + ? h2o[c] + RAD23A[c] =>  XXX-unfold-UBIQP-RAD23A[c] + ? acgam[c] + ? man[c]':
                h2o = str(7*int(number_Nglycans))
                acgam = str(2*int(number_Nglycans))
                man = str(5*int(number_Nglycans))
                NGlyrxns[i] = 'XXX-M5-unfold-UBIQP[c] + ' + h2o +' h2o[c] + RAD23A[c] =>  XXX-unfold-UBIQP[c]-RAD23A[c] + ' + acgam + ' acgam[c] + ' + man +' man[c]'
            if '?' in NGlyrxns[i]:
                    if number_Nglycans == '1':
                        NGlyrxns[i] = NGlyrxns[i].replace('? ', '')
                    else:
                        NGlyrxns[i] = NGlyrxns[i].replace('?', number_Nglycans)

        for reaction in NGlyrxns:
            rxns.append(reaction)
        for reactionName in NGlyrxnNames:
            rxnNames.append(reactionName)
        
        connector = 'XXX-M8B[r]'
   
#------------------------------------------------------------------------------
    #Add  COPII reactions
    if connector == 'XXX-M8B[r]':
        copii_rxns = []
        copii_names = []
        [copii_rxns, copii_names] =  addPathway('COPII_NG',copii_rxns,copii_names)       
        for i in range(len(copii_rxns)):
            copii_rxns[i] = copii_rxns[i].replace("!",str(copii_coeff))
            rxns.append(copii_rxns[i])
            rxnNames.append(copii_names[i])
        
        #[rxns,rxnNames] = addPathway('COPII_NG',rxns,rxnNames)
        connector = 'XXX-M3-GN2[g]'
        
        #-------------- ADD IF PROTEIN = EPO (HUMAN)--------------------------#
        # IMPORTANT! -> Reactions in this pathway assume EPO has NG=3 and Og=1 (see PSIM)
        if protName == 'P01588':
            rxns.append('P01588[c] -> EPO_Human[c]')
            rxnNames.append('make_EPO')
            [rxns,rxnNames] = addPathway('Golgi processing (EPO specific)',rxns,rxnNames)
            connector = 'XXX-M3-GN4-GL4-NA4-F[g]'
            protName = 'EPO_Human'
        
    elif connector == 'XXX-dgpi_cho[r]':
        copii_rxns = []
        copii_names = []
        [copii_rxns, copii_names] =  addPathway('COPII_GPI',copii_rxns,copii_names)       
        for i in range(len(copii_rxns)):
            copii_rxns[i] = copii_rxns[i].replace("!",str(copii_coeff))
            rxns.append(copii_rxns[i])
            rxnNames.append(copii_names[i])
        
        #[rxns,rxnNames] = addPathway('COPII_GPI',rxns,rxnNames)
        connector = 'XXX-dgpi_cho[g]'
        
    elif connector == 'XXX_DSB[r]':
        [rxns,rxnNames] = addPathway('COPII_DSB',rxns,rxnNames)
        
        copii_rxns = []
        copii_names = []
        [copii_rxns, copii_names] =  addPathway('COPII-canonical',copii_rxns,copii_names)       
        for i in range(len(copii_rxns)):
            copii_rxns[i] = copii_rxns[i].replace("!",str(copii_coeff))
            rxns.append(copii_rxns[i])
            rxnNames.append(copii_names[i])
        
        #[rxns,rxnNames] = addPathway('COPII-canonical',rxns,rxnNames)
        connector = 'XXX[g]'
    
    elif connector == 'XXX[r]':
        [rxns,rxnNames] = addPathway('COPII-normal',rxns,rxnNames)
        
        copii_rxns = []
        copii_names = []
        [copii_rxns, copii_names] =  addPathway('COPII-canonical',copii_rxns,copii_names)       
        for i in range(len(copii_rxns)):
            copii_rxns[i] = copii_rxns[i].replace("!",str(copii_coeff))
            rxns.append(copii_rxns[i])
            rxnNames.append(copii_names[i])        
        
        #[rxns,rxnNames] = addPathway('COPII-canonical',rxns,rxnNames)
        connector = 'XXX[g]'

#------------------------------------------------------------------------------
    #Add O-glycosylation reactions
    if PSI[5] != '0': #Has O-glycans?    
        rxns.append(connector + ' -> XXX_preOG[g]')
        rxnNames.append('Start_OG')
        number_Oglycans = PSI[5] #Get number of O-Glycans
        OGlyrxns = []
        OGlyrxnNames = []
        
        [OGlyrxns,OGlyrxnNames] = addPathwayFromCondition('OG>0',OGlyrxns,OGlyrxnNames)
        for i in range(len(OGlyrxns)): #Change the '?' for the number of O-glycans
            if '?' in OGlyrxns[i]:
                    if number_Oglycans == '1':
                        OGlyrxns[i] = OGlyrxns[i].replace('? ', '')
                    else:
                        OGlyrxns[i] = OGlyrxns[i].replace('?', number_Oglycans)

        for reaction in OGlyrxns:
            rxns.append(reaction)
        for reactionName in OGlyrxnNames:
            rxnNames.append(reactionName)
        
        connector = 'XXX-Core2[g]'
        
#------------------------------------------------------------------------------
    #Stay in Golgi if protein is localized there
    if PSI[7] == '[g]' or PSI[7] == '[gm]':
        location = PSI[7]
        #rxns.append(connector + ' -> XXX_mature' + location)
        #rxnNames.append('Final_location_' + location)
        if 'SP_degradation' in rxnNames:
            SPaas = count_AAs(sequence[0:22]) #Amino acids in signal peptide assuming length is 22 on average
            rxns[rxnNames.index('SP_degradation')] = substitute_AAs_count(rxns[rxnNames.index('SP_degradation')],SPaas)
        if 'Ubiquitination_degradation' in rxnNames:
            new_aas = count_AAs(sequence[22:])
            rxns[rxnNames.index('Ubiquitination_degradation')] = rxns[rxnNames.index('Ubiquitination_degradation')].replace("!", "?")
            rxns[rxnNames.index('Ubiquitination_degradation')] = substitute_AAs_count(rxns[rxnNames.index('Ubiquitination_degradation')], new_aas)
            rxns[rxnNames.index('Ubiquitination_degradation')] = rxns[rxnNames.index('Ubiquitination_degradation')].replace("?",str(L-22))        
        # Add GPRs
        GPRs = getGPRsFromRxnNames(rxnNames)
        #Change the reaction names and their formulas to include the UniProt Identifier
        for i in range(len(rxns)):
            rxns[i] = insert_prot_name_in_rxnFormula(rxns[i],protName)
        for i in range(len(rxnNames)):
            rxnNames[i] = insert_prot_name_in_rxnName(rxnNames[i],protName)
        rxns.append(protName + location + ' -> ')
        rxnNames.append(protName + '_Final_demand')
        GPRs.append('')    
        [rxns,rxnNames,GPRs] = addCanonicalRxns(rxns,rxnNames,GPRs)  
        return rxns,rxnNames,GPRs
#------------------------------------------------------------------------------    
    #Add COPI
    elif PSI[7] == '[r]' or PSI[7] == '[rm]':
        rxns.append(connector + ' -> XXX_preCOPI[g]')
        rxnNames.append('Start_COPI')
        
        copi_rxns = []
        copi_names = []
        [copi_rxns, copi_names] =  addPathway('COPI',copi_rxns,copi_names)       
        for i in range(len(copi_rxns)):
            copi_rxns[i] = copi_rxns[i].replace("!",str(copi_coeff))
            rxns.append(copi_rxns[i])
            rxnNames.append(copi_names[i])

        connector = 'XXX_mature[r]'
        location = PSI[7]
        if location == '[r]':
            rxns.append(connector + ' -> ')
            rxnNames.append(protName + '_Final_demand')
        elif location == '[rm]':
            rxns.append(connector + ' -> XXX_mature' + location)
            rxnNames.append('Final_location_' + location)
            rxns.append('XXX_mature' + location + ' -> ')
            rxnNames.append(protName + '_Final_demand')
        
#------------------------------------------------------------------------------
    #Add Clathrin vesicles
    elif PSI[7]=='[x]' or PSI[7]=='[l]' or PSI[7]=='[d]':
        rxns.append(connector + ' -> XXX-preClathrin[g]')
        rxnNames.append('Start_Clathrin_vesicle')
        
        clath_rxns = []
        clath_names = []
        [clath_rxns, clath_names] =  addPathway('Clathrin vesicles',clath_rxns,clath_names)       
        for i in range(len(clath_rxns)):
            clath_rxns[i] = clath_rxns[i].replace("!",str(clathrin_coeff))
            rxns.append(clath_rxns[i])
            rxnNames.append(clath_names[i])        
        
        connector = 'XXX_mature[cv]'
        location = PSI[7]
        rxns.append(connector + ' -> XXX_mature' + location)
        rxnNames.append('Final_location_' + location)
        rxns.append('XXX_mature' + location + ' -> ')
        rxnNames.append(protName + '_Final_demand')
        
#------------------------------------------------------------------------------
    #Send to corresponding location
    else:
        location = PSI[7]
        rxns.append(connector + ' -> XXX-preSV[g]')
        rxnNames.append('Start_Secretion')
        
        sv_rxns = []
        sv_names = []
        [sv_rxns, sv_names] =  addPathway('SV',sv_rxns, sv_names)       
        for i in range(len(sv_rxns)):
            sv_rxns[i] = sv_rxns[i].replace("!",str(clathrin_coeff))
            rxns.append(sv_rxns[i])
            rxnNames.append(sv_names[i])        
        
        #[rxns,rxnNames] = addPathway("SV",rxns,rxnNames)
        if location == '':
            location = '[e]'
        rxns.append('XXX_mature[sv]' + ' -> XXX_mature' + location)
        rxnNames.append('Final_location_' + location)
        rxns.append('XXX_mature' + location + ' -> ')
        rxnNames.append('Final_demand')
#------------------------------------------------------------------------------
    #Add coeffcients to SP_degradation and Ubiquitin_degradation reactions (if applicable)
    if 'SP_degradation' in rxnNames:
        SPaas = count_AAs(sequence[0:22]) #Amino acids in signal peptide assuming length is 22 on average
        rxns[rxnNames.index('SP_degradation')] = substitute_AAs_count(rxns[rxnNames.index('SP_degradation')],SPaas)
    if 'Ubiquitination_degradation' in rxnNames:
        new_aas = count_AAs(sequence[22:])
        rxns[rxnNames.index('Ubiquitination_degradation')] = rxns[rxnNames.index('Ubiquitination_degradation')].replace("!", "?")
        rxns[rxnNames.index('Ubiquitination_degradation')] = substitute_AAs_count(rxns[rxnNames.index('Ubiquitination_degradation')], new_aas)
        rxns[rxnNames.index('Ubiquitination_degradation')] = rxns[rxnNames.index('Ubiquitination_degradation')].replace("?",str(L-22))
#------------------------------------------------------------------------------
    GPRs = getGPRsFromRxnNames(rxnNames)
#------------------------------------------------------------------------------
    #Change the reaction names and their formulas to include the UniProt Identifier
    for i in range(len(rxns)):
        rxns[i] = insert_prot_name_in_rxnFormula(rxns[i],protName)
    for i in range(len(rxnNames)):
        rxnNames[i] = insert_prot_name_in_rxnName(rxnNames[i],protName)
    [rxns,rxnNames,GPRs] = addCanonicalRxns(rxns,rxnNames,GPRs)  
    return rxns,rxnNames,GPRs
  
  
  
  
#------------------------------------------------------------------------------
#PSI = ['IgG','1','17','0','2','0','0','[e]']
#PSI_row = [Name, Sequence, MW, SP, DSB, GPI, NG, OG, Location]
def generateProteinSpecificRxns_B(PSI_row): # Use for non-CHO proteins
    PSI = [PSI_row[0], PSI_row[3],PSI_row[4],PSI_row[5],PSI_row[6],PSI_row[7],'0',PSI_row[8]]
    Kv = 0.7
    connector = ''
    sequence = PSI_row[1]
    L = len(sequence)
    protName = PSI[0]
    MW = PSI_row[2]
    V = MW * 1.21 / 1000.0 # Protein Volume in nm^3
    clathrin_coeff = int(round(29880.01 * Kv / V)) # Number of proteins per clathrin vesicle  
    copi_coeff = int(round(143793.19 * Kv / V))
    copii_coeff = int(round(268082.35 * Kv / V))
    
    #Prepare vectors that will store reactions and components
    rxns = []
    rxnNames = []

    #Add translation reaction
    AAcounts = count_AAs(sequence)
    N = len(sequence) #Number to replace atp's, gtp's, ppi's, pi's, h2o's, amp's, and gpd's
    templateFormula = "? h2o[c] + ? atp[c] + ? gtp[c] + ? gly[c] + ? ala_L[c] + ? val_L[c] + ? leu_L[c] + ? ile_L[c] + ? met_L[c] + ? trp_L[c] + ? phe_L[c] + ? pro_L[c] + ? ser_L[c] + ? thr_L[c] + ? cys_L[c] + ? tyr_L[c] + ? asn_L[c] + ? gln_L[c] + ? glu_L[c] + ? asp_L[c] + ? lys_L[c] + ? arg_L[c] + ? his_L[c] -> ? h[c] + ? amp[c] + adp[c] + ? pi[c] + ? gdp[c] + ? ppi[c] + XXX[c]"
    translationFormula = substitute_AAs_count(templateFormula, AAcounts)
    translationFormula = translationFormula.replace("? h2o[c]",str(2*N-1)+" h2o[c]")
    translationFormula = translationFormula.replace("? h[c]",str(2*N-1)+" h[c]")
    translationFormula = translationFormula.replace("? ppi[c]",str(N)+" ppi[c]")
    translationFormula = translationFormula.replace("? pi[c]",str(2*N-1)+" pi[c]")
    translationFormula = translationFormula.replace("? gdp[c]",str(2*N-2)+" gdp[c]")
    translationFormula = translationFormula.replace("? amp[c]",str(N)+" amp[c]")
    translationFormula = translationFormula.replace("? gtp[c]",str(2*N-2)+" gtp[c]")
    translationFormula = translationFormula.replace("? atp[c]",str(N+1)+" atp[c]")
    translationFormula = translationFormula.replace("XXX[c]",str(PSI[0])+"[c]")
    rxns.append(translationFormula)
    rxnNames.append("TRANSLATION_protein")    
#------------------------------------------------------------------------------ 
    #Add translocation reactions PSI[1] = "Signal Peptide"
    
    if PSI[1] == '0': #If it doesn't have signal peptide then ignore
        GPRs = getGPRsFromRxnNames(rxnNames)
        #Change the reaction names and their formulas to include the UniProt Identifier
        for i in range(len(rxns)):
            rxns[i] = insert_prot_name_in_rxnFormula(rxns[i],protName)
        for i in range(len(rxnNames)):
            rxnNames[i] = insert_prot_name_in_rxnName(rxnNames[i],protName)
        rxns.append(protName + '[c] -> ')
        rxnNames.append(protName + '_Final_demand')
        GPRs.append('')
        return rxns,rxnNames,GPRs
        
    elif PSI[1] == '1': #Translocate protein if it has signal peptide
        if L <= 160:
            [rxns,rxnNames] = addPathway("Post-translational Translocation",rxns,rxnNames)
            if PSI[7] == '[e]' or PSI[7] == '':
                [rxns,rxnNames] = addPathway("Post-translational Translocation (Secretory protein)",rxns,rxnNames)
            if PSI[7] == "[pm]" or PSI[6] != '0':
                [rxns,rxnNames] = addPathway("Post-translational Translocation (Tail anchored membrane protein)",rxns,rxnNames)
        else:
            [rxns,rxnNames] = addPathway("Translocation",rxns,rxnNames)
        
        number_BiP = L/40 #Number of BiPs depends on protein length http://www.cshperspectives.com/content/5/5/a013201.full
        for i in range(len(rxns)):
            rxns[i] = rxns[i].replace("!",str(number_BiP))
        connector = 'XXX[r]'
        
#------------------------------------------------------------------------------
    #Add Disulphide Bond Reactions
    if PSI[2] != '0':
        number_DSB = PSI[2]
        DSBrxns = []
        DSBrxnNames = []
        DSBrxns.append(connector + ' -> XXX_preDSB[r]')
        DSBrxnNames.append('Start_DSB')
        [DSBrxns, DSBrxnNames] = addPathwayFromCondition('DSB>0',DSBrxns,DSBrxnNames)
        for i in range(len(DSBrxns)):
            if '?' in DSBrxns[i]:
                if number_DSB == '1':
                    DSBrxns[i] = DSBrxns[i].replace('? ','')
                else:
                    DSBrxns[i] = DSBrxns[i].replace('?',number_DSB)
        
        connector = 'XXX_DSB[r]'

        for reaction in DSBrxns:
                rxns.append(reaction)
        for reactionName in DSBrxnNames:
                rxnNames.append(reactionName)
#------------------------------------------------------------------------------              
    #Add GPI reactions
    if PSI[3] == '1':
        rxns.append(connector + ' -> XXX_preGPI[r]')
        rxnNames.append('Start_GPI')
        [rxns,rxnNames] = addPathwayFromCondition('GPI=1',rxns,rxnNames)
        connector = 'XXX-dgpi_cho[r]'
#------------------------------------------------------------------------------    
    #Add N-glycosylation reactions
    if PSI[4] != '0': #Has N-glycans?
        rxns.append(connector + ' -> XXX_preNG[r]')
        rxnNames.append('Start_NG')
        number_Nglycans = PSI[4] #Get number of N-Glycans
        NGlyrxns = []
        NGlyrxnNames = []
        
        [NGlyrxns,NGlyrxnNames] = addPathwayFromCondition('NG>0',NGlyrxns,NGlyrxnNames)
        [NGlyrxns,NGlyrxnNames] = addPathway('Golgi processing N',NGlyrxns,NGlyrxnNames)
        for i in range(len(NGlyrxns)): #Change the '?' for the number of N-glycans
            if NGlyrxns[i] == 'XXX-M5-unfold-UBIQP[c] + ? h2o[c] + RAD23A[c] =>  XXX-unfold-UBIQP-RAD23A[c] + ? acgam[c] + ? man[c]':
                h2o = str(7*int(number_Nglycans))
                acgam = str(2*int(number_Nglycans))
                man = str(5*int(number_Nglycans))
                NGlyrxns[i] = 'XXX-M5-unfold-UBIQP[c] + ' + h2o +' h2o[c] + RAD23A[c] =>  XXX-unfold-UBIQP[c]-RAD23A[c] + ' + acgam + ' acgam[c] + ' + man +' man[c]'
            if '?' in NGlyrxns[i]:
                    if number_Nglycans == '1':
                        NGlyrxns[i] = NGlyrxns[i].replace('? ', '')
                    else:
                        NGlyrxns[i] = NGlyrxns[i].replace('?', number_Nglycans)

        for reaction in NGlyrxns:
            rxns.append(reaction)
        for reactionName in NGlyrxnNames:
            rxnNames.append(reactionName)
        
        connector = 'XXX-M8B[r]'
   
#------------------------------------------------------------------------------
    #Add  COPII reactions
    if connector == 'XXX-M8B[r]':
        copii_rxns = []
        copii_names = []
        [copii_rxns, copii_names] =  addPathway('COPII_NG',copii_rxns,copii_names)       
        for i in range(len(copii_rxns)):
            copii_rxns[i] = copii_rxns[i].replace("!",str(copii_coeff))
            rxns.append(copii_rxns[i])
            rxnNames.append(copii_names[i])
        
        #[rxns,rxnNames] = addPathway('COPII_NG',rxns,rxnNames)
        connector = 'XXX-M3-GN2[g]'
        
    elif connector == 'XXX-dgpi_cho[r]':
        copii_rxns = []
        copii_names = []
        [copii_rxns, copii_names] =  addPathway('COPII_GPI',copii_rxns,copii_names)       
        for i in range(len(copii_rxns)):
            copii_rxns[i] = copii_rxns[i].replace("!",str(copii_coeff))
            rxns.append(copii_rxns[i])
            rxnNames.append(copii_names[i])
        
        #[rxns,rxnNames] = addPathway('COPII_GPI',rxns,rxnNames)
        connector = 'XXX-dgpi_cho[g]'
        
    elif connector == 'XXX_DSB[r]':
        [rxns,rxnNames] = addPathway('COPII_DSB',rxns,rxnNames)
        
        copii_rxns = []
        copii_names = []
        [copii_rxns, copii_names] =  addPathway('COPII-canonical',copii_rxns,copii_names)       
        for i in range(len(copii_rxns)):
            copii_rxns[i] = copii_rxns[i].replace("!",str(copii_coeff))
            rxns.append(copii_rxns[i])
            rxnNames.append(copii_names[i])
        
        #[rxns,rxnNames] = addPathway('COPII-canonical',rxns,rxnNames)
        connector = 'XXX[g]'
    
    elif connector == 'XXX[r]':
        [rxns,rxnNames] = addPathway('COPII-normal',rxns,rxnNames)
        
        copii_rxns = []
        copii_names = []
        [copii_rxns, copii_names] =  addPathway('COPII-canonical',copii_rxns,copii_names)       
        for i in range(len(copii_rxns)):
            copii_rxns[i] = copii_rxns[i].replace("!",str(copii_coeff))
            rxns.append(copii_rxns[i])
            rxnNames.append(copii_names[i])        
        
        #[rxns,rxnNames] = addPathway('COPII-canonical',rxns,rxnNames)
        connector = 'XXX[g]'

#------------------------------------------------------------------------------
    #Add O-glycosylation reactions
    if PSI[5] != '0': #Has O-glycans?    
        rxns.append(connector + ' -> XXX_preOG[g]')
        rxnNames.append('Start_OG')
        number_Oglycans = PSI[5] #Get number of O-Glycans
        OGlyrxns = []
        OGlyrxnNames = []
        
        [OGlyrxns,OGlyrxnNames] = addPathwayFromCondition('OG>0',OGlyrxns,OGlyrxnNames)
        for i in range(len(OGlyrxns)): #Change the '?' for the number of O-glycans
            if '?' in OGlyrxns[i]:
                    if number_Oglycans == '1':
                        OGlyrxns[i] = OGlyrxns[i].replace('? ', '')
                    else:
                        OGlyrxns[i] = OGlyrxns[i].replace('?', number_Oglycans)

        for reaction in OGlyrxns:
            rxns.append(reaction)
        for reactionName in OGlyrxnNames:
            rxnNames.append(reactionName)
        
        connector = 'XXX-Core2[g]'
        
#------------------------------------------------------------------------------
    #Stay in Golgi if protein is localized there
    if PSI[7] == '[g]' or PSI[7] == '[gm]':
        location = PSI[7]
        #rxns.append(connector + ' -> XXX_mature' + location)
        #rxnNames.append('Final_location_' + location)
        if 'SP_degradation' in rxnNames:
            SPaas = count_AAs(sequence[0:22]) #Amino acids in signal peptide assuming length is 22 on average
            rxns[rxnNames.index('SP_degradation')] = substitute_AAs_count(rxns[rxnNames.index('SP_degradation')],SPaas)
        if 'Ubiquitination_degradation' in rxnNames:
            new_aas = count_AAs(sequence[22:])
            rxns[rxnNames.index('Ubiquitination_degradation')] = rxns[rxnNames.index('Ubiquitination_degradation')].replace("!", "?")
            rxns[rxnNames.index('Ubiquitination_degradation')] = substitute_AAs_count(rxns[rxnNames.index('Ubiquitination_degradation')], new_aas)
            rxns[rxnNames.index('Ubiquitination_degradation')] = rxns[rxnNames.index('Ubiquitination_degradation')].replace("?",str(L-22))        
        # Add GPRs
        GPRs = getGPRsFromRxnNames(rxnNames)
        #Change the reaction names and their formulas to include the UniProt Identifier
        for i in range(len(rxns)):
            rxns[i] = insert_prot_name_in_rxnFormula(rxns[i],PSI[0])
        for i in range(len(rxnNames)):
            rxnNames[i] = insert_prot_name_in_rxnName(rxnNames[i],PSI[0])
        rxns.append(protName + location + ' -> ')
        rxnNames.append(protName + '_Final_demand')
        GPRs.append('')
        [rxns,rxnNames,GPRs] = addCanonicalRxns(rxns,rxnNames,GPRs)  		
        return rxns,rxnNames,GPRs
#------------------------------------------------------------------------------    
    #Add COPI
    elif PSI[7] == '[r]' or PSI[7] == '[rm]':
        rxns.append(connector + ' -> XXX_preCOPI[g]')
        rxnNames.append('Start_COPI')
        
        copi_rxns = []
        copi_names = []
        [copi_rxns, copi_names] =  addPathway('COPI',copi_rxns,copi_names)       
        for i in range(len(copi_rxns)):
            copi_rxns[i] = copi_rxns[i].replace("!",str(copi_coeff))
            rxns.append(copi_rxns[i])
            rxnNames.append(copi_names[i])

        connector = 'XXX_mature[r]'
        location = PSI[7]
        if location == '[r]':
            rxns.append(connector + ' -> ')
            rxnNames.append(protName + '_Final_demand')
        elif location == '[rm]':
            rxns.append(connector + ' -> XXX_mature' + location)
            rxnNames.append('Final_location_' + location)
            rxns.append('XXX_mature' + location + ' -> ')
            rxnNames.append(protName + '_Final_demand')
        
#------------------------------------------------------------------------------
    #Add Clathrin vesicles
    elif PSI[7]=='[x]' or PSI[7]=='[l]' or PSI[7]=='[d]':
        rxns.append(connector + ' -> XXX-preClathrin[g]')
        rxnNames.append('Start_Clathrin_vesicle')
        
        clath_rxns = []
        clath_names = []
        [clath_rxns, clath_names] =  addPathway('Clathrin vesicles',clath_rxns,clath_names)       
        for i in range(len(clath_rxns)):
            clath_rxns[i] = clath_rxns[i].replace("!",str(clathrin_coeff))
            rxns.append(clath_rxns[i])
            rxnNames.append(clath_names[i])        
        
        connector = 'XXX_mature[cv]'
        location = PSI[7]
        rxns.append(connector + ' -> XXX_mature' + location)
        rxnNames.append('Final_location_' + location)
        rxns.append('XXX_mature' + location + ' -> ')
        rxnNames.append(protName + '_Final_demand')
        
#------------------------------------------------------------------------------
    #Send to corresponding location
    else:
        location = PSI[7]
        rxns.append(connector + ' -> XXX-preSV[g]')
        rxnNames.append('Start_Secretion')
        
        sv_rxns = []
        sv_names = []
        [sv_rxns, sv_names] =  addPathway('SV',sv_rxns, sv_names)       
        for i in range(len(sv_rxns)):
            sv_rxns[i] = sv_rxns[i].replace("!",str(clathrin_coeff))
            rxns.append(sv_rxns[i])
            rxnNames.append(sv_names[i])        
        
        #[rxns,rxnNames] = addPathway("SV",rxns,rxnNames)
        if location == '':
            location = '[e]'
        rxns.append('XXX_mature[sv]' + ' -> XXX_mature' + location)
        rxnNames.append('Final_location_' + location)
        rxns.append('XXX_mature' + location + ' -> ')
        rxnNames.append('Final_demand')
#------------------------------------------------------------------------------
    #Add coeffcients to SP_degradation and Ubiquitin_degradation reactions (if applicable)
    if 'SP_degradation' in rxnNames:
        SPaas = count_AAs(sequence[0:22]) #Amino acids in signal peptide assuming length is 22 on average
        rxns[rxnNames.index('SP_degradation')] = substitute_AAs_count(rxns[rxnNames.index('SP_degradation')],SPaas)
    if 'Ubiquitination_degradation' in rxnNames:
        new_aas = count_AAs(sequence[22:])
        rxns[rxnNames.index('Ubiquitination_degradation')] = rxns[rxnNames.index('Ubiquitination_degradation')].replace("!", "?")
        rxns[rxnNames.index('Ubiquitination_degradation')] = substitute_AAs_count(rxns[rxnNames.index('Ubiquitination_degradation')], new_aas)
        rxns[rxnNames.index('Ubiquitination_degradation')] = rxns[rxnNames.index('Ubiquitination_degradation')].replace("?",str(L-22))
#------------------------------------------------------------------------------
    GPRs = getGPRsFromRxnNames(rxnNames)
#------------------------------------------------------------------------------
    #Change the reaction names and their formulas to include the UniProt Identifier
    for i in range(len(rxns)):
        rxns[i] = insert_prot_name_in_rxnFormula(rxns[i],PSI[0])
    for i in range(len(rxnNames)):
        rxnNames[i] = insert_prot_name_in_rxnName(rxnNames[i],PSI[0])
    
    [rxns,rxnNames,GPRs] = addCanonicalRxns(rxns,rxnNames,GPRs)    
    return rxns,rxnNames,GPRs

#------------------------------------------------------------------------------
