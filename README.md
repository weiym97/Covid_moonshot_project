# Covid_moonshot_project

## Background
SARS-CoV-2 is a single stranded RNA coronavirus. The main protease, $M^{pro}$, also o known as “3-chymotrypsin-like protease". PDB ID 6LU7 contains the $M^{pro}$ with a inhibitor (N3). Molecules that bind specifically to $M^{pro}$ will likely avoid off-target binding to human proteases and thus $M^{pro}$ is an ideal drug target. 
Machine learning models that identify new compounds predicted to bind $M^{pro}$ will accelerate the discovery of new anti-COVID-19 treatments. 

## Binding affinity
Binding affinity can be quantified in different ways, including:  
● free energy change upon binding, ΔG (in kcal/mol)  
● dissociation constant, Kd, (in M, molar)  
● inhibition constant, Ki, (in M, molar)  
● half-maximal inhibitory concentration, IC50, (in M, molar)  

Note: IC50 values can't be compared unless they were measured at the same substrate concentration.

The IC50 is the half maximal inhibitory concentration and indicates the potency of a substance in inhibiting. The lower the IC50 value the less substance is needed to inhibit. To make a drug 

## Molecular descriptors 
Molecular features of the compounds are studied here by the Lipinski's rule of five. This is a rule of thumb to indicate the bioavailablity of the small molecule. According to the 'rule of five'a bioavailable molecule should not violate more than one of the following:  

● 5 or fewer hydrogen bond donors;  
● 10 or fewer hydrogen bond acceptors;  
● A molecular weight (MW) of less than 500 Daltons;  
● An octanol-water partition coefficient (log Po/w) of less than 5. Larger log Po/w means
more lipophilic (i.e., less water soluble).  

Here, RDkit (the rdkit.Chem.Descriptors package) is used to compute the molecular descriptors.  

```python
def calculate_descriptors(smile):
  molecule = Chem.MolFromSmiles(smile)
  if molecule:
    hbd, hba, mw, pow = [x(molecule) for x in [NumHDonors, NumHAcceptors, MolWt, MolLogP]]
    res = [hbd, hba, mw, pow]
  else:
    res = [None] * 4
  return res
```

## ECFP fingerprint
ECFP are topological fingerproints for molecular characterization. They are circular fingerprints. Qualities of these fingerprints are: they can be very rapidly calculated; they are not predefined and can represent an essentially infinite number of different molecular features (including stereochemical information); their features represent the presence of particular substructures, allowing easier interpretation of analysis results; and the ECFP algorithm can be tailored to generate different types of circular fingerprints, optimized for different uses.

```python
#compute ECFP fingerprint with radius=2 and bit vector length=2048
from rdkit.Chem import AllChem
molecule = Chem.MolFromSmiles(smiles)
fp1 = AllChem.GetMorganFingerprintAsBitVect(molecule,2,nBits=1024)
```

The compounds are clustered based on the Tanimoto index:

```python
DataStructs.FingerprintSimilarity(fp1,fp2, metric=DataStructs.DiceSimilarity)
```
By visualising the ...
