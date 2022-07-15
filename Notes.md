 # Notes:

## Domain

- Pab1 - domain is [(125, 200)]
- Bgl3 - no domain because PBD starts as same place as protein sequence (Boolean of PDB lines up with protein without change)
- Ube4b - domain is [(1071, 1173)], but position in Ube4b
- avGFP -
- GB1 -


## STRIDE

Documentation:

       H	    Alpha helix
	   G	    3-10 helix (variation of alpha helix)
	   I	    PI-helix
	   E	    Extended conformation (or strand - beta strands?)
	   B or	b   Isolated bridge (one hydrogen bond of the type found in a beta-sheet)
	   T	    Turn
	   C	    Coil (none of the above) (absence of regular secondary structure)


 - currently not including **C** or **T**

## Secondary Structure Fractions

| Protein | Source | Protein Length (domain) | Dataset Size | Scaled (Ideal) | Sec. Str. % | Alpha Hel. % | Beta Sheet % |
| ------- | ------ | ----------------------- | ------------ | ----------- | ------------ | -------------| ------ |
| Pab1 | Gelman | 577 (75) | 37,710 | 40 (53.6)| 69 | TBD | TBD |
| Bgl3 | Gelman | (501) | 91,031 | 360 (357.9)|54 | TBD | TBD |
| Ube4b | Gelman | 1173 (102) |25,737 | 80 (72.9)| 52
| Thermonuclease | Curated | 231 | 1,068 | 160 (165) | 62 |
| Endolysin | Curated | 164 | 1,376 | 120 (117.1) | 75 |
| Immunoglobulin G-binding protein G | Curated| 448 | 1,221 | 320  | 53
| avGFP | Gelman | 237 | 51,714 | 160 (168)| 64 |
| **GB1** | Gelman | (56) | 536,084 | **40** | 70 |
| Human glucokinase | Biorex | | |  |
| GAL4 | MaveDB | 881 | 1,196 | 640 (629.3)| 47 |
| Small ubiquitin-related modifier 1 | MaveDb | 101 | 1,919 | 80 (72.1)  | 46
| TAR DNA-binding protein 43 | MaveDB | 414 | 1,342 | 280 (295.7)| 36

**CHECK OFFSET\****





Meeting Notes (July 12th, July 21st)
- ~~find the fraction of secondary structure for the proteins we have now (even distribution, if not look for more)~~, and fraction alpha/beta
- train on differing amounts of data and find protein in middle of distribution that has Pearson's R of ~0.4 (middle is Immunoglobulin, close to 55%)
- Scale training value and run other proteins (plot?)
- (Pearson's R on y, sec. str. fraction on x)
- train on alpha/beta subsets (equal alpha and beta fractions)
- maybe check that test size doesn't impact pearson's r?

Meeting Notes (July 12th)
- Keep the test set the same for the same protein in diff. architectures
- Find number of samples in and out of secondary structure and keep that ratio for the train and test sheet
- 3 trials
- Use remaining dataset for test datasets
- Use Ube4b to scale _test_ datasets
- change input layer to architecture


Tuesday/Wednesday - fix code in order to run and get datasets
Thursday/Friday - run code
Monday - extra and plots


Questions


## Choosing Other proteins
- 3 from curated database:
   - nuclease - (ddG 1282 vals)
   - lysozyme - (ddG 1376 vals)
   - protein G - (ddg 1221 vals)
- 2 from Gelman et al.
  - GB1 -
  - avGFP -

## Additional Potential Data Sources
https://www.mavedb.org/scoreset/urn:mavedb:00000001-a-1/
https://fowlerlab.gs.washington.edu/deep-mutational-scanning-data
https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0261829
 - https://github.com/CrisSotomayor/perturbation-networks/tree/main/data

Gal4 - https://www.mavedb.org/scoreset/urn:mavedb:00000012-a-3/
https://www.mavedb.org/scoreset/urn:mavedb:00000051-c-1/
Alpha-synuclein - https://www.mavedb.org/scoreset/urn:mavedb:00000045-a-1/
P63165 (Small ubiquitin-related modifier 1) - https://www.mavedb.org/scoreset/urn:mavedb:00000001-b-1/

Linking protein structural and functional change to mutation using amino acid networks:

https://github.com/CrisSotomayor/perturbation-networks/tree/main/data

- 1be9: 3 (https://www.uniprot.org/uniprotkb/P31016/entry#structure)
- 1d5r:
- 1nd4:
- 3dqw:
- 4bz3:

MaveDB:

Homo Sapiens:
- ~~UBE2I:~~ 3
- ~~SUMO1:~~ around 50%
- ~~CALM1:~~ 3.5
- ~~TPK1:~~ 4
- **hYAP65 WW domain: come back**
- ~~BRCA1 RING domain: (maybe)~~
- ~~CBS:~~ 4
- ~~HMGCR:~~
- ~~LDLRAP1:~~
- ~~alpha-synuclein:~~
- ~~CCR5:~~
- ~~CXCR4:~~
- ~~MTHFR:~~ 4
- ~~MSH2:~~
- ~~ErbB2:~~
- ~~Glycophorin A:~~
- ~~PSD95 PDZ3:~~
- ~~NUDT15:~~
- Aβ42:
- p53:
- **TARDBP:**
- ~~RAF:~~
- ~~CYP2C19:~~
- ~~NCS1:~~
- ~~TP53 (P72R):~~
- ~~IGHG1:~~
- ~~BRCA1:~~
- ~~GCK:~~

Other:
- ~~SARS-CoV-2 receptor binding domain:~~
- ~~human L-Selectin:~~
- ~~DHFR:~~
- ~~LamB:~~
- ~~Dlg4 (PSD95_PDZ3):~~
- ~~HA (H1N1):~~

## Associated Files
- third run for limiting mutations (v3)


## Misc.
- Thermonuclease file is from alpha fold
