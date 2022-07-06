 # Notes:

## Domain

- Pab1 - domain is [(125, 200)]
- Bgl3 - no domain because PBD starts as same place as protein sequence (Boolean of PDB lines up with protein without change)
- Ube4b - domain is [(1071, 1173)], but position in Ube4b


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

| Protein | Protein Length (domain) | Dataset Size | Sec. Str. % | Alpha Hel. % | Beta Sheet % |
| ------- | -------------- | ----------- | ------------ | ------------ | --------|
| Pab1 | 577 (75) | | 69 | TBD | TBD |
| Bgl3 | (501) | |54 | TBD | TBD |
| Ube4b | 1173 (102) || 52
| Thermonuclease | 231 | | 62 |
| Endolysin | 164 | | 75 |
| Immunoglobulin G-binding protein G | 448 | | 53
| avGFP | 237 | | 64 |
| GB1 | (56) | | 70 |
| Human glucokinase | |
| GAL4 | | | 47 |










Meeting Notes (July 12th, July 21st)
- find the fraction of secondary structure for the proteins we have now (even distribution, if not look for more), and fraction alpha/beta
- train on differing amounts of data and find protein in middle of distribution that has Pearson's R of ~0.4
- Scale training value and run other proteins (plot?)
- (Pearson's R on y, sec. str. of x)
- train on alpha/beta subsets (equal alpha and beta fractions)


Questions
- why is pab1 consistently higher in terms of secondary structure?


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




## Associated Files
- third run for limiting mutations (v3)


## Misc.
- Thermonuclease file is from alpha fold
