Methods:

Part 0.5: Confirming Results of Gelman et al.
- adjusted the model to use one-hot encoding
- used models from gelman et al.
- tested full datasets for pab1, bgl3, ube4b to see that results match

Part 1: Testing in and out of secondary Structure
- process dataset
  - pdb files from alphafold or uniprot
  - stride file
  - determine if in sec str. excluding coils and turns
  - assign to each position
- create in sec str. and not in sec. str. datasets
  - match region shape (e.g. secondary structure length if there are fewer ss.
  residues, extra is random) and number of different positions to choose from
  - same dataset size (different sizes)
  - compare pearson's r
  - 3 trials, random datasets each time

Part 2: Comparing model accuracy and secondary structure fraction
  - train on all values
  - scale number of training values according to protein length
  - use ube4b as standard bc. it has around 50% sec str (middleground)
  - num of training values in sec. str. same ratio as percent of sec. str. in protein
  - majority of rest is used as test data (with same ratio)
  - 3 trials each
