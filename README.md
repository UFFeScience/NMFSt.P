# NMFSt.P

<a href="https://sol.sbc.org.br/index.php/bresci/article/view/25492">NMFSt.P: A Notebook for Parallel Identification of Frequent Subtrees in Phylogenetic Tree Ensembles.</a>

## Prerequisites

- Python 3.10.12
- Clustalw (versão 2.1) 
- Arquivo FASTA com sequências de proteínas



## Installation of Dependencies

Before running the project, you must install the Python dependencies specified in the "requirements.txt" file. To do this, run the following command in the terminal:

```
pip install -r requirements.txt
```

To install Clustalw on Linux (Ubuntu):

```
sudo apt update
sudo apt-get install clustalw
```
## File Organization
Ensure that all required files, including protein sequences in FASTA format, are in the directory specified in 'input_path'.

## Workflow Steps

## 1. Construction of Phylogenetic Trees

The first step of the workflow is the construction of phylogenetic trees from the provided protein sequences. To do this, run the "Constructor.ipynb" script in the terminal:

```
python Constructor.ipynb
```
This script performs multiple sequence alignment using ClustalW and then builds the phylogenetic tree using the Neighbor-Joining (NJ) method.

## 2. Construction of Subtrees and MAF Analysis

After constructing the phylogenetic trees, the next step is to generate subtrees from the main trees and perform the MAF (subtree pair frequency matrix) analysis.

Run the "sub_find.ipynb" script in the terminal:

```
python sub_find.ipynb
```
This script will generate all subtrees from the phylogenetic trees and then calculate the subtree pairwise frequency (MAF) matrix. The result will be displayed on the terminal.

## Outputs

The generated subtrees will be saved in the "out/Subtrees" directory. Additionally, the subtree pair frequency (MAF) matrix will be displayed in the terminal while running the "sub_find.ipynb" script.

## Temporary File Cleanup

The "Constructor.ipynb" and "sub_find.ipynb" script will automatically clean up the temporary files generated during the process. Temporary files will be deleted from the "out/tmp/" directory.

## Final considerations

This guide provides an overview of the workflow for building phylogenetic trees and analyzing subtrees. Make sure that the input files are correctly organized in the indicated directories and run the scripts according to the steps described.

## License

This project is licensed under the [MIT License](https://opensource.org/licenses/MIT).
