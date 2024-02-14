## VS-ligPrep-CBDD-filters
This tool assits in applied modified Lipinski's Rule of 5 and our CBDD stuctural alerts before virtural screening.

This is a modifed version of [rd_filters.py](https://github.com/PatWalters/rd_filters/tree/master/rd_filters). We added additional properties and CBDD structural alerts into `alert_collection.csv` file to generate new filters.

# Workflow:
Install RDKit toolkit on Linux using conda:

```
conda create --name rdkit-tools python=3.8
```

Activate RDKit envirionment:
```
conda activate rdkit-tools
```

Run the script:
```
python3 new_filters.py filter --in input.smi --prefix output --rules myrules.json --alert alert_collection.csv --np 4
```

The file `myrules.json` contains the property parameters, adjust the parameters accorrding to your need. If you prefer to use CBDD structural alert, then set the option`"Rule_CBDD"` to true.
