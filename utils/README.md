# Utils
`utils` contains several scripts to interact with the filesystems generated (or being generated) by DiscobaMultimer.

## annotate_AF2_models.sh
Utility to change the names of the filesystem from IDs to any name you want. It requires an `annotation.txt` TSV file with the correspondance between IDs and names (take a look at the example file).

## check_already_computed_AF2_models.sh
As its name suggests, it allows to check if all the IDs in the `IDs_table.txt` file were sucessfully predicted as AF2 models.

## check_if_IDs_are_in_proteome.sh
As its name suggests, it allows to check if all the IDs in the `IDs_table.txt` file are present in the `database.fasta` file. If all are present, you will get no output.
