# DECODE_PQTL_FILE_NAME

DECODE_PQTL_FILE_NAME

## Usage

``` r
decode_pqtl_file_name(unique_id, decode_linker_file, decode_dir)
```

## Arguments

- unique_id:

  String, unique sequence associated with each protein assayed

- decode_linker_file:

  String or Dataframe, the file path to the file containing the deCDOE
  metadata

- decode_dir:

  String, the file path to the directory where the deCODE summary
  statistics are present

## Value

a list of file paths for the deCODE proteomic summary statistics

## Examples

``` r
if (FALSE) { # \dontrun{
# Assuming you have a linker file and the deCODE data directory
linker_file <- "path/to/your/decode_linker.csv"
decode_data_dir <- "path/to/decode_summary_stats/"
# Get the first unique ID from the linker for demonstration
unique_protein_id <- data.table::fread(linker_file)$seqID[1]
file_paths <- decode_pqtl_file_name(unique_id = unique_protein_id,
                                  decode_linker_file = linker_file,
                                  decode_dir = decode_data_dir)
print(file_paths)
} # }
```
