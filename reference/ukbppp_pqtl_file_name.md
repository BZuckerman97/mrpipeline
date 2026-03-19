# UKBPPP_PQTL_FILE_NAME

UKBPPP_PQTL_FILE_NAME

## Usage

``` r
ukbppp_pqtl_file_name(synapse_id, olink_linker_file, olink_dir, olink_rsid_dir)
```

## Arguments

- synapse_id:

  String, synapse id to access from olink linker file

- olink_linker_file:

  String or Dataframe, file containing the olink linker file, or the
  dataframe of the linker file. This function is designed to handle the
  full linker file and will perform the necessary filtering and
  validation internally. It can handle cases where a single `synapse_id`
  maps to multiple rows (e.g., for protein complexes), as long as the
  underlying file paths are consistent.

- olink_dir:

  String, directory of the olink files

- olink_rsid_dir:

  String, directory of the olink rsid files

## Value

a list of the 2 filepaths, one for the ukbppp_pqtl data and the other is
the corresponding rsID metadata file, as well as the name of the assay

## Examples

``` r
if (FALSE) { # \dontrun{
# Assuming you have a linker file and the deCODE data directory
synapse_id <- "syn12345678"
olink_linker_file <- "path/to/your/olink_linker_file.csv"
olink_dir <- "path/to/olink_dir"
olink_rsid_dir <- "path/to/olink_rsid_dir"
file_paths <- ukbppp_pqtl_file_name(synapse_id = synapse_id,
                                  olink_linker_file = olink_linker_file,
                                  olink_dir = olink_dir,
                                  olink_rsid_dir = olink_rsid_dir)
print(file_paths)
} # }
```
