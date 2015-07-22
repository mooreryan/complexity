# Sequencing complexity metrics #

Implementing metrics from

    Orlov, Y. L., & Potapov, V. N. (2004). Complexity: An internet resource for analysis of DNA sequence complexity. Nucleic Acids Research, 32(WEB SERVER ISS.), 628–633. http://doi.org/10.1093/nar/gkh466

## Requirements ##

Gems

- `parse_fasta`
- `ffi`

C libraries

- `libmpfr` (The GNU MPFR Library http://www.mpfr.org/)
- `libgmp` (The GNU Multiple Precision Arithmetic Library https://gmplib.org/)

## Get it working ##

Run `bin/compile_for_ffi` to set up the `.o` and `.so` files.

Assuming you have the `FFI` gem installed already, you should be able
to run `./complexity.rb --help` with no problems.
