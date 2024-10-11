# Transform barcode reads

This is a small program that parses barcodes out of a barcode read. It infers barcode location (and orientation) in the read by comparing sequences to a whitelist, and outputs a new fastq file containing the parsed barcodes as well as a file with the barcode counts.

## Installation

After you've installed the boost C++ libraries, you can install with:

make install PREFIX=. BOOST_ROOT=/path/to/boost_install