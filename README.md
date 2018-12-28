# Insertion Sequence Annotation Tool for Sniffles

A simple insertion sequence annotation tool for sniffles SV result. Provide mobile element (Alu, Line, SVA) annotation and virus genome sequence (HBV, PV) annotation.

**Sniffles provide INS sequence since version v1.0.10. No need to retrive INS sequence through `sv_ins_seq`.**

## Dependencies

**Linux**

~~gcc 4.8+ cmake 3.2+~~

python3

blastn should be in the `PATH`

## Installation

~~shell
git clone --recursive git@gitlab.com:archieyoung/sv_ins_seq.git
cd sv_ins_seq
mkdir build
cd build
cmake ..
make~~


~~A executable `sv_ins_seq` will appear in the build/bin, `sv_ins_seq` gets consensus insert sequnce from bam file.~~ `ins_seq_annot.py` in scripts is the end-user script.

## Work Flow

~~1. Get consensus insert sequences from bam file for `INS` in sniffles vcf.~~
1. Get consensus insert sequences from sniffles vcf.

2. Blast consensus insert sequences against `Mobile element and virus genome` DB using default blastn parametes.

3. Get the best blast hit by choosing the hit which make coverage on query and coverage on target maximum.

## Usage

```shell
python3 ins_seq_annot.py 
usage: ins_seq_annot.py [options]

Insertion sequence annotation for sniffles vcf

optional arguments:
  -h, --help    show this help message and exit
  --vcf FILE    sniffles vcf file [default: None]
  --bam FILE    bam file which must match with sniffles vcf file [default:
                None]
  --prefix STR  Output file prefix [default: None]
```

Output files: prefix.ins.annot.txt is the annotation result.
prefix.ins.fasta is the consensus insert sequences.
prefix.ins.blast.txt is the blast report for the consensus insert sequences.


## Blast Database

Consensus mobile element sequences were downloaded from [A Comprehensive Map of Mobile Element Insertion Polymorphisms in Humans -- Table S11](https://journals.plos.org/plosgenetics/article/file?id=10.1371/journal.pgen.1002236.s029&type=supplementary)

HBV and PV sequence NCBI GIs were downloaded from [ViFi: accurate detection of viral integration and mRNA fusion reveals indiscriminate and unregulated transcription in proximal genomic regions in cervical cancer](https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/nar/46/7/10.1093_nar_gky180/3/gky180_supplemental_files.zip?Expires=1545904569&Signature=pExGC~WRSNARUk6kohIoT4k3VH~Tx6k0NqTnPyrrJWtDQHwL9fxfy0pG39fQfO4YBdljNIND2NYWv4b21QMmPcPiUZCP0Gb8MqDBcm~UeVKFJjxYrRk0sqYtd6KiPmVTOexoHrxRkscS9srHZ9BC2DOUYd1nlzpvM2RGDhR5p60xpsDR1wQUjrIBbaylhfTw6Ik6WkXqy4rFgidPp5H8a4G0cyPNGNJZEYePoNA7Lp5IfDZ-ea9UqJ9fgTeExNIrFRArQVHsAUG4-gpF7wK3tvBQ-eqvQac56hrTUJ8RrOHZclEvSFYLu6jrcDTDQ-R8hOtGesS2u0sXL2tXvw3acw__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA). 410 Sequences were downloaded from [Batch Entrez](https://www.ncbi.nlm.nih.gov/sites/batchentrez) by GIs

Blast database is named `Home_sapiens.mei_virus.db.fasta` and is in database directory.

## Future Plans

1. Retrive insert sequences from split mapping reads.

## Reference

Stewart C, Kural D, Strömberg MP, Walker JA, Konkel MK, et al. (2011) A Comprehensive Map of Mobile Element Insertion Polymorphisms in Humans. PLOS Genetics 7(8): e1002236. https://doi.org/10.1371/journal.pgen.1002236

Nam-phuong D Nguyen, Viraj Deshpande, Jens Luebeck, Paul S Mischel, Vineet Bafna; ViFi: accurate detection of viral integration and mRNA fusion reveals indiscriminate and unregulated transcription in proximal genomic regions in cervical cancer, Nucleic Acids Research, Volume 46, Issue 7, 20 April 2018, Pages 3309–3325, https://doi.org/10.1093/nar/gky180





