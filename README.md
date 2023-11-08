<p align="center">
  <img src="aminos.png" alt="aminos.png"/>
</p>

## Overview

Convert large genetic VCF files into FASTA files corresponding to individual's protein sequences. This repo can handle complex combinations of coding variants and implements the  sequence intermediate representation (SIR) algorithm from VCF2Prot ([paper](https://www.biorxiv.org/content/10.1101/2022.01.21.477084v1.full.pdf), [code](https://github.com/ikmb/vcf2prot)) to enable biobank-scale mapping.

## Currently supported variants

- frameshift
- frameshift&start_lost
- frameshift&stop_retained
- inframe_deletion
- inframe_deletion&stop_retained
- inframe_insertion
- inframe_insertion&stop_retained
- missense
- missense&inframe_altering
- start_lost
- start_lost&splice_region
- stop_gained
- stop_gained&inframe_altering
- stop_lost
- stop_lost&frameshift
