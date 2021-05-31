# BioMSA

> Multiple Sequence Alignment in JavaScript

BioMSA is a JavaScript library for computing alignments of biological
sequences (DNA or Protein) directly in the browser.
It performs progressive alignments using a weighting scheme similarly to
what programs like ClustalW or MUSCLE do.
To the best of my knowledge, it is the only multiple sequence alignment library
written in JavaScript... and it's also the fastest!

## Installing / Getting started

TODO

```shell
npm install biomsa
```

## Use BioMSA in a web page

TODO

## Features

-   Multiple sequence alignment of nucleic and proteic sequences.
-   Tree guided progressive alignment using a weighing scheme, substitution matrices,
    alignment score optimization by dynamic programming.
-   Approximative fast alignment of large DNA sequences (e.g. x5 16kbases mitochondrial DNA sequences in 100ms), using minimizers and diagonals extension.

## Licensing

"The code in this project is licensed under MIT license."
