# BioMSA

> Multiple Sequence Alignment in JavaScript

BioMSA is a JavaScript library for computing alignments of biological
sequences (DNA or Protein) locally in the browser.
It performs progressive alignments using a weighting scheme similarly to
what programs like ClustalW or MUSCLE do.
To the best of my knowledge, it is the only multiple sequence alignment library
written in JavaScript... and it's also the fastest!

## Installing / Getting started

BioMSA is available as a single minified javascript file. Insert the following
snippet in your HTML page to start using it.

```html
<script src="https://cdn.jsdelivr.net/npm/biomsa/dist/biomsa.js"></script>
```

Once loaded, the global `biomsa` object is available and
its `align()` method can be called to align sequences.

```javascript
biomsa.align(["ACTGGGGAGGTGTA", "ACTGAGGTGTA"]).then((result) => {
  console.log(result);
});

// Array [ "ACTGGGGAGGTGTA", "ACT---GAGGTGTA" ]
```

Note: `align()` returns a promise.

### NPM package

BioMSA is also available as an NPM package and can be used as an ECMAScript compatible module too. It is shipped with type declarations too.

```shell
npm install biomsa
```

## Library options

The `align()` method has 2 parameters:

- an array of sequences to align
- an optional configuration object

```javascript
biomsa.align(
    ['SEQVENCE...', 'SEQWANCE...', 'CEQWANSE...'],
    {
        gapopen: -11,
        gapextend: -2,
        matrix: [[.....], [.....], ....],
        method: 'auto',
        type: 'auto',
        gapchar: '-',
        debug: false
    }).then(result => console.log(result))
```

### gapopen

Gap open penalty (a negative number). If not provided, it is set based on the sequence type.

### gapextend

Gap extend penalty (a negative number). If not provided, it is set based on the sequence type.

### matrix

Substitution score matrix as an array of number array. Cells are sorted by amino acid 1 letter code rank (A, C, D, E,...). If not provided, it is set based on the sequence type.

### method

Alignment method. By default, the method is set based on the sequences length.
For sizes greater than 1600 residues, the diagonal based heuristics is used. For shorter sequences a complete Needleman-Wunsch alignment is performed.

- `auto` default value
- `complete` computes an optimal alignment using Needleman-Wunsch algorithm. This can be slow and take a lot of memory for long sequences.
- `diag` computes an alignment by first finding common segments between sequences (called diagonals) and then
  computing the missing segments using NW algorithm.

### type

Sequence type. Can be `amino`, `nucleic` or `auto` when auto-detected. BioMSA encodes non-canonical residues randomly. For example 'B' in a protein sequence which could be "Asn" or "Asp" will be encoded randomly as one of these amino-acids.

### gapchar

Character to use in the aligned sequences to represent a gap. By default `-` (hyphen) is used.

### debug

Boolean. Set to true to report some debugging information to the javascript console.

## Features

- Multiple sequence alignment of nucleic and proteic sequences.
- Tree guided progressive alignment using a weighing scheme, substitution matrices,
  alignment score optimization by dynamic programming.
- Approximative fast alignment of large DNA sequences (e.g. x5 16kbases mitochondrial DNA sequences in 100ms), using minimizers and diagonals extension.

## Licensing

"The code in this project is licensed under MIT license."
