# bioseq-js

A Fast, Lightweight Smith-Waterman Sequence Aligner

## Use

### Demo

[Click here to try it out.](http://cdcgov.github.io/bioseq-js/demo/)

### Installation

`bioseq` is packaged as a module for Node environments, but automatically
attaches itself to `window` if included in a browser. To get it in Node, from
your project directory,

```bash
npm install bioseq --save
```

From there, you can require it wherever is convenient:

```javascript
const bioseq = require("bioseq");
```

Alternately, if you just want it in the browser, [download it](https://github.com/CDCgov/bioseq-js/archive/master.zip)
and include it with a script tag, like so:

```html
<script src="bioseq.js"></script>
```

### Quick Start Example

```javascript
var target = "ATAGCTAGCTAGCATAAGC";
var query = "AGCTAcCGCAT";

var rst = bioseq.align(target, query);

console.log("Score: ", rst.score);
console.log("Starting position: ", rst.position);
console.log("CIGAR: ", bioseq.cigar2str(rst.CIGAR));

var fmt = bioseq.cigar2gaps(target, query, rst.position, rst.CIGAR);

console.log(fmt[0]);
console.log(fmt[1]);
```

For a more comprehensive example, view the source of [demo/index.html](https://github.com/CDCgov/bioseq-js/blob/master/demo/index.html#L67)

### API

bioseq is an object with the following functions: align, cigar2gaps, cigar2str,
bsg_enc_seq, gen_query_profile, gen_score_matrix, makeAlphabetMap, and
makeIntArray.

In general, the only functions of general interest are `align` and `cigar2*`,
so let's talk about those.

#### bioseq.align

`align` is the fundamental function of `bioseq`. It accepts seven parameters:

- _reference_ - the reference sequence
- _query_ - the query sequence or profile

And, optionally:

- is_local - perform local alignment
- matrix - square score matrix or [match,mismatch] array
- gapsc - [gap_open,gap_ext] array; k-length gap costs gap_open+gap_ext
- w - bandwidth, disabled by default
- table - encoding table. It defaults to bioseq.nt5.

Returns an object containing a `score`, `position`, and `CIGAR`. CIGAR is
encoded in the BAM way, where higher 28 bits keeps the length and lower 4 bits
the operation in order of 'MIDNSH'. See `bioseq.cigar2str()` for converting
CIGAR to string.

#### bioseq.cigar2str

Before I get in to what this function does, let's talk about what it needs.
CIGAR is an encoding scheme for sequence alignment instructions. Basically,
it's a list of instructions to make two sequences line up. These instructions
are assembled from a very simple set:

- M for match
- I for insertion
- D for deletion
- N for gap
- S for substitution
- H for hard clipping

Each of these can be preceded by an integer, indicating the number of base pairs
for which this instruction is correct between the two sequences. However, CIGARs
are stored as simpler machine instructions, rather than strings. So, if you
want to show the CIGAR (for some reason), you're going to need this.

#### bioseq.cigar2gaps

Takes a pair of sequences and a CIGAR and returns an array of sequences which
have been aligned according to the instructions in the CIGAR.

- _reference_ - The same reference sequence you passed to `bioseq.align`
- _query_ - The same query you passed to `bioseq.align`
- _start_ - The `position` value of the object returned by `bioseq.align`
- _cigar_ - The `CIGAR` value of the object retruned by `bioseq.align`

And optionally,

- fixed - A boolean indicating whether the reference is fixed or modifiable.
  Basically, if you're aligning many sequences to the same reference, you should
  set fixed to `true`. Otherwise you can [just leave it out](https://developer.mozilla.org/en-US/docs/Glossary/Falsy).

## Theory

`bioseq` implements local and global pairwise alignment with affine gap
penalties. There are two formulations: the Durbin formulation as is
described in his book and the Green formulation as is implemented in [phrap](http://www.phrap.org/phredphrap/phrap.html).
The Durbin formulation is easier to understand, while the Green formulation
is simpler to code and probably faster in practice.

The Durbin formulation is:

    M(i,j) = max{M(i-1,j-1)+S(i,j), E(i-1,j-1), F(i-1,j-1)}

    E(i,j) = max{M(i-1,j)-q-r, F(i-1,j)-q-r, E(i-1,j)-r}

    F(i,j) = max{M(i,j-1)-q-r, F(i,j-1)-r, E(i,j-1)-q-r}

where q is the gap open penalty, r the gap extension penalty and S(i,j) is
the score between the i-th residue in the row sequence and the j-th residue
in the column sequence. Note that the original Durbin formulation disallows
transitions between between E and F states, but we allow them here.

In the Green formulation, we introduce:

    H(i,j) = max{M(i,j), E(i,j), F(i,j)}

The recursion becomes:

    H(i,j) = max{H(i-1,j-1)+S(i,j), E(i,j), F(i,j)}

    E(i,j) = max{H(i-1,j)-q, E(i-1,j)} - r

    F(i,j) = max{H(i,j-1)-q, F(i,j-1)} - r

It is in fact equivalent to the Durbin formulation. In implementation, we
calculate the scores in a different order:

    H(i,j)   = max{H(i-1,j-1)+S(i,j), E(i,j), F(i,j)}

    E(i+1,j) = max{H(i,j)-q, E(i,j)} - r

    F(i,j+1) = max{H(i,j)-q, F(i,j)} - r

i.e. at cell (i,j), we compute E for the next row and F for the next column.

## History

In an online conversation, Istvan Albert was complaining he could not find a
good Smith-Waterman implementation in Javascript.
[Li Heng](https://github.com/lh3) thought he could write one over night by
porting [ksw.c](https://github.com/lh3/bwa/blob/master/ksw.c) to javascript. It
took longer than he planned because he found a couple of subtle bugs in ksw.c.
And while he was porting the C code to javascript, he realized that it is not
that difficult to merge local and banded global alignments in one function.
Achieving that also took extra time.

The end product is a fast and lightweight javascript library for affine-gap
local and banded global pairwise alignment. With a modern Javascript engine, it
is not much slower than a non-SSE C implementation.

This repository adapted Li's original algorithm for
[release on NPM](https://www.npmjs.com/package/bioseq). Li did all the hard
work, and we owe him a large debt of gratitude for it.

## Public Domain

This repository constitutes a work of the United States Government and is not
subject to domestic copyright protection under 17 USC § 105. This repository is in
the public domain within the United States, and copyright and related rights in
the work worldwide are waived through the [CC0 1.0 Universal public domain dedication](https://creativecommons.org/publicdomain/zero/1.0/).
All contributions to this repository will be released under the CC0 dedication. By
submitting a pull request you are agreeing to comply with this waiver of
copyright interest.

## License

The repository utilizes code licensed under the terms of the Apache Software
License and therefore is licensed under ASL v2 or later.

This source code in this repository is free: you can redistribute it and/or modify it under
the terms of the Apache Software License version 2, or (at your option) any
later version.

This source code in this repository is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the Apache Software License for more details.

You should have received a copy of the Apache Software License along with this
program. If not, see http://www.apache.org/licenses/LICENSE-2.0.html

The source code forked from other open source projects will inherit its license.

## Privacy

This repository contains only non-sensitive, publicly available data and
information. All material and community participation is covered by the
Surveillance Platform [Disclaimer](https://github.com/CDCgov/template/blob/master/DISCLAIMER.md)
and [Code of Conduct](https://github.com/CDCgov/template/blob/master/code-of-conduct.md).
For more information about CDC's privacy policy, please visit [http://www.cdc.gov/privacy.html](http://www.cdc.gov/privacy.html).

## Contributing

Anyone is encouraged to contribute to the repository by [forking](https://help.github.com/articles/fork-a-repo)
and submitting a pull request. (If you are new to GitHub, you might start with a
[basic tutorial](https://help.github.com/articles/set-up-git).) By contributing
to this project, you grant a world-wide, royalty-free, perpetual, irrevocable,
non-exclusive, transferable license to all users under the terms of the
[Apache Software License v2](http://www.apache.org/licenses/LICENSE-2.0.html) or
later.

All comments, messages, pull requests, and other submissions received through
CDC including this GitHub page are subject to the [Presidential Records Act](http://www.archives.gov/about/laws/presidential-records.html)
and may be archived. Learn more at [http://www.cdc.gov/other/privacy.html](http://www.cdc.gov/other/privacy.html).

## Records

This repository is not a source of government records, but is a copy to increase
collaboration and collaborative potential. All government records will be
published through the [CDC web site](http://www.cdc.gov).

## Notices

Please refer to [CDC's Template Repository](https://github.com/CDCgov/template)
for more information about [contributing to this repository](https://github.com/CDCgov/template/blob/master/CONTRIBUTING.md),
[public domain notices and disclaimers](https://github.com/CDCgov/template/blob/master/DISCLAIMER.md),
and [code of conduct](https://github.com/CDCgov/template/blob/master/code-of-conduct.md).
