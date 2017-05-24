## Use

### Demo

[Click here to try it out.](http://aaboyles.github.io/bioseq-js/demo/)

### Installation

`bioseq` is packaged as a module for Node environments, but automatically
attaches itself to `window` if included in a browser. To get it in Node, from
your project directory,

```bash
npm install bioseq --save
```

From there, you can require it wherever is convenient:

```javascript
const bioseq = require('bioseq');
```

Alternately, if you just want it in the browser, [download it](https://github.com/AABoyles/bioseq-js/archive/master.zip)
and include it with a script tag, like so:

```html
<script src="bioseq.js"></script>
```
### Quick Start Example

```javascript
var target = 'ATAGCTAGCTAGCATAAGC';
var query  = 'AGCTAcCGCAT';

var rst = bioseq.align(target, query);

console.log('Score: ', rst.score);
console.log('Starting position: ', rst.position);
console.log('CIGAR: ', bioseq.cigar2str(rst.CIGAR));

var fmt = bioseq.cigar2gaps(target, query, rst.position, rst.CIGAR);

console.log(fmt[0]);
console.log(fmt[1]);
```

For a more comprehensive example, view the source of [demo/index.html](https://github.com/AABoyles/bioseq-js/blob/master/demo/index.html#L67)

### API

bioseq is an object with the following functions: align, cigar2gaps, cigar2str,
bsg_enc_seq, gen_query_profile, gen_score_matrix, makeAlphabetMap, and
makeIntArray.

In general, the only functions of general interest are `align` and `cigar2*`,
so let's talk about those.

#### bioseq.align

`align` is the fundamental function of `bioseq`. It accepts seven parameters:

* *target* - target string
* *query*  - query string or query profile

And, optionally:

* is_local - perform local alignment
* matrix   - square score matrix or [match,mismatch] array
* gapsc    - [gap_open,gap_ext] array; k-length gap costs gap_open+gap_ext
* w        - bandwidth, disabled by default
* table    - encoding table. It defaults to bioseq.nt5.

Returns an object containing a `score`, `position`, and `CIGAR`. CIGAR is
encoded in the BAM way, where higher 28 bits keeps the length and lower 4 bits
the operation in order of 'MIDNSH'. See `bioseq.cigar2str()` for converting
CIGAR to string.

#### bioseq.cigar2str

Before I get in to what this function does, let's talk about what it needs.
CIGAR is an encoding scheme for sequence alignment instructions. Basically,
it's a list of instructions to make two sequences line up. These instructions
are assembled from a very simple set:

* M for match
* I for insertion
* D for deletion
* N for gap
* S for substitution
* H for hard clipping

Each of these can be preceded by an integer, indicating the number of base pairs
for which this instruction is correct between the two sequences. However, CIGARs
are stored as simpler machine instructions, rather than strings. So, if you
want to show the CIGAR (for some reason), you're going to need this.

#### bioseq.cigar2gaps

Takes a pair of sequences and a CIGAR and returns an array of sequences which
have been aligned according to the instructions in the CIGAR.

## Theory

`bioseq` implements local and global pairwise alignment with affine gap
penalties. There are two formulations: the Durbin formulation as is
described in his book and the Green formulation as is implemented in phrap.
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
good Smith-Waterman implementation in Javascript. [Heng Li](https://github.com/lh3)
thought he could write one over night by porting [ksw.c](https://github.com/lh3/bwa/blob/master/ksw.c)
to javascript. It took longer than he planned because he found a couple of
subtle bugs in ksw.c. And while he was porting the C code to javascript, he
realized that it is not that difficult to merge local and banded global
alignments in one function. Achieving that also took extra time.

The end product is a fast and lightweight javascript library for affine-gap
local and banded global pairwise alignment. With a modern Javascript engine, it
is not much slower than a non-SSE C implementation.
