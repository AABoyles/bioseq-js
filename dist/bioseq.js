(function(root, factory) {
  if (typeof define === "function" && define.amd) {
    define([], function() {
      return (root.bioseq = factory());
    });
  } else if (typeof module === "object" && module.exports) {
    module.exports = factory();
  } else {
    root.bioseq = factory();
  }
})(typeof self !== "undefined" ? self : this, function() {
  "use strict";
  var bioseq = {};

  bioseq.makeIntArray = function(length, bitSize, fill) {
    var a = bitSize <= 16 ? new Int16Array(length) : new Int32Array(length);
    var b = bitSize <= 8 ? new Int8Array(length) : a;
    var cani8 = typeof Int8Array !== "undefined";
    var arr = cani8 ? b : [];
    if (!cani8 || (fill !== 0 && fill !== undefined)) {
      for (var n = 0; n < length; n++) {
        arr[n] = fill;
      }
    }
    return arr;
  };

  bioseq.makeAlphabetMap = function(str, defaultVal) {
    var aMap = bioseq.makeIntArray(256, 8, defaultVal);
    var lstr = str.toLowerCase();
    for (var n = 0; n < str.length; n++) {
      aMap[str.charCodeAt(n)] = n;
      aMap[lstr.charCodeAt(n)] = n;
    }
    return aMap;
  };

  bioseq.nt5 = bioseq.makeAlphabetMap("ACGT", 4);
  bioseq.nac = bioseq.makeAlphabetMap("ACGTURYMKWSBDHVN", 16);
  bioseq.aac = bioseq.makeAlphabetMap("ARNDCQEGHILKMFPSTWYVBZX", 32);

  bioseq.align = function(reference, query, is_local, matrix, gapsc, w, table) {
    // Set some defaults:
    if (typeof is_local === "undefined") is_local = true;
    if (typeof matrix === "undefined") matrix = [1, -1];
    if (typeof gapsc === "undefined") gapsc = [-1, -1];
    //w is defaulted below...
    if (table == null) table = bioseq.nt5;

    // convert bases to integers
    var t = bioseq.bsg_enc_seq(reference, table);
    var qp = bioseq.gen_query_profile(query, matrix, table);
    var qlen = qp[0].length;

    // adjust band width
    var max_len = qlen > t.length ? qlen : t.length;
    w = w == null || w < 0 ? max_len : w;
    var len_diff = t.reference > qlen ? t.reference - qlen : qlen - t.reference;
    w = w > len_diff ? w : len_diff;

    // set gap score
    // these are penalties which should be non-negative
    var gapo, gape;
    if (typeof gapsc == "number") {
      gapo = 0;
      gape = gapsc > 0 ? gapsc : -gapsc;
    } else {
      gapo = gapsc[0] > 0 ? gapsc[0] : -gapsc[0];
      gape = gapsc[1] > 0 ? gapsc[1] : -gapsc[1];
    }
    var gapoe = gapo + gape; // penalty for opening the first gap

    // initial values
    var H = [],
      E = [],
      z = [],
      score,
      max = 0,
      end_i = -1,
      end_j = -1;
    if (is_local) {
      for (var j = 0; j <= qlen; ++j) {
        H[j] = E[j] = 0;
      }
    } else {
      H[0] = 0;
      E[0] = -gapoe - gapo;
      for (var j = 1; j <= qlen; ++j) {
        if (j >= w) {
          // everything is -inf outside the band
          H[j] = E[j] = Number.NEGATIVE_INFINITY;
        } else {
          H[j] = -(gapoe + gape * (j - 1));
          E[j] = -(gapoe + gapoe + gape * j);
        }
      }
    }

    // the DP loop
    for (var i = 0; i < t.length; ++i) {
      var h1 = 0,
        f = 0,
        m = 0,
        mj = -1;
      var zi,
        qpi = qp[t[i]];
      zi = z[i] = [];
      // only loop through [beg,end) of the query sequence
      var beg = i > w ? i - w : 0;
      var end = i + w + 1 < qlen ? i + w + 1 : qlen;
      if (!is_local) {
        h1 = beg > 0 ? Number.NEGATIVE_INFINITY : -(gapoe + gape * i);
        f = beg > 0 ? Number.NEGATIVE_INFINITY : -(gapoe + gapoe + gape * i);
      }
      for (var j = beg; j < end; ++j) {
        // At the beginning of the loop: h=H[j]=H(i-1,j-1), e=E[j]=E(i,j), f=F(i,j) and h1=H(i,j-1)
        // If we only want to compute the max score, delete all lines involving direction 'd'.
        var e = E[j],
          h = H[j],
          d;
        // set H(i,j-1) for the next row
        H[j] = h1;
        // h = H(i-1,j-1) + S(i,j)
        h += qpi[j];
        d = h > e ? 0 : 1;
        h = h > e ? h : e;
        d = h > f ? d : 2;
        // h = H(i,j) = max{H(i-1,j-1)+S(i,j), E(i,j), F(i,j)}
        h = h > f ? h : f;
        d = !is_local || h > 0 ? d : 1 << 6;
        // save H(i,j) to h1 for the next column
        h1 = h;
        mj = m > h ? mj : j;
        // update the max score in this row
        m = m > h ? m : h;
        h -= gapoe;
        h = !is_local || h > 0 ? h : 0;
        e -= gape;
        d |= e > h ? 1 << 2 : 0;
        // e = E(i+1,j)
        e = e > h ? e : h;
        // save E(i+1,j) for the next row
        E[j] = e;
        f -= gape;
        d |= f > h ? 2 << 4 : 0;
        // f = F(i,j+1)
        f = f > h ? f : h;
        // z[i,j] keeps h for the current cell and e/f for the next cell
        zi[j] = d;
      }
      (H[end] = h1), (E[end] = is_local ? 0 : Number.NEGATIVE_INFINITY);
      if (m > max) {
        max = m;
        end_i = i;
        end_j = mj;
      }
    }
    if (is_local && max == 0) {
      return null;
    }
    score = is_local ? max : H[qlen];

    // backtrack to recover the alignment/cigar
    function push_cigar(ci, op, len) {
      if (ci.length == 0 || op != (ci[ci.length - 1] & 0xf))
        ci.push((len << 4) | op);
      else ci[ci.length - 1] += len << 4;
    }

    var cigar = [],
      tmp,
      which = 0,
      i,
      k,
      start_i = 0;

    if (is_local) {
      (i = end_i), (k = end_j);
      if (end_j != qlen - 1) {
        // then add soft cliping
        push_cigar(cigar, 4, qlen - 1 - end_j);
      }
    } else {
      // (i,k) points to the last cell
      i = t.length - 1;
      k = (i + w + 1 < qlen ? i + w + 1 : qlen) - 1;
    }

    while (i >= 0 && k >= 0) {
      tmp = z[i][k - (i > w ? i - w : 0)];
      which = (tmp >> (which << 1)) & 3;
      if (which == 0 && tmp >> 6) {
        break;
      }
      if (which == 0) {
        which = tmp & 3;
      }
      if (which == 0) {
        // match
        push_cigar(cigar, 0, 1);
        --i, --k;
      } else if (which == 1) {
        // deletion
        push_cigar(cigar, 2, 1);
        --i;
      } else {
        // insertion
        push_cigar(cigar, 1, 1), --k;
      }
    }
    if (is_local) {
      if (k >= 0) {
        // add soft clipping
        push_cigar(cigar, 4, k + 1);
      }
      start_i = i + 1;
    } else {
      // add the first insertion or deletion
      if (i >= 0) {
        push_cigar(cigar, 2, i + 1);
      }
      if (k >= 0) {
        push_cigar(cigar, 1, k + 1);
      }
    }
    for (var i = 0; i < cigar.length >> 1; ++i) {
      // reverse CIGAR
      tmp = cigar[i];
      (cigar[i] = cigar[cigar.length - 1 - i]),
        (cigar[cigar.length - 1 - i] = tmp);
    }
    return { score: score, position: start_i, CIGAR: cigar };
  };

  //## Generate scoring matrix from match/mismatch score
  //
  // * n    - size of the alphabet
  // * a    - match score, positive
  // * b    - mismatch score, negative
  //
  // Returns a sqaure scoring matrix. The last row and column are zero, for
  // matching an ambiguous residue.
  bioseq.gen_score_matrix = function(n, a, b) {
    var m = [];
    // mismatch score b should be non-positive
    if (b > 0) {
      b = -b;
    }
    for (var i = 0; i < n - 1; ++i) {
      m[i] = [];
      for (var j = 0; j < n - 1; ++j) {
        m[i][j] = i == j ? a : b;
      }
      m[i][j] = 0;
    }
    m[n - 1] = [];
    for (var j = 0; j < n; ++j) {
      m[n - 1][j] = 0;
    }
    return m;
  };

  // Generate query profile (a preprocessing step)
  //
  // * _s     - sequence in string or post bsg_enc_seq()
  // * _m     - score matrix or [match,mismatch] array
  // * table  - encoding table; must be consistent with _s and _m
  //
  // Returns a query profile. It is a two-dimensional integer matrix.
  bioseq.gen_query_profile = function(_s, _m, table) {
    var s = typeof _s == "string" ? bioseq.bsg_enc_seq(_s, table) : _s;
    var qp = [],
      matrix;
    if (
      _m.length >= 2 &&
      typeof _m[0] == "number" &&
      typeof _m[1] == "number"
    ) {
      // match/mismatch score
      if (table == null) return null;
      var n = typeof table == "number" ? table : table[table.length - 1] + 1;
      matrix = bioseq.gen_score_matrix(n, _m[0], _m[1]);
    } else {
      //TODO: check if _m is already a matrix!
      //TODO: check if it is really a square matrix!
      matrix = _m;
    }
    for (var j = 0; j < matrix.length; ++j) {
      var qpj,
        mj = matrix[j];
      qpj = qp[j] = [];
      for (var i = 0; i < s.length; ++i) {
        qpj[i] = mj[s[i]];
      }
    }
    return qp;
  };

  //## Generic routines
  // Encode a sequence string with table
  //
  // * seq   - sequence
  // * table - encoding table; must be of size 256
  //
  // Returns an integer array
  bioseq.bsg_enc_seq = function(seq, table) {
    if (table == null) table = bioseq.nt5;
    var s = new Array(seq.length);
    for (var i = 0; i < seq.length; ++i) {
      s[i] = table[seq.charCodeAt(i)];
    }
    return s;
  };

  bioseq.cigar2str = function(cigar) {
    var s = [];
    for (var k = 0; k < cigar.length; ++k) {
      s.push((cigar[k] >> 4).toString() + "MIDNSHP=XB".charAt(cigar[k] & 0xf));
    }
    return s.join();
  };

  bioseq.cigar2gaps = function(reference, query, start, cigar, fixed) {
    var oq = "",
      ot = "",
      lq = 0,
      lt = start;
    for (var k = 0; k < cigar.length; ++k) {
      var op = cigar[k] & 0xf,
        len = cigar[k] >> 4;
      if (op == 0) {
        // match
        oq += query.substr(lq, len);
        ot += reference.substr(lt, len);
        (lq += len), (lt += len);
      } else if (op == 1) {
        // insertion
        if (!fixed) {
          oq += query.substr(lq, len);
          ot += "-".repeat(len);
        }
        lq += len;
      } else if (op == 2) {
        // deletion
        oq += "-".repeat(len);
        ot += reference.substr(lt, len);
        lt += len;
      } else if (op == 4) {
        // soft clip
        lq += len;
      }
    }
    return [ot, oq];
  };

  return bioseq;
});
