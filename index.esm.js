/**
 * Compute a matrix to check if the rectangle between two terminals is empty.
 * 
 * On the matrix, m[i][j] is true if the rectangle between terminals i, j is empty.
 * As this matrix is symmetrical, only the lower half is filled.
 * @param {*} terminals terminals points to consider 
 * @param {*} successors a neighboring table; see rfst.js:getSuccessors()
 * @returns the lower triangular matrix with empty rect info
 */
function getEmptyRects(terminals, successors) {
  const ret = [];
  for (let i = 0; i < terminals.length; i++) {
    ret[i] = [];
    for (let j = 0; j < i; j++) {
      ret[i][j] = false;
    }
  }
  const set = (i, j) => {
    if (i > j) {
      ret[i][j] = true;
    } else if (i < j) {
      ret[j][i] = true;
    }
  };

  for (let i = 0; i < terminals.length; i++) {
    const p = terminals[i];
    const [x, y] = p;

    let topDist = Infinity;
		let bottomDist = Infinity;
		let oldTopDist = Infinity;
		let oldBottomDist = Infinity;
		let topX = x;
		let bottomX = x;
    
    for (let j = successors.east[i]; j >= 0; j = successors.east[j]) {
      const q = terminals[j];
      const dx = q[0] - x;
      let dy = q[1] - y;
      // Are they on the same horizontal/vertical line
      if (dx === 0 || dy === 0) {
        set(i, j);
        continue
      }
      if (dy > 0) {
        // Q is above P
        if (dy <= topDist) {
          set(i, j);
          if (q[0] > topX) {
            oldTopDist = topDist;
            topX = q[0];
          }
          topDist = dy;
        } else if (q[0] === topX && dy <= oldTopDist) {
          set(i, j);
        }
      } else {
        // Q is below P
        dy = -dy;
        if (dy <= bottomDist) {
          set(i, j);
          if (q[0] > bottomX) {
            oldBottomDist = bottomDist;
            bottomX = q[0];
          }
          bottomDist = dy;
        } else if (q[0] === bottomX && dy <= oldBottomDist) {
          set(i, j);
        }
      }
    }
  }
  return ret
}

function isEmptyRect(emptyRects, i, j) {
  if (i > j) {
    return emptyRects[i][j]
  } else if (i < j) {
    return emptyRects[j][i]
  } else {
    return true
  }
}

/**
 * Disjoin Set Union-Find data structure
 * Determines if two terminals are already connected, with a fast method to
 * join sets when they get linked.
 */
class DSUF {

  constructor() {
    this.map = new Map();
  }

  areConnected(i, j) {
    return this.map.has(i) && this.map.get(i).has(j)
  }

  connect(i, j) {
    const seti = this.map.get(i);
    const setj = this.map.get(j);
    if (!seti && !setj) {
      const union = new Set([i, j]);
      this.map.set(i, union);
      this.map.set(j, union);
    } else if (!seti) {
      setj.add(i);
      this.map.set(i, setj);
    } else if (!setj) {
      seti.add(j);
      this.map.set(j, seti);
    } else if (seti !== setj) {
      const [bigger, smaller] = seti.size > setj.size ? [seti, setj] : [setj, seti];
      smaller.forEach(t => {
        bigger.add(t);
        this.map.set(t, bigger);
      });
      bigger.add(i);
      bigger.add(j);
    }
  }

  /**
   * Given a member of a set, returns a canonical member of that set.
   * The result won't change until the set is united with other set.
   * @param {*} i The member to search
   */
  find(i) {
    const set = this.map.get(i);
    return set ? set.values().next().value : i
  }

}

/**
 * Compute a minimum spanning tree for a set of terminals (Kruskal's Algorithm)
 * @param {*} edges list of all edges on the graph, with their lengths
 * @returns the list of edges on the mst
 */
function mst(edges) {
  edges.sort((a, b) => a.len - b.len);
  const dsuf = new DSUF();
  const ret = [];
  edges.forEach(e => {
    if (!dsuf.areConnected(e.p1, e.p2)) {
      ret.push(e);
      dsuf.connect(e.p1, e.p2);
    }
  });
  return ret
}

/**
 * Computes the rectilineal (L1/Manhattan) distance between two terminals
 * @return distance between the two terminals
 */
function RDIST(t1, t2) {
  return Math.abs(t1[0] - t2[0]) + Math.abs(t1[1] - t2[1])
}

/**
 * Returns the distance between two points on a specified direction
 * @param {*} dir axis to measure
 * @param {*} a first point
 * @param {*} b second point
 * @returns distance between the points
 */
function DSTDIR(dir, a, b) {
  const axis = dir % 2;
  return Math.abs(a[axis] - b[axis])
}

/**
 * Returns the distance between two points on a direction perpendicular to
 * the specified one
 * @param {*} dir axis from which get the perpendicular
 * @param {*} a first point
 * @param {*} b second point
 * @returns distance between the points
 */
function DSTDIRP(dir, a, b) {
  const axis = 1 - (dir % 2);
  return Math.abs(a[axis] - b[axis])
}

/**
 * Generates a Steiner point between two terminals
 * @param {*} dir direction on which to generate the Steiner (starting from a)
 * @param {*} a first point
 * @param {*} b second point
 */
function SPOINT(dir, a, b) {
  return dir % 2 ? [a[0], b[1]] : [b[0], a[1]]
}

/**
 * Enumerates rectilineal edges for a terminal set.
 * Optionally can use empty rect info to reduce the number of edges
 * @param {*} terminals terminal points to consider
 * @param {*} emptyRects optional empty rects info to limit edges
 * @returns list of edges between choosen terminals
 */
function getEdges(terminals, emptyRects) {
  // TODO: empty rects not implemented yet
  const edges = [];
  for (let i = 0; i < terminals.length; i++) {
    for (let j = i + 1; j < terminals.length; j++) {
      edges.push({
        p1: i,
        p2: j,
        len: RDIST(terminals[i], terminals[j])
      });
    }
  }
  return edges
}

function getRmst(terminals, emptyRects) {
  const edges = getEdges(terminals);
  const theMst = mst(edges);
  return theMst
}

/**
 * Compute a table of Bottleneck Steiner Distances (BSDs) between terminals
 * @param {*} mst the minimum spanning tree edges
 * @returns the BSD table
 */
function getBSDs(mst) {
  const edges = [{ p1: 0, p2: 0, len: 0 }, ...mst];
  const adjacency = getAdjacency(mst);
  const bsds = { parent: [], edge: [], adjacency, edges };
  let next = mst.length + 1;
  bsds.parent[next] = bsds.edge[next] = null;
  for (let i = 1; i < edges.length; i++) {
    let { p1: u, p2: v } = edges[i];
    let pu = bsds.parent[u];
    let pv = bsds.parent[v];

    while (u !== v && pu && pv) {
      u = pu;
      v = pv;
      pu = bsds.parent[u];
      pv = bsds.parent[v];
    }
    if (!pu && !pv) {
      next++;
      bsds.parent[u] = next;
      bsds.parent[v] = next;
      bsds.edge[u] = i;
      bsds.edge[v] = i;
    } else if (!pu && pv) {
      bsds.parent[u] = pv;
      bsds.edge[u] = i;
    } else if (pu && !pv) {
      bsds.parent[v] = pu;
      bsds.edge[v] = i;
    }
  }
  // Make sure there are no gaps
  for (let i = 0; i <= next; i++) {
    if (!bsds.parent[i]) bsds.parent[i] = null;
    if (!bsds.edge[i]) bsds.edge[i] = null;
  }
  return bsds
}

/**
 * Finds the Bottleneck Steiner Distance in a previously computed BSD table
 * @param {*} bsds previously computed BSD table
 * @param {*} i first terminal
 * @param {*} j second terminal
 * @returns the BSD between the terminals
 */
function bsd(bsds, i, j) {
  if (i === j) return 0
  let index = -1;
  while (i !== j) {
    let ei = bsds.edge[i];
    let ej = bsds.edge[j];
    if (ei > index) index = ei;
    if (ej > index) index = ej;
    i = bsds.parent[i];
    j = bsds.parent[j];
  }
  return bsds.edges[index].len
}

/**
 * Generates adjacency lists for each node of a graph
 * @param {*} edges all edges in the graph
 * @returns array of adjacency lists, indexed by node
 */
function getAdjacency(edges) {
  const ret = [];
  edges.forEach((edge, i) => {
    const { p1, p2 } = edge;
    if (!ret[p1]) ret[p1] = [];
    if (!ret[p2]) ret[p2] = [];
    ret[p1].push({ edge: i, node: p2 });
    ret[p2].push({ edge: i, node: p1 });
  });
  return ret
}

/**
 * Computes the MST measured according to the Bottleneck Steiner Distance
 * @param {*} terms terminal indices to consider (relative to the bsds terms)
 * @param {*} bsds bsd structure
 */
function getBmst(terms, bsds) {
  const edges = [];
  for (let i = 0; i < terms.length; i++) {
    for (let j = i + 1; j < terms.length; j++) {
      edges.push({
        p1: i,
        p2: j,
        len: bsd(bsds, terms[i], terms[j])
      });
    }
  }
  const theMst = mst(edges);
  return theMst
}

/**
 * Computes the MST length according to the BSD
 * @param {*} terms terminal indices to consider (relative to the bsds terms)
 * @param {*} bsds bsd structure
 */
function getBmstLength(terms, bsds) {
  const theMst = getBmst(terms, bsds);
  let total = 0;
  theMst.forEach(e => {
    total += e.len;
  });
  return total
}

/**
 * Compute rectilinear full Steiner trees, or RFSTs, for a given set of terminals.
 * 
 * @param {*} terminals terminals points to consider
 *   Eg:
 *   [[2, 1], [3, 7], ...]
 * 
 * @returns an object with the RFSTs
 *   Eg:
 *   {
 *     terminals: [[2, 1], ...]
 *     fsts: [
 *       {
 *         terminalIndices: [3, 8, 21, 33],
 *         length: 123.45,
 *         steinerPoints: [[2, 7], ...],
 *         edges: [[3, 8], ...],
 *         status: 0,
 *         incompatibleFstIndices: [2, 4, 6, ...] // really needed?
 *       },
 *       ...
 *     ]
 *   }
 */
function rfst(terminals) {
  // Preprocessing
  const successors = getSuccessors(terminals);
  const emptyRects = getEmptyRects(terminals, successors);
  const ub0 = getUb0(terminals, successors);
  const mst = getRmst(terminals);
  // const mstLen = mst.map(e => e.len).reduce((s, l) => s + l)
  const bsds = getBSDs(mst);
  const zt = getZt(terminals, successors, ub0, emptyRects, bsds);
  const ub1 = getUb1(terminals, successors, zt);
  
  // Growing some FSTs
  const ctx = {
    terminals, successors, emptyRects, ub0, bsds, zt, ub1,
    fsp: [],
    fsphash: {}
  };
  for (let dir = 0; dir < 2; dir++) {
    for (let i = 0; i < terminals.length; i++) {
      const fts = {
        terms: [i],
        longterms: [i, -1],
        maxedges: [],
        shortterm: [],
        lrindex: [],
        term_check: [],
        hash: [],
        length: 0,
        dir: dir,
        ubLength: 0,
        ubShortleg: [Infinity, Infinity],
        longindex: 0
      };
      fts.maxedges[i] = 0;
      for (let j = 0; j < terminals.length; j++) {
        fts.lrindex[j] = 0;
        fts.term_check[j] = false;
        fts.hash[j] = null;
      }
      growFST(ctx, fts);
    }
  }

  // Add MST edges too
  mst.forEach(({ p1, p2, len}) => {
    testAndSaveFst(ctx, {
      terms: [p1, p2],
      length: len,
      dir: 0,
      type: 1
    });
  });

  return { terminals, fsts: ctx.fsp }
}

/**
 * Generate a map from each element to its successors in each cardinal direction
 * @param {*} terminals terminal points to consider
 * @returns an object with arrays for each direction which map indices to successors
 */
function getSuccessors(terminals) {
  const { x: xOrdered, y: yOrdered } = getOrderedIndices(terminals);

  const west = [];
  const east = [];
  const north = [];
  const south = [];

  for (let i = 0; i < terminals.length; i++) {
    west[i] = east[i] = north[i] = south[i] = -1;
  }

  for (let i = 1; i < terminals.length; i++) {
    west[xOrdered[i]] = xOrdered[i-1];
    east[xOrdered[i-1]] = xOrdered[i];
    north[yOrdered[i]] = yOrdered[i-1];
    south[yOrdered[i-1]] = yOrdered[i];
  }
  
  return { west, east, north, south }
}

/**
 * Orders terminals by each axis
 * @param {*} terminals terminal points to consider
 * @returns object with arrays of ordered indices for each axis
 */
function getOrderedIndices(terminals) {
  const indicesX = Array.from(Array(terminals.length).keys());
  const indicesY = indicesX.slice();
  indicesX.sort((a, b) => {
    const [tax, tay] = terminals[a];
    const [tbx, tby] = terminals[b];
    if (tax > tbx) return 1
    if (tax < tbx) return -1
    if (tay > tby) return 1
    if (tay < tby) return -1
    return a - b
  });
  indicesY.sort((a, b) => {
    const [tax, tay] = terminals[a];
    const [tbx, tby] = terminals[b];
    if (tay > tby) return 1
    if (tay < tby) return -1
    if (tax > tbx) return 1
    if (tax < tbx) return -1
    return a - b
  });
  return { x: indicesX, y: indicesY }
}

/**
 * Compute UB0 bounds
 * @param {*} terminals terminal points to consider
 * @param {*} successors successors for each axis
 * @returns UB0 object
 */
function getUb0(terminals, successors) {
  const dirs = ['east', 'south', 'west', 'north'];
  const ub0 = {};
  dirs.forEach((d, di) => {
    const succ = successors[d];
    const arr = ub0[d] = [];
    terminals.forEach((p, i) => {
      let bound = Infinity;
      for (let j = succ[i]; j >= 0; j = succ[j]) {
        const p2 = terminals[j];
        const d1 = DSTDIR(di, p, p2);
        if (d1 > bound) break
        const d2 = DSTDIRP(di, p, p2);
        if (d1 > d2) {
          const d3 = d1 + d2;
          if (d3 < bound) {
            bound = d3;
          }
        }
      }
      arr[i] = bound;
    });
  });
  return ub0
}

/**
 * Find short leg candidates for each terminal and direction
 * @returns zt object
 */
function getZt(terminals, successors, ub0s, emptyRects, bsds) {
  const dirs = ['east', 'south', 'west', 'north'];
  const zt = { };

  dirs.forEach((d, di) => {
    zt[d] = [];
    const dirp = dirs[[3, 2, 3, 2][di]];
    const succ = successors[d];
    const ub0 = ub0s[d];
    const ub0p = ub0s[dirp];

    for (let i = 0; i < terminals.length; i++) {
      zt[d][i] = [];
      const p1 = terminals[i];
      let limit = ub0[i];
      for (let j = succ[i]; j >= 0; j = succ[j]) {
        const p2 = terminals[j];
        const d1 = DSTDIR(di, p1, p2);
        if (d1 === 0) continue
        if (d1 > limit) break
        const d2 = DSTDIRP(di, p1, p2);
        if (d2 === 0) break
        const lr = isLeft(di, p1, p2) === [0, 1, 1, 0][di];
        if (!lr) continue
        if (d2 > ub0p[j]) continue

        const b = bsd(bsds, i, j);
        if (d1 > b) continue
        if (d2 > b) continue
        if (isEmptyRect(emptyRects, i, j)) {
          // Candidate found
          zt[d][i].push(j);
        }
      }
    }
  });
  return zt
}

/**
 * Determines if p2 is at the left of the ray from p1 in direction dir
 * @param {*} dir the direction (east=0, south=1, west=2, north=3)
 * @param {*} p1 the first point
 * @param {*} p2 the second point
 * @returns 1 if p2 is strictly to the left, 0 otherwise
 */
function isLeft(dir, p1, p2) {
  const isRTL = [0, 1, 1, 0][dir];
  const axis = [1, 0, 1, 0][dir];
  const ret = isRTL ? p2[axis] >= p1[axis] : p2[axis] <= p1[axis];
  return ret ? 1 : 0
}

/**
 * Compute UB1 bounds
 * @returns UB1 object
 */
function getUb1(terminals, successors, zt) {
  const dirs = ['east', 'south', 'west', 'north'];
  const ub1s = {};
  dirs.forEach((d, di) => {
    const succ = successors[d];
    const ub1 = ub1s[d] = [];
    const dzt = zt[d];
    terminals.forEach((p1, i) => {
      const shortLegCandidates = dzt[i];
      if (!shortLegCandidates.length) {
        ub1[i] = 0;
        return
      }
      const last = shortLegCandidates[shortLegCandidates.length - 1];
      let bound = Infinity;
      const p3 = terminals[last];
      const steiner = SPOINT(di, p1, p3);
      const d3 = DSTDIRP(di, p1, p3);
      for (let j = succ[last]; j >= 0; j = succ[j]) {
        const p2 = terminals[j];
        const d1 = DSTDIR(di, steiner, p2);
        if (d1 > bound) break
        const d2 = DSTDIRP(di, steiner, p2);
        const lr = isLeft(di, p1, p2) === [0, 1, 1, 0][di];
        if (lr && d3 > d2) {
          bound = d1;
          break
        }
        if (d1 > d2) {
          const d4 = d1 + d2;
          if (d4 < bound) {
            bound = d4;
          }
        }
      }
      ub1[i] = DSTDIR(di, p1, p3) + bound;
    });
  });
  return ub1s
}

function growFST(ctx, fst) {
  const r = fst.terms[0];
  const l = fst.terms[fst.terms.length - 1];
  const lastlr = fst.lrindex[l];
  const dirName = ['east', 'south', 'west', 'north'][fst.dir];
  const succ = ctx.successors[dirName];

  const root = ctx.terminals[r];
  const last = ctx.terminals[l];
  let maxBackbone = Infinity;

  const lastDstDirp = DSTDIRP(fst.dir, root, last);

  let needsRestore = false;

  for (;;) {
    let i = fst.longterms[++fst.longindex];
    let lr;
    let p;
    let dstdirp;

    if (i < -1) {
      // No more candidates, and no more can be found
      break
    }

    if (i === -1) {
      // Dynamically add next candidate to longterms
      for (i = succ[fst.longterms[fst.longindex - 1]]; i >= 0; i = succ[i]) {
        p = ctx.terminals[i];
        dstdirp = DSTDIRP(fst.dir, root, p);
        if (dstdirp === 0) {
          lr = 2;
          fst.lrindex[i] = 2;
          fst.shortterm[i] = 0;
          fst.longterms[fst.longindex] = i;
          fst.longterms[fst.longindex + 1] = -2;
          break
        }

        lr = isLeft(fst.dir, root, p);
        const dirp = (fst.dir + (lr ? 1 : -1)) & 0x03;
        const dirpName = ['east', 'south', 'west', 'north'][dirp];

        // Find short leg candidate (if it exists)
        const candidates = ctx.zt[dirpName][i].slice().reverse();
        const j = candidates.find(k => isLeft(fst.dir, root, ctx.terminals[k]) === lr);
        fst.shortterm[i] = j !== undefined ? j : -1;

        // Check upper bounds
        let ub1 = 0;
        if (j >= 0) {
          ub1 = ctx.ub1[dirpName][i];
        }
        let d1 = ctx.ub0[dirpName][i];
        if (d1 < ub1) {
          d1 = ub1;
        }
        if (dstdirp > d1) continue

        fst.lrindex[i] = lr;
        fst.longterms[fst.longindex] = i;
        fst.longterms[fst.longindex + 1] = -1;
        break
      }
      if (i < 0) {
        // No further candidates
        fst.longterms[fst.longindex] = -2;
        break
      }
    } else {
      // next long leg candidate available in longterms
      p = ctx.terminals[i];
      lr = fst.lrindex[i];
      dstdirp = DSTDIRP(fst.dir, root, p);
    }

    const dstdir = DSTDIR(fst.dir, last, p);

    // Check if consecutive terminals share Steiner point
    if (fst.terms.length >= 3 && dstdir === 0) continue

    // Upd. max backbone length using empty diamond property
    if (dstdirp < dstdir) {
      const d1 = dstdir + dstdirp;
      if (d1 < maxBackbone) {
        maxBackbone = d1;
      }
    }

    // Upd. max backbone length using empty rect property
    if (fst.terms.length >= 2 && lr === lastlr && dstdirp < lastDstDirp) {
      if (dstdir < maxBackbone) {
        maxBackbone = dstdir;
      }
    }

    // Check length of new backbone segment
    if (dstdir > maxBackbone) break

    if (lr === 2) {
      // Terminal is on the backbone. Save as type 1 and break.
      if (fst.terms.length >= 2) {
        fst.terms.push(i);
        testAndSaveFst(ctx, {
          ...fst,
          length: fst.length + dstdir + dstdirp,
          type: 1
        });
        fst.terms.pop();
      }
      break
    }

    // Terminal on wrong side of long leg?
    if (fst.terms.length >= 2 && lr === lastlr) continue

    // Empty rect with last term?
    if (!isEmptyRect(ctx.emptyRects, l, i)) continue

    // Check if new backbone segment is longer than any BSD
    let passBsd = true;
    let minBsd = Infinity;
    for (let j = 0; j < fst.terms.length; j++) {
      const k = fst.terms[j];
      let d1 = fst.maxedges[k];
      if (dstdir > d1) {
        d1 = dstdir;
        fst.maxedges[k] = d1;
        needsRestore = true;
      }
      const b = bsd(ctx.bsds, i, k);
      if (d1 > b) {
        passBsd = false;
        break
      }
      if (b < minBsd) {
        minBsd = b;
      }
    }
    if (!passBsd) continue
    let newUbLength = fst.ubLength + minBsd;

    const dirp = (fst.dir + (lr + lr - 1)) & 3;
    const dirpName = ['east', 'south', 'west', 'north'][dirp];

    // Try to make a type 2 FST
    let tryType1 = false;
    let tryGrowing = false;

    let j = fst.shortterm[i];
    if (j < 0) tryType1 = true;

    // Is backbone rect empty?
    if (!isEmptyRect(ctx.emptyRects, r, j)) tryType1 = true;

    // Check UB1
    if (dstdirp > ctx.ub1[dirpName][i]) tryType1 = true;

    // Check BSD for each terminal in current tree
    let q;
    if (!tryType1) {
      q = ctx.terminals[j];
      for (let j2 = 0; j2 < fst.terms.length; j2++) {
        const k = fst.terms[j2];
        let d1 = fst.maxedges[k];
        let d2 = DSTDIRP(fst.dir, root, q);
        let d3 = DSTDIRP(fst.dir, p, q);
        if (d2 > d1) d1 = d2;
        if (d1 > d3) d3 = d1;
        if (d3 > bsd(ctx.bsds, i, k)) tryType1 = true;
        d3 = DSTDIR(fst.dir, p, q);
        if (d1 > d3) d3 = d1;
        if (d3 > bsd(ctx.bsds, j, k)) tryType1 = true;
      }
    }

    // Check short leg upper bound
    if (!tryType1) {
      if (DSTDIRP(fst.dir, root, q) > fst.ubShortleg[lr ? 1 : 0]) tryGrowing = true;

      // Perform FST tests and save as type 2
      if (!tryGrowing) {
        fst.terms.push(i, j);
        testAndSaveFst(ctx, {
          ...fst,
          length: fst.length + DSTDIR(fst.dir, last, q)
            + DSTDIRP(fst.dir, root, p),
          type: 2
        });
        fst.terms.pop();
        fst.terms.pop();
      }
    }

    // Try to make a type 1 FST
    // try_type1:
    if (!tryGrowing) {

      // Check UB0
      if (dstdirp > ctx.ub0[dirpName][i]) continue

      // Check BSD to each terminal in prev tree
      let passBsd = true;
      for (let j2 = 0; j2 < fst.terms.length; j2++) {
        const k = fst.terms[j2];
        let d1 = fst.maxedges[k];
        if (dstdirp > d1) d1 = dstdirp;
        if (d1 > bsd(ctx.bsds, k, i)) {
          passBsd = false;
          break
        }
      }
      if (!passBsd) continue

      // Check if BSD upper bound is shorter than partial tree w/steiner
      if (fst.length + dstdir > newUbLength) continue

      // Check if BSD upper bound is shorter than partial tree
      if (fst.length + dstdir + dstdirp > newUbLength) tryGrowing = true;

      // Do not make 2-terminal FSTs
      if (fst.terms.length <= 1) tryGrowing = true;

      // Is backbone rectangle empty?
      if (!isEmptyRect(ctx.emptyRects, r, i)) tryGrowing = true;

      // Check short leg upper bound
      if (dstdirp > fst.ubShortleg[lr]) tryGrowing = true;

      // Perform FST tests and save type 1 tree
      if (!tryGrowing) {
        fst.terms.push(i);
        newUbLength = testAndSaveFst(ctx, {
          ...fst,
          length: fst.length + dstdir + dstdirp,
          type: 1
        });
        fst.terms.pop();
      }
    }

    // Try to grow the current tree
    // try_growing:

    // Should we generate larger FSTs?
    if (fst.terms.length >= ctx.terminals.length) continue // maxFstSize??

    // Upper bound (A)
    let d1 = fst.ubShortleg[lr];
    if (dstdirp < d1) d1 = dstdirp;
    const newUbShortleg = [];
    newUbShortleg[lr] = d1;
    d1 = fst.ubShortleg[1-lr];
    if (fst.terms.length >= 2) {
      // Upper bound (B)
      let d2 = ctx.ub0[dirpName][i];
      if (minBsd < d2) d2 = minBsd;
      d2 -= dstdirp;
      if (d2 < d1) d1 = d2;

      // Upper bound (C)
      if (dstdir < d1) d1 = dstdir;

      if (fst.terms.length >= 3) {
        // Upper bound (D)
        let lp = fst.terms[fst.terms.length - 2];
        d2 = DSTDIRP(fst.dir, root, ctx.terminals[lp]);
        if (dstdirp < d2) d2 = dstdirp;
        d2 = DSTDIR(fst.dir, ctx.terminals[lp], p) - d2;
        if (d2 < d1) d1 = d2;
      }
    }
    newUbShortleg[1 - lr] = d1;

    fst.terms.push(i);
    fst.maxedges[i] = dstdirp;
    growFST(ctx, {
      ...fst,
      length: fst.length + dstdir + dstdirp,
      ubLength: newUbLength,
      ubShortleg: newUbShortleg
    });
    fst.terms.pop();
  }

  if (needsRestore) {
    // Restore caller's maxedges by recomputing them from scratch...
    let longLegMax = 0;
    for (let i = fst.terms.length - 1; i > 0; i--) {
      const k = fst.terms[i];
      const p = ctx.terminals[k];
      let d1 = DSTDIRP(fst.dir, root, p);
      if (longLegMax > d1) d1 = longLegMax;
      fst.maxedges[k] = d1;
      const j = fst.terms[i - 1];
      const q = ctx.terminals[j];
      d1 = DSTDIR(fst.dir, q, p);
      if (d1 > longLegMax) longLegMax = d1;
    }
    fst.maxedges[r] = longLegMax;
  }
}

function testAndSaveFst(ctx, fst) {
  const size = fst.terms.length;
  const dir = fst.dir;

  let type = fst.type;

  // Is this FST too large?
  if (size >= 2**31 - 1) return fst.length // Doesn't really make sense here

  if (size > 2) {
    const b = getBmstLength(fst.terms, ctx.bsds);
    if (fst.length >= b) return b
  }

  // Simple duplicate tests for size 3 and 4
  if (dir === 1) {
    if (size === 3) return fst.length
    if (size === 4 && fst.type !== 1) return fst.length
  }

  if (size > 4) {
    // No 2 terms on the long leg should share a steiner point
    for (let i = 1; i < size; i++) {
      const p = ctx.terminals[fst.terms[i]];
      const q = ctx.terminals[fst.terms[i-1]];
      const d1 = DSTDIR(dir, p, q);
      if (d1 === 0) return fst.length
    }
  } else if (size === 4) {
    const p1 = ctx.terminals[fst.terms[0]];
    const p2 = ctx.terminals[fst.terms[1]];
    if (DSTDIR(dir, p1, p2) === 0) return fst.length
    const p3 = ctx.terminals[fst.terms[2]];
    const p4 = ctx.terminals[fst.terms[3]];
    if (DSTDIR(dir, p3, p4) === 0) return fst.length

    if (DSTDIR(dir, p2, p3) === 0) {
      if (DSTDIRP(dir, p1, p4) !== 0) return fst.length
      type = 3; // Cross
    }
  } else if (size === 3) {
    // Make sure that 3-terminal FST is not degenerate
    const p1 = ctx.terminals[fst.terms[0]];
    const p2 = ctx.terminals[fst.terms[1]];
    const p3 = ctx.terminals[fst.terms[2]];
    const z12 = (DSTDIR(dir, p1, p2) === 0 ? 1 : 0) + (DSTDIRP(dir, p1, p2) === 0 ? 1 : 0);
    const z13 = (DSTDIR(dir, p1, p3) === 0 ? 1 : 0) + (DSTDIRP(dir, p1, p3) === 0 ? 1 : 0);
    const z23 = (DSTDIR(dir, p2, p3) === 0 ? 1 : 0) + (DSTDIRP(dir, p2, p3) === 0 ? 1 : 0);
    if (z12 + z13 > 1) return fst.length
    if (z12 + z23 > 1) return fst.length
    if (z13 + z23 > 1) return fst.length
  }

  // Check empty diamond property for transformed FST
  let i = 0;
  let last = size - 1;
  if (type === 2) {
    last = size - 2;
    if ((size & 1) === 0) i = 1; // type 2, even
  } else if ((size & 1) !== 0) i = 1; // type 1, odd
  const p1 = ctx.terminals[fst.terms[0]];
  const p2 = ctx.terminals[fst.terms[size - 1]];
  while (i < last) {
    // Check skew diamond
    let s1 = SPOINT(dir, p1, ctx.terminals[fst.terms[i]]);
    let s2 = SPOINT(dir, p2, ctx.terminals[fst.terms[i+1]]);
    if (!diamondEmpty(ctx, s1, s2, fst.terms[i], dir)) return fst.length
    i++;
    if (i >= last) break

    // Check for flat segment
    s1 = SPOINT(dir, p2, ctx.terminals[fst.terms[i]]);
    s2 = SPOINT(dir, p2, ctx.terminals[fst.terms[i+1]]);
    if (!diamondEmpty(ctx, s1, s2, fst.terms[i], dir)) return fst.length
    i++;
  }

  // Test that corner flipped FST does not split into 2+ FST
  const d1 = DSTDIRP(ctx.terminals[fst.terms[size - 1]], ctx.terminals[fst.terms[0]], dir);
  i = size - (fst.type === 1 ? 3 : 4);
  while (i > 0) {
    if (DSTDIRP(ctx.terminals[fst.terms[i]], ctx.terminals[fst.terms[0]], dir) <= d1) {
      return fst.length
    }
    i -= 2;
  }

  // Check for duplicates
  const hash = fst.terms.slice().sort().join(',');
  const dupe = ctx.fsphash[hash];
  if (dupe) {
    if (dupe.length <= fst.length) {
      return dupe.length
    }
    // The new one is shorter! Delete the old one
    const idx = ctx.fsp.indexOf(dupe);
    if (idx >= 0) {
      ctx.fsp.splice(idx, 1);
    }
  }

  // Build FST graph in edge list form
  fst.pterms = fst.terms.map(t => ctx.terminals[t]);
  buildRFSTGraph(ctx, fst);
  const info = {
    terminalIndices: fst.terms.map(t => t + 1), // One-indexed, and must make a copy
    steinerPoints: fst.steins,
    edges: fst.edges.map(e => [
      e.p1 < fst.terms.length ? e.p1 + 1 : (e.p1 - fst.terms.length + 1) * -1,
      e.p2 < fst.terms.length ? e.p2 + 1 : (e.p2 - fst.terms.length + 1) * -1
    ]),
    length: fst.length
  };
  ctx.fsp.push(info);
  ctx.fsphash[hash] = info;

  return fst.length
}

function buildRFSTGraph(ctx, fst) {
  const p1 = fst.pterms[0];
  const p2 = fst.pterms[1];
  fst.steins = [];
  fst.edges = [];
  const size = fst.terms.length;

  if (size <= 2) {
    if (p1[0] !== p2[0] && p1[1] !== p2[1]) {
      const p3 = p1[0] < p2[0] ? SPOINT(1, p1, p2) : SPOINT(1, p2, p1);
      fst.steins.push(p3);
      fst.edges.push({ p1: 0, p2: 2, len: RDIST(p1, p3) });
      fst.edges.push({ p1: 1, p2: 2, len: RDIST(p2, p3) });
      return
    } else {
      fst.edges.push({ p1: 0, p2: 1, len: RDIST(p1, p2) });
      return
    }
  }

  if (fst.type === 3) { // Cross
    const p3 = SPOINT(fst.dir, p1, p2);
    fst.steins.push(p3);
    for (let i = 0; i < 4; i++) {
      fst.edges.push({ p1: i, p2: 4, len: RDIST(fst.pterms[i], fst.steins[0]) });
    }
    return
  }

  // Decide whether to build a normal or corner-flipped topology
  let k = size - 1;
  if (fst.type === 2) k--;
  k = fst.terms[k];
  if (fst.dir === 1 && fst.lrindex[k] === 1) {
    // Build normal (unflipped) topology
    let p1 = fst.type === 1 ? fst.pterms[0] : fst.pterms[size - 1];
    let p2 = fst.pterms[size - 2];
    for (let i = size - 3; i >= 0; i--) {
      fst.steins[i] = SPOINT(fst.dir, p1, p2);
      p1 = fst.pterms[0];
      p2 = fst.pterms[i];
    }
  } else {
    // Build corner-flipped topology
    if (fst.type === 1) {
      k = (size & 1) === 0 ? size - 1 : 0;
    } else {
      k = (size & 1) === 0 ? 0 : size - 1;
    }
    let p1 = fst.pterms[k];
    let p2 = fst.pterms[1];
    for (let i = 0; i < size - 2; i++) {
      fst.steins[i] = SPOINT(fst.dir, p1, p2);
      p1 = fst.pterms[size - 1];
      p2 = fst.pterms[i + 2];
    }
  }

  // Now that Steiner points are in their places, build edges
  let j = 0;
  for (let i = 0; i < size - 2; i++) {
    const nj = size + i;
    fst.edges.push({ p1: j, p2: nj });
    fst.edges.push({ p1: nj, p2: i + 1 });
    j = nj;
  }
  fst.edges.push({ p1: j, p2: size - 1});

  const nedges = (size <= 2 ? 1 : fst.type === 3 ? 4 : 2 * size - 3) + 1;
  let includeCorner = true;
  for (let i = 0; i < nedges - 1; i++) {
    const { p1: j, p2: k } = fst.edges[i];
    const p1 = j < size ? fst.pterms[j] : fst.steins[j - size];
    const p2 = k < size ? fst.pterms[k] : fst.steins[k - size];
    if (includeCorner && p1[0] !== p2[0] && p1[1] !== p2[1]) {
      // Can only happen once
      includeCorner = false;
      const p3 = p1[0] < p2[0] ? SPOINT(1, p1, p2) : SPOINT(1, p2, p1);
      fst.steins.push(p3);

      // Introduce the corner point by replacing one edge by two
      // (one at the end)
      fst.edges[i] = {
        p1: j,
        p2: size + fst.steins.length - 1,
        len: RDIST(p1, p3)
      };
      fst.edges.push({
        p1: k,
        p2: size + fst.steins.length - 1,
        len: RDIST(p2, p3)
      });
    } else {
      fst.edges[i].len = RDIST(p1, p2);
    }
  }
}

function diamondEmpty(ctx, p, q, i, dir) {
  const dirName = ['east', 'south', 'west', 'north'][dir];
  const oppositeDirName = ['east', 'south', 'west', 'north'][dir + 2];
  
  let succ = ctx.successors[dirName];
  const d = RDIST(p, q);

  for (let j = succ[i]; j >= 0; j = succ[j]) {
    const r = ctx.terminals[j];
    const dstdir = DSTDIR(dir, r, q);
    if (dstdir > d) break
    const dstdirp = DSTDIRP(dir, r, q);
    if (RDIST(r, p) < d && dstdir + dstdirp < d) {
      return false
    }
  }

  succ = ctx.successors[oppositeDirName];
  for (let j = succ[i]; j >= 0; j = succ[j]) {
    const r = ctx.terminals[j];
    const dstdir = DSTDIR(dir, r, p);
    if (dstdir > d) break
    const dstdirp = DSTDIRP(dir, r, p);
    if (RDIST(r, q) < d && dstdir + dstdirp < d) {
      return false
    }
  }

  return true
}

/**
 * Intermediate module to use glpk.js without async issues
 */

 const constants = {
  /* direction: */
  'GLP_MIN': 1,  /* minimization */
  'GLP_MAX': 2,  /* maximization */

  /* kind of structural variable: */
  'GLP_CV': 1,  /* continuous variable */
  'GLP_IV': 2,  /* integer variable */
  'GLP_BV': 3,  /* binary variable */

  /* type of auxiliary/structural variable: */
  'GLP_FR': 1,  /* free (unbounded) variable */
  'GLP_LO': 2,  /* variable with lower bound */
  'GLP_UP': 3,  /* variable with upper bound */
  'GLP_DB': 4,  /* double-bounded variable */
  'GLP_FX': 5,  /* fixed variable */

  /* message level: */
  'GLP_MSG_OFF': 0,  /* no output */
  'GLP_MSG_ERR': 1,  /* warning and error messages only */
  'GLP_MSG_ON': 2,  /* normal output */
  'GLP_MSG_ALL': 3,  /* full output */
  'GLP_MSG_DBG': 4,  /* debug output */

  /* solution status: */
  'GLP_UNDEF': 1,  /* solution is undefined */
  'GLP_FEAS': 2,  /* solution is feasible */
  'GLP_INFEAS': 3,  /* solution is infeasible */
  'GLP_NOFEAS': 4,  /* no feasible solution exists */
  'GLP_OPT': 5,  /* solution is optimal */
  'GLP_UNBND': 6,  /* solution is unbounded */
 };

 const solve = async (lp, msgLevel) =>
    import('glpk.js')
      .then(glpk => {
        return glpk.default
      })
      .then(glpk => {
        return glpk(lp, msgLevel)
      });

function updateBestSolutionSet(edges, bbip) {
  // TODO: original function accepts either a LP solution,
  // a list of edge numbers, or a edge set (aka a bitmap).
  // We're only dealing with the list for now.

  const nedges = bbip.cip.edges.length;
  const edgeMask = new Set();

  // We are given an explicit set of edges
  for (let i = 0; i < edges.length; i++) {
    edgeMask.add(edges[i]);
  }

  // Compute length of this solution, in edge order
  let length = 0;
  for (let i = 0; i < nedges; i++) {
    if (edgeMask.has(i)) {
      length += bbip.cip.edges[i].length;
    }
  }

  // Note: Original code deals with keeping track of the K best solutions,
  // but we only need to keep the best one, so we simplified the code.

  if (bbip.solution && length >= bbip.solution.length) {
    // This solution is not the best
    return false
  }

  bbip.solution = {
    length,
    edges: []
  };
  for (let i = 0; i < nedges; i++) {
    if (edgeMask.has(i)) {
      bbip.solution.edges.push(i);
    }
  }

  bbip.upperbound = length;
  return true
}

/**
 * Pre-compute some info for ub heuristic
 * @param {*} cip the hyperedge info from phase 1
 * @returns an object with...
 */
function startupHeuristicUpperBound(cip) {
  const { terminals, edges } = cip;
  const n = edges.length;
  const fstLen = edges.map(e => e.length);
  const mstLen = [];
  const rankings = [];
  let nranks = 0;

  // Compute MST length of each FST
  for (let i = 0; i < n; i++) {
    const terms = edges[i].terminalIndices.map(ti => terminals[ti - 1]);
    const mst = getRmst(terms);
    mstLen[i] = mst.map(e => e.len).reduce((s, l) => s + l);
  }
  let ranking = computeFstRanking(fstLen, mstLen);
  rankings[nranks++] = ranking;

  // Pretend each edge in the MST has length 1
  for (let i = 0; i < n; i++) {
    mstLen[i] = edges[i].terminalIndices.length - 1;
  }
  ranking = computeFstRanking(fstLen, mstLen);
  rankings[nranks++] = ranking;

  const mstEdges = sortedMstEdges(edges);

  return {
    rankings,
    mstEdges
  }
}

/**
 * Produces a ranking of FSTs ordered by increasing (numerator / denominator)
 * @param {*} num list of numerators
 * @param {*} den list of denominators
 */
function computeFstRanking(num, den) {
  // To ensure stability, we round ratios to the 12th decimal place
  const factor = 1000000000000;
  const indices = num.map((_, i) => i);
  const ratios = indices.map(i => Math.round(num[i] / den[i] * factor) / factor);
  indices.sort((a, b) => ratios[a] - ratios[b]);

  const ranking = [];
  for (let i = 0; i < indices.length; i++) ranking[indices[i]] = i;
  return ranking
}

/**
 * Finds all of the MST edges and sorts them by increasing length
 * @param {*} edges the fsts we got from phase 1
 */
function sortedMstEdges(edges) {
  const mstEdges = [];
  for (let i = 0; i < edges.length; i++) {
    if (edges[i].edges.length === 2) {
      mstEdges.push(i);
    }
  }
  
  mstEdges.sort((a, b) => edges[a].length - edges[b].length);
  
  // Remove all non-MST edges from the list
  // TODO

  return mstEdges
}

/**
 * Takes an existing LP solution (a fractional lower bound for an
 * optimal Steiner tree) and attempts to perturb it into a feasible
 * integer solution. If this succeeds and we discover a shorter
 * solution than previously known, it becomes the current best feasible
 * integer solution.
 */
function computeHeuristicUpperBound(solution, bbip) {
  if (!bbip.ubip) {
    bbip.ubip = startupHeuristicUpperBound(bbip.cip);
  }

  const oldUb = bbip.upperbound || Infinity;
  bbip.ubip.bestZ = oldUb;

  // Classify edges in 3 buckets by weight: 1s, fractions, and 0s
  const edgesIntegral = [];
  const edgesFractional = [];
  const edgesZero = [];

  for (let i = 0; i < bbip.cip.edges.length; i++) {
    const weight = solution.vars['e' + i];
    if (weight <= 0.00001) {
      edgesZero.push(i);
    } else if (weight + 0.00001 >= 1) {
      edgesIntegral.push(i);
    } else {
      edgesFractional.push(i);
    }
  }

  // Repeat the greedy heuristic once for each FST rank ordering
  bbip.ubip.rankings.forEach((ranking, i) => {
    // Sort the integral FSTs by rank only
    edgesIntegral.sort(sortByRank(ranking));

    // Sort the fractional FSTs by LP weight, then rank
    edgesFractional.sort(sortByLpAndRank(solution.vars, ranking));

    // Sort the zero-weight FSTs by rank only
    edgesZero.sort(sortByRank(ranking));

    // Try several greedy trees using this ordering of FSTs
    const edgeList = [...edgesIntegral, ...edgesFractional, ...edgesZero];
    // console.log('*** Try Trees [1] for i =', i, 'got edge list:', edgeList)
    tryTrees(edgeList, bbip);

    // Create a second ordering by placing MST edges last
    edgeList.sort((a, b) => {
      const la = bbip.cip.edges[a].edges.length;
      const lb = bbip.cip.edges[b].edges.length;
      if (la === 2 && lb !== 2) return 1
      if (lb === 2 && la !== 2) return -1
      return 0
    });

    // console.log('*** Try Trees [2] for i =', i, 'got edge list:', edgeList)

    // Try a few greedy trees using this ordering
    tryTrees(edgeList, bbip);
  });

  return bbip.ubip.bestZ < oldUb
}

function sortByRank(ranking) {
  return (a, b) => ranking[a] - ranking[b]
}

function sortByLpAndRank(x, ranking) {
  return (a, b) => {
    const w1 = x['e' + a];
    const w2 = x['e' + b];
    if (w1 !== w2) return w1 - w2
    return ranking[a] - ranking[b]
  }
}

/**
 * Construct several trees according to the given sorted list of FSTs
 */
function tryTrees(edgeList, bbip) {
  const used = new Set();

  // Construct the initial tree. Note which edges were used.
  let l = ubKruskal(edgeList, used, bbip);

  if (l === Infinity) return // No initial tree found

  // Prepare a list to be able to quickly add an edge at the front
  const tempEdges = [null, ...edgeList];

  // Determine the limit of edges to try
  let limit;
  for (limit = 2; (1 << limit) < edgeList.length; limit++) {}

  // Compute small absolute gap value
  let smallGap = 0.0001 * Math.abs(bbip.ubip.bestZ || Infinity);
  if (smallGap < 0.0001) smallGap = 0.0001; // Warn: scaled!

  // Greedy local search
  let haveReset = false;
  let orgLimit = limit;
  let oldL = Infinity;
  let k = edgeList.length + 1;

  for (let i = 0; i < edgeList.length; i++) {
    let e = edgeList[i];
    if (used.has(e)) continue

    // Force edge e to be chosen first
    tempEdges[0] = e;

    // Clear old used map
    used.clear();

    l = ubKruskal(tempEdges, used, bbip);

    // If improved solution then replace tempEdges
    if (l < oldL) {
      k = 1;
      for (let j = 0; j < edgeList.length; j++) {
        let e = edgeList[j];
        if (used.has(e) || bbip.cip.edges[e].terminalIndices.length === 2) {
          tempEdges[k++] = e;
        }
      }
      oldL = l;

      // Compute small absolute gap value
      smallGap = 0.0001 * Math.abs(bbip.ubip.bestZ);
      if (smallGap < 0.0001) smallGap = 0.0001; // Warn: scaled!
    }

    // If we are really close to optimum then restart and intensify search
    let currGap = l - bbip.ubip.bestZ;
    if (!haveReset && currGap <= smallGap) {
      // console.log('*** Intensify!! currGap =', currGap, '<= smallGap =', smallGap)
      // Let the new limit be double, plus an a*x**2 term whose coefficient
      // is linearly dependent on the gap
      const fraction = (smallGap - currGap) / smallGap;
      const quadratic = Math.floor(fraction * orgLimit**2);
      limit = 2 * orgLimit + quadratic;

      i = 0;
      haveReset = true;
    }

    if (--limit <= 0) break
  }
}

/**
 * Use Kruskal's algorithm (extended for hypergraphs) to greedily construct
 * a tree from the given sequence of FSTs. If this yields a better solution
 * than previously known, record it and update the upper bound.
 */
function ubKruskal(edgeList, used, bbip) {
  const nverts = bbip.cip.terminals.length;
  const nedges = edgeList.length;

  const dsuf = new DSUF();

  let components = nverts;
  let length = 0;
  const treeEdges = [];
  // let ep1 = 0 // Not needed: just do treeEdges.push()
  let ep2 = 0;
  let ep3 = nedges;
  while (components > 1) {
    if (ep2 >= ep3) {
      // FSTs ran out before tree constructed
      length = Infinity;
      break
    }
    const e = edgeList[ep2++];
    // console.log('% -- New iter: e =', e)
    const mark = new Set();
    let vp1 = 0;
    const vp2 = bbip.cip.edges[e].edges.length * 2;
    while (true) {
      if (vp1 >= vp2) {
        // No cycle! Include e in solution
        // console.log('% No cycle, including', e)
        treeEdges.push(e);
        length += bbip.cip.edges[e].length;
        used.add(e);
        // Unite all subtrees joined
        const [i, ...roots] = Array.from(mark);
        roots.forEach(j => dsuf.connect(i, j));
        components -= mark.size - 1;
        break
      }
      const tij = bbip.cip.edges[e].edges[vp1 >> 1][vp1 & 1] - 1;
      vp1++;
      if (tij < 0) continue
      const oj = bbip.cip.edges[e].terminalIndices[tij] - 1;
      const j = dsuf.find(oj);

      if (mark.has(j)) {
        // console.log('% Cycle with', j, '/ was:', oj)
        break
      }
      // console.log('% Mark', j)
      mark.add(j);
    }
  }

  // console.log('% ****** Got: [', treeEdges.join(', '), '], and', components === 1 ? 'is a' : 'no', 'solution')

  if (components === 1) {
    // A solution was found
    if (updateBestSolutionSet(treeEdges, bbip)) {
      bbip.ubip.bestZ = bbip.upperbound;
    }
  }

  return length
}

/**
 * Gets the constraint pool with the initial set of constraints
 * @param {*} data Problem data, including terminals, masks, and edges
 */
function getConstraintPool(data) {
  const nterms = data.terminals.length;
  const vertMask = new Set(data.initialVertMask);
  const edgeMask = new Set(data.initialEdgeMask);
  const pool = [];

  // Note: we don't need to count sizes first, because
  // we don't deal with memory management ourselves

  // Add constraint for spanning
  pool.push({
    name: "spanning",
    vars: data.edges
      .filter((_, i) => edgeMask.has(i))
      .map((e, i) => ({
        name: 'e' + i,
        coef: e.terminalIndices.length - 1
      })),
    selected: true,
    bnds: {
      type: constants.GLP_FX, // =
      lb: vertMask.size - 1
    }
  });

  // Generate one cutset constraint per terminal
  data.terminals.forEach((_, i) => {
    if (!vertMask.has(i)) return
    pool.push({
      name: 'cutset-' + i,
      vars: data.termTrees[i]
        .filter(j => edgeMask.has(j))
        .map(j => ({
          name: 'e' + j,
          coef: 1
        })),
      selected: true,
      bnds: {
        type: constants.GLP_LO, // >=
        lb: 1,
      }
    });
  });

  // Generate one constraint per incompatible pair
  // TODO ...

  // Generate one constraint for each 2-SEC (Subtour Elimination Constraints)
  const fsmask = new Set();
  for (let i = 0; i < nterms; i++) {
    const tlist = [];
    const counts = data.terminals.map(_ => 0);
    const tmask = new Set();

    if (!vertMask.has(i)) continue
    data.termTrees[i].forEach(fs => {
      if (!edgeMask.has(fs)) return
      fsmask.add(fs);
      const fst = data.edges[fs];
      fst.edges.forEach(edge => {
        edge.forEach(vtx => {
          if (vtx > 0) {
            const j = fst.terminalIndices[vtx - 1] - 1;
            if (j <= i) return
            if (!vertMask.has(j)) return
            counts[j]++;
            if (tmask.has(j)) return
            tmask.add(j);
            tlist.push(j);
          }
        });
      });
    });
    tlist.forEach(j => {
      if (counts[j] < 2) return
      // Generate 2SEC {i, j}
      pool.push({
        name: '2sec-' + i + ',' + j,
        vars: data.termTrees[j]
          .filter(fs => fsmask.has(fs))
          .map(fs => ({
            name: 'e' + fs,
            coef: 1
          })),
        selected: false,
        bnds: {
          type: constants.GLP_UP, // <=
          ub: 1,
        }
      });
    });
    data.termTrees[i].forEach(fs => fsmask.delete(fs));
  }

  // Remove duplicate rows (we do it here instead of row-by-row)
  const poolDupes = new Set();
  const poolUniq = pool.filter(r => {
    const s = JSON.stringify({ ...r, name: null });
    const ret = !poolDupes.has(s);
    poolDupes.add(s);
    return ret
  });

  return poolUniq
}

/**
 * Set up the LP problem instance for the initial constraints of the LP relaxation
 * @param {*} cip Problem data
 * @param {*} pool Constraint pool
 */
function getInitialFormulation(cip, pool) {
  const objective = {
    direction: constants.GLP_MIN,
    name: 'obj',
    vars: cip.edges.map((e, i) => ({
      name: 'e' + i,
      coef: e.length
    }))
  };

  const binaries = cip.edges.map((_, i) => 'e' + i);

  const bounds = cip.edges.map((_, i) => ({
    name: 'e' + i,
    type: constants.GLP_DB,
    lb: 0,
    ub: 1
  }));
 
  return {
    name: 'LP',
    objective,
    // binaries,
    bounds,
    subjectTo: pool,
  }
}

/**
 * Solve the current LP relaxation.
 * First we solve it using only the selected rows, then we check the
 * non-selected ones to determine if they fit the solution.
 * Each row that doesn't fit, gets selected for the next iteration.
 * This process repeats until all constraints are met, or there is no solution.
 */
async function solveLpOverConstraintPool(bbip) {
  const pool = bbip.cpool;
  let solution;
  while (true) {
    const tableaux = pool.filter(r => r.selected);
    const lp = { ...bbip.lp, subjectTo: tableaux };
    solution = await solve(lp, constants.GLP_MSG_OFF);

    console.log({solution})

    if (solution.result.status !== constants.GLP_OPT) break
    
    let anyViolations = false;
    for (let i = 0; i < pool.length; i++) {
      const row = pool[i];
      const slack = computeSlackValue(solution, row);
      
      if (slack > 0.00001) continue // Not binding, much less violated
      // Binding
      if (row.selected) continue // Already in tableaux
      if (slack < -0.00001) {
        // Not in the tableaux, and violated. Add to tableaux
        row.selected = true;
        anyViolations = true;
      }
    }
    
    // Done if no violations found
    if (!anyViolations) break
  }

  // TODO: if (solution.result.status ===  glpk.GLP_OPT) ...

  return {
    status: solution.result.status,
    z: solution.result.z,
    name: solution.name,
    time: solution.time,
    vars: solution.result.vars
  }
}

/**
 * Computes the amount of slack (if any) for a coefficient row
 * with respect to a solution
 * @param {*} solution The solution to check against
 * @param {*} row The target row
 */
function computeSlackValue(solution, row) {
  const sv = solution.result.vars;
  const sum = row.vars.reduce((s, v) => s + v.coef * sv[v.name], 0);
  switch (row.bnds.type) {
    case constants.GLP_UP:
      return row.bnds.ub - sum
    case constants.GLP_FX:
      return Math.abs(sum - row.bnds.lb)
    case constants.GLP_LO:
      return sum - row.bnds.lb
  }
  return 0
}

function createBBTree() {
  let serial = 0;
  return {
    nextSerial: () => serial++,
    first: null,
    heapBest: new Heap(nodeIsBetter),
    heapWorst: new Heap(nodeIsWorse)
  }
}

function nodeIsBetter(a, b) {
  if (a.z < b.z) return true
  if (a.z > b.z) return false
  if (a.num >= b.num) return true
  return false
}

function nodeIsWorse(a, b) {
  return a.z >= b.z
}

/**
 * Binary heap class used to efficiently get the next interesting node
 */
class Heap {

  constructor(isParent) {
    this.isParent = isParent;
    this.array = [];
    this.indices = new Map();
  }

  getRoot() {
    return this.array[0]
  }

  insert(node) {
    // Add it at the end, and sift it up
    let i;
    for (i = this.array.length; i > 0;) {
      const j = (i - 1) >> 1;
      const node2 = this.array[j];
      if (this.isParent(node2, node)) break
      this.indices.set(node2, i);
      this.array[i] = node2;
      i = j;
    }
    this.indices.set(node, i);
    this.array[i] = node;
  }

  remove(node) {
    // Find node being deleted
    let i = this.indices.get(node);
    if (i === undefined) return // Not found

    // Deleted node is no longer here
    this.indices.delete(node);

    // Remove last element from heap
    const node2 = this.array.pop();
    if (node === node2) return // Removed last item, nothing to do

    // Assume that node2 will be in position i
    // First, sift it up...
    while (i > 0) {
      const j = (i - 1) >> 1;
      const node3 = this.array[j];
      if (this.isParent(node3, node2)) break
      this.indices.set(node3, i);
      this.array[i] = node3;
      i = j;
    }

    // Later, sift it down...
    while (i < this.array.length) {
      let j = (i << 1) + 1;
      if (j >= this.array.length) break
      let node3 = this.array[j];
      if (j + 1 < this.array.length) {
        const node4 = this.array[j+1];
        if (this.isParent(node4, node3)) {
          j++;
          node3 = node4;
        }
      }
      if (this.isParent(node2, node3)) break
      this.indices.set(node3, i);
      this.array[i] = node3;
      i = j;
    }
    this.indices.set(node2, i);
    this.array[i] = node2;
  }

  size() {
    return this.array.length
  }

}

/**
 * Prepares the phase 1 results for further processing
 * @param {*} terminals the terminal list used in phase 1
 * @param {*} edges the fsts we got from phase 1
 * @returns an object with all required info
 */
function prepare(terminals, edges) {
  const initialVertMask = new Set(terminals.map((_, i) => i));
  const initialEdgeMask = new Set(edges.map((_, i) => i));
  const termTrees = getTermTrees(edges);

  return {
    terminals, edges,
    initialVertMask, initialEdgeMask,
    termTrees
  }
}

/**
 * Create the "term_trees" array that is indexed by point number
 * and gives a list of all tree-numbers involving that point
 * @param {*} edges the fsts we got from phase 1
 */
function getTermTrees(edges) {
  const termTrees = [];
  edges.forEach((fst, i) => {
    fst.edges.forEach(edge => {
      edge.forEach(vtx => {
        if (vtx >= 0) {
          const vtxIdx = fst.terminalIndices[vtx - 1] - 1;
          if (!termTrees[vtxIdx]) termTrees[vtxIdx] = [];
          if (termTrees[vtxIdx].indexOf(i) === -1) {
            termTrees[vtxIdx].push(i);
          }
        }
      });
    });
  });
  return termTrees
}

/**
 * Create the branch-and-bound info object
 */
function getBbInfo(cip) {
  // Get pool of constraints
  const cpool = getConstraintPool(cip);

  // Get initial formulation
  const lp = getInitialFormulation(cip, cpool);

  // Initialize the branch-and-bound tree
  const bbtree = createBBTree();
  
  // Create vectors to describe the current problem
  // TODO: We never have edge masks or required edges, so we're skipping this

  // Create the root node
  const root = {
    z: -Infinity,
    optimal: false,
    num: bbtree.nextSerial(),
    iter: 0,
    parent: -1,
    var: -1,
    dir: 0,
    depth: 0,
    // ...
  };

  // Fill in the branch-and-bound info structure
  const bbip = {
    cpool,
    lp,
    bbtree,
    cip,
    prevlb: -Infinity,
    bestZ: Infinity,
    // ...
  };

  // Make the root node inactive by putting it in the bbtree
  bbtree.first = root;
  bbtree.heapBest.insert(root);
  bbtree.heapWorst.insert(root);

  return bbip
}

/**
 * Main branch-and-cut algorithm
 * @param {*} cip hyperedge data from prepare()
 * @param {*} bbip data from getBbInfo()
 * @returns the minimum spanning tree
 */
async function branchAndCut(cip, bbip) {

  // Get heuristic upper bound
  bbip.ubip = startupHeuristicUpperBound(cip);

  for (;;) {
    // Select the next node to process
    const node = selectNextNode(bbip.bbtree);
    if (!node) break

    // This is perhaps a new lower bound...
    newLowerBound(node.z, bbip);

    if (node.z > -Infinity) ;

    // Restore the LP tableaux and basis for this node
    const lp = node.lp;

    // Determine new preemption value
    const nextBest = bbip.bbtree.heapBest.getRoot();
    bbip.preemptZ = nextBest ? nextBest.z : bbip.bestZ;

    // Mod LP to represent problem from new node (??)
    // TODO ...

    // Set up new node to be processed
    bbip.node = node;

    // Process the current node
    await computeGoodLowerBound(bbip);

    break // TODO: No infinite loop!
  }
}

/**
 * Select the next node to process from a given bb tree
 * @param {*} bbtree 
 */
function selectNextNode(bbtree) {
  return bbtree.heapBest.getRoot()
}

/**
 * Computes the lower-bound for the current node, which consists of solving
 * the LP and generating violated constraints until either:
 * - It becomes infeasible
 * - Objective meets or exceeds cutoff value
 * - Solution is integral
 * - Separation finds no more violated constraints
 */
async function computeGoodLowerBound(bbip) {

  while (true) {
    let result = await solveLpOverConstraintPool(bbip);
    const z = result.z;
    bbip.node.iter++;

    // console.log('  % Node', bbip.node.num, 'LP', bbip.node.iter, 'Solution, length =', z)

    switch (result.status) {
      case constants.GLP_OPT:
        if (z >= bbip.bestZ) {
          bbip.node.z = bbip.bestZ;
          return 'cutoff'
        }
        break
      // case glpk.cutoff: // We don't have cutoff
      case constants.GLP_INFEAS:
        bbip.node.z = bbip.bestZ;
        return 'infeasible'
      default:
        throw new Error('Solve status = ' + result.status)
    }

    // Solution is feasible, check for integer-feasible...
    const { isInt, numFractional } = integerFeasibleSolution(result, bbip.cip);

    // Check if this node's objective value is now high enough to be preempted
    if (bbip.node.z > bbip.preemptZ) {
      // console.log('preempted')
      return 'preempted'
    }

    // Perhaps we have a new lower bound?
    newLowerBound(z, bbip);

    if (computeHeuristicUpperBound(result, bbip)) {
       newUpperBound(bbip.upperbound, bbip);
    }

    // If we have improved the upper bound, it is possible
    // that this node can now be cutoff
    // TODO

    // TODO: cp = do_separations() ...

    {
      // No more violated constraints found
      break
    }

    // TODO: get rid of slack rows
  }
}

/**
 * Checks if we have an integer feasible solution.
 * First we check for integrality, then we check connectedness.
 */
function integerFeasibleSolution(solution, cip) {
  let numFractional = 0;

  for (let i = 0; i < cip.edges.length; i++) {
    if (solution.vars['e' + i] <= 0.00001) continue
    if (solution.vars['e' + i] + 0.00001 >= 1) continue
    numFractional++;
  }

  if (numFractional) {
    // There are fractional variables; solution is NOT integer feasible
    return { isInt: false, numFractional }
  }

  // All solution variables are either 0 or 1 -- integral.
  // Note all edges present in the solution
  let j = 0;
  let startingEdge = -1;
  const integralEdges = new Set();
  for (let i = 0; i < cip.edges.length; i++) {
    // Object.keys(solution.vars).filter(v => solution.vars[v] >= 0.5))
    if (solution.vars['e' + i] >= 0.5) {
      integralEdges.add(i);
      startingEdge = i;
      j += cip.edges[i].terminalIndices.length - 1;
    }
  }

  // console.log('% *** j =', j, ', num_int', numInt)

  if (j !== cip.terminals.length - 1) {
    // Wrong cardinality of edges -- cannot be a tree
    return { isInt: false, numFractional }
  }

  if (startingEdge < 0) {
    // No edges in solution: problem must have one or fewer vertices.
    // This is connected by default.
    return { isInt: true, numFractional }
  }

  // Create temporary mask of vertices we have not yet seen
  const vertsLeft = new Set(cip.terminals.map((_, i) => i));

  // Find connected component containing the starting edge
  integralEdges.delete(startingEdge);
  const stack = [startingEdge];
  while (stack.length) {
    const fs = stack.pop();
    cip.edges[fs].terminalIndices.forEach(tpp => {
      const t = tpp - 1;
      if (!vertsLeft.has(t)) return
      vertsLeft.delete(t);
      cip.termTrees[t].forEach(fs2 => {
        if (!integralEdges.has(fs2)) return
        integralEdges.delete(fs2);
        stack.push(fs2);
      });
    });
  }

  // See if any vertices were not reached
  const notReached = vertsLeft.size !== 0;

  if (notReached) {
    // At least one more connected component -- solution is not connected,
    // and therefore infeasible (and at least one integer cycle)
    return { isInt: false, numFractional }
  }

  // Solution is a Steiner tree! (Not necessarily minimal)
  return { isInt: true, numFractional }
}

/**
 * Prints a new lower bound
 */
function newLowerBound(lb, bbip) {
  let prev = bbip.prevlb;
  if (lb <= prev) {
    return // No improvement
  }

  if (prev <= -Infinity) {
    prev = lb; // Don't jump from initial value
  }

  // Print the old and new lower bounds
  let oldGap, newGap;
  if (bbip.bestZ >= Infinity || bbip.bestZ === 0) {
    oldGap = newGap = 99.9;
  } else {
    oldGap = 100 * (bbip.bestZ - prev) / bbip.bestZ;
    newGap = 100 * (bbip.bestZ - lb) / bbip.bestZ;
  }

  // TODO: Set solver preempt value...

  // console.log(' % @LO', prev, oldGap)
  // console.log(' % @LN', lb, newGap)

  bbip.prevlb = lb;
}

/**
 * Prints a new upper bound
 */
function newUpperBound(ub, bbip) {
  let prev = bbip.bestZ;
  if (prev >= Infinity) { // Don't jump from infinity
    prev = ub;
  }

  let oldGap, newGap;
  if (bbip.prevlb <= -Infinity) {
    oldGap = newGap = 99.9;
  } else {
    oldGap = 100 * (prev - bbip.prevlb) / prev;
    newGap = 100 * (ub - bbip.prevlb) / ub;
  }

  // TODO: Set solver preempt value...

  // console.log(' % @UO', prev, oldGap)
  // console.log(' % @UN', ub, newGap)
}

/**
 * Returns a standard solution format from a previously computed solution
 * @param {*} bbip problem and solution info
 */
function buildSolution(bbip) {
  const terminals = bbip.cip.terminals;
  const solution = bbip.solution;
  const length = solution.length;
  const steiners = [];
  const edges = [];
  const edgeIds = [];
  solution.edges.forEach(e => {
    const fst = bbip.cip.edges[e];
    const steinerOffset = steiners.length;
    fst.steinerPoints.forEach(p => steiners.push(p));
    
    const id2absolute = n =>
      n > 0 ? fst.terminalIndices[n - 1] : -steinerOffset + n;
    const id2coords = n =>
      n > 0 ? terminals[n - 1] : steiners[-n - 1];
    
    fst.edges.forEach(v => {
      edgeIds.push(v.map(id2absolute));
      edges.push(v.map(id2coords));
    });
  });

  return {
    terminals: terminals,
    steiners: steiners,
    edges: edges,
    edgeIds: edgeIds,
    length: length
  }
}

/**
 * Our wrapper function to simplify usage
 * @param {*} terminals the terminal list used in phase 1
 * @param {*} edges the fsts we got from phase 1
 */
async function bb(terminals, edges) {
  const data = prepare(terminals, edges);
  const bbip = getBbInfo(data);
  await branchAndCut(data, bbip);
  return buildSolution(bbip)
}

/**
 * Main entry point
 * @param {*} terminals List of terminals to use
 * @returns The RSMT for the given terminals
 */
function rsmt(terminals) {
  if (terminals.length < 2) {
    return {
      terminals,
      steiners: [],
      edges: [],
      edgeIds: [],
      length: 0
    }
  }
  const { fsts } = rfst(terminals);
  return bb(terminals, fsts)
}

export default rsmt;
