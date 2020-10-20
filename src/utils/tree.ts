/**
 * @file tree
 * @author Paul Pillot <paul.pillot@libmol.org>
 * @license MIT
 * @copyright 2020
 */

import { TSequence } from './sequence';

export type TreeNode = {
    seq: TSequence;
    profile: any[];
    childA: number;
    childB: number;
    distance: number;
    numSeq: number[];
    msa: string[];
};

export type InternalNode = Omit<TreeNode, 'seq'>;

export function isLeafNode(node: TreeNode | InternalNode): node is TreeNode {
    return 'seq' in node;
}

/**
 * Cluster sequences based on a distance Matrix using the UPGMA algorithm
 *
 * @export
 * @param {*} mD
 * @param {*} tSeq
 * @returns
 */
export function makeTree(mD: number[][], tSeq: TSequence[]) {
    let clusters: (TreeNode | InternalNode)[] = [];
    let nbSeq = mD.length;
    let tabIdx: number[] = [];

    for (var i = 0; i < nbSeq; i++) {
        clusters[i] = {
            seq: tSeq[i],
            profile: [],
            childA: i,
            childB: i,
            distance: 0,
            numSeq: [i],
            msa: [],
        };
        tabIdx[i] = i;
    }

    for (var i = 0, nbIterations = nbSeq - 1; i < nbIterations; i++) {
        // make a node from the min value in distance matrix

        const node = findMinInDistanceMatrix(mD, tabIdx);

        // Add references to find back node

        tabIdx.push(nbSeq + i);
        clusters.push(node);

        // Recompute matrix

        mD = recomputeDistMatrix(mD, node['previousA'], node['previousB'], tabIdx);
    }

    return clusters;
}

/**
 * Loop through the distance matrix to find the minimum
 * TODO: this is O(n^2) and can easily be made O(n) by maintaining a vector
 * of minimum values and their position in each row of the matrix.
 * @param matrix
 * @param tabIdx
 */
function findMinInDistanceMatrix(matrix: number[][], tabIdx: number[]) {
    var min = 1,
        minX = 0,
        minY = 0,
        l = matrix.length;

    for (var i = 0; i < l; i++) {
        // *2 speedup by restricting to the triangular matrix
        for (var j = i + 1; j < l; j++) {
            if (matrix[i][j] < min) {
                min = matrix[i][j];
                minX = i;
                minY = j;
            }
        }
    }

    var a = {
        profile: [],
        childA: tabIdx[minX],
        childB: tabIdx[minY],
        distance: min,
        previousA: minX,
        previousB: minY,
        msa: [],
        numSeq: [],
    };
    return a;
}

function recomputeDistMatrix(matrix: number[][], x: number, y: number, tI: number[]) {
    const xyMax = Math.max(x, y);
    const xyMin = Math.min(x, y);
    const averages = [];
    const l = matrix.length;

    // loop matrix rows
    for (var i = 0; i < l; i++) {
        if (i == x || i == y) continue;

        // avg computation from MAFFT
        // When the distances between merged branches and the current one are
        // dissimilar, this corrects the value towards the distance to the
        // closest branch.

        const avg = (0.1 * (matrix[i][x] + matrix[i][y])) / 2 + 0.9 * Math.min(matrix[i][x], matrix[i][y]);

        // Add avg at the end of the row(i)

        matrix[i].push(avg);

        // collect the avg values to build the column

        averages.push(avg);

        // remove the cells at positions row(x) and row(y)

        matrix[i].splice(xyMax, 1);
        matrix[i].splice(xyMin, 1);
    }

    // the bottom of the future last column has the identity value (dist=0)

    averages.push(0);

    // Append column to matrix

    matrix.push(averages);

    // remove the values that were merged in averages

    matrix.splice(xyMax, 1);
    matrix.splice(xyMin, 1);
    tI.splice(xyMax, 1);
    tI.splice(xyMin, 1);

    return matrix;
}
