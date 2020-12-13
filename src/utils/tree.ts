/**
 * @file tree
 * @author Paul Pillot <paul.pillot@libmol.org>
 * @license MIT
 * @copyright 2020
 */

import { TSequence } from './sequence';

export enum NODE_TYPE {
    LEAF,
    NODE,
    ROOT
}

export type LeafNode = {
    seq: TSequence;
    profile: any[];
    childA: number;
    childB: number;
    distance: number;
    numSeq: number[];
    msa: string[];
    id: string,
    weight: number,
    parent: number,
    depth: number,
    type: NODE_TYPE
};

export type InternalNode = Omit<LeafNode, 'seq'> & {
    profile: any[],
    tabWeight: number[]
};

export function isLeafNode(node: LeafNode | InternalNode): node is LeafNode {
    return 'seq' in node;
}

export type Tree = (LeafNode|InternalNode)[];

/**
 * Cluster sequences based on a distance Matrix using the UPGMA algorithm
 *
 * @export
 * @param {*} mD
 * @param {*} tSeq
 * @returns
 */
export function makeTree(mD: number[][], tSeq: TSequence[]) {
    let clusters: Tree = [];
    let nbSeq = mD.length;
    let tabIdx: number[] = [];

    for (var i = 0; i < nbSeq; i++) {
        clusters[i] = {
            type: NODE_TYPE.LEAF,
            seq: tSeq[i],
            profile: [],
            childA: i,
            childB: i,
            distance: 0,
            numSeq: [i],
            msa: [],
            id: i.toString(),
            weight: 0,
            parent: -1,
            depth: 0
        };
        tabIdx[i] = i;
    }

    let lNode: InternalNode;
    let lMinX: number;  // index in matrix of column with the min value
    let lMinY: number;  // index in matrix of row with the min value.
    for (var i = 0, nbIterations = nbSeq - 1; i < nbIterations; i++) {
        // make a node from the min value in distance matrix

        [lNode, lMinX, lMinY] = findMinInDistanceMatrix(mD, tabIdx);

        // finalize node data with child properties
        const lChildA = clusters[lNode.childA];
        const lChildB = clusters[lNode.childB];

        // TODO: lexicographic order here. Is it stable?
        lNode.id = (lChildA.id < lChildB.id) ?
              `|${ lChildA.id },${ lChildB.id }|`
            : `|${ lChildB.id },${ lChildA.id }|`;

        lNode.depth = Math.max(lChildA.depth, lChildB.depth) + 1;

        lChildB.parent = lChildA.parent = clusters.length;

        // Add references to find back node

        tabIdx.push(nbSeq + i);
        clusters.push(lNode);

        // Recompute matrix

        mD = recomputeDistMatrix(mD, lMinX, lMinY, tabIdx);
    }

    clusters[clusters.length - 1].type = NODE_TYPE.ROOT;

    return clusters;
}

/**
 * Loop through the distance matrix to find the minimum
 * TODO: this is O(n^2) and can easily be made O(n) by maintaining a vector
 * of minimum values and their position in each row of the matrix.
 * @param matrix
 * @param tabIdx
 */
function findMinInDistanceMatrix (
    matrix: number[][],
    tabIdx: number[]
): [InternalNode, number, number] {
    let min = Infinity,
        minX = 0,
        minY = 0,
        l = matrix.length;

    for (let i = 0; i < l; i++) {
        // *2 speedup by restricting to the triangular matrix
        for (let j = i + 1; j < l; j++) {
            if (matrix[i][j] < min) {
                min = matrix[i][j];
                minX = i;
                minY = j;
            }
        }
    }

    const a: InternalNode = {
        profile: [],
        childA: tabIdx[minX],
        childB: tabIdx[minY],
        distance: min,
        msa: [],
        numSeq: [],
        type: NODE_TYPE.NODE,
        depth: 0,
        id: '',
        tabWeight: [],
        parent: -1,
        weight: 0
    };

    return [a, minX, minY];
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

        const avg = (
            0.1 * (matrix[i][x] + matrix[i][y])) / 2
            + 0.9 * Math.min(matrix[i][x], matrix[i][y]
        );

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

/**
 * Provides sequence weight following clustalw weighting scheme
 * From clwwt.cpp in MUSCLE (Edgar) :
 * Compute weights by the CLUSTALW method.
    Thompson, Higgins and Gibson (1994), CABIOS (10) 19-29;
    see also CLUSTALW paper.

    Weights are computed from the edge lengths of a rooted tree.

    Define the strength of an edge to be its length divided by the number
    of leaves under that edge. The weight of a sequence is then the sum
    of edge strengths on the path from the root to the leaf.

    Example.

            0.2
           |-----A     0.1
         --|x       |------- B     0.7
           |--------|y          |----------- C
            0.3     |-----------|z
                        0.4     |-------------- D
                                        0.8

    Edge	Length	Leaves	Strength
    ----	-----	------	--------
    xy		0.3		3		0.1
    xA		0.2		1		0.2
    yz		0.4		2		0.2
    yB		0.1		1		0.1
    zC		0.7		1		0.7
    zD		0.8		1		0.8

    Leaf	Path		Strengths			Weight
    ----	----		---------			------
    A		xA			0.2					0.2
    B		xy-yB		0.1 + 0.1			0.2
    C		xy-yz-zC	0.1 + 0.2 + 0.7		1.0
    D		xy-yz-zD	0.1 + 0.2 + 0.8		1.1
    * @param   {object} cluster cluster tree computed
    * @returns {array} array containing weights for each node in the cluster tree
    */
 export function clustalWeights (cluster: Tree) {
    const tabWeights: number[] = [];
    let totalWeight = 0;
    const root = cluster[cluster.length - 1];

    function setWeight (node: LeafNode|InternalNode) {
        var edgeLength = 0,
            nbLeaves = 0;

        switch (node.type) {
        case NODE_TYPE.ROOT:
            node.weight = 0;
            setWeight(cluster[node.childA]);
            setWeight(cluster[node.childB]);
            break;
        case NODE_TYPE.NODE:
            edgeLength = cluster[node.parent].distance - node.distance;
            nbLeaves = node.id.split(',').length;
            node.weight = cluster[node.parent].weight + edgeLength / nbLeaves;

            setWeight(cluster[node.childA]);
            setWeight(cluster[node.childB]);
            break;
        case NODE_TYPE.LEAF:
            edgeLength = cluster[node.parent].distance;
            node.weight = cluster[node.parent].weight + edgeLength;
            totalWeight += node.weight;
            break;
        }
    };

    setWeight(root);

    //normalize weights to 1
    // TODO: does is matter that internal nodes weights are not normalized?
    var i = 0;
    while (cluster[i].type == NODE_TYPE.LEAF) {
        tabWeights[i] = cluster[i].weight /= totalWeight;
        i++;
    }

    return tabWeights;
};
