/**
 * @file progressive.alignment.ts
 * @author Paul Pillot <paul.pillot@libmol.org>
 * @license MIT
 * @copyright 2021
 * @description Controller for the optimal multiple alignment procedure
 * Most of this algorithm is inspired by R.C. Edgar MUSCLE library (a).
 *
 * (a) Edgar, R.C. MUSCLE: a multiple sequence alignment method with reduced
 * time and space complexity. BMC Bioinformatics 5, 113 (2004).
 * https://doi.org/10.1186/1471-2105-5-113
 *
 */

import { DEBUG, TAlignmentParam } from "./params";
import { TSequence, distanceMatrix, sortMSA, distanceKimura } from "../sequence/sequence";
import { InternalNode, isLeafNode, makeTree, clustalWeights,
    compareTrees,
    setProfiles} from "../sequence/tree";
import Log from '../utils/logger';
import { pairwiseAlignment, MSASeqAlignment, MSAMSAAlignment } from "./align";
import { estringProduct, estringTransform } from "../utils/estring";
import { mergeProfiles } from "../sequence/profile";



export function progressiveAlignment(seq: TSequence[], pParam: TAlignmentParam) {

    // Recursive depth first tree based alignment
    const treeAlign = (node: InternalNode, tabWeight: number[]) => {

        const nodeA = tree[node.childA],
              nodeB = tree[node.childB];
        let result = {} as { score: number; };

        // note: msa property has already been computed when the tree
        // has had subtrees from a previous alignment, copied over.
        // Then, they can be skipped to avoid costly computations.
        if (!isLeafNode(nodeA) && (nodeA.msa.length === 0)) {
            treeAlign(nodeA, tabWeight);
        }

        if (!isLeafNode(nodeB) && (nodeB.msa.length === 0)) {
            treeAlign(nodeB, tabWeight);
        }

        // Now nodeA and nodeB are either leaf or msa, we can align them
        if (isLeafNode(nodeA)) { //A is a single sequence
            if (isLeafNode(nodeB)) { //B is a single sequence
                const lR = pairwiseAlignment(nodeA.seq, nodeB.seq, pParam);
                result.score = lR.score;

                node.numSeq = [nodeA.numSeq[0], nodeB.numSeq[0]];
                node.estring = [lR.estringA, lR.estringB];
                node.profile = mergeProfiles(
                    nodeA.profile,
                    nodeB.profile,
                    lR.estringA,
                    lR.estringB,
                    pParam
                )

                if (DEBUG)
                    Log.add(`seq ${nodeA.numSeq} - seq ${nodeB.numSeq}`);

            } else { //B is a MSA
                nodeB.tabWeight = tabWeight.filter((_, idx) => nodeB.numSeq.includes(idx));

                    // !!! B and A order are inverted

                const lR = MSASeqAlignment(nodeB, nodeA, pParam);
                result.score = lR.score;

                    // !!! estringA applies to nodeB

                node.numSeq = [nodeA.numSeq[0], ...nodeB.numSeq];
                node.estring = [
                    lR.estringB,
                    ...nodeB.estring.map(es => estringProduct(lR.estringA, es))
                ];
                node.profile = mergeProfiles(
                    nodeA.profile,
                    nodeB.profile,
                    lR.estringB,
                    lR.estringA,
                    pParam
                );

                if (DEBUG)
                    Log.add(`seq ${nodeA.numSeq} - MSA ${nodeB.numSeq}`);
            }
        } else if (!isLeafNode(nodeB)) { // A & B are both MSA

            nodeB.tabWeight = tabWeight.filter((_, idx) => nodeB.numSeq.includes(idx));
            nodeA.tabWeight = tabWeight.filter((_, idx) => nodeA.numSeq.includes(idx));

            const lR = MSAMSAAlignment(nodeA, nodeB, pParam);
            result.score = lR.score;

            node.estring = [
                ...nodeA.estring.map(es => estringProduct(lR.estringA, es)),
                ...nodeB.estring.map(es => estringProduct(lR.estringB, es))
            ];
            node.profile = mergeProfiles(
                nodeA.profile,
                nodeB.profile,
                lR.estringA,
                lR.estringB,
                pParam
            );

            if (DEBUG)
                Log.add(`MSA ${nodeA.numSeq} - MSA ${nodeB.numSeq}`);
        }
        node.numSeq = [...nodeA.numSeq, ...nodeB.numSeq];

        return result.score;
    };

    if (DEBUG)
        Log.add('Start Progressive Alignment');

    // Compute fast pair similarity (see k-mers binary Muscle)
    const lDistMatrixKMers = distanceMatrix(seq);

    if (DEBUG)
        Log.add('K-mer distance matrix');

    // Compute tree from distance matrix
    const tree = makeTree(lDistMatrixKMers, seq);
    const root = tree[tree.length - 1] as InternalNode;

    if (DEBUG)
        Log.add('Build Tree');

    // Compute weights (this is used to prevent close sequences to skew
    // the computation toward them, which would result in not disfavouring
    // gaps in their already aligned sequences).
    const weights = clustalWeights(tree);

    if (DEBUG)
        Log.add('Compute Weights - Start MSA');

    setProfiles(tree, pParam);

    // First alignment following guide tree
    const score1 = treeAlign(root, weights);
    root.msa = computeMSA(root, seq);

    let msa = sortMSA(root.msa, root.numSeq); // sort sequences in the order they came in

    if (DEBUG)
        Log.add('End MSA computation');

    return msa;

//================== For now, don't refine ====================

    // Refinement.
    // Now that we got a first alignment of all sequences, let see if some
    // sequences have a better likeliness now than what was found in the
    // first place with the fast k-mer based distance.
    // First step: build a distance matrix based on identities, pondered
    // using Kimura method for long branches (effective mutations are >
    // to measured mutations)
    const lDistMatrixKimura = distanceKimura(msa);

    if (DEBUG)
        Log.add('Distance Kimura');

    // compute new tree
    var tree2 = makeTree(lDistMatrixKimura, seq);
    const root2 = tree2[tree2.length - 1] as InternalNode;

    if (DEBUG)
        Log.add('Build tree 2');

    // compare trees.
    // Note: Impure function, it will modify tree2 if differences are found
    // to take advantage of all the alignments already compared.
    var treeIdentity = compareTrees(tree, tree2);

    if (treeIdentity === true) { // No possible improvement
        return msa;
    } else {
        //Repeat progressive alignment along this new tree
        var tabWeight2 = clustalWeights(tree2);

        //$log.debug(tabWeight);
        const score2 = treeAlign(root2, tabWeight2);

        //$log.debug('scores:', score1, score2);
        //var matriceDistancesKimura2 = _utils.distancesKimura(tree2[tree2.length - 1].msa);
        if (score2 > score1) {
            msa = sortMSA(root2.msa, root2.numSeq);
        }
        return msa;
    }

}

function computeMSA (root: InternalNode, seq: TSequence[]) {
    let lMSA: string[] = [];
    root.numSeq.forEach((seqIdx, alIdx) => {
        lMSA.push(estringTransform(seq[seqIdx].rawSeq, root.estring[alIdx]));
    });
    return lMSA;
}
