import { pairwiseAlignment } from '../align';
import { makeSequence } from '../../sequence/sequence'
import { getAlignmentParameters, SEQUENCE_TYPE } from '../params';
import { estringTransform } from '../../utils/estring';

test('aligns ends', () => {
    const lSeqA = makeSequence('CCCAATTCTATACCAACAC');
    const lSeqB = makeSequence('CCCAATTCTA');
    const lParam = getAlignmentParameters(SEQUENCE_TYPE.NUCLEIC);
    const lResult = pairwiseAlignment(lSeqA, lSeqB, lParam);
    const lRB = estringTransform(lSeqB.rawSeq, lResult.estrings[1], { gapchar: '-' });
    expect(lRB).toBe('CCCAATTCTA---------');
});

test('add gaps in sequence', () => {
    const lSeqA = makeSequence(  'AGCCCTTA');
    const lSeqB = makeSequence('ACTGCCCTCA');
    const lParam = getAlignmentParameters(SEQUENCE_TYPE.NUCLEIC);
    const lResult = pairwiseAlignment(lSeqA, lSeqB, lParam);
    const lRA = estringTransform(lSeqA.rawSeq, lResult.estrings[0], { gapchar: '-' });
    expect(lRA).toBe('A--GCCCTTA');
});

test('aligns starts', () => {
    const lSeqA = makeSequence('ACTGGGGAGGTGTA');
    const lSeqB = makeSequence('GAGGTGTA');
    const lParam = getAlignmentParameters(SEQUENCE_TYPE.NUCLEIC);
    const lResult = pairwiseAlignment(lSeqA, lSeqB, lParam);
    const lRA = estringTransform(lSeqB.rawSeq, lResult.estrings[1], { gapchar: '-' });
    expect(lRA).toBe('------GAGGTGTA');
});
