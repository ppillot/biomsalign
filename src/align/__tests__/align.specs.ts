import { pairwiseAlignment } from '../align';
import { makeSequence } from '../../sequence/sequence'
import { getAlignmentParameters, SEQUENCE_TYPE } from '../params';
import { estringTransform } from '../../utils/estring';

test('aligns ends', () => {
    const lSeqA = makeSequence('CCCAATTCTATACCAACAC');
    const lSeqB = makeSequence('CCCAATTCTA');
    const lParam = getAlignmentParameters(SEQUENCE_TYPE.NUCLEIC);
    const lResult = pairwiseAlignment(lSeqA, lSeqB, lParam);
    const lRB = estringTransform(lSeqB.rawSeq, lResult.estrings[1]);
    expect(lRB).toBe('CCCAATTCTA---------');
});

test('aligns start', () => {
    const lSeqA = makeSequence(  'AGCCCTTA');
    const lSeqB = makeSequence('ACTGCCCTCA');
    const lParam = getAlignmentParameters(SEQUENCE_TYPE.NUCLEIC);
    const lResult = pairwiseAlignment(lSeqA, lSeqB, lParam);
    const lRA = estringTransform(lSeqA.rawSeq, lResult.estrings[0]);
    expect(lRA).toBe('--AGCCCTTA');
});
