import { pairwiseAlignment } from '../align';
import { makeSequence } from '../../sequence/sequence'
import { getAlignmentParameters, SEQUENCE_TYPE } from '../params';

test('aligns ends', () => {
    const lSeqA = makeSequence('CCCAATTCTATACCAACAC');
    const lSeqB = makeSequence('CCCAATTCTA');
    const lParam = getAlignmentParameters(SEQUENCE_TYPE.NUCLEIC);
    const lResult = pairwiseAlignment(lSeqA, lSeqB, lParam);

    expect(lResult.alignment[1]).toBe('CCCAATTCTA---------');
});

test('aligns start', () => {
    const lSeqA = makeSequence(  'AGCCCTTA');
    const lSeqB = makeSequence('ACTGCCCTCA');
    const lParam = getAlignmentParameters(SEQUENCE_TYPE.NUCLEIC);
    const lResult = pairwiseAlignment(lSeqA, lSeqB, lParam);

    expect(lResult.alignment[0]).toBe('--AGCCCTTA');
});
