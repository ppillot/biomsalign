import { pairwiseAlignment } from '../align';
import { makeSequence } from '../../sequence/sequence'
import { setAlignmentParameters } from '../params';

test('aligns ends', () => {
    const lSeqA = makeSequence('CCCAATTCTATACCAACAC');
    const lSeqB = makeSequence('CCCAATTCTA');
    setAlignmentParameters();
    const lResult = pairwiseAlignment(lSeqA, lSeqB);

    expect(lResult.alignment[1]).toBe('CCCAATTCTA---------');
});

test('aligns start', () => {
    const lSeqA = makeSequence(  'AGCCCTTA');
    const lSeqB = makeSequence('ACTGCCCTCA');
    setAlignmentParameters();
    const lResult = pairwiseAlignment(lSeqA, lSeqB);

    expect(lResult.alignment[0]).toBe('--AGCCCTTA');
});
