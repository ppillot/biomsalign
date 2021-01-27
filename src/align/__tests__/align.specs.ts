import { pairwiseAlignment } from '../align';
import { makeSequence } from '../../sequence/sequence'
import { setAlignmentParameters } from '../params';

test('aligns ends', () => {
    const lSeqA = makeSequence('CCCAATTCTATACCAACAC');
    const lSeqB = makeSequence('CCCAATTCTA');
    setAlignmentParameters();
    const lResult = pairwiseAlignment(lSeqA, lSeqB, [0], [1]);

    expect(lResult.alignment[1]).toBe('CCCAATTCTA---------');
});
