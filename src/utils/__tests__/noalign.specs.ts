import { extractMinimizers, TMinimizer, dnaHammingDistance } from '../noalign';
import { SEQUENCE_TYPE } from '../params';

test('Extracts Minimizers from sequence', () => {
    const lSeq = {
        rawSeq:'',
        compressedSeq: [],
        encodedSeq: [3,2,1,1,0,1,1,0,3,2,1,0,2,3,2,1,0,0,2,3,2,0,1,2],
        type: SEQUENCE_TYPE.NUCLEIC
    }
    const lMin = (1 << 12) + (1 << 10) + (3 << 6) + (2 << 4) + (1 << 2);
    const [lMinzMap, _] = extractMinimizers(lSeq, 8, 12);

    expect(lMinzMap.has(lMin));
    const lMinz0 = (lMinzMap.get(lMin) as TMinimizer[])[0];
    expect(lMinz0.minimizedSubarray).toEqual([3,2,1,1,0,1,1,0,3,2,1,0,2,3,2,1])
    expect(lMinz0.winPos).toBe(0);
});

test('Get Hamming distance from a pair of Kmers', () => {
    const lKmerA = 0b0111001011001001;
    const lKmerB = 0b0101001100111001;

    expect(dnaHammingDistance(lKmerA, lKmerB)).toBe(4);
});
