import { SEQUENCE_TYPE } from '../params';
import { getSequenceType, encodeSeqToNum } from '../sequence';

test('Detects DNA', () => {
    const seq = 'ATCG';
    expect(getSequenceType(seq)).toBe(SEQUENCE_TYPE.NUCLEIC);
});

test('Detects RNA', () => {
    const seq = 'AUCG';
    expect(getSequenceType(seq)).toBe(SEQUENCE_TYPE.NUCLEIC);
});

test('Detects Proteic', () => {
    const seq = 'ATMVCIWCG';
    expect(getSequenceType(seq)).toBe(SEQUENCE_TYPE.PROTEIN);
});

test('Throws error when cannot detect sequence type', () => {
    const seq = 'BOU';
    expect(() => {
        getSequenceType(seq);
    }).toThrowError();
});

test('Encodes proteic sequence', () => {
    const seq = 'ACDEFWY';
    expect(encodeSeqToNum(seq, SEQUENCE_TYPE.PROTEIN)).toEqual([0, 1, 2, 3, 4, 18, 19]);
});
