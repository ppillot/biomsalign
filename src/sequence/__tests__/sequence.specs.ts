import { SEQUENCE_TYPE } from '../../align/params';
import { getSequenceType, encodeSeqToNum, sortMSA } from '../sequence';

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

test('Detects extended amino acids codes', () => {
    const seq = 'ATMXZJ';
    expect(getSequenceType(seq)).toBe(SEQUENCE_TYPE.PROTEIN);
});

test('Estimates likelyhood of nucleic type', () => {
    const seq = 'ATGCNTGAATGTACCATACGGTA';
    expect(getSequenceType(seq)).toBe(SEQUENCE_TYPE.NUCLEIC);
});

test('Throws error when cannot detect sequence type', () => {
    const seq = 'BOU';
    expect(() => {
        getSequenceType(seq);
    }).toThrowError();
});

test('Encodes proteic sequence', () => {
    const seq = 'ACDEFWY';
    expect(Array.from(encodeSeqToNum(seq, SEQUENCE_TYPE.PROTEIN))).toEqual([0, 1, 2, 3, 4, 18, 19]);
});

test('Sorts MSA', () => {
    const lMSA = ['a', 'b', 'c', 'd'];
    const lOrder = [3, 0, 1, 2];
    expect(sortMSA(lMSA, lOrder)).toEqual(['b', 'c', 'd', 'a']);
})