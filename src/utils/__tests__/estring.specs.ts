import { estringTransform, estringProduct, estringCat, epath2estring, estringMerge, estringDifference } from '../estring';

test('Transforms a string', () => {
    const lTxt = 'MQTIF';
    const lEString = [3,-1,2];
    let lTxt2 = estringTransform(lTxt, lEString);
    expect(lTxt2).toBe('MQT-IF');
});

test('Multiplies estrings', () => {
    const lTxt = 'MQTIFAGH';
    const lEString1 = [3,-1,2,-1,3];
    const lEString2 = [5,-2,2,-2,1];
    const lProduct = estringProduct(lEString2, lEString1);
    const lTxt2 = estringTransform(lTxt, lProduct);
    let lTxt3 = estringTransform(
        estringTransform(lTxt, lEString1),
        lEString2
    );
    expect(lTxt2).toBe(lTxt3);
});

test('Catenates estrings', () => {
    const lEString1 = [3,-1,3];
    const lEString2 = [5,-2,1];
    const lEString3 = [-5,2];
    let lCat = estringCat(lEString1, lEString2);
    expect(lCat).toEqual([3,-1,8,-2,1]);
    lCat = estringCat(lEString1, lEString3);
    expect(lCat).toEqual([3,-1,3,-5,2]);
});

test('Builds estrings from edit paths', () => {
    const lEpath = [-1, -1, -1, 1, 1, 1, -1, 1, 4];
    let lEstring = epath2estring(lEpath);
    expect(lEstring).toEqual([5, -1, 3, -3]);
});

test('Merges estrings', () => {
    const lEString1 = [2,-1,2,-1,1];
    const lEString2 = [2,-1,3,-1];
    const lEString3 = [3,-4,2];
    const lEString4 = [2,-2,3];
    let lMerge = estringMerge(lEString1, lEString2);
    expect(lMerge).toEqual([2,-1,2,-1,1,-1]);
    lMerge = estringMerge(lEString3, lEString4);
    expect(lMerge).toEqual([2,-2,1,-4,2]);
});

test('Diffs estrings', () => {
    const lEString3 = [2,-1,2,-1,1,-1];
    const lEString1 = [2,-1,2,-1,1];
    const lEString2 = [2,-1,3,-1];
    let lDiff = estringDifference(lEString3, lEString1);
    expect(lDiff).toEqual([7,-1]);
    lDiff = estringDifference(lEString3, lEString2);
    expect(lDiff).toEqual([5,-1,2]);
});
