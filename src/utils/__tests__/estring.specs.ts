import { estringTransform, estringProduct, estringCat } from '../estring';

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