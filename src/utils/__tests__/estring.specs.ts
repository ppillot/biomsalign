import { estringTransform } from '../estring';

test('Transforms a string', () => {
    const lTxt = 'MQTIF';
    const lEString = [3,-1,2];
    let lTxt2 = estringTransform(lTxt, lEString);
    expect(lTxt2).toBe('MQT-IF');
});
