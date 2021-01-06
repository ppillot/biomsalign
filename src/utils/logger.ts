/**
 * @file logger
 * @author Paul Pillot <paul.pillot@libmol.org>
 * @license MIT
 * @copyright 2020
 */

let gEvents: Map<string, number>

function start () {
    gEvents = new Map();
    add('START');
}

function add (name: string) {
    if (gEvents.has(name)) {
        add(name + '_');
        return;
    }
    gEvents.set(name, performance.now());
}

function summary () {
    let lSummary: {[k: string]: number} = {
        'start': 0
    };

    add ('END');

    let lPrevTime = 0;
    gEvents.forEach((val, key) => {
        if (key === 'START') {
            lPrevTime = val;
            return;
        }
        lSummary[key] = val - lPrevTime;
        lPrevTime = val;
    });

    console.table(lSummary);
}

export default {
    start,
    add,
    summary,
}
