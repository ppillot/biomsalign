import { DEQueue } from '../queue';

test('Inits a queue', () => {
    const lQueue = new DEQueue(10);
    expect(lQueue.size).toBe(0);
    expect(lQueue.isEmpty).toBe(true);
});

test('Pushes at the head of the queue', () => {
    const lQueue = new DEQueue<number>(10);
    lQueue.pushHead(1);
    const lHead = lQueue.getHead();
    expect(lHead).toBe(1);
});

test('Retrieves items at the tail of the queue', () => {
    const lQueue = new DEQueue<number>(10);
    lQueue.pushHead(1);
    lQueue.pushHead(2);
    const lTail = lQueue.getTail();
    expect(lTail).toBe(1);
});

test('Pushes at the tail of the queue', () => {
    const lQueue = new DEQueue<number>(10);
    lQueue.pushTail(1);
    const lTail = lQueue.getTail();
    expect(lTail).toBe(1);
});

test('Retrieves items pushed at the tail, at the head of the queue', () => {
    const lQueue = new DEQueue<number>(10);
    lQueue.pushTail(1);
    lQueue.pushTail(2);
    const lHead = lQueue.getHead();
    expect(lHead).toBe(1);
});

test('Pops items from the head of the queue', () => {
    const lQueue = new DEQueue<number>(10);
    lQueue.pushHead(1);
    lQueue.pushHead(2);
    lQueue.pushHead(3);
    let lHead = lQueue.popHead();
    expect(lHead).toBe(3);
    expect(lQueue.size).toBe(2);
    expect(lQueue.getHead()).toBe(2);
});

test('Pops items from the tail of the queue', () => {
    const lQueue = new DEQueue<number>(10);
    lQueue.pushHead(1);
    lQueue.pushHead(2);
    lQueue.pushHead(3);
    let lTail = lQueue.popTail();
    expect(lTail).toBe(1);
    expect(lQueue.size).toBe(2);
    expect(lQueue.getTail()).toBe(2);
});
