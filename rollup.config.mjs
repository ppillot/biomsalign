import typescript from '@rollup/plugin-typescript'
import pkg from './package.json' assert { type: 'json' }
import terser from '@rollup/plugin-terser';

export default {
    input: 'src/index.ts',
    output: [
      {
        file: 'build/biomsa.js',
        format: 'umd',
        name: 'biomsa',
        sourcemap: true
      },
      {
        file: 'build/biomsa.esm.js',
        format: 'es',
        sourcemap: true
      },
    ],
    external: [
      ...Object.keys(pkg.dependencies || {}),
      ...Object.keys(pkg.peerDependencies || {}),
    ],plugins: [
        typescript(),
        (process.env.NODE_ENV === 'production'
          && terser({
            toplevel: true,
            format: {comments: false},
            compress: {ecma: 2015, passes: 3, unsafe: true},
            mangle: {
              properties: {
                reserved: ['biomsa', 'gapopen', 'gapextend', 'matrix', 'method', 'type', 'gapchar', 'debug']
              }
            },
            nameCache: {}
          })
        )
      ],
}