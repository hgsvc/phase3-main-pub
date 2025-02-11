// vite config for serverless deployment (i.e., files only)

import { defineConfig } from 'vite'
import vue from '@vitejs/plugin-vue'
import { fileURLToPath, URL } from 'node:url'
import { viteSingleFile } from "vite-plugin-singlefile"

// https://vitejs.dev/config/
export default defineConfig({
  plugins: [vue(), viteSingleFile()],
  build: {
    cssCodeSplit: false,
    assetsInlineLimit: 100000000,
  },
  resolve: {
    alias: [
      { find: '@', replacement: '/src' },
      { find: 'views', replacement: '/src/views' },
      { find: 'components', replacement: '/src/components' },
    ]
  },
  server: {
    open: true,
  }
})
