import { createApp } from 'vue'

import App from './App.vue'

// Vuetify3
import 'vuetify/styles'
import { createVuetify } from 'vuetify'
import * as components from 'vuetify/components'
import * as directives from 'vuetify/directives'
import '@mdi/font/css/materialdesignicons.min.css'

// D3
import { select } from 'd3-selection';
import 'd3-transition';

// EasyDataTable
import Vue3EasyDataTable from 'vue3-easy-data-table'
import 'vue3-easy-data-table/dist/style.css'

const vuetify = createVuetify({
  components,
  directives,
  theme: {
    defaultTheme: 'light'
  },
})

import './assets/main.css'

const app = createApp(App).use(vuetify)

// D3 axes
app.directive('axis', (el, binding) => {
    const axisMethod = binding.value
    select(el)
      .transition()
      .call(axisMethod)
})

// EasyDataTable
app.component('EasyDataTable', Vue3EasyDataTable)

app.mount('#app')
