import { defineConfig } from 'vitepress'
import { tabsMarkdownPlugin } from 'vitepress-plugin-tabs'
import mathjax3 from "markdown-it-mathjax3";
import footnote from "markdown-it-footnote";
import path from 'path'

function getBaseRepository(base: string): string {
  if (!base || base === '/') return '/';
  const parts = base.split('/').filter(Boolean);
  return parts.length > 0 ? `/${parts[0]}/` : '/';
}

const baseTemp = {
  base: '/ADRIA.jl/',// TODO: replace this in makedocs!
}

const navTemp = {
  nav: [
{ text: 'Introduction', link: '/index' },
{ text: 'Concepts', collapsed: false, items: [
{ text: 'Dynamic Multi-Criteria Decision Analysis', link: '/concepts/dMCDA' },
{ text: 'Disturbances', link: '/concepts/disturbances' }]
 },
{ text: 'Usage', collapsed: false, items: [
{ text: 'Getting Started', link: '/usage/getting_started' },
{ text: 'Loading a Domain', link: '/usage/loading_a_domain' },
{ text: 'Loading Results', link: '/usage/loading_results' },
{ text: 'Generating scenarios', link: '/usage/generating_scenarios' },
{ text: 'Running scenarios', link: '/usage/scenario_runs' },
{ text: 'Scenario Discovery', link: '/usage/scenario_discovery' },
{ text: 'Analysis', link: '/usage/analysis' },
{ text: 'Cookbook examples', link: '/usage/cookbook' }]
 },
{ text: 'Architecture', collapsed: false, items: [
{ text: 'Architectural overview', link: '/architecture/architecture' },
{ text: 'Inputs and Outputs', link: '/architecture/domain_and_resultsets' }]
 },
{ text: 'Development', collapsed: false, items: [
{ text: 'Development setup', link: '/development/development_setup' },
{ text: 'Contributing a metric', link: '/development/metrics' },
{ text: 'Release Guide', link: '/development/release_guide' },
{ text: 'Building Documentation', link: '/development/building_docs' }]
 },
{ text: 'ADRIA API', link: '/API' }
]
,
}

const nav = [
  ...navTemp.nav,
  {
    component: 'VersionPicker'
  }
]

// https://vitepress.dev/reference/site-config
export default defineConfig({
  base: '/ADRIA.jl/',// TODO: replace this in makedocs!
  title: 'ADRIA.jl',
  description: 'Documentation for ADRIA.jl',
  lastUpdated: true,
  cleanUrls: true,
  outDir: '../final_site', // This is required for MarkdownVitepress to work correctly...
  head: [
    ['link', { rel: 'icon', href: '/favicon.ico' }],
    ['script', {src: `${getBaseRepository(baseTemp.base)}versions.js`}],
    // ['script', {src: '/versions.js'], for custom domains, I guess if deploy_url is available.
    ['script', {src: `${baseTemp.base}siteinfo.js`}]
  ],
  
  vite: {
    resolve: {
      alias: {
        '@': path.resolve(__dirname, '../components')
      }
    },
    optimizeDeps: {
      exclude: [ 
        '@nolebase/vitepress-plugin-enhanced-readabilities/client',
        'vitepress',
        '@nolebase/ui',
      ], 
    }, 
    ssr: { 
      noExternal: [ 
        // If there are other packages that need to be processed by Vite, you can add them here.
        '@nolebase/vitepress-plugin-enhanced-readabilities',
        '@nolebase/ui',
      ], 
    },
  },
  markdown: {
    math: true,
    config(md) {
      md.use(tabsMarkdownPlugin),
      md.use(mathjax3),
      md.use(footnote)
    },
    theme: {
      light: "github-light",
      dark: "github-dark"}
  },
  themeConfig: {
    outline: 'deep',
    logo: { src: '/logo.png', width: 24, height: 24},
    search: {
      provider: 'local',
      options: {
        detailedView: true
      }
    },
    nav,
    sidebar: [
{ text: 'Introduction', link: '/index' },
{ text: 'Concepts', collapsed: false, items: [
{ text: 'Dynamic Multi-Criteria Decision Analysis', link: '/concepts/dMCDA' },
{ text: 'Disturbances', link: '/concepts/disturbances' }]
 },
{ text: 'Usage', collapsed: false, items: [
{ text: 'Getting Started', link: '/usage/getting_started' },
{ text: 'Loading a Domain', link: '/usage/loading_a_domain' },
{ text: 'Loading Results', link: '/usage/loading_results' },
{ text: 'Generating scenarios', link: '/usage/generating_scenarios' },
{ text: 'Running scenarios', link: '/usage/scenario_runs' },
{ text: 'Scenario Discovery', link: '/usage/scenario_discovery' },
{ text: 'Analysis', link: '/usage/analysis' },
{ text: 'Cookbook examples', link: '/usage/cookbook' }]
 },
{ text: 'Architecture', collapsed: false, items: [
{ text: 'Architectural overview', link: '/architecture/architecture' },
{ text: 'Inputs and Outputs', link: '/architecture/domain_and_resultsets' }]
 },
{ text: 'Development', collapsed: false, items: [
{ text: 'Development setup', link: '/development/development_setup' },
{ text: 'Contributing a metric', link: '/development/metrics' },
{ text: 'Release Guide', link: '/development/release_guide' },
{ text: 'Building Documentation', link: '/development/building_docs' }]
 },
{ text: 'ADRIA API', link: '/API' }
]
,
    editLink: { pattern: "https://github.com/open-AIMS/ADRIA.jl/edit/main/docs/src/:path" },
    socialLinks: [
      { icon: 'github', link: 'https://github.com/open-AIMS/ADRIA.jl' }
    ],
    footer: {
      message: 'Made with <a href="https://luxdl.github.io/DocumenterVitepress.jl/dev/" target="_blank"><strong>DocumenterVitepress.jl</strong></a><br>',
      copyright: `© Copyright ${new Date().getUTCFullYear()}.`
    }
  }
})
