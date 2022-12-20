# Building Documentation

## Building documentation locally

From the ADRIA project directory:

```bash
$ cd docs
$ julia --project=. make.jl
```

Locally generated documentation can be found under `docs/build`. Open the `index.html` file with any web browser.


## Documentation deployment

Documentation is hosted on [GitHub Pages](https://pages.github.com/) via [GitHub Actions](https://github.com/features/actions).

Configuration is found [here](https://github.com/open-AIMS/ADRIA.jl/blob/main/.github/workflows/documentation.yml).

Documentation is automatically built and deployed:

- When a PR targeting `main` is submitted  
  In this case, a preview URL is created: e.g., a URL with `previews/PR###` at the end, where `PR###` refers to the PR ID.
- On commit/merge to `main`  
  In this case the main documentation website is updated