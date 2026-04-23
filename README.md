# Alpha-beta-methylome-app

A web application tool for exploration and visualization of DNA methylation- and gene expression data presented in [Cell-specific DNA methylation in human
alpha and beta cells regulates gene
expression in type 2 diabetes](https://doi.org/10.1038/s42255-026-01498-9).
  
Only the code, Dockerfiles, and r-environment requirments to run the app is provided in this repo. 


#### File tree structure

    ├── Dockerfile
    ├── renv.lock
    └ app
        ├── Data
        │   ├── Annotations
        │   ├── Expression
        │   │   ├── alpha
        │   │   │   ├── Age
        │   │   │   ├── Sex
        │   │   │   └── ctrl_vs_T2D
        │   │   ├── alpha_vs_beta
        │   │   └── beta
        │   │       ├── Age
        │   │       ├── Sex
        │   │       └── ctrl_vs_T2D
        │   ├── Results_10
        │   │   ├── alpha
        │   │   │   ├── Age
        │   │   │   ├── Sex
        │   │   │   └── ctrl_vs_T2D
        │   │   ├── alpha_vs_beta
        │   │   └── beta
        │   │       ├── Age
        │   │       ├── Sex
        │   │       └── ctrl_vs_T2D
        │   └── Results_5
        │       ├── alpha
        │       │   ├── Age
        │       │   ├── Sex
        │       │   └── ctrl_vs_T2D
        │       ├── alpha_vs_beta
        │       └── beta
        │           ├── Age
        │           ├── Sex
        │           └── ctrl_vs_T2D
        └── Scripts
