# Actions README

## On pull request build/merge/sync

```mermaid
  stateDiagram-v2
  %% Steps
  cbf: Checkout blastFOAM
  tar: Build Tar for Cache
  cc: Check Cache for Match
  make: Build from Makefile
  art: Create Artifact
  da: Download Artifact
  inst: Install
  cbf2: Checkout blastFOAM
  run: Run Tests
  pl: Upload Validation Plots
  uninst: Uninstall

  %% Jobs
  val: run-validation-tests
  state build {
    Direction LR
    [*] --> cbf
    cbf --> tar
    tar --> cc

    cc --> make: false
    make --> art
    art --> [*]

    cc --> [*]: true
  }
  state val {
    Direction LR
    [*] --> da
    da --> inst
    inst --> cbf2
    cbf2 --> run
    run --> pl
    pl --> uninst
    uninst --> [*]

  }

  %% Workflows
  state "pull-request-check.yaml" as prc
  state prc {
  [*] --> build
  build --> val
  val --> [*]
  }
```

## Continuous Distribution

_Run on push to staging, master, dev_

```mermaid
  stateDiagram-v2

%% Steps
dep: Install dependencies
vers: Version (actions/composites/semantic-versioning)

cbf: Checkout blastFOAM
tar: Build Tar for Cache
cc: Check Cache for Match
make: Build from Makefile
art: Create Artifact

cbf2: Checkout blastFOAM
da: Download Artifact
rel: Add Artifact to Versioned Release

cbf3: Checkout blastFOAM
da2: Download Artifact
dl: Docker Login to GHCR
bpi: Docker Build and Push Image


%% Jobs
state version{
  [*] --> dep
  dep --> vers
  vers --> [*]
}
state build {
    [*] --> cbf
    cbf --> tar
    tar --> cc

    cc --> make: false
    make --> art
    art --> [*]

    cc --> [*]: true
  }
state release_deb {
  [*] --> cbf2
  cbf2 --> da
  da --> rel
  rel --> [*]
}
state build_docker {
  [*] --> cbf3
  cbf3 --> da2
  da2 --> dl
  dl --> bpi
  bpi --> [*]
}
%% Workflows
state deployment {
  direction LR
  [*] --> version
  version --> build: New Version
  build --> release_deb
  release_deb --> build_docker
  build_docker --> [*]
  version --> [*]: No New Version
}

state public_deployment {
  direction LR
  state code {
    cbf4: Checkout blastFOAM
    pcp: Push Code Publically
    crp: Push Release Publically
    [*] --> cbf4
    cbf4 --> pcp
    pcp --> crp
    crp --> [*]
  }
  state docker {
    ghcr: Login to GitHub Registry
    dh: Login to DockerHub
    dt: Construct Dev Tag
    pdi: Pull Docker Image
    tags: Create Tags
    ti: Tag Images
    pi: Push Images

    [*]--> ghcr
    ghcr --> dh
    dh --> dt
    dt --> pdi
    pdi --> tags
    tags --> ti
    ti --> pi
    pi --> [*]
  }
  [*] --> code
  code --> docker
  docker --> [*]
}
[*] --> deployment
deployment --> public_deployment
public_deployment --> [*]
```
