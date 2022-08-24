# Actions README

Note: Each step requires the previous step(s) to complete (unless otherwise noted).
_This is deprecated and replaced with the below_

```mermaid
  flowchart TB
    subgraph build [Build, Push, Test]
    Push -->  sem[Semantic Version]
    sem --> deb[Build Debian Package]
    deb -->
    docBuild[Build Docker Image] & vt[Validation Testing]
    vt & docBuild --> prerel[Pre-release based on branch]
    end

    prerel -- if main branch --> pubdep
    prerel -- else --> pridep

    subgraph pubdep [Public Deployment]

    bgh[Release on Blastfoam GitHub]
    dh[Upload to DockerHub]
    dp[Upload Debian Package]
    end

    subgraph pridep [Private Deployment]
    dgh[Upload Docker Image to Private Repo]
    dd[Upload Debian Package to Private Repo]

    end
```

## On pull request build/merge/sync

```mermaid
  stateDiagram-v2
  state "cache hit" as ch
  [*] --> ch
  ch --> build: false
  build --> validation
  ch --> validation: true
  validation --> [*]
```

## Continuous Distribution

_Run on push to staging, master, dev_

```mermaid
  stateDiagram-v2

  state "build debian package" as deb
  state "build docker image" as doc
  state "private github prerelease" as pre
  state "public github release" as pub
  state "push docker image" as docPush
  state "push debian package" as debPush

  [*] --> version
  version --> deb
  deb --> doc
  doc --> pre: staging
  pre --> [*]

  doc --> pub: main
  pub --> debPush
  debPush --> docPush
  docPush --> [*]




```
