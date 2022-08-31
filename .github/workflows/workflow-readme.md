# Actions README

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
  state new_release <<choice>>
  [*] --> version
  version --> new_release
  new_release --> deb: new release
  new_release --> [*]: no new release
  deb --> doc
  doc --> pre: staging
  pre --> [*]

  doc --> pub: main
  pub --> debPush
  debPush --> docPush
  docPush --> [*]




```
