# Repository development instructions

Before changing any repository script or layout:

1. read `REPOSITORY_DESIGN.md` in full, including **Scientific-script coding standard**;
2. read the README of the module being changed;
3. inspect the current entry point and follow its complete data flow before proposing abstractions;
4. check `git status` and preserve unrelated work;
5. state which scientific behavior must remain unchanged and how it will be verified.

This repository contains publication data-analysis scripts, not a general-purpose
application framework. Prefer direct, sequential analysis code. Do not introduce
runner, adapter, registry, provider, manager, workflow, service, class hierarchy,
or generic configuration layers unless a present, concrete scientific or file-format
variation requires them.

When refactoring, count the final files, functions, and lines. The default outcome
must be the same size or smaller. Any increase must be justified by a new scientific
capability, a real external-data boundary, or a regression check that protects a
publication claim.
