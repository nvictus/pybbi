# Upgrade and patch the UCSC source

To upgrade the vendored UCSC tools C source code, we need to apply two small patches to expose the private functions `bbiChromId`, `bbiSummarySlice`, and `bbiIntervalSlice` to [Cython](https://github.com/nvictus/pybbi/blob/master/bbi/cbbi.pxd) so we can implement querying routines to populate Python/NumPy data structures.

1. Select the latest stable release from https://github.com/ucscGenomeBrowser/kent/releases.
2. Replace the value of `VERSION` in `patch_source.sh`.
3. Run `patch_source.sh` to download and extract the library code as `src/` and header files as `include/`, and apply the patches.
4. If all succeeds, replace the repo's `src/` and `include/` directories with these new ones.

If nothing dramatically changed, the patch should work, our `cbbi.pxd` declaration file should remain valid, and everything should compile.
