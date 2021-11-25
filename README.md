# Malaria (Plasmodium Vivax) example pipeline

This is an example [BioNix](https://github.com/PapenfussLab/BioNix) pipeline for
variant calling in a target region over 5k publicly available Plasmodium Vivax
sequencing samples. The data processing aligns each sample against the reference
genome, sorts the reads, then calls variants on the population using
[GATK](https://github.com/broadinstitute/gatk).

Notable design choices:

- there are a lot of samples, so `fetchurl` is overridden in nixpkgs to prevent
  substitutions, which avoids querying the cache for each sequence;
- similarly, BioNix's `stage` is overridden to prevent substitutions.

## Building/Executing

This repository is configured as a Nix
[flake](https://nixos.wiki/wiki/Flakes#Installing_nix_flakes) and can be built
with:

``` sh
nix build github:jbedo/malaria-variant-calling
```

## Slurm cluster execution

The repository is also setup to exercise experimental [Slurm support
patches](https://github.com/jbedo/static-nix). If building on slurm then using
nix built with the patch should be sufficient to submit the jobs to the queue.
