{
  description = "Variant calling for Plasmodium Vivax";
  inputs = {
    bionix.url = "github:papenfusslab/bionix";
    flake-utils.url = "github:numtide/flake-utils";
    nixpkgs.url = "nixpkgs";
  };

  outputs = { self, nixpkgs, bionix, flake-utils }:
    flake-utils.lib.eachDefaultSystem
      (system: with bionix.lib
        {
          nixpkgs = import nixpkgs {
            inherit system; overlays = [
            (_: super: with super; {
              fetchurl = x: (fetchurl x).overrideAttrs (_: {
                MEMORY = "100M";
                PPN = "1";
                WALLTIME = "10:00:00";
                allowSubstitutes = false;
              });
            })
          ];
          };
          overlays = [
            (self: super: {
              exec = f: x@{ ppn ? 1, mem ? 1, walltime ? "2:00:00", ... }: y:
                (f (removeAttrs x [ "ppn" "mem" "walltime" ]) y).overrideAttrs (attrs: {
                  PPN = if attrs.passthru.multicore or false then ppn else 1;
                  MEMORY = toString mem + "G";
                  WALLTIME = walltime;
                  preferLocalBuild = true;
                  allowSubstitutes = false;
                });
            })

            ./resources.nix
          ];
        }; {
        defaultPackage = callBionix ./. { };
      });
}
