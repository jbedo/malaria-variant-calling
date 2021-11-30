{ bionix }:

with bionix;
with lib;

let
  genmap = pkgs.callPackage ./genmap.nix { };

in
rec {
  index = exec' (ref: stage {
    name = "genmap-index";
    buildInputs = [ genmap ];
    buildCommand = ''
      genmap index -F ${ref} -I $out
    '';
  });

  calcmap = exec
    ({ k ? 30, e ? 2 }: ref: stage {
      name = "genmap-mappability";
      buildInputs = [ genmap ];
      buildCommand = ''
        mkdir $out
        genmap map -K ${toString k} -E ${toString e} --wig -I ${index ref} -O $out -T $NIX_BUILD_CORES
      '';
      passthru.multicore = true;
    });

  toBW = exec
    (_: wig: stage {
      name = "wig2bigwig.bw";
      buildInputs = [ pkgs.kent ];
      buildCommand = ''
        wigToBigWig ${wig}/*.{wig,chrom.sizes} $out
      '';
    });
}
