{ bionix ? import <bionix> { } }:

with bionix;
with lib;
with types;

let
  ref = fetchFastA {
    url = "https://plasmodb.org/common/downloads/release-54/PvivaxP01/fasta/data/PlasmoDB-54_PvivaxP01_Genome.fasta";
    sha256 = "sha256-1ZahI54i1iZ3kvJGzEz00Gv3l9BgzHL+1gJ1I+1wA3k=";
  };
  PMIX = "PvP01_13_v2:876959-881753";
  PMX = "PvP01_01_v2:555499-559411";

  ena = pkgs.callPackage ./ena.nix { inherit fetchFastQGZ; };
  samples = ena.importFromJSON ./PRJEB2140.json;

  QDNAseq = callBionix ./QDNAseq.nix { inherit ref; targets = [ PMIX PMX ]; };

  preprocess = id: flip pipe [
    (bwa.align { inherit ref; flags = "-R'@RG\\tID:${id}\\tSM:${id}'"; })
    (samtools.sort { })
  ];

  callSNP = flip pipe [
    attrValues
    (map (gatk.callHaplotype { targets = [ PMIX PMX ]; }))
    (gatk.merge { })
    (gatk.callGenotypes { })
  ];

in
linkOutputs {
  cnv = pipe samples [ (mapAttrs preprocess) (mapAttrs (_: QDNAseq.call { })) (mapAttrs' (n: nameValuePair "${n}.rds")) ];
  "variants.vcf.gz" = pipe samples [ (mapAttrs preprocess) callSNP ];
}
