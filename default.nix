{ bionix ? import <bionix> { }
, small ? false
, chunkSize ? 10
}:

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
  samples' = ena.importFromJSON ./PRJEB2140.json;
  samples = if small then takeAttrs 9 samples' else samples';
  takeAttrs = n: xs: listToAttrs (zipListsWith nameValuePair (take n (attrNames xs)) (attrValues xs));

  QDNAseq = callBionix ./QDNAseq.nix { inherit ref; targets = [ PMIX PMX ]; };

  preprocess = id: flip pipe [
    (bwa.align { inherit ref; RG = { ID = id; SM = id; }; })
    (samtools.sort { })
  ];

  chunkFold1 = n: f: xs:
    if length xs <= n then
      f xs
    else
      chunkFold1 n f (chunkFold1' n f xs);
  chunkFold1' = n: f: xs:
    if length xs <= n then
      xs
    else
      [ (f (take n xs)) ] ++ chunkFold1' n f (drop n xs);

in
linkOutputs {
  cnv = mapAttrs'
    (name: flip pipe [
      (preprocess name)
      (QDNAseq.call { inherit name; })
      (nameValuePair "${name}.rds")
    ])
    samples;
  "variants.vcf.gz" = pipe samples [
    (mapAttrsToList preprocess)
    (map (gatk.callHaplotype { targets = [ PMIX PMX ]; }))
    (chunkFold1 chunkSize (gatk.merge { }))
    (gatk.callGenotypes { })
  ];
}
